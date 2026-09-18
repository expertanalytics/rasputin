#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <terrain/cdt/detria_backend.hpp>
#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>

#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace py = pybind11;

using terrain::Chain;
using terrain::ChainRole;
using terrain::IndexedMesh2;
using terrain::Point2;
using terrain::Point3;
using terrain::Pslg;
using terrain::PslgBuilder;
using terrain::PslgBuildResult;
using terrain::PslgDiagnostic;
using terrain::PslgError;
using terrain::TriangleIndices;
using terrain::cdt::CdtOptions;
using terrain::cdt::CdtOutcome;
using terrain::cdt::CdtStatus;
using terrain::cdt::DetriaBackend;

namespace {

// A zero-copy numpy view over a core buffer: read-only, and anchored to the
// Python object that owns the memory. Both halves are load-bearing. Without the
// base object the view is a use-after-free the moment its owner is dropped
// (06-cdt-viewer.md risk 4); without clearing the writeable flag, Python writes
// straight through into a type whose whole contract is immutability.
template <typename T>
[[nodiscard]] py::array_t<T> readonly_view(const py::object& owner, const T* data,
                                           std::vector<py::ssize_t> shape,
                                           std::vector<py::ssize_t> strides) {
    // An empty std::vector may hand out a null data pointer, and numpy responds
    // to a null pointer by ignoring the base object and copying. Point an empty
    // view at an address it will never read rather than at nothing.
    static constexpr double kEmpty[1]{};
    const T* ptr = data != nullptr ? data : reinterpret_cast<const T*>(&kEmpty[0]);
    py::array_t<T> view{std::move(shape), std::move(strides), ptr, owner};
    py::detail::array_proxy(view.ptr())->flags &= ~py::detail::npy_api::NPY_ARRAY_WRITEABLE_;
    return view;
}

[[nodiscard]] py::array_t<double> point_view(const py::object& owner,
                                             std::span<const Point2> points) {
    static_assert(sizeof(Point2) == 2 * sizeof(double) && std::is_standard_layout_v<Point2>,
                  "the (N, 2) float64 view reinterprets the Point2 buffer as coordinate pairs");
    return readonly_view<double>(owner, reinterpret_cast<const double*>(points.data()),
                                 {static_cast<py::ssize_t>(points.size()), 2},
                                 {static_cast<py::ssize_t>(sizeof(Point2)),
                                  static_cast<py::ssize_t>(sizeof(double))});
}

[[nodiscard]] py::array_t<std::uint32_t> index_view(const py::object& owner,
                                                    std::span<const std::uint32_t> indices) {
    return readonly_view<std::uint32_t>(owner, indices.data(),
                                        {static_cast<py::ssize_t>(indices.size())},
                                        {static_cast<py::ssize_t>(sizeof(std::uint32_t))});
}

// TypeError for a path, ValueError for a mis-shaped array: the first is a
// category error at the I/O boundary -- file decoding is Python's and the core
// never sees a path -- and the second is a well-typed argument with wrong data.
[[nodiscard]] std::vector<Point2> as_points(const py::object& vertices) {
    if (py::isinstance<py::str>(vertices) || py::isinstance<py::bytes>(vertices) ||
        py::hasattr(vertices, "__fspath__")) {
        throw py::type_error("vertices must be coordinates, not a path or a filename: decode "
                             "the file in Python and hand the core an (N, 2) array");
    }
    const auto array =
        py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(vertices);
    if (!array) {
        throw py::type_error("vertices must be convertible to a float64 array of shape (N, 2)");
    }
    if (array.ndim() != 2 || array.shape(1) != 2) {
        throw py::value_error(
            std::format("vertices must have shape (N, 2), got a {}-dimensional array",
                        array.ndim()));
    }
    const auto n = static_cast<std::size_t>(array.shape(0));
    const double* xy = array.data();
    std::vector<Point2> points;
    points.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        points.push_back(Point2{xy[2 * i], xy[2 * i + 1]});
    }
    return points;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    // The GIL is released around exactly one call: triangulate, at the bottom
    // of this file. Everything else exposed here is O(1) arithmetic or an
    // array view that copies nothing, where acquiring and releasing would cost
    // more than the call itself. The triangulation kernel is the case that
    // does not fit -- a parallel refinement pass holding the GIL would
    // serialise every worker -- so it, and anything as long-running that lands
    // later, is wrapped in py::gil_scoped_release.
    m.doc() = "C++20 geometry primitives and the constrained Delaunay\n"
            "triangulator for the rasputin terrain engine.";

    py::class_<Point2>(m, "Point2", R"doc(
A 2D point with double-precision coordinates.

Immutable: the coordinates are read-only because the type is hashable, and a
mutable hashable value silently corrupts any set or dict holding it.

Equality is exact, component-wise floating-point comparison, so NaN coordinates
never compare equal and +0.0 equals -0.0. Hashing is consistent with that.
Coordinate de-duplication that must tolerate floating-point drift should snap
-round upstream rather than rely on this.
)doc")
        .def(py::init([](double x, double y) { return Point2{x, y}; }),
             py::arg("x") = 0.0, py::arg("y") = 0.0,
             "Construct from coordinates, defaulting to the origin.")
        .def_readonly("x", &Point2::x, "The x coordinate (read-only).")
        .def_readonly("y", &Point2::y, "The y coordinate (read-only).")
        .def(py::self == py::self, "Exact component-wise equality.")
        .def(py::self != py::self, "Exact component-wise inequality.")
        .def(py::self + py::self, "Component-wise sum.")
        .def(py::self - py::self, "Component-wise difference.")
        .def(py::self * double(), "Scale by a scalar.")
        .def(double() * py::self, "Scale by a scalar.")
        .def(py::self / double(), "Divide by a scalar. Division by zero yields inf or nan.")
        .def(-py::self, "Component-wise negation.")
        // Mirrors the std::formatter specialization in point.hpp, so the C++
        // and Python renderings of a point cannot drift apart.
        .def("__repr__", [](const Point2& p) { return std::format("{}", p); },
             "Render as Point2(x, y), matching the C++ std::formatter output.")
        .def("__hash__", [](const Point2& p) { return std::hash<Point2>{}(p); },
             "Hash of the coordinate pair, consistent with equality.");

    py::class_<Point3>(m, "Point3", R"doc(
A 3D point with double-precision coordinates.

Immutable and hashable on the same terms as Point2: exact component-wise
equality, read-only coordinates, and a hash consistent with equality.
)doc")
        .def(py::init([](double x, double y, double z) { return Point3{x, y, z}; }),
             py::arg("x") = 0.0, py::arg("y") = 0.0, py::arg("z") = 0.0,
             "Construct from coordinates, defaulting to the origin.")
        .def_readonly("x", &Point3::x, "The x coordinate (read-only).")
        .def_readonly("y", &Point3::y, "The y coordinate (read-only).")
        .def_readonly("z", &Point3::z, "The z coordinate (read-only).")
        .def(py::self == py::self, "Exact component-wise equality.")
        .def(py::self != py::self, "Exact component-wise inequality.")
        .def(py::self + py::self, "Component-wise sum.")
        .def(py::self - py::self, "Component-wise difference.")
        .def(py::self * double(), "Scale by a scalar.")
        .def(double() * py::self, "Scale by a scalar.")
        .def(py::self / double(), "Divide by a scalar. Division by zero yields inf or nan.")
        .def(-py::self, "Component-wise negation.")
        .def("__repr__", [](const Point3& p) { return std::format("{}", p); },
             "Render as Point3(x, y, z), matching the C++ std::formatter output.")
        .def("__hash__", [](const Point3& p) { return std::hash<Point3>{}(p); },
             "Hash of the coordinate triple, consistent with equality.");

    // Overloaded rather than method-bound: dot and cross are symmetric binary
    // operations, and mixing dimensions is a TypeError rather than a coercion.
    m.def("dot", py::overload_cast<const Point2&, const Point2&>(&terrain::dot),
          py::arg("a"), py::arg("b"),
          R"doc(
Dot product of two points of the same dimension.

Accepts two Point2 or two Point3. Mixing dimensions raises TypeError rather
than coercing. Runs in constant time.
)doc");
    m.def("dot", py::overload_cast<const Point3&, const Point3&>(&terrain::dot),
          py::arg("a"), py::arg("b"), "Dot product of two Point3.");

    m.def("cross", py::overload_cast<const Point2&, const Point2&>(&terrain::cross),
          py::arg("a"), py::arg("b"),
          R"doc(
Cross product, whose result type depends on the dimension.

For two Point2 this returns a float: the scalar z-component of the 3D cross
product, i.e. a.x*b.y - a.y*b.x. Its sign gives the orientation of the turn
from a to b. For two Point3 it returns a Point3 orthogonal to both.

Mixing dimensions raises TypeError. Runs in constant time.
)doc");
    m.def("cross", py::overload_cast<const Point3&, const Point3&>(&terrain::cross),
          py::arg("a"), py::arg("b"),
          "Cross product of two Point3, returning a Point3 orthogonal to both.");

    // ----------------------------------------------------------------------
    // The CDT surface. docs/increments/06-cdt-viewer.md, "The binding
    // surface", is exhaustive: an addition here is a design change needing a
    // reason in an increment file, not a line in a PR (risk 3). Nothing below
    // takes or returns a path, a file, a CRS, a raw span, a mutable view, a
    // detria type, an adjacency or a kernel parameter.
    // ----------------------------------------------------------------------

    py::enum_<ChainRole>(m, "ChainRole", R"doc(
The role a constraint chain plays in the domain.

Outer and Hole are closed rings whose winding the validator checked; Breakline
is an open polyline. The integer values are the C++ enumerator order, which a
renderer keying a colour off .value depends on.
)doc")
        .value("Outer", ChainRole::Outer, "A counterclockwise ring bounding the domain.")
        .value("Hole", ChainRole::Hole, "A clockwise ring carving a hole out of the domain.")
        .value("Breakline", ChainRole::Breakline, "An open polyline constraint.");

    // Not in the design's "what crosses out" table, but PslgDiagnostic.error
    // has to be *something* and an enum is what ChainRole and CdtStatus have
    // already established a C++ enumeration crosses as. A string would make
    // every Python consumer compare magic text; an int would lose the names.
    py::enum_<PslgError>(m, "PslgError", R"doc(
Why the validator rejected a proposed PSLG. Carried by PslgDiagnostic.error.

The validator is exhaustive rather than early-returning, so bulk-wrong input
commonly yields several of these at once.
)doc")
        .value("NoOuterChain", PslgError::NoOuterChain)
        .value("ChainTooShort", PslgError::ChainTooShort)
        .value("IndexOutOfRange", PslgError::IndexOutOfRange)
        .value("NonFiniteVertex", PslgError::NonFiniteVertex)
        .value("StoredClosure", PslgError::StoredClosure)
        .value("WrongWinding", PslgError::WrongWinding)
        .value("DegenerateRing", PslgError::DegenerateRing)
        .value("VertexCountOverflow", PslgError::VertexCountOverflow);

    py::enum_<CdtStatus>(m, "CdtStatus", R"doc(
What the triangulation backend did, grouped by what the caller should do about
it rather than by where in the backend it arose. describe() renders each one.
)doc")
        .value("Ok", CdtStatus::Ok)
        .value("NotRun", CdtStatus::NotRun)
        .value("NotNoded", CdtStatus::NotNoded)
        .value("DegenerateGeometry", CdtStatus::DegenerateGeometry)
        .value("InvalidTopology", CdtStatus::InvalidTopology)
        .value("MalformedInput", CdtStatus::MalformedInput)
        .value("BackendFailure", CdtStatus::BackendFailure);

    m.def("describe", [](CdtStatus status) { return std::string{terrain::cdt::describe(status)}; },
          py::arg("status"),
          "One sentence of prose for a CdtStatus, including every failure "
          "status. Raises TypeError for anything that is not a CdtStatus.");

    py::class_<Chain>(m, "Chain", R"doc(
One constraint chain: a half-open run [begin, begin + count) of chain_indices,
its role, and whether it is a river. Frozen; count counts DISTINCT vertices, so
a closed ring does not store its closing index.
)doc")
        .def_readonly("begin", &Chain::begin, "Offset of this chain's first index.")
        .def_readonly("count", &Chain::count, "Number of distinct vertices in this chain.")
        .def_readonly("role", &Chain::role, "The chain's ChainRole.")
        .def_readonly("is_river", &Chain::is_river, "Whether the chain is a river feature.");

    py::class_<PslgDiagnostic>(m, "PslgDiagnostic", R"doc(
One reason a proposed PSLG was rejected: an error, the offending chain, the
offending vertex, and a human-readable message carrying the values that do not
fit the two index fields. Frozen.
)doc")
        .def_readonly("error", &PslgDiagnostic::error, "The PslgError kind.")
        .def_readonly("chain", &PslgDiagnostic::chain,
                      "Index into chains, or 2**32 - 1 when no single chain is at fault.")
        // The header's warning was written for exactly this binding: it is a
        // valid index or the sentinel, never the offending out-of-range value,
        // so indexing vertices with it cannot perform the out-of-bounds read
        // that IndexOutOfRange exists to prevent.
        .def_readonly("vertex", &PslgDiagnostic::vertex,
                      "A VALID index into vertices, or the kNoVertex sentinel 2**32 - 1 -- "
                      "never the offending out-of-range value.")
        .def_readonly("message", &PslgDiagnostic::message, "Human-readable detail; not parsed.");

    py::class_<PslgBuildResult>(m, "PslgBuildResult", R"doc(
The outcome of build_pslg: either a Pslg, or the WHOLE list of reasons why not.

Failure is data rather than an exception, and the list is complete rather than
first-failure, because input at this boundary is wrong in bulk and one
diagnostic per call turns one fix into N round trips.
)doc")
        .def_property_readonly("ok", &PslgBuildResult::ok,
                               "True iff a Pslg was produced and diagnostics is empty.")
        .def_property_readonly(
            "pslg", [](const PslgBuildResult& r) { return r.ok() ? &*r.pslg : nullptr; },
            py::return_value_policy::reference_internal, "The Pslg, or None if ok is False.")
        .def_readonly("diagnostics", &PslgBuildResult::diagnostics,
                      "Every PslgDiagnostic the validator found, not just the first.");

    py::class_<Pslg>(m, "Pslg", R"doc(
A validated planar straight-line graph: the constraint set, checked once.

Opaque and not constructible from Python. There is exactly one producer,
build_pslg, and holding a Pslg is the proof that validation ran -- a Python
constructor would void that proof. Every accessor is a read-only view.
)doc")
        .def_property_readonly(
            "vertices",
            [](const py::object& self) {
                return point_view(self, self.cast<const Pslg&>().vertices());
            },
            "Read-only (N, 2) float64 view of the vertex buffer.")
        .def_property_readonly(
            "chains",
            [](const Pslg& self) {
                return std::vector<Chain>{self.chains().begin(), self.chains().end()};
            },
            "The chains, in order, as frozen Chain records.")
        .def_property_readonly(
            "chain_indices",
            [](const py::object& self) {
                return index_view(self, self.cast<const Pslg&>().chain_indices());
            },
            "Read-only (M,) uint32 view of the flat index buffer.")
        // The C++ precondition is a debug assert, which in a release build is
        // no precondition at all; reached from Python it has to be an
        // exception rather than a crash vector.
        .def(
            "indices_of",
            [](const py::object& self, std::size_t c) {
                const Pslg& pslg = self.cast<const Pslg&>();
                if (c >= pslg.chains().size()) {
                    throw py::index_error(std::format("chain {} is out of range: this PSLG has "
                                                      "{} chains",
                                                      c, pslg.chains().size()));
                }
                return index_view(self, pslg.indices_of(c));
            },
            py::arg("c"),
            "Read-only uint32 view of chain c's indices. Raises IndexError out of range.");

    py::class_<IndexedMesh2>(m, "IndexedMesh2", R"doc(
A flat indexed triangle mesh: vertices, triangles, and one constraint mask per
triangle. Read-only, and not constructible from Python.

Bit e of a mask is set iff the edge (v[e], v[(e + 1) % 3]) is a constraint edge
-- NOT CGAL's "edge e is opposite vertex e", which is a rotation of this and
draws a plausible picture with the constraints on the wrong edges.

The three arrays are zero-copy views that keep this mesh alive, so an array
stays valid after every other reference to the mesh is dropped.
)doc")
        .def_property_readonly(
            "vertices",
            [](const py::object& self) {
                return point_view(self, self.cast<const IndexedMesh2&>().vertices());
            },
            "Read-only (N, 2) float64 view. Begins with the PSLG's vertices, in order.")
        .def_property_readonly(
            "triangles",
            [](const py::object& self) {
                static_assert(sizeof(TriangleIndices) == 3 * sizeof(std::uint32_t));
                const std::span<const TriangleIndices> t =
                    self.cast<const IndexedMesh2&>().triangles();
                return readonly_view<std::uint32_t>(
                    self, reinterpret_cast<const std::uint32_t*>(t.data()),
                    {static_cast<py::ssize_t>(t.size()), 3},
                    {static_cast<py::ssize_t>(sizeof(TriangleIndices)),
                     static_cast<py::ssize_t>(sizeof(std::uint32_t))});
            },
            "Read-only (T, 3) uint32 view of counterclockwise triangles.")
        .def_property_readonly(
            "constrained_edges",
            [](const py::object& self) {
                const std::span<const std::uint8_t> mask =
                    self.cast<const IndexedMesh2&>().constrained_edges();
                return readonly_view<std::uint8_t>(self, mask.data(),
                                                   {static_cast<py::ssize_t>(mask.size())}, {1});
            },
            "Read-only (T,) uint8 view of the per-triangle constraint masks.")
        .def_property_readonly("triangle_count", &IndexedMesh2::triangle_count,
                               "Number of triangles.")
        // Asks about triangles, not vertices: a point set with no outline
        // triangulates successfully to zero interior triangles, which is the
        // "forgot the outline" silent failure.
        .def_property_readonly("empty", &IndexedMesh2::empty, "True iff there are no triangles.");

    py::class_<CdtOutcome>(m, "CdtOutcome", R"doc(
What triangulate returned: a status, a message, and a mesh.

On Ok the message is empty and the mesh is non-empty; on failure the message
says something and the mesh is empty. ok() is a method rather than a property,
mirroring the C++ accessor.
)doc")
        .def_readonly("status", &CdtOutcome::status, "The CdtStatus.")
        .def_readonly("message", &CdtOutcome::message,
                      "The backend's diagnosis, empty on success.")
        .def_property_readonly(
            "mesh", [](const CdtOutcome& self) { return &self.mesh; },
            py::return_value_policy::reference_internal,
            "The mesh, which keeps this outcome alive.")
        .def("ok", &CdtOutcome::ok, "True iff status is Ok.");

    // A free function rather than a bound PslgBuilder: build() is
    // rvalue-ref-qualified, so binding it would leave a moved-from builder that
    // any Python name could still call. Python declares what PSLG it wants and
    // receives one, or the reasons why not.
    m.def(
        "build_pslg",
        [](const py::object& vertices, const py::iterable& chains) {
            PslgBuilder builder{as_points(vertices)};
            for (const py::handle chain : chains) {
                const auto [indices, role, is_river] =
                    chain.cast<std::tuple<std::vector<std::uint32_t>, ChainRole, bool>>();
                builder.add_chain(std::span<const std::uint32_t>{indices}, role, is_river);
            }
            return std::move(builder).build<terrain::pred::DefaultKernel>();
        },
        py::arg("vertices"), py::arg("chains"), R"doc(
Validate a constraint set and return a PslgBuildResult.

vertices is any (N, 2) float64-convertible array-like; chains is a sequence of
(indices, role, is_river). Invalid input is reported as the result's whole
diagnostics list rather than raised. A mis-shaped vertex array is a ValueError;
a path or a filename is a TypeError, because the core never sees a path.
)doc");

    m.def(
        "triangulate",
        [](const Pslg& pslg, bool delaunay) {
            const CdtOptions options{delaunay};
            // The one call in this module long enough to be worth the cost of
            // releasing: a refinement pass running these in parallel would
            // otherwise serialise on the GIL. Released before the backend runs
            // and reacquired before the outcome is converted.
            const py::gil_scoped_release unlocked;
            return terrain::cdt::triangulate<DetriaBackend>(pslg, options);
        },
        py::arg("pslg"), py::arg("delaunay") = true, R"doc(
Triangulate a validated PSLG, releasing the GIL for the duration.

The input must already be noded: crossing constraints are a failure status, not
a repair. With delaunay=False the backend skips the Delaunay flips, so the same
input can be seen both ways.
)doc");
}
