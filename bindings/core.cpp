#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <terrain/cdt/detria_backend.hpp>
#include <terrain/cdt/result.hpp>
#include <terrain/cdt/triangulate.hpp>
#include <terrain/core/edge_properties.hpp>
#include <terrain/core/indexed_mesh.hpp>
#include <terrain/core/noded_pslg.hpp>
#include <terrain/core/point.hpp>
#include <terrain/core/pslg.hpp>
#include <terrain/core/pslg_builder.hpp>
#include <terrain/noding/node.hpp>
#include <terrain/noding/noded_pslg_builder.hpp>
#include <terrain/predicates/default_kernel.hpp>
#include <terrain/raster/geometry.hpp>
#include <terrain/raster/sample.hpp>
#include <terrain/raster/view.hpp>
#include <terrain/refinement/refine.hpp>

#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace py = pybind11;

using terrain::Chain;
using terrain::ChainRole;
using terrain::EdgeProperties;
using terrain::IndexedMesh2;
using terrain::NodedPslg;
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
using terrain::noding::NodeOptions;
using terrain::noding::NodeOutcome;
using terrain::noding::NodeStatus;
using terrain::refinement::RefineOutcome;
using terrain::refinement::RefineStatus;

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

// The dense per-edge property masks, as bare uint32. EdgeProperties is
// deliberately not bound -- Python already has an integer with |, & and
// bit_count(), and Chain.properties already crosses as bits() -- so this
// reinterprets the buffer rather than wrapping each element. The static_assert
// is what stops that being a silent lie if the type ever grows a member.
[[nodiscard]] py::array_t<std::uint32_t> properties_view(const py::object& owner,
                                                         std::span<const EdgeProperties> props) {
    static_assert(sizeof(EdgeProperties) == sizeof(std::uint32_t) &&
                      std::is_standard_layout_v<EdgeProperties>,
                  "the (E,) uint32 view reinterprets the EdgeProperties buffer as bare masks");
    return readonly_view<std::uint32_t>(owner,
                                        reinterpret_cast<const std::uint32_t*>(props.data()),
                                        {static_cast<py::ssize_t>(props.size())},
                                        {static_cast<py::ssize_t>(sizeof(EdgeProperties))});
}

// The four accessors of viz/protocols.py's PslgLike, bound once for both graph
// types. NOT ~25 lines saved: PslgLike is a structural contract with three
// implementations, and a member added to one hand-written block and forgotten
// on the other is a drift the protocol suite catches only for whichever type it
// happens to name. One template makes the two bindings the same object.
//
// The return is deliberately NOT [[nodiscard]], where 05c-noder-wiring.md's
// sketch of this signature had it: Pslg chains nothing onto it, so the
// attribute would buy one -Wunused-result cast at the only call that correctly
// ignores the value. The registration is the effect; the handle is a
// convenience for the caller that adds accessors.
template <typename T>
py::class_<T> bind_pslg_like(py::module_& m, const char* name, const char* doc) {
    return py::class_<T>(m, name, doc)
        .def_property_readonly(
            "vertices",
            [](const py::object& self) {
                return point_view(self, self.template cast<const T&>().vertices());
            },
            "Read-only (N, 2) float64 view of the vertex buffer.")
        .def_property_readonly(
            "chains",
            [](const T& self) {
                return std::vector<Chain>{self.chains().begin(), self.chains().end()};
            },
            "The chains, in order, as frozen Chain records. Unlike the array\n"
            "accessors on this type, this COPIES: a fresh list of C chains is\n"
            "built on every read, so bind it once rather than re-reading it\n"
            "inside a per-edge loop.")
        .def_property_readonly(
            "chain_indices",
            [](const py::object& self) {
                return index_view(self, self.template cast<const T&>().chain_indices());
            },
            "Read-only (M,) uint32 view of the flat index buffer.")
        // The C++ precondition is a debug assert, which in a release build is
        // no precondition at all; reached from Python it has to be an
        // exception rather than a crash vector.
        .def(
            "indices_of",
            [](const py::object& self, std::size_t c) {
                const T& graph = self.template cast<const T&>();
                if (c >= graph.chains().size()) {
                    throw py::index_error(std::format("chain {} is out of range: this PSLG has "
                                                      "{} chains",
                                                      c, graph.chains().size()));
                }
                return index_view(self, graph.indices_of(c));
            },
            py::arg("c"),
            "Read-only uint32 view of chain c's indices. Raises IndexError out of range.");
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
        // The whole shape, not just the rank: the predicate rejects (3, 3) on
        // its second dimension, and "a 2-dimensional array" both reads as a
        // contradiction and never names the dimension that is wrong.
        std::string shape{"("};
        for (py::ssize_t i = 0; i < array.ndim(); ++i) {
            shape += std::format("{}{}", array.shape(i), i + 1 < array.ndim() ? ", " : "");
        }
        shape += array.ndim() == 1 ? ",)" : ")";
        throw py::value_error(std::format("vertices must have shape (N, 2), got shape {}", shape));
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

// The one admission check on an untrusted property mask, and a ValueError
// rather than a PslgDiagnostic: PslgError enumerates structural defects of a
// constraint set, every one of which a C++ caller can also commit, and this one
// no C++ caller can -- EdgeProperties::bit is the only route to a set bit and
// i >= kMaxProperties is a precondition violation rather than a datum. It joins
// the mis-shaped vertex array above as a marshalling failure.
//
// Compared AS PYTHON INTEGERS rather than cast first, because a Python int is
// unbounded: casting 1 << 64 to any C++ integer type raises before the check
// could run, and `mask >> 32` on a negative int is -1 in Python -- truthy, and
// not a bit position.
[[nodiscard]] EdgeProperties as_properties(const py::object& mask, std::size_t chain) {
    // bool first, and separately, because Python's bool IS an int: it satisfies
    // the isinstance check below and 0 and 1 are both in range, so the
    // pre-migration spelling (idx, role, is_river) would pass silently. True
    // would even mean "river" -- but only because bit 0 happens to be river in
    // DEFAULT_VOCABULARY, which lives in Python and can be renumbered without
    // any C++ suite noticing. The C++ half refuses the same spelling outright
    // (EdgeProperties has no conversion from bool in either direction), and the
    // boundary must not be the one place it survives.
    if (py::isinstance<py::bool_>(mask)) {
        throw py::value_error(std::format(
            "chain {}: the property mask is a bool ({}); pass a mask of property bits, "
            "not the old is_river flag",
            chain, std::string{py::str(mask)}));
    }
    // The type arm is separate from the range arm below because the two
    // defects are different objects and the message has to name the one that
    // was rejected. isinstance<py::int_> is a strict PyLong_Check, so a
    // numpy.int64 -- what a caller holds after indexing any integer array --
    // is refused here however legal its value; folded into the range arm it
    // was reported as "the property mask 1 is not a set of bits below 32",
    // telling the caller to fix the one thing that was not wrong. The bool arm
    // above already names the type for exactly this reason.
    if (!py::isinstance<py::int_>(mask)) {
        throw py::value_error(
            std::format("chain {}: the property mask has type {}; pass a built-in int of "
                        "property bits",
                        chain, std::string{py::str(py::type::of(mask))}));
    }
    const py::int_ ceiling{std::uint64_t{1} << EdgeProperties::kMaxProperties};
    if (mask < py::int_{0} || mask >= ceiling) {
        throw py::value_error(
            std::format("chain {}: the property mask {} is not a set of bits below {}", chain,
                        std::string{py::str(mask)}, EdgeProperties::kMaxProperties));
    }
    // Bit by bit, because EdgeProperties has no conversion from a word in
    // either direction and must never grow one: that absence is what makes
    // add_chain(idx, role, true) a compile error rather than a changed meaning.
    const auto word = mask.cast<std::uint32_t>();
    EdgeProperties properties;
    for (unsigned i = 0; i < EdgeProperties::kMaxProperties; ++i) {
        if (((word >> i) & 1u) != 0u) {
            properties = properties | EdgeProperties::bit(i);
        }
    }
    return properties;
}

// The bound RasterView. The variant picks the cell type once, at construction;
// `array` is the numpy object the view points into, held so the buffer lives as
// long as the view does. Destroyed by pybind11 with the GIL held.
struct BoundRasterView {
    py::object array;
    std::variant<terrain::raster::RasterView<float>, terrain::raster::RasterView<double>> view;
};

// Exact dtype and C order, or nothing: a conversion would be a silent copy the
// view then outlives, so anything else is a TypeError and the caller converts.
template <typename T>
[[nodiscard]] std::optional<BoundRasterView> try_view(const py::object& array,
                                                      const terrain::raster::RasterGeometry* g,
                                                      std::optional<double> nodata) {
    using Exact = py::array_t<T, py::array::c_style>;
    if (!py::isinstance<Exact>(array))
        return std::nullopt;
    const auto a = array.cast<Exact>();
    if (a.ndim() != 2)
        throw py::value_error(std::format("raster_view: array must be 2-D, got {} dimensions",
                                          a.ndim()));
    // Rows and columns come from the array's own shape, never from an argument.
    const terrain::raster::RasterGeometry geometry{
        g->x_min(), g->y_max(), g->delta_x(), g->delta_y(),
        static_cast<std::size_t>(a.shape(1)), static_cast<std::size_t>(a.shape(0))};
    // One conversion of the sentinel to T. decode_dem has already refused a
    // sentinel the cell type cannot hold, so it is exact.
    const std::optional<T> sentinel =
        nodata ? std::optional<T>{static_cast<T>(*nodata)} : std::nullopt;
    return BoundRasterView{array, terrain::raster::RasterView<T>{geometry, a.data(), sentinel}};
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
than coercing.
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

Mixing dimensions raises TypeError.
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
        .value("NoOuterChain", PslgError::NoOuterChain,
               "No chain has role Outer, so there is no bounded domain to mesh.")
        .value("ChainTooShort", PslgError::ChainTooShort,
               "A chain has fewer vertices than its role needs: 3 for a ring, 2 for a "
               "breakline.")
        .value("IndexOutOfRange", PslgError::IndexOutOfRange,
               "A chain names a vertex the buffer does not hold. The offending value is in "
               "the message, never in .vertex.")
        .value("NonFiniteVertex", PslgError::NonFiniteVertex,
               "A vertex coordinate is NaN or infinite.")
        .value("StoredClosure", PslgError::StoredClosure,
               "A ring repeats its first vertex as its last. Drop the trailing index: a "
               "chain stores distinct vertices and closes implicitly.")
        .value("WrongWinding", PslgError::WrongWinding,
               "A ring turns the wrong way: Outer must be counterclockwise, Hole clockwise.")
        .value("DegenerateRing", PslgError::DegenerateRing,
               "A ring encloses no area because all its vertices are collinear.")
        .value("VertexCountOverflow", PslgError::VertexCountOverflow,
               "The vertices or the chain indices do not fit Chain's 32-bit fields.");

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

    py::enum_<NodeStatus>(m, "NodeStatus", R"doc(
What the noder did, grouped by what the caller should do about it. describe()
renders each one, and three of the rows point at the snap spacing -- two of them
in opposite directions, which is why the lever is part of the prose rather than
inferred from the name.
)doc")
        .value("Ok", NodeStatus::Ok)
        .value("NotRun", NodeStatus::NotRun)
        .value("InvalidSnapSpacing", NodeStatus::InvalidSnapSpacing)
        .value("CoordinateOutOfRange", NodeStatus::CoordinateOutOfRange)
        .value("RingCollapsed", NodeStatus::RingCollapsed)
        .value("RingDegenerateAfterSnap", NodeStatus::RingDegenerateAfterSnap)
        .value("NonSimpleRing", NodeStatus::NonSimpleRing)
        .value("NotConverged", NodeStatus::NotConverged)
        .value("MalformedOutput", NodeStatus::MalformedOutput);

    // ONE NAME OVER TWO ENUMERATIONS, and the safety question is answered by
    // pybind11's two-pass resolution: the first pass runs every overload with
    // conversions disabled and py::enum_ registers a distinct Python type per
    // enumeration, so describe(CdtStatus.Ok) and describe(NodeStatus.Ok) resolve
    // exactly even though both are integer 0 underneath. Two names would force
    // every caller holding a status to branch on its type to choose a function,
    // which is precisely the branch a status band does not otherwise need.
    m.def("describe", [](CdtStatus status) { return std::string{terrain::cdt::describe(status)}; },
          py::arg("status"),
          "One sentence of prose for a CdtStatus, including every failure "
          "status. Raises TypeError for anything that is neither a CdtStatus "
          "nor a NodeStatus.");

    m.def("describe",
          [](NodeStatus status) { return std::string{terrain::noding::describe(status)}; },
          py::arg("status"),
          "One sentence of prose for a NodeStatus, including the self-checks. "
          "Raises TypeError for anything that is neither a NodeStatus nor a "
          "CdtStatus.");

    py::class_<Chain>(m, "Chain", R"doc(
One constraint chain: a half-open run [begin, begin + count) of chain_indices,
its role, and its property set. Frozen; count counts DISTINCT vertices, so a
closed ring does not store its closing index.
)doc")
        .def_readonly("begin", &Chain::begin, "Offset of this chain's first index.")
        .def_readonly("count", &Chain::count, "Number of distinct vertices in this chain.")
        .def_readonly("role", &Chain::role, "The chain's ChainRole.")
        .def_property_readonly(
            "properties", [](const Chain& self) { return self.properties.bits(); },
            R"doc(The chain's feature property set, as a bare 32-bit mask.

An int, because EdgeProperties is deliberately not bound: Python already has an
integer with |, & and bit_count(), and the meaning of each bit lives in one
Pydantic model, tin_engine.features.EdgeVocabulary. A bound class would grow a
second, competing vocabulary object next to it and carry no meaning the int does
not.
)doc");

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
                      "Every PslgDiagnostic the validator found, not just the first.\n"
                      "Rebuilt on every read, like Pslg.chains and unlike the array\n"
                      "accessors; bind it once if you read it more than once.");

    bind_pslg_like<Pslg>(m, "Pslg", R"doc(
A validated planar straight-line graph: the constraint set, checked once.

Opaque and not constructible from Python. There is exactly one producer,
build_pslg, and holding a Pslg is the proof that validation ran -- a Python
constructor would void that proof. Every accessor is a read-only view.
)doc");

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

    bind_pslg_like<NodedPslg>(m, "NodedPslg", R"doc(
The constraint set after snap rounding: every crossing and every hot-pixel
incidence resolved into a node.

Opaque, not constructible from Python, and NOT a subclass of Pslg -- there is no
conversion either way. That is the point: triangulate takes one of these, so
un-noded input is unrepresentable at the entry point rather than diagnosed
inside it. Holding one is the proof that NodedPslgBuilder's verification ran,
which is why the builder itself does not cross.

It has the same SHAPE as Pslg -- vertices, chains, chain_indices, indices_of --
because the renderer's scene join is structural and must be fed the noded graph
rather than the input one. Node ids are NOT input vertex ids: the node set is
ordered by its lattice point, so use node_of_input_vertex to follow an input
vertex across.
)doc")
        .def_property_readonly(
            "grid_spacing", [](const NodedPslg& self) { return self.grid().spacing(); },
            "The snap spacing this graph was noded at, in the input's units.\n"
            "The SnapGrid itself does not cross: a grid vocabulary in Python\n"
            "next to no consumer, where this is the one value a caller needs\n"
            "and it is the value it handed in.")
        .def_property_readonly(
            "edge_properties",
            [](const py::object& self) {
                return properties_view(self, self.cast<const NodedPslg&>().edge_properties());
            },
            "Read-only (E,) uint32 view of the per-edge property masks, DENSE\n"
            "and index-aligned with the flat edge enumeration: chain c's edge k\n"
            "is at sum(edge_count(j) for j < c) + k. Each entry is the UNION\n"
            "over every input chain that contributed geometry to that edge, so\n"
            "a road noded along a river carries both bits. 0 means\n"
            "unclassified, which is a legal value and not an error.")
        .def_property_readonly(
            "node_of_input_vertex",
            [](const py::object& self) {
                return index_view(self, self.cast<const NodedPslg&>().node_of_input_vertex());
            },
            "Read-only (N,) uint32 view mapping each INPUT vertex position onto\n"
            "its node id. Total over the input's vertex array, unreferenced\n"
            "vertices included. The only thing that answers \"where did vertex 4\n"
            "go\" without a coordinate search.");

    py::class_<NodeOutcome>(m, "NodeOutcome", R"doc(
What node() returned: a NodedPslg, or a status and a message saying why not.

Mirrors PslgBuildResult and CdtOutcome. On Ok the message is empty and pslg is
present; on failure pslg is None and the message says something. ok() is a
method rather than a property, mirroring the C++ accessor -- and note that `if
outcome.ok` is truthy for a bound method, so the parentheses are load-bearing.
)doc")
        .def_readonly("status", &NodeOutcome::status, "The NodeStatus.")
        .def_readonly("message", &NodeOutcome::message,
                      "The noder's diagnosis, empty on success.")
        .def_property_readonly(
            "pslg",
            [](const NodeOutcome& self) { return self.ok() ? &*self.pslg : nullptr; },
            py::return_value_policy::reference_internal,
            "The NodedPslg, which keeps this outcome alive, or None if ok is False.")
        .def("ok", &NodeOutcome::ok, "True iff a NodedPslg was produced.");

    m.def(
        "node",
        [](const Pslg& pslg, double spacing, std::uint32_t max_rounds) {
            const NodeOptions options{spacing, max_rounds};
            // Licensed by node.hpp:5-11 and by nothing else: node<K> is a pure
            // function of (pslg, options) with no statics and no caches, so two
            // threads on one Pslg produce bit-identical output. Nothing
            // Python-owned is touched between release and reacquire -- the
            // input is a C++ object the wrapper owns, and the outcome is
            // converted after the lock returns.
            const py::gil_scoped_release unlocked;
            return terrain::noding::node<terrain::pred::DefaultKernel>(pslg, options);
        },
        py::arg("pslg"), py::arg("spacing"), py::arg("max_rounds") = 4, R"doc(
Snap-round a validated PSLG into a NodedPslg, releasing the GIL for the duration.

spacing has NO DEFAULT, mirroring NodeOptions::spacing: the right value is a
policy question about the data and this layer is not the composition root.
max_rounds caps the settling loop and reproduces the C++ default; raising it can
only ever turn a refusal into a success, never into a different mesh.

A refused noding is a fact about the terrain and about the spacing, so it comes
back as a status rather than as an exception -- build_pslg's rule, applied to the
second engine entry point. Only a mis-shaped argument raises: a NodedPslg, a
path or a filename is a TypeError, because noding twice is a caller error and
the core never sees a path.
)doc");

    // A free function rather than a bound PslgBuilder: build() is
    // rvalue-ref-qualified, so binding it would leave a moved-from builder that
    // any Python name could still call. Python declares what PSLG it wants and
    // receives one, or the reasons why not.
    m.def(
        "build_pslg",
        [](const py::object& vertices, const py::iterable& chains) {
            PslgBuilder builder{as_points(vertices)};
            std::size_t index = 0;
            for (const py::handle chain : chains) {
                const auto [indices, role, mask] =
                    chain.cast<std::tuple<std::vector<std::uint32_t>, ChainRole, py::object>>();
                builder.add_chain(std::span<const std::uint32_t>{indices}, role,
                                  as_properties(mask, index));
                ++index;
            }
            return std::move(builder).build<terrain::pred::DefaultKernel>();
        },
        py::arg("vertices"), py::arg("chains"), R"doc(
Validate a constraint set and return a PslgBuildResult.

vertices is any (N, 2) float64-convertible array-like; chains is a sequence of
(indices, role, properties), where properties is a bare mask whose bits are
named by a tin_engine.features.EdgeVocabulary. The coordinates are copied, so
the returned Pslg neither aliases nor keeps alive the array handed in. Invalid
input is reported as the result's whole diagnostics list rather than raised. A
mis-shaped vertex array is a ValueError, and so is a mask that is negative or
carries a bit at or above 32 -- both are marshalling failures rather than
PslgDiagnostics, because no C++ caller can commit either. A path or a filename
is a TypeError, because the core never sees a path.
)doc");

    m.def(
        "triangulate",
        [](const NodedPslg& pslg, bool delaunay) {
            const CdtOptions options{delaunay};
            // The one call in this module long enough to be worth the cost of
            // releasing: a refinement pass running these in parallel would
            // otherwise serialise on the GIL. Released before the backend runs
            // and reacquired before the outcome is converted.
            const py::gil_scoped_release unlocked;
            return terrain::cdt::triangulate<DetriaBackend>(pslg, options);
        },
        py::arg("pslg"), py::arg("delaunay") = true, R"doc(
Triangulate a noded PSLG, releasing the GIL for the duration.

The parameter is a NodedPslg and there is no Pslg overload: as of increment 5c,
un-noded input is unrepresentable here rather than diagnosed inside, so handing
this a Pslg is a TypeError and not a status. Call node() first. With
delaunay=False the backend skips the Delaunay flips, so the same input can be
seen both ways.
)doc");

    py::class_<BoundRasterView>(m, "RasterView", R"doc(
A zero-copy view over a 2-D float32 or float64 numpy array of DEM nodes.
Built by raster_view(), not constructible from Python. Holds the array alive.
)doc");

    m.def(
        "raster_view",
        [](const py::object& array, double x_min, double y_max, double delta_x, double delta_y,
           std::optional<double> nodata) {
            // Placeholder cols/rows of 1: only the affine half is read from it,
            // and the RasterGeometry constructor validates the deltas here.
            const terrain::raster::RasterGeometry affine{x_min, y_max, delta_x, delta_y, 1, 1};
            if (auto v = try_view<float>(array, &affine, nodata))
                return std::move(*v);
            if (auto v = try_view<double>(array, &affine, nodata))
                return std::move(*v);
            throw py::type_error("raster_view: array must be a C-contiguous float32 or float64 "
                                 "numpy array; convert it explicitly rather than have it copied");
        },
        py::arg("array"), py::kw_only(), py::arg("x_min"), py::arg("y_max"), py::arg("delta_x"),
        py::arg("delta_y"), py::arg("nodata") = py::none(), R"doc(
View a row-major DEM array as a raster without copying it.

Row 0 lies at y_max and y decreases with the row index; both deltas are
positive. Rows and columns are the array's shape. The affine scalars are
keyword-only. Any other dtype or layout is a TypeError, never a silent copy.
)doc");

    m.def(
        "sample",
        [](const BoundRasterView& raster, const py::object& points) {
            const auto xy =
                py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(points);
            if (!xy || xy.ndim() != 2 || xy.shape(1) != 2)
                throw py::value_error("sample: points must be a float64 array of shape (N, 2)");
            const auto n = static_cast<std::size_t>(xy.shape(0));
            py::array_t<double> z(static_cast<py::ssize_t>(n));
            py::array_t<bool> valid(static_cast<py::ssize_t>(n));
            static_assert(sizeof(Point2) == 2 * sizeof(double) && std::is_standard_layout_v<Point2>,
                          "the (N, 2) float64 buffer is read as Point2 pairs");
            const std::span<const Point2> pts{reinterpret_cast<const Point2*>(xy.data()), n};
            const std::span<double> zs{z.mutable_data(), n};
            const std::span<bool> flags{valid.mutable_data(), n};
            {
                // Every buffer touched below is held by a local or by `raster`.
                const py::gil_scoped_release unlocked;
                std::visit([&](const auto& v) { terrain::raster::bilinear_batch(v, pts, zs, flags); },
                           raster.view);
            }
            return py::make_tuple(z, valid);
        },
        py::arg("view"), py::arg("points"), R"doc(
Bilinear z at each of the (N, 2) points, as (z, valid): float64 (N,) and bool (N,).

valid is False outside the grid, for a non-finite point, and where any of the
four surrounding nodes is NoData or NaN. z is 0.0 there, never NaN.
)doc");

    py::enum_<RefineStatus>(m, "RefineStatus", R"doc(
Why refine produced a mesh or did not. Everything but Ok is a refusal of the
input: a start vertex that is not a DEM node, a start triangle that is not
counter-clockwise, or a tolerance that is negative or not finite.
)doc")
        .value("Ok", RefineStatus::Ok)
        .value("OffLattice", RefineStatus::OffLattice)
        .value("NotCounterClockwise", RefineStatus::NotCounterClockwise)
        .value("InvalidTolerance", RefineStatus::InvalidTolerance);

    py::class_<RefineOutcome>(m, "RefineOutcome", R"doc(
What refine returned: a status and message, the refined mesh as read-only
arrays that keep this outcome alive, and four numbers. The arrays are empty
unless ok().
)doc")
        .def_readonly("status", &RefineOutcome::status, "The RefineStatus.")
        .def_readonly("message", &RefineOutcome::message, "Empty on success.")
        .def("ok", &RefineOutcome::ok, "True iff status is Ok.")
        .def_property_readonly(
            "vertices",
            [](const py::object& self) {
                return point_view(self, self.cast<const RefineOutcome&>().vertices);
            },
            "(M, 2) float64 world coordinates; every one is a DEM node.")
        .def_property_readonly(
            "z",
            [](const py::object& self) {
                const auto& z = self.cast<const RefineOutcome&>().z;
                return readonly_view<double>(self, z.data(), {static_cast<py::ssize_t>(z.size())},
                                             {static_cast<py::ssize_t>(sizeof(double))});
            },
            "(M,) float64, the node's DEM value; 0.0 where valid is False.")
        .def_property_readonly(
            "valid",
            [](const py::object& self) {
                static_assert(sizeof(bool) == sizeof(std::uint8_t));
                const auto& v = self.cast<const RefineOutcome&>().valid;
                return readonly_view<bool>(self, reinterpret_cast<const bool*>(v.data()),
                                           {static_cast<py::ssize_t>(v.size())}, {1});
            },
            "(M,) bool, False where the node is NoData.")
        .def_property_readonly(
            "triangles",
            [](const py::object& self) {
                const auto& t = self.cast<const RefineOutcome&>().triangles;
                return readonly_view<std::uint32_t>(
                    self, reinterpret_cast<const std::uint32_t*>(t.data()),
                    {static_cast<py::ssize_t>(t.size()), 3},
                    {static_cast<py::ssize_t>(sizeof(TriangleIndices)),
                     static_cast<py::ssize_t>(sizeof(std::uint32_t))});
            },
            "(K, 3) uint32 counter-clockwise triangles.")
        .def_property_readonly(
            "edges",
            [](const py::object& self) {
                const auto& e = self.cast<const RefineOutcome&>().edges;
                return readonly_view<std::uint32_t>(
                    self, reinterpret_cast<const std::uint32_t*>(e.data()),
                    {static_cast<py::ssize_t>(e.size()), 2},
                    {static_cast<py::ssize_t>(2 * sizeof(std::uint32_t)),
                     static_cast<py::ssize_t>(sizeof(std::uint32_t))});
            },
            "(F, 2) uint32 constraint edges, each once.")
        .def_property_readonly(
            "masks",
            [](const py::object& self) {
                return index_view(self, self.cast<const RefineOutcome&>().masks);
            },
            "(F,) uint32 property masks, one per edge.")
        .def_readonly("rounds", &RefineOutcome::rounds, "Scan rounds run.")
        .def_readonly("inserted", &RefineOutcome::inserted, "Vertices inserted.")
        .def_readonly("max_error", &RefineOutcome::max_error,
                      "Largest |z - plane| over triangles with three valid vertices.")
        .def_readonly("uncovered", &RefineOutcome::uncovered,
                      "Valid DEM nodes left inside triangles with a NoData vertex.");

    m.def(
        "refine",
        [](const BoundRasterView& raster, const IndexedMesh2& mesh, const py::object& edges,
           const py::object& masks, double tolerance, unsigned threads) {
            using U32 = py::array_t<std::uint32_t, py::array::c_style | py::array::forcecast>;
            const auto e = U32::ensure(edges);
            const auto k = U32::ensure(masks);
            if (!e || !k || e.ndim() != 2 || e.shape(1) != 2 || k.ndim() != 1
                || k.shape(0) != e.shape(0))
                throw py::value_error("refine: edges must be (E, 2) and masks (E,), uint32");
            const auto n = static_cast<std::size_t>(k.shape(0));
            const std::span<const std::array<std::uint32_t, 2>> pairs{
                reinterpret_cast<const std::array<std::uint32_t, 2>*>(e.data()), n};
            const std::span<const std::uint32_t> bits{k.data(), n};
            const terrain::refinement::RefineOptions options{tolerance, threads};
            // Every buffer read below is held by a local or by `raster`, and
            // the outcome is converted after the lock returns.
            const py::gil_scoped_release unlocked;
            return std::visit(
                [&](const auto& v) {
                    return terrain::refinement::refine(v, mesh, pairs, bits, options);
                },
                raster.view);
        },
        py::arg("view"), py::arg("mesh"), py::arg("edges"), py::arg("masks"), py::kw_only(),
        py::arg("tolerance"), py::arg("threads") = 0, R"doc(
Refine a start mesh against the DEM until every triangle is within tolerance.

mesh's vertices must all be DEM nodes and its triangles counter-clockwise;
edges (E, 2) and masks (E,) are its constraint edges as the CLI builds them.
tolerance is in the DEM's vertical unit. threads only sets how the scan is
split; the output is identical for every value, and 0 means all cores.
A refused input comes back as a status; a mis-shaped array is a ValueError.
Releases the GIL.
)doc");
}
