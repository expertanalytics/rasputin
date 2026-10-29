---
name: modern-cxx
description: C++20 standards for the terrain core: concepts and templates, ownership, zero-overhead abstractions and Pybind11 bindings.
---

# Modern C++ & Python Bindings

## 1. Swappable components
Core geometric entities (`Point`, `Mesh`, `RefinementPolicy`) use C++20
concepts and templates, not virtual dispatch; the triangulation and
refinement pipelines are templates over those concepts.

## 2. Memory
RAII; `std::unique_ptr` only at resource boundaries; spans and string views
for zero-copy passing. Contiguous storage; mesh connectivity as indices into
arrays, not pointers.

## 3. Bindings
* `pybind11` code lives only in the wrapper layer; core headers know nothing
  of Python.
* Share coordinates through `py::array_t` without copying.
* Release the GIL (`py::gil_scoped_release`) around long triangulation or
  refinement calls.

## 4. Quality
Passes the compiler gate (`CLAUDE.md` §4). No `using namespace std;` in
headers.
