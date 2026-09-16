# Agent Skill: Modern C++ Engineering & Python Bindings

Enforces C++20/C++23 standards, zero-overhead abstractions, component swappability, and Pybind11 integration.

## 1. Component Swappability & Abstraction (C++20)
* **Compile-Time Polymorphism:** Use C++20 **Concepts** and templates instead of runtime virtual tables (`vtable`) for core geometric entities (e.g., `Point`, `Mesh`, `RefinementPolicy`).
* **Interface Definition:** Define strict concepts for standard geometry actions:
  ```cpp
  template<typename T>
  concept GeometryKernel = requires(T k, Point p1, Point p2) {
      { k.orient2d(p1, p2, p3) } -> std::same_as<Orientation>;
  };
  ```
* **Future-Proofing:** Implement the triangulation and refinement pipelines as templates parameterized by these concepts, allowing drop-in replacement of the core engine later.

## 2. Memory, Data Layout & Execution
* **Ownership:** Enforce RAII. Use smart pointers (`std::unique_ptr`) only for resource boundaries; use stack allocation, string views (`std::string_view`), and spans (`std::span`) for zero-copy data passing.
* **Asynchrony & Streaming:** Use C++20 **Coroutines** (`std::generator`, `co_yield`, `co_await`) for lazy terrain sampling, heavy batch processing, or streaming mesh updates to/from Python without massive allocations.
* **Memory Locality:** Prefer contiguous memory layouts (`std::vector`). Avoid pointer-chasing topologies; represent mesh connectivity via indices into contiguous arrays.

## 3. Python Bindings (Pybind11)
* **Binding Isolation:** Keep `pybind11` code strictly inside a dedicated wrapper layer (e.g., `bindings/`). Core C++ headers must remain 100% agnostic of Pybind11.
* **Type Mapping:** Map internal C++ types safely to Python objects or NumPy arrays using `py::array_t` for high-throughput coordinate sharing without copying.
* **GIL Management:** Release the Global Interpreter Lock (`py::gil_scoped_release`) for intensive C++ triangulation or refinement kernels to allow native multi-threading.

## 4. Code Quality & Formatting
* **Safety:** Compiles with `-Wall -Wextra -Werror -Wpedantic`. No `using namespace std;` in headers.
* **Modernity:** Use `auto` for type deduction where it improves readability, `constexpr` for compile-time math, and structured bindings for tuple/struct unpacking.

