# Notices

Rasputin is released under the MIT licence (see `LICENSE`),
`Copyright (c) Expert Analytics AS`.

This file credits the third-party code and data that the repository or its
built package contains. Why each entry is here is recorded in
`docs/increments/release-hygiene.md`, sections 4 and 5.

## Code compiled into the package

The extension module `tin_engine._core` contains compiled code from these two
projects. Their licences require their notice in every copy, binary copies
included, so the texts are reproduced in full.

### detria

Constrained Delaunay triangulation, vendored in `lib/detria/` from
<https://github.com/Kimbatt/detria> at commit
`8aa25f3e0dedf8d37623e7085b69f985e5845ad7`. Upstream offers it under WTFPL or
MIT; rasputin uses it under MIT (`lib/detria/README.md`).

```text
Copyright (c) Kimbatt (https://github.com/Kimbatt)

Permission is hereby granted, free of charge, to any person obtaining a copy of this software
and associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
```

### pybind11

The Python binding layer; its headers compile into `tin_engine._core`.
<https://github.com/pybind/pybind11>, version 3.1.0, BSD-3-Clause.

```text
Copyright (c) 2016 Wenzel Jakob <wenzel.jakob@epfl.ch>, All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## Fetched for tests, not redistributed

- **Catch2** v3.6.0 (<https://github.com/catchorg/Catch2>), BSL-1.0. CMake
  fetches it to build the C++ test suites (`tests/cpp/CMakeLists.txt`); it is
  not part of the package.

## Data in the repository

- **Kartverket DTM10**, © Kartverket, CC BY 4.0: the benchmark tile
  `tests/fixtures/dem_archive/7908_3_10m_z33.tif` and the windows in
  `tests/fixtures/dtm10/` (`extract.py` there records their source release).
- **CORINE Land Cover 2018**, © European Union, Copernicus Land Monitoring
  Service 2018, European Environment Agency (EEA): the files in
  `tests/fixtures/corine/`. The full attribution, the conditions of use and how
  each file was modified are in `tests/fixtures/corine/NOTICE`.
- **Benchmark evidence** in `docs/benchmarks/` (pictures, domains, logs) is
  derived from the same DTM10 and CORINE data and carries the same two credits.

## Runtime dependencies

The package's runtime dependencies (numpy, pydantic, shapely, pyproj, typer,
tifffile and what they pull in) are not redistributed with rasputin: pip
installs each under its own licence. None requires rasputin to use a copyleft
licence. The check, dependency by dependency, is section 4 of
`docs/increments/release-hygiene.md`.
