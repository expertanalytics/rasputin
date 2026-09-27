| measure | quarter | tile |
|---|---:|---:|
| refine loop: incircle calls | 1,615,895 | 1,639,699 |
| refine loop: exact path | 146,962 (9.095 %) | 154,502 (9.423 %) |
| refine loop: exact path, Cocircular | 146,962 | 154,502 |
| refine loop: exact path, four node corners | 146,962 | 154,502 |
| refine loop: exact path, lattice det = 0 | 146,962 | 154,502 |
| refine loop: exact path, lattice det != 0 | 0 | 0 |
| refine loop: exact path, QW2 conditions hold | 146,962 | 154,502 |
| refine loop: filtered path, four node corners | 1,468,493 | 1,485,197 |
| refine loop: filtered path, QW2 conditions hold | 1,468,493 | 1,485,197 |
| refine loop: all calls QW2 would answer | 1,615,455 (99.973 %) | 1,639,699 (100.000 %) |
| refine loop: calls with fewer than four node corners | 440 | 0 |
| refine loop: lattice sign differs from kernel (four nodes) | 0 | 0 |
| legalise_all: incircle calls | 533 | 48,133 |
| legalise_all: exact path | 61 (11.445 %) | 16,129 (33.509 %) |
| legalise_all: exact path, Cocircular | 61 | 16,129 |
| legalise_all: exact path, four node corners | 61 | 16,129 |
| legalise_all: exact path, lattice det = 0 | 61 | 16,129 |
| legalise_all: exact path, lattice det != 0 | 0 | 0 |
| legalise_all: exact path, QW2 conditions hold | 61 | 16,129 |
| legalise_all: filtered path, four node corners | 61 | 32,004 |
| legalise_all: filtered path, QW2 conditions hold | 61 | 32,004 |
| legalise_all: all calls QW2 would answer | 122 (22.889 %) | 48,133 (100.000 %) |
| legalise_all: calls with fewer than four node corners | 411 | 0 |
| legalise_all: lattice sign differs from kernel (four nodes) | 0 | 0 |
| quality pass: incircle calls | 7,677 | 4,017 |
| quality pass: exact path | 22 (0.287 %) | 0 (0.000 %) |
| quality pass: exact path, Cocircular | 22 | 0 |
| quality pass: exact path, four node corners | 22 | 0 |
| quality pass: exact path, lattice det = 0 | 22 | 0 |
| quality pass: exact path, lattice det != 0 | 0 | 0 |
| quality pass: exact path, QW2 conditions hold | 22 | 0 |
| quality pass: filtered path, four node corners | 4,600 | 4,017 |
| quality pass: filtered path, QW2 conditions hold | 4,600 | 4,017 |
| quality pass: all calls QW2 would answer | 4,622 (60.206 %) | 4,017 (100.000 %) |
| quality pass: calls with fewer than four node corners | 3,055 | 0 |
| quality pass: lattice sign differs from kernel (four nodes) | 0 | 0 |
| max spread from d, four-node calls, refine loop (nodes) | 951 | 88 |
| max spread from d, exact-path calls, refine loop (nodes) | 146 | 48 |
| dx, dy; rows x cols | 10, 10; 5051 x 5051 | 10, 10; 5051 x 5051 |
| frame-exact sufficient condition (bits of dx + bit_width) <= 53 | 3 + 13 = 16: True | 3 + 13 = 16: True |
| rounds, inserted, flips | 41, 213,464, 445,657 | 53, 219,837, 445,675 |

Shapes of the exact-path (tie) quads, refine loop:

| shape | quarter | tile |
|---|---:|---:|
| other cyclic quad (no parallel sides) | 47,547 (32.353 %) | 47,269 (30.594 %) |
| axis-aligned square | 39,977 (27.202 %) | 39,911 (25.832 %) |
| rotated square | 17,100 (11.636 %) | 17,203 (11.134 %) |
| axis-aligned rectangle, not square | 16,261 (11.065 %) | 16,567 (10.723 %) |
| isosceles trapezoid, parallel sides on rows or columns | 12,697 (8.640 %) | 20,146 (13.039 %) |
| isosceles trapezoid, other direction | 11,169 (7.600 %) | 11,080 (7.171 %) |
| rotated rectangle, not square | 2,211 (1.504 %) | 2,326 (1.505 %) |

Most frequent axis-aligned sizes (short x long side, nodes), refine loop exact path:

- quarter: square_1x1 37,887, rectangle_1x2 13,119, square_2x2 1,845, rectangle_1x3 1,182, rectangle_2x3 1,107, rectangle_1x4 240
- tile: square_1x1 37,866, rectangle_1x2 13,103, square_2x2 1,794, rectangle_2x3 1,133, rectangle_1x3 1,131, rectangle_1x4 207
