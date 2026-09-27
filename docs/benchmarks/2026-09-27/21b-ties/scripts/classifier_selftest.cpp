#include <terrain/mesh/prof_ties.hpp>
#include <cstdio>
using namespace terrain::mesh;
using namespace terrain::mesh::prof_ties;
static MeshVertex V(double c, double r) { return MeshVertex{c, r}; }
int main() {
    // shape on (col, row); the classifier reads them as a set
    std::printf("%s\n", shape({V(0,0), V(1,0), V(1,1), V(0,1)}).c_str());        // axis_square_1x1
    std::printf("%s\n", shape({V(0,0), V(2,0), V(2,1), V(0,1)}).c_str());        // axis_rectangle_1x2
    std::printf("%s\n", shape({V(0,0), V(1,1), V(0,2), V(-1,1)}).c_str());       // rotated_square
    std::printf("%s\n", shape({V(0,0), V(2,2), V(1,3), V(-1,1)}).c_str());       // rotated_rectangle
    std::printf("%s\n", shape({V(0,0), V(3,0), V(2,1), V(1,1)}).c_str());        // trapezoid axis? (not cyclic in general, shape only)
    std::printf("%s\n", shape({V(5,0), V(3,4), V(-4,3), V(0,-5)}).c_str());      // other_cyclic (r = 5)
    std::printf("%s\n", shape({V(0,0), V(1,0), V(2,0), V(0,1)}).c_str());        // collinear_triple
    // det: unit square cocircular -> 0; d at centre-ish inside -> positive
    std::printf("det sq %d\n", (int)lattice_det({V(0,0), V(1,0), V(1,1), V(0,1)}));
    // a,b,c CCW in (col,-row): (0,0),(2,0) at row 0; (1,-1) is row 1 -> (1,1) in world up? use row -> -y
    std::printf("det in %d\n", (int)(lattice_det({V(0,0), V(2,0), V(1,-2), V(1,-1)}) > 0));  // world (0,0),(2,0),(1,2), query (1,1): inside
    std::printf("fma 0.1*3 exact %d, 10*4999 exact %d\n", product_exact(3, 0.1), product_exact(4999, 10));
    std::printf("bits 10=%d 0.1=%d 1=%d\n", significant_bits(10), significant_bits(0.1), significant_bits(1));
}
