#include <terrain/mesh/prof_ties.hpp>
#include <cstdio>
using namespace terrain::mesh;
using namespace terrain::mesh::prof_ties;
static MeshVertex V(double c, double r) { return MeshVertex{c, r}; }
int main() {
    std::printf("%s\n", shape({V(5,0), V(3,4), V(-5,0), V(0,-5)}).c_str());  // other_cyclic
    phase = 2;
    record(V(0,0), V(1,0), V(1,1), V(0,1), 0.1, 0.1, 0, true);         // frame inexact? 0.1*1 exact; row 1 -> exact
    record(V(0,3), V(7,0), V(7,3), V(0,0), 0.1, 0.1, 0, true);         // 7*0.1, 3*0.1 inexact
    record(V(0,0), V(1,0), V(1,1), V(0,1), 10, 20, 0, true);           // dx != dy
    record(V(0,0), V(20000,0), V(1,1), V(0,1), 10, 10, 0, true);       // spread > 2^14, det != 0
    record(V(0,0), V(1,0), V(1,1), V(0,1), 10, 10, 1, false);          // sign differs
    record(V(0.5,0), V(1,0), V(1,1), V(0,1), 10, 10, 0, true);         // 3 nodes
    for (auto& [k, n] : counts) std::printf("%s %llu\n", k.c_str(), n);
}
