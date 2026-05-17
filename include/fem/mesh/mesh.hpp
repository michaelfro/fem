#pragma once
#include <array>
#include <optional>
#include <cstddef> 

namespace fem::mesh::triangular_2d{
    struct Mesh { // templates for dimension independence?
        std::vector<std::array<double,2>> nodes; // x, y coordinates of the node
        std::vector<std::array<int, 3>> elements; // indices of the 3 corner nodes of each element 
    };
}


namespace fem::mesh::hierarchical {

// Data for the k_th level mapping back to the finest grid
struct LevelData {
    int level;
    int n_nodes; // needed?
    std::vector<int> coarse_to_fine; // maps coarse node index to fine node index
    std::vector<int> fine_to_coarse; // maps fine node index to coarse node index (-1 if not a coarse node)     
    std::vector<double> inv_diag; // inverted diagonal of the stiffness matrix at this level
};

struct MidpointNode { // newly introduced node
    int fine_index; // index in finest mesh
    int coarse_parent_1; // coarse node indices 
    int coarse_parent_2;
};

struct HierarchicalMesh {
    fem::mesh::triangular_2d::Mesh fine_mesh;
    std::vector<LevelData> levels; // One LevelData per level, including finest
    
    // midpoints[k] = nodes introduced when going from level k to k+1
    std::vector<std::vector<MidpointNode>> midpoints;
};

}
