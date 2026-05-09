#pragma once

namespace fem::mesh::triangular_2d{
    struct Mesh { // templates for dimension independence?
        std::vector<std::array<double,2>> nodes; // x, y coordinates of the node
        std::vector<std::array<int, 3>> elements; // indices of the 3 corner nodes of each element 
    };
}
