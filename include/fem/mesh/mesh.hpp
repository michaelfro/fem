#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

/**
 * @file mesh.hpp
 * @brief Mesh data structures used by FEM assembly and multilevel code.
 */

/**
 * @namespace fem::mesh::triangular_2d
 * @brief Types for two-dimensional triangular meshes.
 */
namespace fem::mesh::triangular_2d {

/**
 * @brief Linear triangular mesh in two spatial dimensions.
 *
 * Nodes are stored as Cartesian coordinates `(x, y)`. Elements store the three
 * zero-based node indices of each triangle.
 */
struct Mesh {
    /**
     * @brief Coordinates of all mesh vertices.
     */
    std::vector<std::array<double, 2>> nodes;

    /**
     * @brief Connectivity of triangular elements.
     */
    std::vector<std::array<int, 3>> elements;
};

} // namespace fem::mesh::triangular_2d

/**
 * @namespace fem::mesh::hierarchical
 * @brief Data structures for nested mesh hierarchies.
 */
namespace fem::mesh::hierarchical {

/**
 * @brief Mapping and solver data for one mesh level.
 *
 * The level is represented with respect to the finest mesh so that vectors can
 * be prolongated, restricted, or accumulated on the finest-grid indexing.
 */
struct LevelData {
    /**
     * @brief Zero-based level number in the hierarchy.
     */
    int level;

    /**
     * @brief Number of nodes on this level.
     */
    int n_nodes;

    /**
     * @brief Maps this level's local node index to the corresponding finest-grid node index.
     */
    std::vector<int> coarse_to_fine;

    /**
     * @brief Maps finest-grid node indices back to this level's local node indices.
     *
     * Entries are `-1` for finest-grid nodes that are not present on this level.
     */
    std::vector<int> fine_to_coarse;

    /**
     * @brief Inverse diagonal of the stiffness matrix assembled on this level.
     */
    std::vector<double> inv_diag;
};

/**
 * @brief Node introduced by one refinement step.
 *
 * The node is represented by its finest-grid index and the two parent nodes on
 * the previous coarser level.
 */
struct MidpointNode {
    /**
     * @brief Index of the introduced node in the finest mesh.
     */
    int fine_index;

    /**
     * @brief First parent node index on the coarser level.
     */
    int coarse_parent_1;

    /**
     * @brief Second parent node index on the coarser level.
     */
    int coarse_parent_2;
};

/**
 * @brief Nested mesh hierarchy for multilevel methods.
 */
struct HierarchicalMesh {
    /**
     * @brief Finest mesh in the hierarchy.
     */
    fem::mesh::triangular_2d::Mesh fine_mesh;

    /**
     * @brief Level data for all meshes, including the finest level.
     */
    std::vector<LevelData> levels;

    /**
     * @brief Newly introduced nodes per refinement step.
     *
     * `midpoints[k]` contains the nodes introduced when going from level `k` to
     * level `k + 1`.
     */
    std::vector<std::vector<MidpointNode>> midpoints;
};

} // namespace fem::mesh::hierarchical
