#pragma once

#include <cstddef>
#include <fem/mesh/mesh.hpp>
#include <string>
#include <vector>

/**
 * @file hierarchy.hpp
 * @brief Construction and transfer operators for nested mesh hierarchies.
 */

namespace fem::mesh::hierarchical {

/**
 * @brief Builds a nested mesh hierarchy from progressively refined VTK meshes.
 *
 * The input meshes are expected to be nested refinements where nodes introduced
 * between adjacent levels are edge midpoints of the coarser level.
 *
 * @param vtk_files Mesh files ordered from coarsest to finest.
 * @return Hierarchical mesh data with finest-grid mappings and midpoint parents.
 */
HierarchicalMesh build_hierarchical_mesh(const std::vector<std::string> &vtk_files);

/**
 * @brief Prolongates a vector from level `level` to level `level + 1`.
 *
 * Existing coarse nodes are copied to the finer level. Midpoint nodes receive
 * the arithmetic average of their two coarse parent values.
 *
 * @param hm Hierarchical mesh data.
 * @param level Coarse source level.
 * @param coarse_values Values on `level`.
 * @return Values on `level + 1`.
 */
std::vector<double> prolongate_to_next_level(const HierarchicalMesh &hm, size_t level,
                                             const std::vector<double> &coarse_values);

/**
 * @brief Applies the transpose of the level prolongation operator.
 *
 * This maps values from level `level + 1` back to level `level`. Values at
 * unchanged coarse nodes are injected, and midpoint values contribute one half
 * to each parent.
 *
 * @param hm Hierarchical mesh data.
 * @param level Coarse target level.
 * @param fine_values Values on `level + 1`.
 * @return Values on `level`.
 */
std::vector<double> restrict_to_previous_level(const HierarchicalMesh &hm, size_t level,
                                               const std::vector<double> &fine_values);

/**
 * @brief Prolongates values from an arbitrary level to the finest level.
 *
 * @param hm Hierarchical mesh data.
 * @param level Source level.
 * @param values Values on the source level.
 * @return Values on the finest level.
 */
std::vector<double> prolongate_to_finest(const HierarchicalMesh &hm, size_t level,
                                         const std::vector<double> &values);

} // namespace fem::mesh::hierarchical
