#include <fem/io/input.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/mesh/hierarchy.hpp>
#include <fem/problems/poisson.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace fem::mesh::hierarchical {
namespace {

using CoordinateKey = std::pair<long, long>;

CoordinateKey make_coordinate_key(double x, double y) {
    // The current VTK meshes store coordinates with about six decimal digits.
    // Use a slightly coarser key so computed midpoints match stored nodes.
    constexpr double coordinate_scale = 1e5;
    return {std::llround(x * coordinate_scale), std::llround(y * coordinate_scale)};
}

void check_level(const HierarchicalMesh &hm, size_t level) {
    if (level + 1 >= hm.levels.size()) {
        throw std::out_of_range("Hierarchy level has no finer successor");
    }
}

} // namespace

HierarchicalMesh build_hierarchical_mesh(const std::vector<std::string> &vtk_files) {

    if (vtk_files.empty()) {
        throw std::invalid_argument("Cannot build a hierarchy without mesh files");
    }

    HierarchicalMesh hm;

    std::vector<fem::mesh::triangular_2d::Mesh> meshes;
    meshes.reserve(vtk_files.size());
    for (const auto &file : vtk_files) {
        meshes.push_back(fem::io::read_mesh_2d(file));
    }

    hm.fine_mesh = meshes.back();
    const int n_fine = static_cast<int>(hm.fine_mesh.nodes.size());

    std::map<CoordinateKey, int> fine_lookup;
    for (int i = 0; i < n_fine; ++i) {
        auto [x, y] = hm.fine_mesh.nodes[i];
        fine_lookup[make_coordinate_key(x, y)] = i;
    }

    hm.midpoints.resize(meshes.size() - 1);

    for (size_t k = 0; k < meshes.size(); ++k) {
        const auto &mesh_k = meshes[k];
        const int n_coarse = static_cast<int>(mesh_k.nodes.size());

        LevelData level_data;
        level_data.level = static_cast<int>(k);
        level_data.n_nodes = n_coarse;
        level_data.coarse_to_fine.resize(n_coarse);
        level_data.fine_to_coarse.assign(n_fine, -1);

        const auto coarse_matrix = fem::problems::poisson::assemble_system(mesh_k).stiffness_matrix;
        auto diag = fem::linalg::get_diag(coarse_matrix).values;
        for (auto &entry : diag) {
            entry = 1.0 / entry;
        }
        level_data.inv_diag = std::move(diag);

        for (int i = 0; i < n_coarse; ++i) {
            auto [x, y] = mesh_k.nodes[i];
            const int fine_index = fine_lookup.at(make_coordinate_key(x, y));
            level_data.coarse_to_fine[i] = fine_index;
            level_data.fine_to_coarse[fine_index] = i;
        }

        if (k + 1 < meshes.size()) {
            const auto &mesh_next = meshes[k + 1];
            std::map<CoordinateKey, std::pair<int, int>> midpoint_to_parents;

            for (const auto &tri : mesh_k.elements) {
                const std::array<std::pair<int, int>, 3> edges = {
                    std::make_pair(tri[0], tri[1]),
                    std::make_pair(tri[1], tri[2]),
                    std::make_pair(tri[0], tri[2]),
                };

                for (auto [a, b] : edges) {
                    const double mx = 0.5 * (mesh_k.nodes[a][0] + mesh_k.nodes[b][0]);
                    const double my = 0.5 * (mesh_k.nodes[a][1] + mesh_k.nodes[b][1]);
                    midpoint_to_parents[make_coordinate_key(mx, my)] = {a, b};
                }
            }

            for (size_t i = 0; i < mesh_next.nodes.size(); ++i) {
                auto [x, y] = mesh_next.nodes[i];
                const int fine_index = fine_lookup.at(make_coordinate_key(x, y));

                if (level_data.fine_to_coarse[fine_index] != -1) {
                    continue;
                }

                auto parents = midpoint_to_parents.find(make_coordinate_key(x, y));
                if (parents == midpoint_to_parents.end()) {
                    throw std::runtime_error("Could not find parents for new node at level " +
                                             std::to_string(k + 1));
                }

                auto [a, b] = parents->second;
                hm.midpoints[k].push_back({fine_index, a, b});
            }

            if (hm.midpoints[k].size() + static_cast<size_t>(n_coarse) != mesh_next.nodes.size()) {
                throw std::runtime_error("Unexpected number of midpoint nodes at level " +
                                         std::to_string(k + 1));
            }
        }

        hm.levels.push_back(std::move(level_data));
    }

    return hm;
}

std::vector<double> prolongate_to_next_level(const HierarchicalMesh &hm, size_t level,
                                             const std::vector<double> &coarse_values) {
    check_level(hm, level);

    const auto &coarse_level = hm.levels[level];
    const auto &fine_level = hm.levels[level + 1];
    assert(coarse_values.size() == static_cast<size_t>(coarse_level.n_nodes));

    std::vector<double> fine_values(fine_level.n_nodes, 0.0);

    for (int coarse_index = 0; coarse_index < coarse_level.n_nodes; ++coarse_index) {
        const int finest_index = coarse_level.coarse_to_fine[coarse_index];
        const int fine_index = fine_level.fine_to_coarse[finest_index];
        assert(fine_index >= 0);
        fine_values[fine_index] = coarse_values[coarse_index];
    }

    for (const auto &midpoint : hm.midpoints[level]) {
        const int fine_index = fine_level.fine_to_coarse[midpoint.fine_index];
        assert(fine_index >= 0);
        fine_values[fine_index] = 0.5 * (coarse_values[midpoint.coarse_parent_1] +
                                         coarse_values[midpoint.coarse_parent_2]);
    }

    return fine_values;
}

std::vector<double> restrict_to_previous_level(const HierarchicalMesh &hm, size_t level,
                                               const std::vector<double> &fine_values) {
    check_level(hm, level);

    const auto &coarse_level = hm.levels[level];
    const auto &fine_level = hm.levels[level + 1];
    assert(fine_values.size() == static_cast<size_t>(fine_level.n_nodes));

    std::vector<double> coarse_values(coarse_level.n_nodes, 0.0);

    for (int coarse_index = 0; coarse_index < coarse_level.n_nodes; ++coarse_index) {
        const int finest_index = coarse_level.coarse_to_fine[coarse_index];
        const int fine_index = fine_level.fine_to_coarse[finest_index];
        assert(fine_index >= 0);
        coarse_values[coarse_index] += fine_values[fine_index];
    }

    for (const auto &midpoint : hm.midpoints[level]) {
        const int fine_index = fine_level.fine_to_coarse[midpoint.fine_index];
        assert(fine_index >= 0);
        coarse_values[midpoint.coarse_parent_1] += 0.5 * fine_values[fine_index];
        coarse_values[midpoint.coarse_parent_2] += 0.5 * fine_values[fine_index];
    }

    return coarse_values;
}

std::vector<double> prolongate_to_finest(const HierarchicalMesh &hm, size_t level,
                                         const std::vector<double> &values) {
    if (level >= hm.levels.size()) {
        throw std::out_of_range("Hierarchy level does not exist");
    }

    std::vector<double> result = values;
    for (size_t current_level = level; current_level + 1 < hm.levels.size(); ++current_level) {
        result = prolongate_to_next_level(hm, current_level, result);
    }
    return result;
}

} // namespace fem::mesh::hierarchical
