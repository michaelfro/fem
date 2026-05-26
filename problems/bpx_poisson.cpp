#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <fem/io/input.hpp>
#include <fem/io/output.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/cg.hpp>
#include <fem/mesh/hierarchy.hpp>
#include <fem/mesh/mesh.hpp>
#include <fem/problems/poisson.hpp>
#include <iostream>
#include <vector>

using fem::linalg::csr_matrix;
using fem::mesh::hierarchical::build_hierarchical_mesh;
using fem::mesh::hierarchical::HierarchicalMesh;
using fem::mesh::hierarchical::prolongate_to_next_level;
using fem::mesh::hierarchical::restrict_to_previous_level;
using fem::problems::poisson::assemble_system;
using fem::problems::poisson::lin_system;

/**
 * @brief Creates a BPX preconditioner action on the finest mesh.
 *
 * @param hm Hierarchical mesh data containing coarse-to-fine mappings and inverse diagonals.
 * @return Callable that applies `sum_l P_l D_l^{-1} P_l^T` to a residual vector.
 */
auto create_bpx(HierarchicalMesh &hm) {

    return [&hm](const std::vector<double> &r) -> std::vector<double> {
        assert(r.size() == hm.fine_mesh.nodes.size());

        std::vector<double> result(r.size(), 0.0);
        const size_t finest_level = hm.levels.size() - 1;

        std::vector<std::vector<double>> restricted_residuals(hm.levels.size());
        restricted_residuals[finest_level] = r;
        for (size_t current_level = finest_level; current_level > 0; --current_level) {
            restricted_residuals[current_level - 1] = restrict_to_previous_level(
                hm, current_level - 1, restricted_residuals[current_level]);
        }

        for (size_t level = 0; level < hm.levels.size(); ++level) {
            std::vector<double> level_correction = restricted_residuals[level];
            const auto &level_data = hm.levels[level];
            assert(level_correction.size() == level_data.inv_diag.size());

            for (size_t i = 0; i < level_correction.size(); ++i) {
                level_correction[i] *= level_data.inv_diag[i];
            }

            for (size_t current_level = level; current_level < finest_level; ++current_level) {
                level_correction = prolongate_to_next_level(hm, current_level, level_correction);
            }

            for (size_t i = 0; i < result.size(); ++i) {
                result[i] += level_correction[i];
            }
        }

        return result;
    };
}

int main() {

    std::cout << "Assembling hierachical mesh structure \n";
    std::vector<std::string> meshfiles;
    const int number_meshes = 5;
    for (int i = 1; i <= number_meshes; ++i) {
        meshfiles.push_back("../data/meshes/quader" + std::to_string(i) + ".vtk");
    }
    HierarchicalMesh hm = build_hierarchical_mesh(meshfiles);

    std::cout << "Assembling linear system \n";
    lin_system linear_system = assemble_system(hm.fine_mesh);
    std::vector<double> rhs = linear_system.rhs;
    csr_matrix stiffness_matrix = linear_system.stiffness_matrix;

    std::cout << "Creating bpx_preconditioner \n";
    auto bpx_preconditioner = create_bpx(hm);

    size_t max_iterations = 10000;
    double reduction_factor = 10e-10;
    std::vector<double> x_0(rhs.size(), 0.);

    std::cout << "Solving with PCG \n";
    fem::linalg::cg::solution solution_bpx = fem::linalg::cg::PCG(
        stiffness_matrix, rhs, x_0, bpx_preconditioner, max_iterations, reduction_factor);
    std::cout << "PCG finished in " << solution_bpx.iterations
              << " iterations with a final residuum of " << solution_bpx.final_residuum << "\n";

    std::cout << "Solving without BPX with CG \n";
    fem::linalg::cg::solution solution =
        fem::linalg::cg::CG(stiffness_matrix, rhs, max_iterations, reduction_factor);
    std::cout << "Normal CG finished in " << solution.iterations
              << " iterations with a final residuum of " << solution.final_residuum << "\n";

    std::cout << "Error in nodal l2 norm: "
              << fem::problems::poisson::nodal_l2_error(hm.fine_mesh, solution_bpx.solution)
              << "\n";

    fem::io::write_solution_vtk_2d(solution_bpx, hm.fine_mesh, "poisson_bpx.vtk");
    fem::io::write_error(solution_bpx, hm.fine_mesh, fem::problems::poisson::analytical_solution,
                         "poisson_bpx_error.vtk");

    return 0;
}
