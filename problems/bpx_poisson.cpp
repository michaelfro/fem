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
using fem::problems::poisson::assemble_system;
using fem::problems::poisson::lin_system;

/**
 * @brief Creates a diagonal BPX-style preconditioner action on the finest mesh.
 *
 * @param hm Hierarchical mesh data containing coarse-to-fine mappings and inverse diagonals.
 * @return Callable that applies the current BPX approximation to a residual vector.
 *
 * @note This is currently a componentwise approximation. Full BPX should use the reusable
 * prolongation and restriction functions from `fem::mesh::hierarchical`.
 */
auto create_bpx(HierarchicalMesh &hm) {

    std::vector<double> bpx(hm.fine_mesh.nodes.size(), 0.0);

    for (size_t level = 0; level < hm.levels.size(); ++level) {
        const auto &ld = hm.levels[level];

        for (int i = 0; i < ld.n_nodes; ++i) {
            bpx[ld.coarse_to_fine[i]] += ld.inv_diag[i];
        }

        if (level + 1 < hm.levels.size()) {
            for (const auto &mid : hm.midpoints[level]) {
                bpx[mid.fine_index] += 0.5 * (bpx[ld.coarse_to_fine[mid.coarse_parent_1]] +
                                              bpx[ld.coarse_to_fine[mid.coarse_parent_2]]);
            }
        }
    }

    return [bpx](const std::vector<double> &r) -> std::vector<double> {
        assert(bpx.size() == r.size());
        std::vector<double> result(r.size());
        for (size_t i = 0; i < bpx.size(); ++i) {
            result[i] = bpx[i] * r[i];
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
