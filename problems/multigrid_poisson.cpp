#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fem/io/input.hpp>
#include <fem/io/output.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/cg.hpp>
#include <fem/mesh/hierarchy.hpp>
#include <fem/mesh/mesh.hpp>
#include <fem/problems/poisson.hpp>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct MultigridLevel {
    fem::mesh::triangular_2d::Mesh mesh;
    fem::linalg::csr_matrix stiffness_matrix;
    std::vector<double> rhs;
};

double constant_rhs(double, double) { return 1.0; }

std::vector<double> compute_residual(const fem::linalg::csr_matrix &matrix,
                                     const std::vector<double> &rhs, const std::vector<double> &x) {
    std::vector<double> Ax(x.size(), 0.0);
    std::vector<double> residual(x.size(), 0.0);
    fem::linalg::matrix_vector_product(matrix, x, Ax);
    fem::linalg::vector::vector_sum(1.0, rhs, -1.0, Ax, residual);
    return residual;
}

void gauss_seidel(const fem::linalg::csr_matrix &matrix, const std::vector<double> &rhs,
                  std::vector<double> &x, size_t iterations) {
    assert(matrix.row_index.size() == x.size() + 1);
    assert(rhs.size() == x.size());

    for (size_t sweep = 0; sweep < iterations; ++sweep) {
        for (size_t row = 0; row < x.size(); ++row) {
            double diagonal = 0.0;
            double off_diagonal_sum = 0.0;

            for (size_t entry = matrix.row_index[row]; entry < matrix.row_index[row + 1]; ++entry) {
                const size_t col = matrix.col_index[entry];
                const double value = matrix.values[entry];

                if (col == row) {
                    diagonal = value;
                } else {
                    off_diagonal_sum += value * x[col];
                }
            }

            if (std::abs(diagonal) < 1e-14) {
                throw std::runtime_error("Gauss-Seidel encountered a zero diagonal entry");
            }

            x[row] = (rhs[row] - off_diagonal_sum) / diagonal;
        }
    }
}

fem::linalg::csr_matrix compress_triplets(std::vector<fem::problems::poisson::triplet> triplets,
                                          size_t n_rows) {
    std::sort(triplets.begin(), triplets.end(), [](const auto &a, const auto &b) {
        if (a.global_row == b.global_row) {
            return a.global_col < b.global_col;
        }
        return a.global_row < b.global_row;
    });

    std::vector<fem::problems::poisson::triplet> compressed;
    compressed.reserve(triplets.size());

    for (size_t i = 0; i < triplets.size(); ++i) {
        const size_t row = triplets[i].global_row;
        const size_t col = triplets[i].global_col;
        double sum = 0.0;

        size_t j = i;
        while (j < triplets.size() && triplets[j].global_row == row &&
               triplets[j].global_col == col) {
            sum += triplets[j].contribution;
            ++j;
        }

        if (std::abs(sum) > 1e-12) {
            compressed.push_back({row, col, sum});
        }
        i = j - 1;
    }

    fem::linalg::csr_matrix matrix;
    matrix.row_index.assign(n_rows + 1, 0);
    matrix.values.reserve(compressed.size());
    matrix.col_index.reserve(compressed.size());

    size_t current_row = 0;
    for (const auto &triplet : compressed) {
        while (current_row < triplet.global_row) {
            matrix.row_index[current_row + 1] = matrix.values.size();
            ++current_row;
        }

        matrix.values.push_back(triplet.contribution);
        matrix.col_index.push_back(triplet.global_col);
    }

    while (current_row < n_rows) {
        matrix.row_index[current_row + 1] = matrix.values.size();
        ++current_row;
    }

    return matrix;
}

MultigridLevel assemble_unit_poisson_system(const fem::mesh::triangular_2d::Mesh &mesh) {
    std::vector<double> rhs(mesh.nodes.size(), 0.0);
    std::vector<fem::problems::poisson::triplet> triplets;
    triplets.reserve(9 * mesh.elements.size());

    for (const auto &element : mesh.elements) {
        const std::array<double, 2> p1 = mesh.nodes[element[0]];
        const std::array<double, 2> p2 = mesh.nodes[element[1]];
        const std::array<double, 2> p3 = mesh.nodes[element[2]];

        const double J[2][2] = {
            {p2[0] - p1[0], p3[0] - p1[0]},
            {p2[1] - p1[1], p3[1] - p1[1]},
        };

        const double det_J = J[0][0] * J[1][1] - J[1][0] * J[0][1];
        assert(std::abs(det_J) > 1e-12);

        const double area = 0.5 * std::abs(det_J);

        const double inv_J_transposed[2][2] = {
            {J[1][1] / det_J, -J[1][0] / det_J},
            {-J[0][1] / det_J, J[0][0] / det_J},
        };

        const double grad_ref[3][2] = {{-1.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}};
        double grad_basis[3][2];

        for (int i = 0; i < 3; ++i) {
            grad_basis[i][0] =
                inv_J_transposed[0][0] * grad_ref[i][0] + inv_J_transposed[0][1] * grad_ref[i][1];
            grad_basis[i][1] =
                inv_J_transposed[1][0] * grad_ref[i][0] + inv_J_transposed[1][1] * grad_ref[i][1];
        }

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                const double dot =
                    grad_basis[i][0] * grad_basis[j][0] + grad_basis[i][1] * grad_basis[j][1];
                triplets.push_back(
                    {static_cast<size_t>(element[i]), static_cast<size_t>(element[j]), area * dot});
            }
        }

        const std::array<double, 2> centroid = {
            (p1[0] + p2[0] + p3[0]) / 3.0,
            (p1[1] + p2[1] + p3[1]) / 3.0,
        };
        const double load = area * constant_rhs(centroid[0], centroid[1]) / 3.0;

        rhs[element[0]] += load;
        rhs[element[1]] += load;
        rhs[element[2]] += load;
    }

    fem::linalg::csr_matrix stiffness_matrix =
        compress_triplets(std::move(triplets), mesh.nodes.size());
    fem::problems::poisson::impose_boundary_conditions(rhs, stiffness_matrix, mesh);

    return {mesh, stiffness_matrix, rhs};
}

void add_in_place(std::vector<double> &x, const std::vector<double> &correction) {
    assert(x.size() == correction.size());
    for (size_t i = 0; i < x.size(); ++i) {
        x[i] += correction[i];
    }
}

/**
 * @brief Applies one recursive multigrid V-cycle.
 *
 * @param levels Discretizations from coarsest level `0` to finest level.
 * @param hierarchy Mesh hierarchy providing inter-level transfer operators.
 * @param level Current zero-based level.
 * @param x Current approximation on @p level, modified in place.
 * @param rhs Right-hand side on @p level.
 * @param nu1 Number of pre-smoothing Gauss-Seidel sweeps.
 * @param nu2 Number of post-smoothing Gauss-Seidel sweeps.
 */
void multigrid_v_cycle(const std::vector<MultigridLevel> &levels,
                       const fem::mesh::hierarchical::HierarchicalMesh &hierarchy, size_t level,
                       std::vector<double> &x, const std::vector<double> &rhs, size_t nu1,
                       size_t nu2) {
    assert(level < levels.size());
    assert(x.size() == rhs.size());

    const auto &level_data = levels[level];

    if (level == 0) {
        fem::linalg::cg::solution coarse_solution =
            fem::linalg::cg::CG(level_data.stiffness_matrix, rhs, 1000, 1e-14);
        x = std::move(coarse_solution.solution);
        return;
    }

    gauss_seidel(level_data.stiffness_matrix, rhs, x, nu1);

    const std::vector<double> residual = compute_residual(level_data.stiffness_matrix, rhs, x);
    const std::vector<double> coarse_rhs =
        fem::mesh::hierarchical::restrict_to_previous_level(hierarchy, level - 1, residual);

    std::vector<double> coarse_error(coarse_rhs.size(), 0.0);
    multigrid_v_cycle(levels, hierarchy, level - 1, coarse_error, coarse_rhs, nu1, nu2);

    const std::vector<double> fine_error =
        fem::mesh::hierarchical::prolongate_to_next_level(hierarchy, level - 1, coarse_error);
    add_in_place(x, fine_error);

    gauss_seidel(level_data.stiffness_matrix, rhs, x, nu2);
}

std::vector<MultigridLevel> build_multigrid_levels(const std::vector<std::string> &mesh_files) {
    std::vector<MultigridLevel> levels;
    levels.reserve(mesh_files.size());

    for (const auto &mesh_file : mesh_files) {
        levels.push_back(assemble_unit_poisson_system(fem::io::read_mesh_2d(mesh_file)));
    }

    return levels;
}

} // namespace

int main(int argc, char **argv) {
    const size_t nu1 = (argc > 1) ? std::stoul(argv[1]) : 3;
    const size_t nu2 = (argc > 2) ? std::stoul(argv[2]) : 3;
    const size_t finest_level = (argc > 3) ? std::stoul(argv[3]) : 4;
    const size_t max_cycles = (argc > 4) ? std::stoul(argv[4]) : 50;
    const double tolerance = 1e-10;

    std::vector<std::string> mesh_files;
    const size_t number_meshes = 5;
    for (size_t i = 1; i <= number_meshes; ++i) {
        mesh_files.push_back("../data/meshes/quader" + std::to_string(i) + ".vtk");
    }

    if (finest_level >= mesh_files.size()) {
        throw std::out_of_range("Requested multigrid level does not exist");
    }

    std::vector<std::string> active_mesh_files(mesh_files.begin(),
                                               mesh_files.begin() + finest_level + 1);

    std::cout << "Building multigrid hierarchy with " << active_mesh_files.size() << " levels\n";
    const auto hierarchy = fem::mesh::hierarchical::build_hierarchical_mesh(active_mesh_files);
    const auto levels = build_multigrid_levels(active_mesh_files);

    std::vector<double> x(levels[finest_level].rhs.size(), 0.0);

    double residual_norm = fem::linalg::vector::l2_norm(
        compute_residual(levels[finest_level].stiffness_matrix, levels[finest_level].rhs, x));
    const double initial_residual_norm = residual_norm;

    std::cout << "Initial residual norm: " << residual_norm << "\n";

    size_t cycle = 0;
    for (; cycle < max_cycles && residual_norm > tolerance * initial_residual_norm; ++cycle) {
        multigrid_v_cycle(levels, hierarchy, finest_level, x, levels[finest_level].rhs, nu1, nu2);
        residual_norm = fem::linalg::vector::l2_norm(
            compute_residual(levels[finest_level].stiffness_matrix, levels[finest_level].rhs, x));
        std::cout << "V-cycle " << cycle + 1 << ": residual norm = " << residual_norm << "\n";
    }

    fem::linalg::cg::solution solution;
    solution.solution = x;
    solution.iterations = static_cast<int>(cycle);
    solution.final_residuum = residual_norm;

    fem::io::write_solution_vtk_2d(solution, levels[finest_level].mesh, "multigrid_poisson.vtk");

    std::cout << "Multigrid finished after " << cycle << " V-cycles with residual norm "
              << residual_norm << "\n";

    return 0;
}
