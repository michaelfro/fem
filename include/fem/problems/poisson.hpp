#pragma once

#include <fem/io/input.hpp>
#include <fem/io/output.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/cg.hpp>
#include <fem/mesh/mesh.hpp>
#include <iostream>
#include <vector>

/**
 * @file poisson.hpp
 * @brief Assembly and demonstration helpers for a two-dimensional Poisson problem.
 */

/**
 * @namespace fem::problems::poisson
 * @brief Linear finite element discretization of a model Poisson equation.
 */
namespace fem::problems::poisson {

/**
 * @brief Function pointer type for scalar right-hand sides.
 */
using rhs_function = double (*)(double, double);

/**
 * @brief Sparse assembly entry before compression into CSR format.
 */
struct triplet {
    /**
     * @brief Global matrix row index.
     */
    size_t global_row;

    /**
     * @brief Global matrix column index.
     */
    size_t global_col;

    /**
     * @brief Value contributed to the matrix entry.
     */
    double contribution;
};

/**
 * @brief Element assembly output before global compression.
 */
struct local_result {
    /**
     * @brief Stiffness contributions as row-column-value triplets.
     */
    std::vector<triplet> local_contributions;

    /**
     * @brief Assembled right-hand side vector.
     */
    std::vector<double> rhs;
};

/**
 * @brief Fully assembled linear system.
 */
struct lin_system {
    /**
     * @brief Global stiffness matrix in CSR format.
     */
    fem::linalg::csr_matrix stiffness_matrix;

    /**
     * @brief Global load vector.
     */
    std::vector<double> rhs;
};

/**
 * @brief Right-hand side function of the model Poisson problem.
 *
 * @param x First coordinate.
 * @param y Second coordinate.
 * @return Function value `f(x, y)`.
 */
double f(const double x, const double y);

/**
 * @brief Imposes homogeneous Dirichlet boundary conditions on the square boundary.
 *
 * Boundary nodes are detected on `x = +/-1` and `y = +/-1`. The corresponding
 * rows and columns are eliminated symmetrically, and boundary rows are replaced
 * by identity rows. This keeps the global node numbering intact while
 * preserving matrix symmetry.
 *
 * @param rhs Right-hand side vector to modify in place.
 * @param stiffness_matrix Stiffness matrix to modify in place.
 * @param mesh Mesh whose boundary nodes are constrained.
 */
void impose_boundary_conditions(std::vector<double> &rhs, fem::linalg::csr_matrix &stiffness_matrix,
                                const fem::mesh::triangular_2d::Mesh mesh);

/**
 * @brief Computes local finite element contributions for all mesh elements.
 *
 * Uses linear triangular basis functions and centroid quadrature for the load
 * vector.
 *
 * @param mesh Triangular mesh to assemble on.
 * @return Triplet stiffness contributions and the assembled load vector.
 */
local_result local_contribution(const fem::mesh::triangular_2d::Mesh mesh);

/**
 * @brief Computes local finite element contributions for a custom right-hand side.
 *
 * Uses linear triangular basis functions and centroid quadrature for the load
 * vector.
 *
 * @param mesh Triangular mesh to assemble on.
 * @param rhs Function used for the load vector.
 * @return Triplet stiffness contributions and the assembled load vector.
 */
local_result local_contribution(const fem::mesh::triangular_2d::Mesh mesh, rhs_function rhs);

/**
 * @brief Assembles the global Poisson linear system.
 *
 * The function computes local contributions, compresses duplicate matrix
 * entries into CSR format, and imposes homogeneous Dirichlet boundary
 * conditions.
 *
 * @param mesh Triangular mesh to assemble on.
 * @return Assembled stiffness matrix and right-hand side.
 */
lin_system assemble_system(const fem::mesh::triangular_2d::Mesh mesh);

/**
 * @brief Assembles the global Poisson linear system for a custom right-hand side.
 *
 * @param mesh Triangular mesh to assemble on.
 * @param rhs Function used for the load vector.
 * @return Assembled stiffness matrix and right-hand side.
 */
lin_system assemble_system(const fem::mesh::triangular_2d::Mesh mesh, rhs_function rhs);

/**
 * @brief Analytical reference solution for the model problem.
 *
 * @param x First coordinate.
 * @param y Second coordinate.
 * @return Exact solution value at `(x, y)`.
 */
double analytical_solution(const double x, const double y);

/**
 * @brief Computes a nodal L2-like error against @ref analytical_solution.
 *
 * @param mesh Mesh whose nodes define the sample points.
 * @param approx_sol Approximate nodal solution.
 * @return Scaled root-mean-square nodal error.
 */
double nodal_l2_error(const fem::mesh::triangular_2d::Mesh mesh,
                      const std::vector<double> &approx_sol);

/**
 * @brief Runs the standalone Poisson example.
 *
 * @return Process-style status code.
 */
int run_poisson_problem();

} // namespace fem::problems::poisson
