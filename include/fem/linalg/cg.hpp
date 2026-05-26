#pragma once

#include <cmath>
#include <cstddef>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/Vector.hpp>
#include <iostream>
#include <vector>

/**
 * @file cg.hpp
 * @brief Conjugate gradient solvers for sparse linear systems.
 */

/**
 * @namespace fem::linalg::cg
 * @brief Iterative Krylov solvers based on conjugate gradients.
 */
namespace fem::linalg::cg {

/**
 * @brief Result returned by CG-like solvers.
 */
struct solution {
    /**
     * @brief Approximate solution vector.
     */
    std::vector<double> solution;

    /**
     * @brief Number of iterations performed.
     */
    int iterations;

    /**
     * @brief Final residual measure used by the solver.
     */
    double final_residuum;
};

/**
 * @brief Solves `matrix * x = rhs` with the conjugate gradient method.
 *
 * The initial guess is the zero vector. The method assumes that @p matrix is
 * symmetric positive definite.
 *
 * @param matrix System matrix in CSR format.
 * @param rhs Right-hand side vector.
 * @param max_iterations Maximum number of iterations.
 * @param reduction_factor Absolute stopping tolerance for the residual norm.
 * @return Solver result containing the approximate solution and convergence data.
 */
solution CG(const fem::linalg::csr_matrix &matrix, const std::vector<double> &rhs,
            const size_t max_iterations, double reduction_factor);

/**
 * @brief Solves `matrix * x = rhs` with preconditioned conjugate gradients.
 *
 * The preconditioner must be callable as `apply_preconditioner(r)` and return a
 * vector of the same size as `r`. The method assumes that @p matrix and the
 * effective preconditioner are symmetric positive definite.
 *
 * @tparam Preconditioner Callable type implementing the preconditioner action.
 * @param matrix System matrix in CSR format.
 * @param rhs Right-hand side vector.
 * @param x_0 Initial guess.
 * @param apply_preconditioner Callable applying the preconditioner to a vector.
 * @param max_iterations Maximum number of iterations.
 * @param reduction_factor Absolute stopping tolerance for the preconditioned residual measure.
 * @return Solver result containing the approximate solution and convergence data.
 */
template <typename Preconditioner>
solution PCG(const fem::linalg::csr_matrix &matrix, const std::vector<double> &rhs,
             const std::vector<double> &x_0, Preconditioner &apply_preconditioner,
             const size_t max_iterations, const double reduction_factor) {

    std::vector<double> p;
    std::vector<double> x;
    std::vector<double> r(x_0.size(), 0.);
    std::vector<double> Ax(x_0.size(), 0.);
    std::vector<double> z(x_0.size(), 0.);
    x = x_0;

    fem::linalg::matrix_vector_product(matrix, x_0, Ax);
    fem::linalg::vector::vector_sum(1., rhs, -1., Ax, r);
    z = apply_preconditioner(r);
    p = z;

    double error_norm = fem::linalg::vector::l2_norm(r);

    if (error_norm < reduction_factor) {
        solution my_sol;
        my_sol.solution = x_0;
        my_sol.iterations = 0;
        my_sol.final_residuum = error_norm;
        return my_sol;
    }

    double alpha;
    double beta;
    double rTz = fem::linalg::vector::inner_product(r, z);
    double rTz_new;

    for (size_t k = 0; k < max_iterations; k++) {

        std::vector<double> Ap(x_0.size(), 0.);
        fem::linalg::matrix_vector_product(matrix, p, Ap);

        alpha = rTz / fem::linalg::vector::inner_product(p, Ap);

        fem::linalg::vector::vector_sum(1, x, alpha, p, x);
        fem::linalg::vector::vector_sum(1., r, -1. * alpha, Ap, r);
        z = apply_preconditioner(r);
        rTz_new = fem::linalg::vector::inner_product(r, z);
        error_norm = std::sqrt(rTz_new);
        if (error_norm < reduction_factor) {
            solution my_sol;
            my_sol.solution = x;
            my_sol.iterations = k + 1;
            my_sol.final_residuum = error_norm;
            return my_sol;
        }

        beta = rTz_new / rTz;
        rTz = rTz_new;

        fem::linalg::vector::vector_sum(1, z, beta, p, p);
    }

    solution my_sol;
    my_sol.solution = x;
    my_sol.iterations = max_iterations;
    my_sol.final_residuum = error_norm;
    return my_sol;
}

/**
 * @brief Runs a small demonstration of the unpreconditioned CG solver.
 */
void test_cg();

} // namespace fem::linalg::cg
