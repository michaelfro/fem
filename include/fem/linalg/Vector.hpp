#pragma once

#include <cassert>
#include <vector>

/**
 * @file Vector.hpp
 * @brief Basic dense vector operations used by iterative solvers.
 */

/**
 * @namespace fem::linalg::vector
 * @brief Utility functions for `std::vector<double>` arithmetic.
 */
namespace fem::linalg::vector {

/**
 * @brief Computes the weighted vector sum `res = w1 * v1 + w2 * v2`.
 *
 * @param w1 Weight applied to @p v1.
 * @param v1 First input vector.
 * @param w2 Weight applied to @p v2.
 * @param v2 Second input vector.
 * @param res Output vector. It must have the same size as both inputs.
 */
void vector_sum(double w1, const std::vector<double> &v1, double w2, const std::vector<double> &v2,
                std::vector<double> &res);

/**
 * @brief Computes the Euclidean inner product of two vectors.
 *
 * @param v First vector.
 * @param w Second vector.
 * @return Sum of componentwise products.
 */
double inner_product(const std::vector<double> &v, const std::vector<double> &w);

/**
 * @brief Computes the Euclidean norm of a vector.
 *
 * @param v Input vector.
 * @return `sqrt(inner_product(v, v))`.
 */
double l2_norm(const std::vector<double> &v);

/**
 * @brief Prints vector entries to standard output.
 *
 * @param v Vector to print.
 */
void print(const std::vector<double> &v);

} // namespace fem::linalg::vector
