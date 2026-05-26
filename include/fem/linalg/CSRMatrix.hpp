#pragma once

#include <cstddef>
#include <vector>

/**
 * @file CSRMatrix.hpp
 * @brief Sparse matrix utilities for compressed sparse row storage.
 */

/**
 * @namespace fem::linalg
 * @brief Linear algebra data structures and algorithms.
 */
namespace fem::linalg {

/**
 * @brief Sparse matrix stored in compressed sparse row (CSR) format.
 *
 * For row `i`, nonzero entries are stored in the half-open range
 * `[row_index[i], row_index[i + 1])`. `col_index[k]` gives the column of
 * `values[k]`.
 */
struct csr_matrix {
    /**
     * @brief Nonzero matrix entries in row-major CSR order.
     */
    std::vector<double> values;

    /**
     * @brief Column index for each entry in @ref values.
     */
    std::vector<size_t> col_index;

    /**
     * @brief Row offsets into @ref values and @ref col_index.
     *
     * A matrix with `n` rows has `n + 1` row offsets.
     */
    std::vector<size_t> row_index;
};

/**
 * @brief Computes `result = matrix * x`.
 *
 * @param matrix Sparse matrix in CSR format.
 * @param x Input vector. Its size must equal the number of matrix rows.
 * @param result Output vector. It must already have the correct size.
 */
void matrix_vector_product(const csr_matrix &matrix, const std::vector<double> &x,
                           std::vector<double> &result);

/**
 * @brief Extracts the strict upper triangular part of a CSR matrix.
 *
 * @param matrix Input matrix.
 * @return Matrix containing entries with column index greater than row index.
 */
csr_matrix get_triu(const csr_matrix &matrix);

/**
 * @brief Extracts the strict lower triangular part of a CSR matrix.
 *
 * @param matrix Input matrix.
 * @return Matrix containing entries with column index smaller than row index.
 */
csr_matrix get_tril(const csr_matrix &matrix);

/**
 * @brief Extracts the diagonal entries of a CSR matrix.
 *
 * @param matrix Input matrix.
 * @return Sparse diagonal matrix containing entries with equal row and column index.
 */
csr_matrix get_diag(const csr_matrix &matrix);

/**
 * @brief Prints a dense representation of a square CSR matrix to standard output.
 *
 * @param matrix Matrix to print.
 */
void print_csr_matrix(const csr_matrix &matrix);

/**
 * @brief Checks whether a square CSR matrix is symmetric.
 *
 * @param matrix Matrix to inspect.
 * @return `true` if every stored entry has a matching transposed entry within tolerance.
 */
bool is_symmetric(const csr_matrix &matrix);

/**
 * @brief Builds the standard one-dimensional Poisson tridiagonal matrix.
 *
 * The returned matrix has diagonal entries `2` and off-diagonal entries `-1`.
 *
 * @param n Number of unknowns.
 * @return CSR matrix for the one-dimensional Poisson stencil.
 */
csr_matrix make_poisson_1d(const int n);

} // namespace fem::linalg
