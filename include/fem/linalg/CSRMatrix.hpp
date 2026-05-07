#pragma once

#include <vector>
#include <cstddef>

namespace fem::linalg::matrix{

    struct csr_matrix {
        std::vector<double> values;
        std::vector<size_t> col_index;
        std::vector<size_t> row_index;
        // int size; // do I need this?
    };

    void matrix_vector_product(const csr_matrix &matrix, const std::vector<double> &x, std::vector<double> &result); 
    csr_matrix get_triu(const csr_matrix &matrix);
    csr_matrix get_tril(const csr_matrix &matrix); 
    csr_matrix get_diag(const csr_matrix &matrix);
    void print(const csr_matrix &matrix);
    bool is_symmetric(const csr_matrix &matrix);
    csr_matrix make_poisson_1d(const int n);
    // Do I need norm or condition number?
}

