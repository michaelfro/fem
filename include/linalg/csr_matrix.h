#pragma once

#include <vector>
#include <cstddef>

namespace fem::linalg::crs_matrix{

    struct csr_matrix {
        std::vector<double> values;
        std::vector<size_t> col_index;
        std::vector<size_t> row_index;
        // int size; // do I need this?
    };

    void matrix_vector_product(const crs_matrix &matrix, const std::vector<double> &x, std::vector<double> &result); 
    void get_triu(const crs_matrix &matrix, crs_matrix &U);
    void get_tril(const crs_matrix &matrix, crs_matrix &L);
    // Do I need norm or condition number?
}

