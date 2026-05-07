#include "crs_matrix.h"
#include <cmath>
#include <cassert>


namespace fem::linalg::crs_matrix{

    /// Computes the matrix product result = matrix * x for a sparse matrix.
    void matrix_vector_product(const crs_matrix &matrix, const std::vector<double> &x, std::vector<double> &result) {
        assert(x.size() == result.size());
        assert(x.size() == matrix.row_index.size() - 1);
        
        //loop over rows
        for(std::size_t i = 0; i < x.size(); ++i){
            result[i] = 0.0;

            const std::size_t start_index = matrix.row_index[i];
            const std::size_t end_index = matrix.row_index[i + 1];

            //loop over used non-zero entries of row i
            for(std::size_t j = start_index; j < end_index; ++j){
                result[i] += matrix.values.at(j) * x[matrix.col_index[j]];
            }
        }
    }


    void get_triu(const crs_matrix &matrix, crs_matrix &U){
        return
    }


    void get_tril(const crs_matrix &matrix, crs_matrix &L){
        return
    }
}
