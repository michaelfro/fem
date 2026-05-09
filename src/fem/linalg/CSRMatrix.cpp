#include <fem/linalg/CSRMatrix.hpp>
#include <cmath>
#include <cassert>
#include <iostream>
#include <cstddef>
#include <vector>


namespace fem::linalg{

    /// Computes the matrix product result = matrix * x for a sparse matrix.
    void matrix_vector_product(const csr_matrix &matrix, const std::vector<double> &x, std::vector<double> &result) {
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

    csr_matrix get_tril(const csr_matrix &matrix){
        std::vector<double> values;
        std::vector<size_t> row_index;
        std::vector<size_t> col_index;

        row_index.push_back(0);

        //loop over rows
        for(std::size_t i = 0; i < matrix.row_index.size() - 1; ++i){

            const std::size_t start_index = matrix.row_index[i];
            const std::size_t end_index = matrix.row_index[i + 1];
            int counter = 0;
            //loop over non-zero entries of row i
            for(std::size_t j = start_index; j < end_index; ++j){
               
               if (i > matrix.col_index[j]){
                  values.push_back(matrix.values[j]);
                  col_index.push_back(matrix.col_index[j]);
                  counter += 1; 
               }
            }
            row_index.push_back(row_index[i] + counter);
        }


        csr_matrix U; 
        U.values = values;
        U.row_index = row_index; 
        U.col_index = col_index;

        return U;
    }

    csr_matrix get_triu(const csr_matrix &matrix){
        std::vector<double> values;
        std::vector<size_t> row_index;
        std::vector<size_t> col_index;

        row_index.push_back(0);

        //loop over rows
        for(std::size_t i = 0; i < matrix.row_index.size() - 1; ++i){

            const std::size_t start_index = matrix.row_index[i];
            const std::size_t end_index = matrix.row_index[i + 1];
            int counter = 0;
            //loop over non-zero entries of row i
            for(std::size_t j = start_index; j < end_index; ++j){
               
               if (i < matrix.col_index[j]){
                  values.push_back(matrix.values[j]);
                  col_index.push_back(matrix.col_index[j]);
                  counter += 1; 
               }
            }
            row_index.push_back(row_index[i] + counter);
        }

               // matrix.values.at(j);

        csr_matrix U; 
        U.values = values;
        U.row_index = row_index; 
        U.col_index = col_index;

        return U;
    }

    
    csr_matrix get_diag(const csr_matrix &matrix) {
        std::vector<double> values;
        std::vector<size_t> row_index;
        std::vector<size_t> col_index;

        row_index.push_back(0);

        //loop over rows
        for(std::size_t i = 0; i < matrix.row_index.size() - 1; ++i){

            const std::size_t start_index = matrix.row_index[i];
            const std::size_t end_index = matrix.row_index[i + 1];
            int counter = 0;
            //loop over non-zero entries of row i
            for(std::size_t j = start_index; j < end_index; ++j){
               if (i == matrix.col_index[j]){
                  values.push_back(matrix.values[j]);
                  col_index.push_back(matrix.col_index[j]);
                  counter += 1; 
               }
            }
            row_index.push_back(row_index[i] + counter);
        }


        csr_matrix U; 
        U.values = values;
        U.row_index = row_index; 
        U.col_index = col_index;

        return U;

    }

    void print_csr_matrix(const csr_matrix &matrix){
          // Infer the number of rows and columns assuming it is a square matrix
          size_t rows = matrix.row_index.size() - 1;
          size_t cols = rows;

          // Print the dense matrix directly without storing it
          // std::cout << "Dense Matrix Representation:\n";
          for (size_t i = 0; i < rows; ++i)
          {
              size_t row_start = matrix.row_index[i];
              size_t row_end = matrix.row_index[i + 1];
              size_t current_col = 0;

              for (size_t j = row_start; j < row_end; ++j)
              {
                  size_t col = matrix.col_index[j];
                  while (current_col < col)
                  {
                      std::cout << "0 ";
                      ++current_col;
                  }
                  std::cout << matrix.values[j] << " ";
                  ++current_col;
              }

              while (current_col < cols)
              {
                  std::cout << "0 ";
                  ++current_col;
              }
              std::cout << "\n";
          }
    }


    bool is_symmetric(const csr_matrix &matrix) {
        size_t rows = matrix.row_index.size() - 1;

        for (size_t i = 0; i < rows; ++i) {
            size_t row_start = matrix.row_index[i];
            size_t row_end = matrix.row_index[i + 1];

            for (size_t j = row_start; j < row_end; ++j) {
                size_t col = matrix.col_index[j];
                double value = matrix.values[j];

                // Check if the corresponding entry exists and matches
                size_t col_start = matrix.row_index[col];
                size_t col_end = matrix.row_index[col + 1];
                bool found = false;

                for (size_t k = col_start; k < col_end; ++k) {
                    if (matrix.col_index[k] == i) {
                        // Use a relative tolerance for floating-point comparison
                        if (std::abs(matrix.values[k] - value) > 1e-12 * std::max(std::abs(matrix.values[k]), std::abs(value))) {
                            return false; // Values do not match
                        }
                        found = true;
                        break;
                    }
                }

                if (!found) {
                    return false; // Symmetric entry not found
                }
            }
        }

        return true;
    }

    csr_matrix make_poisson_1d(int n) {
        csr_matrix A;
        A.row_index.push_back(0);

        for (int i = 0; i < n; ++i) { 
            if (i > 0) {
                A.col_index.push_back(i - 1);
                A.values.push_back(-1.0);
            }

            A.col_index.push_back(i);
            A.values.push_back(2.0);

            if (i < n - 1) {
                A.col_index.push_back(i + 1);
                A.values.push_back(-1.0);
            }

            A.row_index.push_back(A.col_index.size());
        }

        return A;
    }
}
