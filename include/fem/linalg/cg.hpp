#pragma once

#include <fem/linalg/CSRMatrix.hpp>
#include <vector>
#include <cstddef>

namespace fem::linalg::cg{
    struct solution {
        std::vector<double> solution;
        int iterations;
        double final_residuum;  
    };

    
    solution CG(const fem::linalg::csr_matrix &matrix, const std::vector<double> &rhs,
                            const size_t max_iterations, double reduction_factor);

    // solution PCG(const fem::linalg::csr_matrix &matrix, const std::vector<double> &rhs,
    //                         const std::vector<double> &x_0, const fem::linalg::csr_matrix &precondtioner,
    //                         const size_t max_iterations, const double reduction_factor);
    
    template<typename Preconditioner>
    solution PCG(const fem::linalg::csr_matrix &matrix, const std::vector<double> &rhs,
                            const std::vector<double> &x_0, Preconditioner &apply_precondtioner,
                            const size_t max_iterations, const double reduction_factor);

    void test_cg();
    
    // void test_pcg();
}
