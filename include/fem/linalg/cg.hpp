#pragma once

#include <fem/linalg/CSRMatrix.hpp>
#include <vector>
#include <cstddef>
#include <fem/linalg/cg.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/Vector.hpp>
#include <cstddef> 
#include <cmath>
#include <iostream>

using fem::linalg::csr_matrix;
using fem::linalg::matrix_vector_product;
using fem::linalg::vector::inner_product;

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
    
    /**
     * Solves the linear system Ax = b with the preconditioned CG method.
     */
    template<typename Preconditioner>
    solution PCG(const csr_matrix &matrix, const std::vector<double> &rhs,
                 const std::vector<double> &x_0, Preconditioner & apply_preconditioner, 
                     const size_t max_iterations, const double reduction_factor) {

        std::vector<double> p;
        std::vector<double> x;
        std::vector<double> r(x_0.size(),0.);
        std::vector<double> Ax(x_0.size(),0.);
        std::vector<double> z(x_0.size(),0.);
        x = x_0;

        matrix_vector_product(matrix,x_0,Ax);
        fem::linalg::vector::vector_sum(1.,rhs,-1.,Ax,r);
        z = apply_preconditioner(r);
        p = z;

        double error_norm = fem::linalg::vector::l2_norm(r);

        if (error_norm < reduction_factor){
            solution my_sol;
            my_sol.solution = x_0;
            my_sol.iterations = 0;
            my_sol.final_residuum = error_norm;
            return my_sol;
        }

        double alpha;
        double beta;
        double rTz = inner_product(r,z); 
        double rTz_new;

        for(size_t k = 0; k < max_iterations; k++){

            std::vector<double> Ap(x_0.size(),0.); // = Ap
            matrix_vector_product(matrix,p,Ap);

            alpha = rTz / inner_product(p,Ap);

            fem::linalg::vector::vector_sum(1,x,alpha,p,x); 
            fem::linalg::vector::vector_sum(1.,r,-1. * alpha,Ap,r);
            z = apply_preconditioner(r);
            rTz_new = inner_product(r,z);
            error_norm = std::sqrt(rTz_new);
            if (error_norm < reduction_factor){ // covergence check:
                solution my_sol;
                my_sol.solution = x;
                my_sol.iterations = k + 1;
                my_sol.final_residuum = error_norm;
                return my_sol;
            }

            beta = rTz_new / rTz;
            rTz = rTz_new;

            fem::linalg::vector::vector_sum(1,z,beta,p,p);
        }

        solution my_sol;
        my_sol.solution = x;
        my_sol.iterations = max_iterations;
        my_sol.final_residuum = error_norm;
        return my_sol;

    }


    void test_cg();
    
    // void test_pcg();
}
