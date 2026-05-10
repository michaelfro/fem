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
    /**
     * Solves the linear system Ax = b with the CG method.
     */
    solution CG(const csr_matrix &matrix, const std::vector<double> &rhs,
                                const size_t max_iterations, double reduction_factor) {

        std::vector<double> x_0;
        std::vector<double> r_0;
        std::vector<double> p;
        std::vector<double> x;
        std::vector<double> r;

        // I chose as an inital guess is the zero vector in the correct dimension:
        // this implies that r_0 is just b:

        for (size_t i = 0; i < rhs.size(); i++){
            x_0.push_back(0.);
            r_0.push_back(rhs[i]); // r_0 = b = rhs
            p.push_back(rhs[i]); // p_0 = r_0
        }

        x = x_0;
        r = r_0;

        double error_norm = std::sqrt(inner_product(r_0,r_0));

        if (error_norm < reduction_factor){
            solution my_sol;
            my_sol.solution = x_0;
            my_sol.iterations = 0;
            my_sol.final_residuum = error_norm;
            return my_sol;
        }

        double alpha;
        double beta;
        double rTr = inner_product(r,r); // = r.Tr
        double rTr_new;

        for(size_t k = 0; k < max_iterations; k++){

            std::vector<double> temp(x_0.size(),0.); // = Ap
            matrix_vector_product(matrix,p,temp);

            double pos_def_check = inner_product(p,temp);

            // if (pos_def_check < 0) {
            //     std::cout << "CG: Matrix is not positive definite\n";
            // }

            alpha = rTr / pos_def_check;

            fem::linalg::vector::vector_sum(1,x,alpha,p,x); 
            fem::linalg::vector::vector_sum(1.,r,-1. * alpha,temp,r);

            rTr_new = inner_product(r,r);
            error_norm = std::sqrt(rTr_new);
            if (error_norm < reduction_factor){ // covergence check:
                solution my_sol;
                my_sol.solution = x;
                my_sol.iterations = k + 1;
                my_sol.final_residuum = error_norm;
                return my_sol;
            }

            beta = rTr_new / rTr;
            rTr = rTr_new;

            fem::linalg::vector::vector_sum(1,r,beta,p,p);
        }

        solution my_sol;
        my_sol.solution = x;
        my_sol.iterations = max_iterations;
        my_sol.final_residuum = error_norm;
        return my_sol;

    }


    /**
     * Solves the linear system Ax = b with the CG method.
     */
    solution PCG(const csr_matrix &matrix, const std::vector<double> &rhs,
                 const std::vector<double> &x_0, const csr_matrix &preconditioner, 
                     const size_t max_iterations, const double reduction_factor) {

        std::vector<double> p;
        std::vector<double> x;
        std::vector<double> r(x_0.size(),0.);
        std::vector<double> Ax(x_0.size(),0.);
        std::vector<double> z(x_0.size(),0.);
        x = x_0;

        matrix_vector_product(matrix,x_0,Ax);
        fem::linalg::vector::vector_sum(1.,rhs,-1.,Ax,r);
        matrix_vector_product(preconditioner,r,z);
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
            matrix_vector_product(preconditioner,r,z);
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


    void test_cg() {
        // Define the 4x4 matrix A in CSR format
        csr_matrix A;
        A.values = {4, -1, -1, 4, -1, -1, 4, -1, -1, 4};
        A.col_index = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
        A.row_index = {0, 2, 5, 8, 10};
        std::cout << "Matrix A:" << std::endl;
        fem::linalg::print_csr_matrix(A);
        // Define the right-hand side vector b
        std::vector<double> b = {1, 2, 2, 1};

        // Call the CG solver
        solution sol = CG(A, b, 100, 10e-12);

        // Output the solution
        std::cout << "Solution vector x:" << std::endl;
        fem::linalg::vector::print(sol.solution);

        // Output the number of iterations and final residual
        std::cout << "CG iterations: " << sol.iterations << std::endl;
        std::cout << "Final residual: " << sol.final_residuum << std::endl;

    }


    void test_pcg() {
        int n = 100;
        csr_matrix A = fem::linalg::make_poisson_1d(n);
        // print(A);
        std::vector<double> x_true(n, 1.0);

        // Define the right-hand side vector b
        std::vector<double> b(n, 1.0);
        matrix_vector_product(A,x_true,b);

        // Define the preconditioner M^-1 in CSR format:
        csr_matrix M;

        M.row_index.resize(n+1);
        M.col_index.resize(n);
        M.values.resize(n);

        for (int i = 0; i < n; ++i) {
            M.row_index[i] = i;
            M.col_index[i] = i;
            M.values[i] = 0.5; // Jacobi preconditioner, Diagonal is 2, so 1/2
            // M.values[i] = 1.; // identity preconditioner
        }

        M.row_index[n] = n;

        std::vector<double> x_0(n, 0.0);

        // Call the PCG solver
        solution sol = PCG(A, b, x_0, M, 100, 10e-12);

        // Output the solution
        std::cout << "Solution vector x:" << std::endl;
        fem::linalg::vector::print(sol.solution);

        // Output the number of iterations and final residual
        std::cout << "PCG iterations: " << sol.iterations << std::endl;
        std::cout << "Final residual: " << sol.final_residuum << std::endl;

    }

}
