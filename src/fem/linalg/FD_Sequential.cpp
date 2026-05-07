#include <iostream>
#include <iomanip>
#include <cmath>
#include <math.h> 
#include <vector>
#include <chrono>
#include <cassert>

#define _USE_MATH_DEFINES

/**
 * @brief Structure holding a mtrix in CRS format
 */
struct crs_matrix
{
    std::vector<double> values;
    std::vector<size_t> col_index;
    std::vector<size_t> row_index;
};

/**
 * @brief A structure holding the information about the result from the solver
 */
struct solution
{
    std::vector<double> solution;
    size_t iterations;
    double final_residuum, initial_residuum;
};

/**
 * @brief A function that compute the value of f at the point (x,y)
 * @return double
 */
double rhs_function(const double x, const double y)
{
    return 2.0 * M_PI*M_PI*std::sin(M_PI * x) * std::cos(M_PI * y);
}

/**
 * @brief A function that compute the value of g at the point (x,y)
 * @return double
 */
double dirichlet_bc(const double x, const double y)
{

    const double tol = 1e-12; // Tolerance for floating-point comparisons

    if (std::abs(x - 1.0) <= tol) {
        return 0.0; // u(1, y) = 0
    }
    if (std::abs(x - 0.0) <= tol) {
        return 0.0; // u(0, y) = 0
    }
    if (std::abs(y - 1.0) <= tol) {
        return -std::sin(M_PI * x); // u(x, 1) = -sin(pi * x)
    }
    if (std::abs(y - 0.0) <= tol) {
        return std::sin(M_PI * x); // u(x, 0) = sin(pi * x)
    }

    throw std::logic_error("Used dirichlet_bc for non-boundary point");
    
}

/**
 * @brief A function that compute the analytical solution at the point (x,y)
 * @return double
 */
double analytical_solution(const double x, const double y)
{
    return std::sin(M_PI * x) * std::cos(M_PI * y);
}


/**
 * @brief A function that assembles the stiffness matrix for a system with a uniform grid
 * and same number of nodes in each dimension
 *
 * @param N number of elements in one dimension
 * @return crs_matrix
 */
crs_matrix assemble_FD_stiffness_matrix(const size_t N) {
    const std::size_t system_size = (N+1)*(N+1);
    std::vector<double> values;
    std::vector<size_t> col_index;
    std::vector<size_t> row_index(system_size + 1, 0);

    for (size_t i = 0; i <= N; ++i) {
        for (size_t j = 0; j <= N; ++j) {
            size_t node_idx = i * (N+1) + j;

            // Check if the node is on the boundary
            if (i == 0 || i == N || j == 0 || j == N) { 
                // Boundary node: Set diagonal to 1
                values.push_back(1.0);
                col_index.push_back(node_idx);
                row_index[node_idx + 1] = values.size();
                continue;
            }

            if (i > 0) {
                values.push_back(-1);
                col_index.push_back(node_idx - (N+1));
            }
            if (j > 0) {
                values.push_back(-1);
                col_index.push_back(node_idx - 1);
            }

            values.push_back(4);
            col_index.push_back(node_idx);

            if (j < N) {
                values.push_back(-1);
                col_index.push_back(node_idx + 1);
            }
            if (i < N) {
                values.push_back(-1);
                col_index.push_back(node_idx + (N+1));
            }

            row_index[node_idx + 1] = values.size();
        }
    }

    return {values, col_index, row_index};
}

crs_matrix assemble_reduced_stiffness_matrix(const size_t N) {
    const std::size_t system_size = (N-1)*(N-1);
    std::vector<double> values;
    std::vector<size_t> col_index;
    std::vector<size_t> row_index(system_size + 1, 0);

    for (size_t i = 0; i < N - 1; ++i) {
        for (size_t j = 0; j < N - 1; ++j) {
            size_t node_idx = i * (N-1) + j;

            if (i > 0) {
                values.push_back(-1);
                col_index.push_back(node_idx - (N-1));
            }

            if (j > 0) {
                values.push_back(-1);
                col_index.push_back(node_idx - 1);
            }

            values.push_back(4);
            col_index.push_back(node_idx);

            if (j < N-2) {
                values.push_back(-1);
                col_index.push_back(node_idx + 1);
            }

            if (i < N-2) {
                values.push_back(-1);
                col_index.push_back(node_idx + (N-1));
            }

            row_index[node_idx + 1] = values.size(); 

        }
    }

    return {values, col_index, row_index};
}

/**
 * @brief A function that assembles the right-hand side vector for a system with a uniform grid
 * and same number of nodes in each dimension
 *
 * @param N number of elements in one dimension
 * @return std::vector<double>
 */
std::vector<double> assemble_rhs(const size_t N)
{
    const double h = 1.0 / N;
    const std::size_t system_size = (N+1)*(N+1);
    std::vector<double> rhs(system_size, 0.0);

    for (size_t i = 0; i <= N; ++i) {
        for (size_t j = 0; j <= N; ++j) {
            size_t node_idx = i * (N+1) + j;

            if (i == 0 || i == N || j == 0 || j == N) {
                rhs[node_idx] = dirichlet_bc(j * h, i * h);
            } else {
                rhs[node_idx] = h * h * rhs_function(j * h, i * h);
            }
        }
    }
    return rhs;
}

std::vector<double> assemble_reduced_rhs(const size_t N)
{
    const double h = 1.0 / N;
    const double h_squared = h*h;

    const std::size_t system_size = (N-1)*(N-1);

    std::vector<double> rhs(system_size, 0.0);

    for (size_t i = 1; i < N; ++i) { // Loop over inner rows
        for (size_t j = 1; j < N; ++j) { // Loop over inner columns
            size_t node_idx = (i - 1) * (N - 1) + (j - 1); // Map to reduced system index

            // Compute the contribution from the right-hand side function
            rhs[node_idx] = h_squared * rhs_function(j * h, i * h);

            // Adjust for neighbors on the boundary using Dirichlet boundary conditions
            if (i == 1) { // Bottom boundary
                rhs[node_idx] += dirichlet_bc(j * h, 0.0);
            }
            if (i == N - 1) { // Top boundary
                rhs[node_idx] += dirichlet_bc(j * h, 1.0);
            }
            if (j == 1) { // Left boundary
                rhs[node_idx] += dirichlet_bc(0.0, i * h);
            }
            if (j == N - 1) { // Right boundary
                rhs[node_idx] += dirichlet_bc(1.0, i * h);
            }
        }
    }

    return rhs;
}

/**
 * @brief A function that computes a matrix vector multiplication between a crs matrix and a vector
 *
 * @param matrix
 * @param x
 * @param result
 */
void sparse_matrix_vector_mult(const crs_matrix &matrix, const std::vector<double> &x, std::vector<double> &result)
{

    for (size_t i = 0; i < x.size(); ++i) {
        size_t start_index = matrix.row_index[i], end_index = matrix.row_index[i + 1];
        result[i] = 0.0;
        for (size_t j = start_index; j < end_index; ++j) {
            result[i] += matrix.values[j] * x[matrix.col_index[j]];
        }
    }
}

void print_sparse_matrix_as_dense(const crs_matrix &matrix)
{
    // Infer the number of rows and columns assuming it is a square matrix
    if(not (matrix.row_index.size() > 1)){
        return;
    }
    size_t rows = matrix.row_index.size() - 1;
    size_t cols = rows;

    // Print the dense matrix directly without storing it
    std::cout << "Dense Matrix Representation:\n";
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


/**
 * @brief A function that computes the dot product of two vectors
 *
 * @param vec_1
 * @param vec_2
 * @return double
 */
double dot_product(const std::vector<double> &vec_1, const std::vector<double> &vec_2)
{
    // assert vec_1.size() == vec_2.size() // might be slowing down the code
    double res = 0; 
    for(std::size_t i = 0; i < vec_1.size(); ++i) {
        res += vec_1[i] * vec_2[i];
    }
    return res;
}

/**
 * @brief A function that computes the Euclidean norm of a vector
 *
 * @param x
 * @return double
 */
double euclidean_norm(const std::vector<double> &x)
{
    double sum = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        sum += x[i] * x[i];
    }
    return std::pow(sum,1./2.);
}

/**
 * @brief A function that does a scaled vector addition between two vectors and returns the result
 * vec1 = alpha_1*vec1 + alpha_2*vec2
 *
 * @param vec_1
 * @param vec_2
 * @param alpha_1
 * @param alpha_2
 */

std::vector<double> scaled_vector_addition (std::vector<double> &vec_1, const std::vector<double> &vec_2, const double alpha_1, const double alpha_2)
{
    //  asserts ?
    std::vector<double> res(vec_1.size(),0);
    for (size_t i = 0; i < vec_1.size(); ++i) {
        res[i] = alpha_1 * vec_1[i] + alpha_2 * vec_2[i];
    }
    return res;
}

/**
 * @brief A function that does an inplace scaled vector addition between two vectors
 * vec1 = alpha_1*vec1 + alpha_2*vec2
 *
 * @param vec_1
 * @param vec_2
 * @param alpha_1
 * @param alpha_2
 */
void scaled_vector_addition_inplace(std::vector<double> &vec_1, std::vector<double> &vec_2, const double alpha_1, const double alpha_2)
{
    // asserts?
    for (size_t i = 0; i < vec_1.size(); ++i) {
        vec_1[i] = alpha_1 * vec_1[i] + alpha_2 * vec_2[i];
    }

}

/**
 * @brief Sparse CG Solver
 *
 * @param A crs matrix
 * @param b right hand side
 * @param reduction_factor for the stopping condition
 * @param max_iterations
 * @return solution
 */
solution CG(crs_matrix &A, std::vector<double> &b, const double reduction_factor, const size_t max_iterations)
{    
    std::vector<double> p = b;
    std::vector<double> x (b.size(),0.0);
    std::vector<double> r = b;

    std::vector<double> Ap(b.size(), 0.0); 

    double bb = dot_product(b,b);

    double rTr, rTr_new, alpha, beta;
    rTr = dot_product(r,r);

    size_t counter = 1;
    bool converged = (rTr <= reduction_factor * reduction_factor * bb || counter >= max_iterations);
    
    while(!converged){
        
        sparse_matrix_vector_mult(A,p,Ap);
        double pAp = dot_product(Ap,p); 

        alpha = rTr / pAp;
        scaled_vector_addition_inplace(x,p,1,alpha); 
        scaled_vector_addition_inplace(r,Ap,1,-1. * alpha);

        rTr_new = dot_product(r,r);

        beta = rTr_new / rTr;
        rTr = rTr_new;
        scaled_vector_addition_inplace(p,r,beta,1);
        counter++;
        converged = (rTr_new <= reduction_factor * reduction_factor * bb || counter >= max_iterations);
    }
    
    return {x,counter,rTr_new};
}

/**
 * @brief A function that compute the discrete L2 norm of the error
 * between the approximate solution and the analytical solution
 *
 * @param N
 * @param approx_sol
 * @return double
 */
double discrete_l2_norm(const size_t N, const std::vector<double> &approx_sol)
{
    double h = 1.0 / N;
    double sum = 0;
    for (size_t i = 0; i < N+1; ++i) {
        for (size_t j = 0; j < N+1; ++j) {
            double x = j * h;
            double y = i * h;
            double analytical = analytical_solution(x, y);
            double numerical = approx_sol[i * (N+1) + j];
            double err = numerical - analytical;
            sum += err * err;
        }
    }
    return std::pow(sum * h * h, 1./2.);
}

double discrete_reduced_l2_norm(const size_t N, const std::vector<double> &approx_sol)
{
    double h = 1.0 / N;
    double sum = 0;
    for (size_t i = 1; i < N; ++i) {
        for (size_t j = 1; j < N; ++j) {
            double x = j * h;
            double y = i * h;
            double analytical = analytical_solution(x, y);
            double numerical = approx_sol[(i-1) * (N-1) + (j-1)];
            double err = numerical - analytical;
            // std::cout << "The error is " << err << "\n";
            sum += err * err;
        }
    }
    return std::pow(sum * h * h, 1./2.);
}

void test_cg_solver() {
    // Define the 4x4 matrix A in CRS format
    crs_matrix A;
    A.values = {4, -1, -1, 4, -1, -1, 4, -1, -1, 4};
    A.col_index = {0, 1, 0, 1, 2, 1, 2, 3, 2, 3};
    A.row_index = {0, 2, 5, 8, 10};
    std::cout << "Matrix A:" << std::endl;
    print_sparse_matrix_as_dense(A);
    // Define the right-hand side vector b
    std::vector<double> b = {1, 2, 2, 1};

    // Call the CG solver
    solution sol = CG(A, b, 1e-12, 100);

    // Output the solution
    std::cout << "Solution vector x:" << std::endl;
    for (size_t i = 0; i < sol.solution.size(); ++i) {
        std::cout << sol.solution[i] << " ";
    }
    std::cout << std::endl;

    // Output the number of iterations and final residual
    std::cout << "CG iterations: " << sol.iterations << std::endl;
    std::cout << "Final residual: " << sol.final_residuum << std::endl;
}

bool is_symmetric(const crs_matrix &matrix) {
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

int main(int argc, char **argv)
{
    // Choose N
    std::size_t N = atoi(argv[1]);
    
    // Assemble Stiffness Matrix
    crs_matrix A = assemble_reduced_stiffness_matrix(N);
    
    // debug prints:
    // print_sparse_matrix_as_dense(A);

    // if (is_symmetric(A)) {
    //     std::cout << "The matrix is symmetric." << std::endl;
    // } else {
    //     std::cout << "The matrix is not symmetric." << std::endl;
    // }

    // Assemble Rhs (with Dirichlet Boundary Conditions)
    std::vector<double> b = assemble_reduced_rhs(N);

    auto start_time = std::chrono::high_resolution_clock::now();

    // Solve the System

    solution sol = CG(A, b, 1e-12, N*N);

    // Stop timing the CG solver
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast <std::chrono::duration<double>> (end_time - start_time).count();

    // Compute the L2 error
    double l2_err = discrete_reduced_l2_norm(N,sol.solution);

    // Output relevent information
    std::cout << "After " << sol.iterations 
              << " iterations we are left with a L2 error of " 
              << l2_err << " for N = " << N  
              << " within "<< duration << " seconds." << std::endl;

    return 0;

}