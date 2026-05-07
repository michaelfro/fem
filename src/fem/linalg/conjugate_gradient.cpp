#include "conjugate_gradient.h"
#include <cmath>
#include <cassert>

/// Computes the vector sum out = v1 + alpha * v2. out, v1 and v2 must be the same size.
void vector_sum(std::vector<double> &out, double w1, const std::vector<double> &v1, 
                double w2, const std::vector<double> &v2) {

    if (v1.size() != v2.size()){
        return;
    }
    const int N = v1.size();
    for (int i = 0; i < N; i++){

        out[i] = w1 * v1[i] + w2 * v2[i];

    }

}

/// Computes the inner product v'w. v and w must be the same size.
double inner_product(const std::vector<double> &v, const std::vector<double> &w) {

    double res = 0;

    assert (v.size() == w.size());

    for (int i = 0; i < v.size(); i++){
        res += v[i] * w[i];
    }

    return res;
}

/// Computes the matrix product result = matrix * x for a sparse matrix.
void matrix_vector_product(const sparse_matrix &matrix, const std::vector<double> &x, std::vector<double> &result) {

    const int dim = x.size();

    // Initialize result with correct size
    result.clear();
    for (int i = 0; i < dim; i++){
        result.push_back(0.);
    }

    for (int i = 0; i < dim; i++){
        
        const int nmbr_entries = matrix.row_index[i+1] - matrix.row_index[i]; // Number of entries in each row of the matrix
        for(int j = 0; j < nmbr_entries; j++){
            int my_index = matrix.row_index[i] + j; // I'm sure there is a prettier way, sorry
            result[i] += matrix.values[my_index] * x[matrix.col_index[my_index]];

        }
        
    }

    return;

}

/**
 * Solves the linear system Ax = b with the CG method.
 */
solution conjugate_gradient(const sparse_matrix &matrix, const std::vector<double> &rhs,
                            const size_t max_iterations, double reduction_factor) {

    std::vector<double> x_0;
    std::vector<double> r_0;
    std::vector<double> p;
    std::vector<double> x;
    std::vector<double> r;

    // I chose as an inital guess is the zero vector in the correct dimsension:
    // this implies that r_0 is just b:

    for (int i = 0; i < rhs.size(); i++){
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

    for(int k = 0; k < max_iterations; k++){

        std::vector<double> temp; // = Ap
        matrix_vector_product(matrix,p,temp);

        alpha = rTr / inner_product(p,temp);

        vector_sum(x,1,x,alpha,p); 
        vector_sum(r,1.,r,-1. * alpha,temp);

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

        vector_sum(p,1,r,beta,p);
    }

    solution my_sol;
    my_sol.solution = x;
    my_sol.iterations = max_iterations;
    my_sol.final_residuum = error_norm;
    return my_sol;

}

