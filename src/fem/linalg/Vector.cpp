#include <fem/linalg/Vector.hpp>
#include <cassert>
#include <cstddef>
#include <vector>
#include <iostream>
#include <cmath>

namespace fem::linalg::vector{

    // Computes the vector sum res = w1 * v1 + w2 * v2. res, v1 and v2 must be the same size.
    void vector_sum(double w1, const std::vector<double> &v1, double w2, const std::vector<double> &v2, std::vector<double> &res) {
        assert(res.size() == v1.size());
        assert(res.size() == v2.size());

        for(std::size_t i = 0; i < res.size(); ++i){
            res[i] = w1*v1[i] + w2*v2[i]; 
        }
    }

    // Computes the inner product v'w. v and w must be the same size.
    double inner_product(const std::vector<double> &v, const std::vector<double> &w) {
        assert(v.size() == w.size());

        double result = 0.0;
        for(std::size_t i = 0; i < v.size(); ++i){
            result += v[i]*w[i];
        }
        return result;
    }

    double l2_norm(const std::vector<double> &v) {
        return std::sqrt(inner_product(v,v));
    }

    void print(const std::vector<double> &v) {
        for (double d : v)
            std::cout << d << " ";
        std::cout << "\n";
    }
}
