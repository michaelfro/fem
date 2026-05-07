#include"linalg/vector_ops.h"

namespace fem::linalg{


    // Computes the vector sum out = w1*v1 + w2 * v2. out, v1 and v2 must be the same size.
    void vector_sum(std::vector<double> &out, double w1, const std::vector<double> &v1, double w2, const std::vector<double> &v2) {
        assert(out.size() == v1.size());
        assert(out.size() == v2.size());

        for(std::size_t i = 0; i < out.size(); ++i){
            out[i] = w1*v1[i] + w2*v2[i]; 
        }
    }

    /// Computes the inner product v'w. v and w must be the same size.
    double inner_product(const std::vector<double> &v, const std::vector<double> &w) {
        assert(v.size() == w.size());

        double result = 0.0;
        for(std::size_t i = 0; i < v.size(); ++i){
            result += v[i]*w[i];
        }
        return result;
    }

    double l2_norm(const std::vector<double> &v) {
        return inner_product(v,v);
    }
}
