#pragma once
#include<vector>
#include<cassert>

namespace fem::linalg{ 
    void vector_sum(std::vector<double> &out, double w1, const std::vector<double> &v1, double w2, const std::vector<double> &v2); 
    double inner_product(const std::vector<double> &v, const std::vector<double> &w); 
    double l2_norm(const std::vector<double> &v); 
}
