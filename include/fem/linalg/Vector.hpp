#pragma once
#include<vector>
#include<cassert>

namespace fem::linalg::vector{ 
    void vector_sum(double w1, const std::vector<double> &v1, double w2, const std::vector<double> &v2,std::vector<double> &res); 
    double inner_product(const std::vector<double> &v, const std::vector<double> &w); 
    double l2_norm(const std::vector<double> &v); 
    void print(const std::vector<double> &v);
}
