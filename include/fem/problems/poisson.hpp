# pragma once

#include <iostream>
#include <vector>
#include <fem/io/input.hpp>
#include <fem/io/output.hpp>
#include <fem/mesh/mesh.hpp>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/cg.hpp>

namespace fem::problems::poisson {

    struct triplet {
        size_t global_row;
        size_t global_col; 
        double contribution;
    };

    struct local_result {
        std::vector<triplet> local_contributions;
        std::vector<double> rhs;
    };

    struct lin_system {
        fem::linalg::csr_matrix stiffness_matrix;
        std::vector<double> rhs;
    };


    double f(const double x, const double y);
    void impose_boundary_conditions(std::vector<double>& rhs, fem::linalg::csr_matrix& stiffness_matrix, const fem::mesh::triangular_2d::Mesh mesh);
    local_result local_contribution(const fem::mesh::triangular_2d::Mesh mesh);
    lin_system assemble_system(const fem::mesh::triangular_2d::Mesh mesh);
    double analytical_solution(const double x, const double y);
    double nodal_l2_error(const fem::mesh::triangular_2d::Mesh mesh, const std::vector<double> &approx_sol);
    int run_poisson_problem();
}
