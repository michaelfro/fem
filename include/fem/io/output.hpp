#pragma once
#include <fem/linalg/cg.hpp>
#include <fem/mesh/mesh.hpp>
#include <string>

using Function_2D = double (*) (double, double);

namespace fem::io{
    void write_solution_csv_2d(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh, const std::string& filename);
    void print_solution(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh);
    void write_error(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh, Function_2D analytical_solution);

}
