#pragma once
#include <fem/linalg/cg.hpp>
#include <fem/mesh/mesh.hpp>

namespace fem::io{
    void write_solution(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh);
    void print_solution(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh);
}
