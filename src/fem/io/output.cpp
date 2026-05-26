#include <fem/io/output.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fem::io {
namespace {

std::filesystem::path make_output_path(const std::string &filename) {
    std::filesystem::path output_dir = "../output";

    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
    }

    return output_dir / filename;
}

void check_solution_size(const fem::linalg::cg::solution &solution,
                         const fem::mesh::triangular_2d::Mesh &mesh) {
    if (solution.solution.size() != mesh.nodes.size()) {
        throw std::runtime_error("Solution size does not match mesh node count");
    }
}

void write_vtk_header(std::ofstream &file, const fem::mesh::triangular_2d::Mesh &mesh) {
    file << "# vtk DataFile Version 3.0\n";
    file << "FEM solution\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS " << mesh.nodes.size() << " double\n";
    for (const auto &node : mesh.nodes) {
        file << node[0] << " " << node[1] << " 0\n";
    }

    file << "\nPOLYGONS " << mesh.elements.size() << " " << 4 * mesh.elements.size() << "\n";
    for (const auto &element : mesh.elements) {
        file << "3 " << element[0] << " " << element[1] << " " << element[2] << "\n";
    }
}

void write_scalar_field(std::ofstream &file, const std::string &name,
                        const std::vector<double> &values) {
    file << "SCALARS " << name << " double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (double value : values) {
        file << value << "\n";
    }
}

} // namespace

void write_solution_csv_2d(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh,
                           const std::string &filename) {
    check_solution_size(solution, mesh);

    std::ofstream file(make_output_path(filename));

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    file << "index,x,y,u\n";
    for (size_t i = 0; i < solution.solution.size(); ++i) {
        file << i << "," << mesh.nodes[i][0] << "," << mesh.nodes[i][1] << ","
             << solution.solution[i] << "\n";
    }
}

void write_solution_vtk_2d(const fem::linalg::cg::solution &solution,
                           const fem::mesh::triangular_2d::Mesh &mesh,
                           const std::string &filename) {
    check_solution_size(solution, mesh);

    std::ofstream file(make_output_path(filename));
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    write_vtk_header(file, mesh);
    file << "\nPOINT_DATA " << mesh.nodes.size() << "\n";
    write_scalar_field(file, "solution", solution.solution);
}

void write_error(const fem::linalg::cg::solution &solution,
                 const fem::mesh::triangular_2d::Mesh &mesh, Function_2D analytical_solution,
                 const std::string &filename) {
    check_solution_size(solution, mesh);

    std::vector<double> exact(mesh.nodes.size());
    std::vector<double> error(mesh.nodes.size());
    std::vector<double> abs_error(mesh.nodes.size());

    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
        exact[i] = analytical_solution(mesh.nodes[i][0], mesh.nodes[i][1]);
        error[i] = solution.solution[i] - exact[i];
        abs_error[i] = std::abs(error[i]);
    }

    std::ofstream file(make_output_path(filename));
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    write_vtk_header(file, mesh);
    file << "\nPOINT_DATA " << mesh.nodes.size() << "\n";
    write_scalar_field(file, "solution", solution.solution);
    write_scalar_field(file, "analytical_solution", exact);
    write_scalar_field(file, "error", error);
    write_scalar_field(file, "absolute_error", abs_error);
}

void write_error(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh,
                 Function_2D analytical_solution) {
    write_error(solution, mesh, analytical_solution, "solution_error.vtk");
}

void print_solution(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh) {
    check_solution_size(solution, mesh);

    for (size_t i = 0; i < solution.solution.size(); ++i) {
        std::cout << i << ": (" << mesh.nodes[i][0] << ", " << mesh.nodes[i][1] << ") -> "
                  << solution.solution[i] << "\n";
    }
}

} // namespace fem::io
