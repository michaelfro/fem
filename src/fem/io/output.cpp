#include <fem/io/output.hpp>
#include <fem/linalg/cg.hpp>
#include <fem/mesh/mesh.hpp>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <filesystem>

using Function_2D = double (*) (double, double);

namespace fem::io{
    
    void write_solution_csv_2d(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh, const std::string& filename) { // this all depends on dimension - how to generalize this
        
        std::filesystem::path output_dir = "../output";

        if (!std::filesystem::exists(output_dir)) {
            std::filesystem::create_directories(output_dir);
        }

        std::filesystem::path file_path = output_dir / filename;

        std::ofstream file(file_path);

        if (!file.is_open()) {
           throw std::runtime_error("Could not open file");
        }
        file << "index,x,y,u\n"; 
        for (size_t i = 0; i < solution.solution.size(); ++i) {
            file << i << "," << mesh.nodes[i][0] << "," << mesh.nodes[i][1] << "," << solution.solution[i] << "\n";
        }

        file.close();
    }
    
    void write_error(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh, Function_2D analytical_solution) { // might need a fstream argument?
        return; 
    } 
    
    void print_solution(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh) {
        return;
    }
}
