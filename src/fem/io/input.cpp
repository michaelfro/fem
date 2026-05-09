#include <fem/io/input.hpp>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>  
#include <stdexcept>
#include <iostream>

using fem::mesh::triangular_2d::Mesh;

namespace fem::io{

    Mesh read_mesh_2d (const std::string& filename) {
        
        std::vector<std::array<double,2>> nodes;
        std::vector<std::array<int,3>> elements;

        std::ifstream file(filename);
        if (!file) {
            throw std::runtime_error("Could not open file");
        }

        std::string line;

        while(std::getline(file,line)){
            if (line.find("POINTS") == 0){
                std::stringstream ss(line);

                std::string temp;
                int number_points;

                ss >> temp >> number_points >> temp; 

                for (int i = 0; i < number_points; ++i) {
                    double x, y, z;
                    file >> x >> y >> z; // z values are in the vtk file but we ignore them for our mesh
                    nodes.push_back({x,y});
                }
            }

            if (line.find("POLYGONS") == 0) {
                std::stringstream ss(line);

                std::string temp; 
                int number_polys;

                ss >> temp >> number_polys >> temp; 

                for (int i = 0; i < number_polys; ++i) {
                    int a, b, c;
                    file >> a >> b >> c;
                    elements.push_back({a,b,c});
                }

            }
        }


        Mesh mesh;
        mesh.nodes = nodes;
        mesh.elements = elements;
        return mesh;
    }

}
