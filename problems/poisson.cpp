#include <iostream>
#include <vector>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/io/MeshImporter.hpp>
#include <fem/linalg/CG.hpp>

int main(){ 

  fem::mesh::triangular_2d::Mesh mesh = fem::mesh::triangular_2d::read("../data/meshes/square.vtk"); // assumes you run the code from /build
  std::cout << "Number of nodes in the mesh: " << mesh.nodes.size() << std::endl;

  return 0;
}
