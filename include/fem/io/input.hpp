#include <vector>
#include <array>
#include <string>
#include <fem/mesh/mesh.hpp>

namespace fem::io {

  fem::mesh::triangular_2d::Mesh read_mesh_2d(const std::string &filename);

}
