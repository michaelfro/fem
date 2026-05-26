#pragma once

#include <array>
#include <fem/mesh/mesh.hpp>
#include <string>
#include <vector>

/**
 * @file input.hpp
 * @brief Mesh input utilities.
 */

/**
 * @namespace fem::io
 * @brief File input and output helpers.
 */
namespace fem::io {

/**
 * @brief Reads a two-dimensional triangular mesh from an ASCII VTK PolyData file.
 *
 * The reader consumes the `POINTS` section and triangular entries from the
 * `POLYGONS` section. The third point coordinate in the VTK file is ignored.
 *
 * @param filename Path to the VTK file.
 * @return Mesh containing node coordinates and triangle connectivity.
 * @throws std::runtime_error If the file cannot be opened.
 */
fem::mesh::triangular_2d::Mesh read_mesh_2d(const std::string &filename);

} // namespace fem::io
