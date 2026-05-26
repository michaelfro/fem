#pragma once

#include <fem/linalg/cg.hpp>
#include <fem/mesh/mesh.hpp>
#include <string>

/**
 * @file output.hpp
 * @brief Output helpers for FEM solutions and diagnostics.
 */

/**
 * @brief Function pointer type for scalar functions on two-dimensional coordinates.
 */
using Function_2D = double (*)(double, double);

namespace fem::io {

/**
 * @brief Writes a two-dimensional nodal solution to a CSV file.
 *
 * The output file contains the columns `index,x,y,u` and is written below the
 * `../output` directory relative to the process working directory.
 *
 * @param solution Solver result containing nodal values.
 * @param mesh Mesh whose node coordinates correspond to the solution entries.
 * @param filename Name of the CSV file to create.
 * @throws std::runtime_error If the output file cannot be opened.
 */
void write_solution_csv_2d(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh,
                           const std::string &filename);

/**
 * @brief Prints a nodal solution.
 *
 * @param solution Solver result containing nodal values.
 * @param mesh Mesh whose node coordinates correspond to the solution entries.
 *
 * @note This function is currently a placeholder.
 */
void print_solution(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh);

/**
 * @brief Writes or reports solution error against an analytical reference.
 *
 * @param solution Solver result containing nodal values.
 * @param mesh Mesh whose node coordinates correspond to the solution entries.
 * @param analytical_solution Function used as the exact solution.
 *
 * @note This function is currently a placeholder.
 */
void write_error(fem::linalg::cg::solution solution, fem::mesh::triangular_2d::Mesh mesh,
                 Function_2D analytical_solution);

} // namespace fem::io
