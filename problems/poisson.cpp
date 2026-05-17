#include <iostream>
#include <vector>
#include <fem/io/input.hpp>
#include <fem/io/output.hpp>
#include <fem/mesh/mesh.hpp>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/cg.hpp>
#include <cmath>
#include <array>
#include <algorithm>
#include <cassert>

double f(const double x, const double y) {
    return x*x + y*y - 2;
}

void impose_boundary_conditions(std::vector<double>& rhs, fem::linalg::csr_matrix& stiffness_matrix, const fem::mesh::triangular_2d::Mesh mesh){ // homogenous dirichlet boundary conditions
    
    // it would be better to get rid of the dirichlet indices for performance..

    for (size_t i = 0; i < rhs.size(); ++i) {
       double eps = 1e-12;
       if  (std::abs(mesh.nodes[i][0] - 1) <= eps || std::abs(mesh.nodes[i][0] + 1) <= eps 
            || std::abs(mesh.nodes[i][1] - 1) <= eps || std::abs(mesh.nodes[i][1] +1) <= eps)
       {
          rhs[i] = 0.;

          // Change column entries to preserve symmetry, very expensive
          // for (size_t idx = 0; idx < stiffness_matrix.values.size(); ++ idx) {
          //     if (stiffness_matrix.col_index[k] == i) {
          //         stiffness_matrix.values[k] = 0.;
          //     }
          // }
          
          
          // for a boundary node change rhs and row in the stiffness_matrix - this breaks symmetry
          for (size_t k = stiffness_matrix.row_index[i]; k < stiffness_matrix.row_index[i+1]; ++k) {
              stiffness_matrix.values[k] = (stiffness_matrix.col_index[k] == i) ? 1. : 0.;
          }
       }
    }
    return;
}

local_result local_contribution(const fem::mesh::triangular_2d::Mesh mesh) {

    std::vector<double> rhs(mesh.nodes.size(),0.);
    std::vector<triplet> local_contributions;
    local_contributions.reserve(9 * mesh.elements.size());

    for (auto ele : mesh.elements) {
          // three corner nodes of the current element:
          std::array <double,2> p1 = mesh.nodes[ele[0]];
          std::array <double,2> p2 = mesh.nodes[ele[1]];
          std::array <double,2> p3 = mesh.nodes[ele[2]];
          
          double x1 = p1[0], y1 = p1[1];
          double x2 = p2[0], y2 = p2[1];
          double x3 = p3[0], y3 = p3[1];

          // Build Jacobian:
          double J[2][2] = {
            {x2 - x1, x3 - x1},
            {y2 - y1, y3 - y1}
          };

          const double det_J = J[0][0] * J[1][1] - J[1][0] * J[0][1];

          assert (std::abs(det_J) > 1e-12);

          const double area = std::abs(det_J) * 0.5;

          // Inverted transposed Jacobian:
          double inv_J_transposed[2][2];

          inv_J_transposed[0][0] = J[1][1] / det_J;
          inv_J_transposed[0][1] = -J[1][0] / det_J;
          inv_J_transposed[1][0] = -J[0][1] / det_J;
          inv_J_transposed[1][1] = J[0][0] / det_J;
          
          // reference Gradients:
          double grad_ref[3][2] = {
              {-1, -1},
              {1, 0},
              {0, 1}
          };

          double grad_L[3][2]; // = inv_J_transposed.dot(grad_ref)

          for (int i = 0; i < 3; ++i) {
              grad_L[i][0] = inv_J_transposed[0][0] * grad_ref[i][0] + 
                           inv_J_transposed[0][1] * grad_ref[i][1];
              grad_L[i][1] = inv_J_transposed[1][0] * grad_ref[i][0] + 
                           inv_J_transposed[1][1] * grad_ref[i][1];
          }
          
          // local stiffness matrix

          // double K_e[3][3];

          for (int k = 0; k < 3; ++k){
              for(int k_dash = 0; k_dash < 3; ++k_dash){
                  // dot = grad_L_k dot grad_L_k_dash
                  double dot = grad_L[k][0] * grad_L[k_dash][0] + 
                               grad_L[k][1] * grad_L[k_dash][1];

                  // K_e[k][k_dash] = area * dot; // local stiffness_matrix
                  triplet triplet;
                  triplet.contribution = area*dot;
                  triplet.global_row = ele[k];
                  triplet.global_col = ele[k_dash];
                  local_contributions.push_back(triplet);
              }
          }

          // debug print the local stiffness_matrices
          // for (int k = 0; k < 3; ++k){
          //    for(int k_dash = 0; k_dash < 3; ++k_dash){
          //      std::cout << K_e[k][k_dash] << ", ";
          //    }
          //    std::cout << "\n";
          // }
          // std::cout << "\n";
          

          // rhs contribution
          // compute centroid of element:
          double c [2];
          c[0] = (p1[0] + p2[0] + p3[0]) / 3;
          c[1] = (p1[1] + p2[1] + p3[1]) / 3;
          double f_c = f(c[0],c[1]);
          double val = (area / 3.0) * f_c;
            
          rhs[ele[0]] += val;
          rhs[ele[1]] += val;
          rhs[ele[2]] += val;
    }
    return {local_contributions,rhs};
 }

lin_system assemble_system(const fem::mesh::triangular_2d::Mesh mesh) {
    
    local_result local_res = local_contribution(mesh);

    std::vector<triplet> local_contributions = local_res.local_contributions;
    std::vector<double> rhs = local_res.rhs;

    // std::sort magic to sort first by row then column
    std::sort(local_contributions.begin(),local_contributions.end(),[](const triplet& a, const triplet& b)
        {
            if (a.global_row == b.global_row)
                return a.global_col < b.global_col;

            return a.global_row < b.global_row;
        });

    // debug print of sorted local_contributions:
    // std::cout << "Local contribution triplets:\n";
    // for (auto i : local_contributions) {
    //     std::cout << i.global_row << " " << i.global_col << " " << i.contribution << "\n";
    // }

    // compress duplicate indice pairs of sorted local_contributions
    std::vector<triplet> compressed; 
    compressed.reserve(local_contributions.size()); // can be improved right?

    for (size_t i = 0; i < local_contributions.size(); ++i) {
        size_t r = local_contributions[i].global_row;
        size_t c = local_contributions[i].global_col;
        double sum = 0.;

        size_t j = i;
        while (j < local_contributions.size() 
                  && local_contributions[j].global_row == r 
                  && local_contributions[j].global_col == c) 
        {
            sum += local_contributions[j].contribution;
            ++j;
        }
        if(std::abs(sum) > 1e-12) // save storage in the case the sum is zero
            compressed.push_back({r, c, sum});
        i = j - 1; // avoid counting the ones we looped over in the while loop twice
    }

    //debug print of compressed local_contributions:
    // std::cout << "compressed contribution triplets:\n";
    // for (auto i : compressed) {
    //     std::cout << i.global_row << " " << i.global_col << " " << i.contribution << "\n";
    // }

    // assemble global stiffness_matrix from compressed triplets
    std::vector<size_t> row_index;
    std::vector<size_t> col_index;
    std::vector<double> values;
    values.reserve(compressed.size());
    col_index.reserve(compressed.size());
    row_index.reserve(compressed.size()+1);

    row_index.push_back(0);

    size_t r = 0;
    size_t counter = 0;

    for (auto triplet : compressed) {

        if (triplet.global_row > r) { 
            row_index.push_back(row_index[r] + counter);
            r = triplet.global_row;
            counter = 0;
        }

        col_index.push_back(triplet.global_col);
        values.push_back(triplet.contribution);
        ++counter;
    }
    row_index.push_back(row_index[r] + counter);

    fem::linalg::csr_matrix stiffness_matrix = {values, col_index, row_index};

    impose_boundary_conditions(rhs,stiffness_matrix,mesh);

    return {stiffness_matrix, rhs};
}

double analytical_solution(const double x, const double y) {
    return -0.5 * (x * x - 1) * (y * y - 1);
}

double nodal_l2_error(const fem::mesh::triangular_2d::Mesh mesh, const std::vector<double> &approx_sol) //todo move this somewhere else - to the vector stuff?
{
    double sum = 0;
    for (size_t i = 0; i < mesh.nodes.size(); ++i) {
            double x = mesh.nodes[i][0]; 
            double y = mesh.nodes[i][1];
            double analytical = analytical_solution(x, y);
            double numerical = approx_sol[i];
            double err = numerical - analytical;
            sum += err * err;
    }
    double domain_volume = 4.; // todo make this depend on the mesh?
    
    return std::pow((sum * domain_volume)/ mesh.nodes.size(), 1./2.);
}


int run_poisson_problem(){ 

  fem::mesh::triangular_2d::Mesh mesh = fem::io::read_mesh_2d("../data/meshes/square.vtk"); // assumes you run the code from /build
  std::cout << "Number of nodes in the mesh: " << mesh.nodes.size() << "\n";
  std::cout << "Number of elements in the mesh " << mesh.elements.size() << "\n";

  lin_system linear_system = assemble_system(mesh);
  std::vector<double> rhs = linear_system.rhs;
  fem::linalg::csr_matrix stiffness_matrix = linear_system.stiffness_matrix;

  size_t max_iterations = 10000;
  double reduction_factor = 10e-10;

  std::cout << "Number of entries in stiffness matrix " << stiffness_matrix.values.size() << "\n" 
            << "Number of entries in row_index: " << stiffness_matrix.row_index.size() << "\n" 
            << "Number of entries in col_index " << stiffness_matrix.col_index.size() << "\n"
            << "Size of rhs: " << rhs.size() << "\n";

  // fem::linalg::print_csr_matrix(stiffness_matrix);
  
  std::cout << "Is the stiffness_matrix symmetric? " << fem::linalg::is_symmetric(stiffness_matrix) << "\n";
  
  fem::linalg::cg::solution solution = fem::linalg::cg::CG(stiffness_matrix, rhs, max_iterations, reduction_factor);
  std::cout << "Cg finished in " << solution.iterations << " iterations with a final residuum of " << solution.final_residuum << "\n";

  std::cout << "Error in nodal l2 norm: " << nodal_l2_error(mesh, solution.solution) << "\n"; 

  fem::io::write_solution_csv_2d(solution, mesh,"poisson_dirchlet_2d.csv");

  return 0;
}
