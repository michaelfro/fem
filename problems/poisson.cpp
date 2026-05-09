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

struct triplet {
    size_t global_row;
    size_t global_col; 
    double contribution;
};

struct local_result {
    std::vector<triplet> local_contributions;
    std::vector<double> rhs;
};

struct lin_system {
    fem::linalg::csr_matrix stiffness_matrix;
    std::vector<double> rhs;
};

double f(const double x, const double y) {
    return x*x + y*y - 2;
}

void employ_boundary_conditions(std::vector<double>& rhs, const fem::mesh::triangular_2d::Mesh mesh){ // homogenous dirichlet boundary conditions
    
    for (size_t i = 0; i < rhs.size(); ++i) {
       double eps = 1e-12;
       if  (std::abs(mesh.nodes[i][0] - 1) <= eps || std::abs(mesh.nodes[i][0] + 1) <= eps 
            || std::abs(mesh.nodes[i][1] - 1) <= eps || std::abs(mesh.nodes[i][1] +1)<= eps)
       {
          rhs[i] = 0.;
       }
    }
    return;
}

local_result local_contribution(const fem::mesh::triangular_2d::Mesh mesh) {

    std::vector<double> rhs(mesh.elements.size(),0.);
    std::vector<triplet> local_contributions;

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

          double det_J = J[0][0] * J[1][1]- J[1][0] * J[0][1];
          double area = std::abs(det_J) * 0.5;

          // Inverted transposed Jacobian:
          double inv_J_transposed[2][2];

          inv_J_transposed[0][0] = J[0][0] / det_J;
          inv_J_transposed[0][1] = -J[0][0] / det_J;
          inv_J_transposed[1][0] = -J[1][0] / det_J;
          inv_J_transposed[1][1] = J[1][1] / det_J;
          
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
          for (int k = 0; k < 3; ++k){
              for(int k_dash = 0; k_dash < 3; ++k_dash){
                  // dot = grad_L_k dot grad_L_k_dash
                  double dot = grad_L[k][0] * grad_L[k_dash][0] + 
                               grad_L[k][1] * grad_L[k_dash][1];

                  // K_e [k][k'] = area * dot; local stiffness_matrix
                  triplet triplet;
                  triplet.contribution = area*dot;
                  triplet.global_row = ele[k];
                  triplet.global_col = ele[k_dash];
                  local_contributions.push_back(triplet);
              }
          }
          // rhs contribution
          // compute centroid of element:
          double c [2];
          c[0] = (p1[0] + p2[0] + p3[0]) / 3;
          c[1] = (p1[1] + p2[1] + p3[1]) / 3;
          double f_c = f(c[0],c[1]);
          double val = (area / 3.0) * f_c;
            
          // todo check if this is correct
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
            if (a.global_row != b.global_row)
                return a.global_row < b.global_row;

            return a.global_col < b.global_col;
        });
    
    // compress duplicate indice pairs
    std::vector<triplet> compressed; 

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
        compressed.push_back({r, c, sum});
        i = j;
    }

    std::vector<size_t> row_index;
    std::vector<size_t> col_index;
    std::vector<double> values;
    values.reserve(compressed.size());
    col_index.reserve(compressed.size());
    row_index.reserve(compressed.size()+1);

    size_t r = 0;
    size_t counter = 0;

    for (auto triplet : compressed) {
        if (triplet.global_row > r) { 
            r = triplet.global_row;
            row_index.push_back(row_index[r] + counter);
            counter = 0;
        }
        col_index.push_back(triplet.global_col);
        values.push_back(triplet.contribution);
        
    }

    fem::linalg::csr_matrix stiffness_matrix = {values, col_index, row_index};
    employ_boundary_conditions(rhs,mesh);

    return {stiffness_matrix, rhs};
}


int main(){ 

  fem::mesh::triangular_2d::Mesh mesh = fem::io::read_mesh_2d("../data/meshes/square.vtk"); // assumes you run the code from /build
  std::cout << "Number of nodes in the mesh: " << mesh.nodes.size() << std::endl;
  std::cout << "Number of elements in the mesh " << mesh.elements.size() << std::endl;

  lin_system linear_system = assemble_system(mesh);
  std::vector<double> rhs = linear_system.rhs;
  fem::linalg::csr_matrix stiffness_matrix = linear_system.stiffness_matrix;

  size_t max_iterations = 10000;
  double reduction_factor = 10e-10;

  fem::linalg::cg::solution solution = fem::linalg::cg::CG(stiffness_matrix, rhs, max_iterations, reduction_factor);
  
  fem::io::write_solution(solution, mesh);

  return 0;
}
