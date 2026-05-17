#include <iostream>
#include <vector>
#include <fem/io/input.hpp>
#include <fem/io/output.hpp>
#include <fem/mesh/mesh.hpp>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/linalg/cg.hpp>
#include <fem/problems/poisson.hpp>
#include <cmath>
#include <array>
#include <algorithm>
#include <cassert>
#include <map>

using fem::problems::poisson::assemble_system;
using fem::problems::poisson::lin_system;
using fem::mesh::hierarchical::HierarchicalMesh;
using fem::mesh::hierarchical::LevelData;
using fem::linalg::csrmatrix;

HierarchicalMesh build_hierarchical_mesh(const std::vector<std::string>& vtk_files) { 

    HierarchicalMesh hm;
    
    // Load all meshes
    std::vector<fem::mesh::triangular_2d::Mesh> meshes;
    for (auto& f : vtk_files)
        meshes.push_back(fem::io::read_mesh_2d(f));
    
    hm.fine_mesh = meshes.back();
    int n_fine = hm.fine_mesh.nodes.size();
    
    // Build coordinate lookup for finest mesh:
    auto make_key = [](double x, double y) -> std::pair<long,long> {
        return { std::llround(x * 1e10), 
                 std::llround(y * 1e10) };
    };
    
    std::map<std::pair<long,long>, int> fine_lookup;
    for (int i = 0; i < n_fine; ++i) {
        auto [x, y] = hm.fine_mesh.nodes[i];
        fine_lookup[make_key(x, y)] = i;
    }
    
    hm.midpoints.resize(meshes.size() - 1);
    
    // loop over the meshes:
    for (int k = 0; k < meshes.size(); ++k) {
        auto& mesh_k = meshes[k];
        int n_coarse = mesh_k.nodes.size();
        
        LevelData ld;
        ld.level = k;
        ld.n_nodes = n_coarse;
        ld.coarse_to_fine.resize(n_coarse);
        ld.fine_to_coarse.assign(n_fine, -1);

        // Compute inverse diagonal of coarse stiffness matrix
        const auto coarse_matrix = assemble_system(coarse_mesh).stiffness_matrix; 
        auto diag = fem::linalg::get_diag(coarse_matrix);
        for (auto ele in diag) {
            ele = 1 / ele;
        }
        ld.inv_diag = diag;

        
        // Map coarse nodes to fine nodes 
        for (int i = 0; i < n_coarse; ++i) {
            auto [x, y] = mesh_k.nodes[i];
            auto key = make_key(x, y);
            int fi = fine_lookup.at(key); 
            ld.coarse_to_fine[i] = fi;
            ld.fine_to_coarse[fi] = i;
        }
        
        // Find midpoints introduced at this level transition
        if (k + 1 < meshes.size()) {
            auto& mesh_next = meshes[k+1];
            
            // Build edge midpoint lookup from level k
            std::map<std::pair<long,long>, std::pair<int,int>> midpoint_lookup;
            for (auto& tri : mesh_k.elements) {
                std::array<std::pair<int,int>,3> edges = {
                    std::make_pair(tri[0], tri[1]),
                    std::make_pair(tri[1], tri[2]),
                    std::make_pair(tri[0], tri[2])
                };
                for (auto [a, b] : edges) {
                    double mx = (mesh_k.nodes[a][0] + mesh_k.nodes[b][0]) / 2.0;
                    double my = (mesh_k.nodes[a][1] + mesh_k.nodes[b][1]) / 2.0;
                    midpoint_lookup[make_key(mx, my)] = {a, b};
                }
            }
            
            // Check which next-level nodes are midpoints
            for (int i = 0; i < mesh_next.nodes.size(); ++i) {
                auto [x, y] = mesh_next.nodes[i];
                auto key = make_key(x, y);
                if (midpoint_lookup.count(key)) {
                    auto [p, q] = midpoint_lookup[key];
                    int fi = fine_lookup.at(key);
                    hm.midpoints[k].push_back({fi, p, q});
                }
            }
        }
        
        hm.levels.push_back(std::move(ld)); // some nice optimization
    }
    
    return hm;
}

 auto create_bpx(HierarchicalMesh& hm) {

    std::vector<double> bpx(hm.fine_mesh.nodes.size(), 0.0);

    for(size_t level = 0; level < hm.levels.size(); ++level) {
        
        auto& ld = hm.levels[level]; // LevelData

        // Coarse contribution
        for (size_t i = 0; i < ld.n_nodes; ++i) {
            bpx[ld.coarses_to_fine[i]] += ld.inv_diag[i];
        }

        // midpoint contribution
        if(k + 1 < hm.levels.size()) { // if we have midpoints --  is this correct?
            for (auto& mid: hm.midpoints[level]) {
                bpx[mid.fine_index] += 0.5 * (ld.inv_diag[mid.coarse_parent_1] + ld.inv_diag[mid.coarse_parent_2]);
            } 
        }

    } 

    // now create a function that returns componentwise multiplication with this vector
    return [bpx](std::vector<double>& r) {
        assert (bpx.size() == r.size());
        for (size_t i = 0; i < bpx.size(); ++i) {
            r[i] = bpx[i] * r[i];  
        }
    } 
}

int main(){ 

  // fem::mesh::triangular_2d::Mesh mesh = fem::io::read_mesh_2d("../data/meshes/quader5.vtk"); // assumes you run the code from /build
  // std::cout << "Number of nodes in the mesh: " << mesh.nodes.size() << "\n";
  // std::cout << "Number of elements in the mesh " << mesh.elements.size() << "\n";

  lin_system linear_system = assemble_system(mesh);
  std::vector<double> rhs = linear_system.rhs;
  csr_matrix stiffness_matrix = linear_system.stiffness_matrix;

  auto bpx_preconditioner = create_bpx(stiffness_matrix, "../data/meshes/quader1.vtk", 4);

  size_t max_iterations = 10000;
  double reduction_factor = 10e-10;

  // std::cout << "Number of entries in stiffness matrix " << stiffness_matrix.values.size() << "\n" 
  //           << "Number of entries in row_index: " << stiffness_matrix.row_index.size() << "\n" 
  //           << "Number of entries in col_index " << stiffness_matrix.col_index.size() << "\n"
  //           << "Size of rhs: " << rhs.size() << "\n";

  // initial guess:
  std::vector<double> x_0(0.,rhs.size());

  fem::linalg::cg::solution solution = fem::linalg::cg::PCG(stiffness_matrix, rhs, x_0, bpx_preconditioner, max_iterations, reduction_factor);
  std::cout << "PCG finished in " << solution.iterations << " iterations with a final residuum of " << solution.final_residuum << "\n";

  fem::linalg::cg::solution solution = fem::linalg::cg::CG(stiffness_matrix, rhs, max_iterations, reduction_factor);
  std::cout << "Normal CG finished in " << solution.iterations << " iterations with a final residuum of " << solution.final_residuum << "\n";

  std::cout << "Error in nodal l2 norm: " << nodal_l2_error(mesh, solution.solution) << "\n"; 

  fem::io::write_solution_csv_2d(solution, mesh,"poisson_bpx.csv");

  return 0;
}
