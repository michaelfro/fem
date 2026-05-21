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
using fem::linalg::csr_matrix;

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
        const auto coarse_matrix = assemble_system(mesh_k).stiffness_matrix; 
        auto diag = fem::linalg::get_diag(coarse_matrix).values;
        for (auto& ele : diag) {
            ele = 1 / ele;
        }
        ld.inv_diag = std::move(diag);

        
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

            // Build edge to parents lookup
            std::map<std::pair<int,int>, std::pair<int,int>> edge_to_parents;

            for (auto& tri : mesh_k.elements) {
                std::array<std::pair<int,int>,3> edges = {
                    std::make_pair(tri[0], tri[1]),
                    std::make_pair(tri[1], tri[2]),
                    std::make_pair(tri[0], tri[2])
                };
                for (auto [a, b] : edges) {
                    auto key = std::make_pair(std::min(a,b), std::max(a,b));
                    edge_to_parents[key] = {a, b};
                }
            }

            // For each new node at level k+1, find its two coarse parents by checking
            // all edges — the midpoint lies between the two coarse nodes whose
            // fine indices bracket it

            for (int i = 0; i < (int)mesh_next.nodes.size(); ++i) { // for all nodes in the next finer mesh
                auto [x, y] = mesh_next.nodes[i];
                int fi = fine_lookup.at(make_key(x, y)); //find fine_index for that node

                if (ld.fine_to_coarse[fi] == -1) { // if that node is not in the coarse mesh

                    // Find the edge whose midpoint this is by checking all coarse edges
                    bool found = false;

                    for (auto& [edge_key, parents] : edge_to_parents) {
                        auto [a, b] = parents;
                        int fi_a = ld.coarse_to_fine[a];
                        int fi_b = ld.coarse_to_fine[b];
                        double mx = (hm.fine_mesh.nodes[fi_a][0] + hm.fine_mesh.nodes[fi_b][0]) / 2.0;
                        double my = (hm.fine_mesh.nodes[fi_a][1] + hm.fine_mesh.nodes[fi_b][1]) / 2.0;
                        if (make_key(mx, my) == make_key(x, y)) {
                            hm.midpoints[k].push_back({fi, a, b});
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                        throw std::runtime_error("Could not find parents for new node at level " 
                            + std::to_string(k+1));
                }
            }

            // assert(ld.n_nodes + hm.midpoints[k].size() == mesh_next.nodes.size());
        
            if (hm.midpoints[k].size() + n_coarse < mesh_next.nodes.size()) {
                std::cout << "Problem with midpoint size.\n" << "We have " << hm.midpoints[k].size() 
                            << " midpoints and " << n_coarse << " coarse nodes but the next level has " 
                            << meshes[k+1].nodes.size() << " nodes.\n";

            }
            

        }
        
        hm.levels.push_back(std::move(ld)); 
    }
    
    return hm;
}

 auto create_bpx(HierarchicalMesh& hm) {

    std::vector<double> bpx(hm.fine_mesh.nodes.size(), 0.0);

    for (size_t level = 0; level < hm.levels.size(); ++level) {
        
        auto& ld = hm.levels[level]; // LevelData

        // Coarse contribution
        for (size_t i = 0; i < ld.n_nodes; ++i) {
            bpx[ld.coarse_to_fine[i]] += ld.inv_diag[i];
        }

        // Midpoint contribution in the next level - todo fix this
        if (level + 1 < hm.levels.size()) { // if we have midpoints - finest level has no children 
            for (auto& mid: hm.midpoints[level]) {
                // bpx[mid.fine_index] += 0.5 * (ld.inv_diag[mid.coarse_parent_1] + ld.inv_diag[mid.coarse_parent_2]); // parents are at level k, midpoint at k + 1
                bpx[mid.fine_index] += 0.5 * (bpx[ld.coarse_to_fine[mid.coarse_parent_1]] + bpx[ld.coarse_to_fine[mid.coarse_parent_2]}); // parents are at level k, midpoint at k + 1                                                                                                                            
            } 
        }
    } 

    // debug check if all bpx values are positive - spd precond
    for (auto val : bpx) {
        if (val < 0) {
            std::cout << "Oh no, bad values in bpx\n";
        }
    }

    // return a function that returns componentwise multiplication with this vector - Max uses the restriction transposed to prolongate
    return [bpx](const std::vector<double>& r) -> std::vector<double> {
        assert (bpx.size() == r.size());
        std::vector<double> result(r.size());
        for (size_t i = 0; i < bpx.size(); ++i) {
            result[i] = bpx[i] * r[i];  
        }
        return result;
    };
}

int main(){ 

  std::cout << "Assembling hierachical mesh structure \n";
  std::vector<std::string> meshfiles; 
  const int number_meshes = 5;
  for (int i = 1; i <= number_meshes; ++i) {
      meshfiles.push_back("../data/meshes/quader" + std::to_string(i) + ".vtk");
  }
  HierarchicalMesh hm = build_hierarchical_mesh(meshfiles);

  std::cout << "Assembling linear system \n";
  lin_system linear_system = assemble_system(hm.fine_mesh);
  std::vector<double> rhs = linear_system.rhs;
  csr_matrix stiffness_matrix = linear_system.stiffness_matrix;

  std::cout << "Creating bpx_preconditioner \n";
  auto bpx_preconditioner = create_bpx(hm);

  size_t max_iterations = 10000;
  double reduction_factor = 10e-10;
  std::vector<double> x_0(rhs.size(),0.);

  std::cout << "Solving with PCG \n";
  fem::linalg::cg::solution solution_bpx = fem::linalg::cg::PCG(stiffness_matrix, rhs, x_0, bpx_preconditioner, max_iterations, reduction_factor);
  std::cout << "PCG finished in " << solution_bpx.iterations << " iterations with a final residuum of " << solution_bpx.final_residuum << "\n";

  std::cout << "Solving without BPX with CG \n";
  fem::linalg::cg::solution solution = fem::linalg::cg::CG(stiffness_matrix, rhs, max_iterations, reduction_factor);
  std::cout << "Normal CG finished in " << solution.iterations << " iterations with a final residuum of " << solution.final_residuum << "\n";

  std::cout << "Error in nodal l2 norm: " << fem::problems::poisson::nodal_l2_error(hm.fine_mesh, solution_bpx.solution) << "\n"; 

  fem::io::write_solution_csv_2d(solution_bpx, hm.fine_mesh,"poisson_bpx.csv");

  return 0;
}
