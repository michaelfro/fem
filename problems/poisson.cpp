#include <iostream>
#include <vector>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/io/VTKWriter.hpp>
#include <fem/linalg/CG.hpp>

int main(){
  // fem::linalg::cg::test_cg();
  // fem::linalg::cg::test_pcg(); 
  //
  auto A = fem::linalg::matrix::make_poisson_1d(10);

  auto L = fem::linalg::matrix::get_triu(A);
  fem::linalg::matrix::print(L);

  auto U = fem::linalg::matrix::get_tril(A);
  fem::linalg::matrix::print(U);

  auto D = fem::linalg::matrix::get_diag(A);
  fem::linalg::matrix::print(D);
  return 0;
}
