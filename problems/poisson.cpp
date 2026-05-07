#include <iostream>
#include <vector>
#include <fem/linalg/Vector.hpp>
#include <fem/linalg/CSRMatrix.hpp>
#include <fem/io/VTKWriter.hpp>
#include <fem/linalg/CG.hpp>

int main(){
  fem::linalg::cg::test_cg();
  fem::linalg::cg::test_pcg(); 
  return 0;
}
