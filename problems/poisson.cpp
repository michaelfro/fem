#include<iostream>
#include<vector>
#include"fem/linalg/Vector.hpp"
#include"fem/linalg/CSRMatrix.hpp"
#include"fem/io/VTKWriter.hpp"

int main(){
  std::vector<double> x = {1.,1.,1.};
  double norm = fem::linalg::vector::l2_norm(x);
  std::cout <<"Testing compilation, the l2_norm of x is " << norm << std::endl;
  return 0;
}
