

#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/interpolator.hpp>
#include <iostream>
#include <stdio.h>
int main() {

  gsplines::basis::BasisLegendre basis(6);
  std::size_t dim = 2;
  std::size_t intervals = 1;
  gsplines::Interpolator inter(dim, intervals, basis);
  Eigen::VectorXd tauv(1);
  Eigen::MatrixXd wp(2, dim);
  wp = Eigen::MatrixXd::Random(2, dim);
  std::cout << wp << '\n';
  fflush(stdout);
  tauv(0) = 1.0;
  inter.interpolate(tauv, wp);
}
