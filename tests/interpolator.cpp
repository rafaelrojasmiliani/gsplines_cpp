

#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Interpolator.hpp>
#include <iostream>
#include <stdio.h>

#include <random>

int main() {

  for (int i = 1; i < 10; i++) {
    gsplines::basis::BasisLegendre basis(2 * i);
    std::size_t dim = 2;
    std::size_t intervals = 4;
    Eigen::MatrixXd wp = Eigen::MatrixXd::Random(intervals + 1, dim);
    gsplines::Interpolator inter(dim, intervals, basis);
    Eigen::VectorXd tauv(intervals);
    inter.fill_interpolating_matrix(tauv);
    printf("--------- fillinf vector ------------\n");
    inter.fill_interpolating_vector(wp);
    printf("-----------------------------\n");
  }
}
