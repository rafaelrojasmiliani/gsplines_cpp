

#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Interpolator.hpp>
#include <iostream>
#include <stdio.h>

#include <random>

int main() {

  for (std::size_t i = 1; i < 3; i++) {
    gsplines::basis::BasisLegendre basis(4);
    std::size_t dim = 2;
    std::size_t intervals = 4;
    Eigen::MatrixXd wp = Eigen::MatrixXd::Random(intervals + 1, dim);
    gsplines::Interpolator inter(dim, intervals, basis);
    Eigen::VectorXd tauv(intervals);
    inter.fill_interpolating_matrix(tauv);
    printf("--------- fillinf vector ------------ %zu %zu\n", 2u * i,
           basis.get_dim());
    inter.fill_interpolating_vector(wp);
    printf("-----------------------------\n");
  }
}
