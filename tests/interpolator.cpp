

#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Interpolator.hpp>
#include <iostream>
#include <stdio.h>

#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <random>

TEST(Interpolator, Value) {

  for (std::size_t i = 1; i < 3; i++) {
    gsplines::basis::BasisLegendre basis(4);
    std::size_t dim = 2;
    std::size_t intervals = 4;
    Eigen::VectorXd tau = Eigen::VectorXd::Random(4);
    Eigen::MatrixXd wp = Eigen::MatrixXd::Random(intervals + 1, dim);
    gsplines::Interpolator inter(dim, intervals, basis);
    Eigen::VectorXd tauv(intervals);
    inter.fill_interpolating_vector(wp);

    gsplines::GSpline res = inter.interpolate(tau, wp);

    Eigen::VectorXd bp = res.get_domain_breakpoints();

    /*
    for (std::size_t m = 0; m < bp.size(); m++) {
      bp(m);
    }*/
  }
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
