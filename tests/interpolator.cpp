

#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Interpolator.hpp>
#include <iostream>
#include <stdio.h>

#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <random>
using namespace gsplines;
TEST(Interpolator, Value) {

  for (std::size_t i = 1; i < 3; i++) {
    basis::BasisLegendre basis(4);
    std::size_t dim = 2;
    std::size_t intervals = 4;
    Eigen::VectorXd tau = Eigen::VectorXd::Random(4).array() + 1.2;
    Eigen::MatrixXd wp = Eigen::MatrixXd::Random(intervals + 1, dim);
    Interpolator inter(dim, intervals, basis);
    inter.fill_interpolating_vector(wp);
    GSpline res = inter.interpolate(tau, wp);

    Eigen::MatrixXd bp = res.get_waypoints();
    EXPECT_TRUE(tools::approx_equal(bp, wp, 1.0e-9))
        << "actual break points \n"
        << bp << "\n ---\n desired break points \n"
        << wp;
    /*
    for (std::size_t m = 0; m < bp.size(); m++) {
      bp(m);
    }*/
  }
  EXPECT_TRUE(true);
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
