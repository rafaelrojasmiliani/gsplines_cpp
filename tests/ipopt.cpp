#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

TEST(Iport, TestRun) {
  long number_of_wp = 3;
  long codom_dim = 7;

  double exec_time = (double)number_of_wp - 1.0;

  Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

  EXPECT_TRUE(gsplines::optimization::optimal_sobolev_norm(
      wp, gsplines::basis::BasisLegendre(6), {{1.0, 3}}, exec_time));

  EXPECT_TRUE(gsplines::optimization::broken_lines_path(wp));

  EXPECT_TRUE(gsplines::optimization::minimum_acceleration_path(wp));
  EXPECT_TRUE(gsplines::optimization::minimum_jerk_path(wp));
  EXPECT_TRUE(gsplines::optimization::minimum_snap_path(wp));
}
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
