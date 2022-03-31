#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <iostream>

TEST(Iport, TestRun) {
  std::size_t number_of_wp = 3;
  std::size_t codom_dim = 7;

  double exec_time = (double)number_of_wp - 1.0;

  Eigen::MatrixXd wp(Eigen::MatrixXd::Random(number_of_wp, codom_dim));

  gsplines::GSpline trj = gsplines::optimization::optimal_sobolev_norm(
      wp, gsplines::basis::BasisLegendre(6), {{1.0, 3}}, exec_time);

  gsplines::GSpline trj_2 = gsplines::optimization::broken_lines_path(wp);
  gsplines::GSpline trj_3 =
      gsplines::optimization::minimum_acceleration_path(wp);
  gsplines::GSpline trj_4 = gsplines::optimization::minimum_jerk_path(wp);
  gsplines::GSpline trj_5 = gsplines::optimization::minimum_snap_path(wp);
  gsplines::GSpline trj_6 = gsplines::optimization::minimum_crackle_path(wp);
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
