

#include <eigen3/Eigen/Core>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

using namespace gsplines;
TEST(GLLSpline, Derivative_Operator) {

  std::size_t dim = 8;
  std::size_t nglp = 10;
  std::size_t n_inter = 10;
  std::size_t wpn = 5;
  collocation::GLLSpline q1 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::GLLSpline::Derivative dmat(q1);

  collocation::GLLSpline q_d_test = q1.derivate();
  collocation::GLLSpline q_d_nom = dmat * q1;

  EXPECT_TRUE(tools::approx_equal(q_d_nom, q_d_test, 1.0e-9));
}

TEST(GLLSpline, Transpose_Left_Multiplication) {

  std::size_t dim = 8;
  std::size_t nglp = 10;
  std::size_t n_inter = 10;
  std::size_t wpn = 5;
  collocation::GLLSpline q1 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::GLLSpline q2 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::GLLSpline q_nom =
      collocation::GaussLobattoLagrangeSpline::approximate(
          q1.dot(q2.derivate()), nglp, n_inter);

  collocation::GLLSpline::TransposeLeftMultiplication q1_t(q1);

  collocation::GLLSpline::Derivative dmat(q1);

  collocation::GLLSpline q_test = q1_t * dmat * q1;

  EXPECT_TRUE(tools::approx_equal(q_nom, q_test, 1.0e-9));
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
