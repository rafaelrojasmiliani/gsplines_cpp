

#include <eigen3/Eigen/Core>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

#include <random>

using namespace gsplines;
TEST(GLLSpline, Derivative_Operator) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  std::size_t dim = uint_dist(mt);
  std::size_t nglp = 2 * uint_dist(mt);
  std::size_t n_inter = uint_dist(mt);
  std::size_t wpn = uint_dist(mt);

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

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  std::size_t dim = uint_dist(mt);
  std::size_t nglp = 2 * uint_dist(mt);
  std::size_t n_inter = uint_dist(mt);
  std::size_t wpn = uint_dist(mt);

  collocation::GLLSpline q1 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::GLLSpline q2 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::GLLSpline q_nom =
      collocation::GaussLobattoLagrangeSpline::approximate(q1.dot(q2), nglp,
                                                           n_inter);

  collocation::GLLSpline::TransposeLeftMultiplication q1_t(q1);

  collocation::GLLSpline q_test = q1_t * q2;

  EXPECT_TRUE(tools::approx_equal(q_nom, q_test, 1.0e-9))
      << "\n Nom:\n"
      << q_nom.get_coefficients().transpose() << "\n Test:\n"
      << q_test.get_coefficients().transpose() << "\n"
      << "Error: "
      << (q_nom.get_coefficients() - q_test.get_coefficients())
             .array()
             .abs()
             .maxCoeff();
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
