

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <random>
using namespace gsplines;
/** Test the following properties
 *
 * P1. That the sum of all the lagrange polynomials is 1.0
 *
 * P2. That the j-th lagrange polynomial is 1 at the j-th point and zero at the
 * other points.
 * */
TEST(BasisLagrange, Value_Properties) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  for (int i = 0; i < 100; i++) {
    std::size_t dim = uint_dist(mt);

    Eigen::VectorXd points = Eigen::VectorXd::Random(dim);
    gsplines::basis::BasisLagrange basis_lagrange(points);
    Eigen::VectorXd buff(dim);
    // Test P1
    for (int j = 0; j < 30; j++) {
      basis_lagrange.eval_on_window(real_dist(mt), 1.0, buff);
      double res = buff.sum();
      EXPECT_TRUE(tools::approx_equal(res, 1.0, 1.0e-7))
          << "value " << std::to_string(res);
    }
    // Test P2
    for (std::size_t j = 0; j < dim; j++) {
      basis_lagrange.eval_on_window(points(j), 1.0, buff);
      EXPECT_TRUE(tools::approx_equal(buff(j), 1.0, 1.0e-7))
          << "value " << std::to_string(buff(j));

      for (std::size_t k = 0; k < dim; k++)
        if (k != j)
          EXPECT_TRUE(tools::approx_equal(buff(k), 0.0, 1.0e-7))
              << "value " << std::to_string(buff(k));
    }
  }
}
/** Here we test that the interpolator gives the same curve
 * for two differente basis of the vector space of Polynomials.
 * In fact, the interpolation using the Legendre of the Lagrange basis
 * must be equal.
 * */
TEST(BasisLagrange, GSpline_Value) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> uint_dist(3, 10);

  std::size_t dim = uint_dist(mt);
  std::size_t codom_dim = uint_dist(mt);
  std::size_t n_intervals = uint_dist(mt);

  Eigen::VectorXd glp =
      gsplines::collocation::legendre_gauss_lobatto_points(dim);

  Eigen::VectorXd points = glp;
  gsplines::basis::BasisLagrange basis_lagrange(points);
  gsplines::basis::BasisLegendre basis_legendre(dim);
  Eigen::MatrixXd waypoints(
      Eigen::MatrixXd::Random(n_intervals + 1, codom_dim));
  Eigen::VectorXd tau(Eigen::VectorXd::Random(n_intervals).array() + 1.5);
  Eigen::VectorXd time_spam =
      Eigen::VectorXd::LinSpaced(250, 0.0, tau.array().sum());

  gsplines::GSpline curve_1 =
      gsplines::interpolate(tau, waypoints, basis_legendre);

  Eigen::MatrixXd val_1 = curve_1(time_spam);
  for (int i = 0; i < 100; i++) {
    gsplines::GSpline curve_2 =
        gsplines::interpolate(tau, waypoints, basis_lagrange);
    Eigen::MatrixXd val_2 = curve_2(time_spam);
    EXPECT_TRUE(gsplines::tools::approx_equal(val_1, val_2, 1.0e-6));
  }
  /*
    {
      gsplines::GSpline curve_1 = gsplines::optimization::optimal_sobolev_norm(
          waypoints, basis_lagrange, {{3, 1.0}}, 5);

      gsplines::GSpline curve_2 = gsplines::optimization::optimal_sobolev_norm(
          waypoints, basis_legendre, {{3, 1.0}}, 5);

      Eigen::MatrixXd val_1 = curve_1(time_spam);
      Eigen::MatrixXd val_2 = curve_2(time_spam);

      double err = Eigen::abs((val_1 - val_2).array()).maxCoeff();
      double max_1 = Eigen::abs((val_1).array()).maxCoeff();
      double max_2 = Eigen::abs((val_1).array()).maxCoeff();

      // assert(err / max_1 < 1.0e-2 and err / max_2 < 1.0e-2);
    }
    */
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
