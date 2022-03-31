

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
TEST(BasisLagrange, Value) {

  Eigen::VectorXd glp = gsplines::collocation::legendre_gauss_lobatto_points(4);
  gsplines::basis::BasisLagrange basis_lagrange(glp);
  gsplines::basis::BasisLegendre basis_legendre(4);

  Eigen::MatrixXd waypoints(Eigen::MatrixXd::Random(5, 8));
  Eigen::VectorXd tau(Eigen::VectorXd::Random(4));
  Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(30, 0.0, 5.0);

  gsplines::GSpline curve_1 =
      gsplines::interpolate(tau, waypoints, basis_lagrange);

  Eigen::MatrixXd val_1 = curve_1(time_spam);
  for (int i = 0; i < 10; i++) {
    gsplines::GSpline curve_2 =
        gsplines::interpolate(tau, waypoints, basis_legendre);
    Eigen::MatrixXd val_2 = curve_2(time_spam);
    gsplines::tools::approx_equal(val_1, val_2, 1.0e-9);
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
