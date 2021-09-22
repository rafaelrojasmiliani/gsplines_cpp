#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>

int main() {

  std::size_t number_of_intervals = 12;
  std::size_t codom_dim = 8;
  double exect_time = number_of_intervals;

  std::size_t basis_dim = 4;

  Eigen::VectorXd glp =
      gsplines::collocation::legendre_gauss_lobatto_points(basis_dim);

  gsplines::basis::BasisLagrange basis_lagrange(glp);

  gsplines::basis::BasisLegendre basis_legendre(basis_dim);

  Eigen::MatrixXd waypoints(
      Eigen::MatrixXd::Random(number_of_intervals + 1, codom_dim));

  Eigen::VectorXd tau(
      (Eigen::VectorXd::Random(number_of_intervals).array() + 1.0) / 2.0);

  Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(30, 0.0, exect_time);

  gsplines::GSpline curve_1 =
      gsplines::interpolate(tau, waypoints, basis_lagrange);
  gsplines::GSpline curve_2 =
      gsplines::interpolate(tau, waypoints, basis_legendre);

  Eigen::MatrixXd val_1 = curve_1(time_spam);
  Eigen::MatrixXd val_2 = curve_2(time_spam);

  Eigen::VectorXd res =
      basis_lagrange.continuity_matrix(number_of_intervals, codom_dim, 1, tau) *
      curve_1.get_coefficients();
  double val = Eigen::abs(res.array()).maxCoeff();
  assert(val < 1.0e-9);

  return 0;
}
