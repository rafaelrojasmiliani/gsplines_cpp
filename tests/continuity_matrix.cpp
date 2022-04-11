#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <random>

TEST(ContinuityMatrix, Value) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> uint_dist(2, 4);

  for (int _ = 0; _ < 100; _++) {
    std::size_t number_of_intervals = uint_dist(mt);
    std::size_t codom_dim = uint_dist(mt);

    std::size_t basis_dim = 6;

    Eigen::VectorXd glp =
        gsplines::collocation::legendre_gauss_lobatto_points(basis_dim);

    gsplines::basis::BasisLagrange basis_lagrange(glp);

    gsplines::basis::BasisLegendre basis_legendre(basis_dim);

    Eigen::MatrixXd waypoints(
        6 * Eigen::MatrixXd::Random(number_of_intervals + 1, codom_dim));

    Eigen::VectorXd tau(
        (Eigen::VectorXd::Random(number_of_intervals).array() + 7.0) / 2.0);

    gsplines::GSpline curve_1 =
        gsplines::interpolate(tau, waypoints, basis_lagrange);

    Eigen::VectorXd res = basis_lagrange.continuity_matrix(number_of_intervals,
                                                           codom_dim, 2, tau) *
                          curve_1.get_coefficients();
    double val = Eigen::abs(res.array()).maxCoeff();
    EXPECT_TRUE(gsplines::tools::approx_equal(val, 0.0, 0.1))
        << res.transpose() << "\n"
        << "codo dim = " << codom_dim << "\n number of interval "
        << number_of_intervals << "\n"
        << basis_lagrange.continuity_matrix(number_of_intervals, codom_dim, 2,
                                            tau);
  }
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
