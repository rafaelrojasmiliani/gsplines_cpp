

#include <eigen3/Eigen/Core>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
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
TEST(BasisLagrange_GL, Value_Properties) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> real_dist(0.0, 1.0);
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  for (int i = 0; i < 100; i++) {
    std::size_t dim = uint_dist(mt);

    Eigen::VectorXd points = collocation::legendre_gauss_lobatto_points(dim);
    gsplines::basis::BasisLagrangeGaussLobatto basis_lagrange(dim);
    Eigen::VectorXd buff(dim);
    // Test P1
    for (int j = 0; j < 30; j++) {
      basis_lagrange.eval_on_window(real_dist(mt), 1.0, buff);
      double res = buff.sum();
      EXPECT_TRUE(tools::approx_equal(res, 1.0, 5.0e-5))
          << "value " << std::to_string(res);
    }
    // Test P2
    for (std::size_t j = 0; j < dim; j++) {
      basis_lagrange.eval_on_window(points(j), 1.0, buff);
      EXPECT_TRUE(tools::approx_equal(buff(j), 1.0, 5.0e-5))
          << "value " << std::to_string(buff(j));

      for (std::size_t k = 0; k < dim; k++)
        if (k != j)
          EXPECT_TRUE(tools::approx_equal(buff(k), 0.0, 5.0e-5))
              << "value " << std::to_string(buff(k));
    }

    Eigen::MatrixXd mat = basis_lagrange.get_derivative_matrix_block();
    Eigen::VectorXd buff_2(dim);
    for (std::size_t k = 1; k < 2; k++) {

      double s = real_dist(mt);
      basis_lagrange.eval_derivative_on_window(s, 2.0, k, buff);
      basis_lagrange.eval_on_window(s, 2.0, buff_2);

      EXPECT_TRUE(tools::approx_equal(buff, mat.transpose() * buff_2, 1.0e-8))
          << buff << "\n--\n"
          << mat.transpose() * buff_2;

      mat *= mat;
    }
  }
}
TEST(GLP, Derivative_Operator) {

  std::size_t dim = 2;
  std::size_t nglp = 4;
  std::size_t n_inter = 2;
  std::size_t wpn = 5;

  Eigen::VectorXd glp =
      collocation::legendre_gauss_lobatto_points({0, 1}, nglp, n_inter);

  GSpline curve =
      optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim));
  collocation::GLLSpline q1 =
      collocation::GaussLobattoLagrangeSpline::approximate(curve, nglp,
                                                           n_inter);

  std::cout << glp.transpose() << "\n\n";
  std::cout << curve(glp) << "\n\n";
  std::cout << q1.get_coefficients().transpose() << "\n\n";
  std::cout << q1(glp) << "\n\n";
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
