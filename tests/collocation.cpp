

#include "gsplines/Collocation/GaussLobattoPointsWeights.hpp"
#include <eigen3/Eigen/Core>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoLagrangeFunctionals.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <iostream>
#include <random>

using namespace gsplines;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<double> real_dist(0.0, 1.0);
std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

std::size_t dim = uint_dist(mt);
std::size_t nglp = 2 * uint_dist(mt) + 2;
std::size_t n_inter = uint_dist(mt);
std::size_t wpn = uint_dist(mt);
/** Test the following properties
 *
 * P1. That the sum of all the lagrange polynomials is 1.0
 *
 * P2. That the j-th lagrange polynomial is 1 at the j-th point and zero at the
 * other points.
 * */
TEST(Collocation, Derivative_Operator) {

  collocation::GLLSpline q1 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::Derivative dmat(q1);

  collocation::GLLSpline q_d_test = q1.derivate();
  collocation::GLLSpline q_d_nom = dmat * q1;

  EXPECT_TRUE(tools::approx_equal(q_d_nom, q_d_test, 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(2 * (q1.derivate() + q1.derivate()),
                                  2 * (dmat + dmat) * q1, 1.0e-9));
}

TEST(GLLSpline, Transpose_Left_Multiplication) {

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

  collocation::TransposeLeftMultiplication q1_t(q1);

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

  EXPECT_TRUE(tools::approx_equal(4 * q_nom, 2 * (q1_t + q1_t) * q2, 1.0e-9));
}

/* Test the evaluation of the continuity functions
 * Confirm that the ouput of the interpolaiton is contiuous up the the first
 * deriviate  */
TEST(Collocation, ContinuityError) {

  gsplines::basis::BasisLagrangeGaussLobatto bgl(6);

  Eigen::MatrixXd waypoints(Eigen::MatrixXd::Random(n_inter + 1, dim));

  Eigen::VectorXd tauv((Eigen::VectorXd::Random(n_inter).array() + 3.0) / 2.0);

  gsplines::GSpline curve_1 = gsplines::interpolate(tauv, waypoints, bgl);

  collocation::ContinuityError cerr(curve_1, 3);

  Eigen::VectorXd res = cerr(curve_1);

  EXPECT_TRUE(tools::approx_zero(cerr(curve_1), 1.0e-7))
      << cerr(curve_1).array().abs().maxCoeff() << std::endl;
}
TEST(Collocation, Operations) {

  collocation::GLLSpline q1 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::GLLSpline q2 =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)),
          nglp, n_inter);

  collocation::Derivative der(q1);
  collocation::Integral integral(q1);
  collocation::TransposeLeftMultiplication tr1_(q1 - q2);
  collocation::TransposeLeftMultiplication tr2_(q1 - der * q2);

  EXPECT_TRUE(
      tools::approx_equal(tr1_ * (q1 - q2), tr1_ * q1 - tr1_ * q2, 1.0e-9));

  EXPECT_TRUE(
      tools::approx_equal(tr2_ * (-der * q2), -tr2_ * der * q2, 1.0e-9));

  EXPECT_TRUE(
      tools::approx_equal(tr2_ * (q1 - q2), tr2_ * q1 - tr2_ * q2, 1.0e-9));

  EXPECT_TRUE(tools::approx_equal(
      (der * (q1 - der * q2)).get_coefficients(),
      der.to_matrix() * (q1 - der * q2).get_coefficients(), 1.0e-9));

  EXPECT_TRUE(tools::approx_equal(
      (der * (q1 - der * q2)).get_coefficients(),
      der.to_matrix() *
          (q1.get_coefficients() - der.to_matrix() * q2.get_coefficients()),
      1.0e-9));

  GSpline b = der * q1 - der * der * q2;
  GSpline a = der * (q1 - der * q2);
  EXPECT_TRUE(tools::approx_equal(a, b, 1.0e-9));

  EXPECT_TRUE(tools::approx_equal(tr2_ * (q1 - der * q2),
                                  tr2_ * q1 - tr2_ * der * q2, 1.0e-9));
}

TEST(Collocation, IfOpt) {

  n_inter = 2;
  dim = 1;
  nglp = 4;
  collocation::GLLSpline _in =
      collocation::GaussLobattoLagrangeSpline::approximate(
          optimization::minimum_jerk_path(
              Eigen::MatrixXd::Random(n_inter + 1, dim)),
          nglp, n_inter);

  Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
      n_inter + 1, _in.get_domain().first, _in.get_domain().second);

  Eigen::VectorXd tauv =
      Eigen::VectorXd::Ones(n_inter) * _in.get_domain_length() / n_inter;

  std::cout << "------------------------\n"
            << tauv.transpose() << "\n"
            << time_spam << "\n------------------------ " << n_inter << "\n"
            << _in.get_waypoints();

  collocation::GLLSpline first_guess =
      interpolate(tauv, _in(time_spam), basis::BasisLagrangeGaussLobatto(nglp));

  std::shared_ptr<collocation::GLLSplineVariable> variable =
      std::make_shared<collocation::GLLSplineVariable>(first_guess);
  // 1.2 Constraints
  std::shared_ptr<collocation::ConstraintWrapper<collocation::ContinuityError>>
      continuity = std::make_shared<
          collocation::ConstraintWrapper<collocation::ContinuityError>>(
          0.0, 0.0, first_guess, 1);

  std::shared_ptr<collocation::CostWrapper<collocation::SobolevDistance>> cost =
      std::make_shared<collocation::CostWrapper<collocation::SobolevDistance>>(
          first_guess, nglp, n_inter, 1);

  ifopt::Problem nlp;

  nlp.AddVariableSet(variable);
  nlp.AddConstraintSet(continuity);
  nlp.AddCostSet(cost);
  // nlp.PrintCurrent();

  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  ipopt.SetOption("derivative_test", "first-order");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("fast_step_computation", "yes");
  ipopt.SetOption("hessian_approximation", "limited-memory");
  ipopt.SetOption("jac_c_constant", "yes");
  ipopt.SetOption("print_level", 5);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
}

collocation::GLLSpline vector =
    collocation::GaussLobattoLagrangeSpline::approximate(
        optimization::minimum_jerk_path(
            Eigen::MatrixXd::Random(n_inter + 1, dim)),
        nglp, n_inter);

collocation::GLLSpline scalar =
    collocation::GaussLobattoLagrangeSpline::approximate(
        optimization::minimum_jerk_path(
            Eigen::MatrixXd::Random(n_inter + 1, 1)),
        nglp, n_inter);

TEST(Collocation, MultiplicationConstConst) {

  Eigen::VectorXd time_spam =
      collocation::legendre_gauss_lobatto_points({0, 1}, nglp, n_inter);
  {
    collocation::GLLSpline res = vector * scalar;
    collocation::GLLSpline res1 = scalar * vector;
    EXPECT_TRUE(tools::approx_equal(res, res1, 1.0e-9));
    EXPECT_TRUE(tools::approx_equal(res(time_spam), res1(time_spam), 1.0e-9));
    EXPECT_TRUE(tools::approx_equal(vector(time_spam).array().colwise() *
                                        scalar(time_spam).col(0).array(),
                                    res1(time_spam), 1.0e-7));
    {
      collocation::GLLSpline res = vector * (scalar * scalar);
      collocation::GLLSpline res1 = (vector * scalar) * scalar;
      EXPECT_TRUE(tools::approx_equal(res, res1, 1.0e-9));
      EXPECT_TRUE(tools::approx_equal(res(time_spam), res1(time_spam), 1.0e-9));
      EXPECT_TRUE(tools::approx_equal(vector(time_spam).array().colwise() *
                                          (scalar(time_spam).col(0).array() *
                                           scalar(time_spam).col(0).array()),
                                      res1(time_spam), 1.0e-7));
    }
  }
}
TEST(Collocation, MultiplicationConstNonConst) {

  Eigen::VectorXd time_spam =
      collocation::legendre_gauss_lobatto_points({0, 1}, nglp, n_inter);

  collocation::GLLSpline res1 = ((vector + vector) * scalar) * scalar;
  collocation::GLLSpline res2 = ((vector + vector) * scalar) * scalar;
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
