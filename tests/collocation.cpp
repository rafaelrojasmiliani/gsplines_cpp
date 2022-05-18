

#include "gsplines/Collocation/GaussLobattoPointsWeights.hpp"
#include <eigen3/Eigen/Core>
#include <fenv.h>
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
using namespace gsplines::collocation;

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
/*
TEST(Collocation, Derivative_Operator) {

  GLLSpline q1 = GaussLobattoLagrangeSpline::approximate(
      optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)), nglp,
      n_inter);

  Derivative dmat(q1);

  GLLSpline q_d_test = q1.derivate();
  GLLSpline q_d_nom = dmat * q1;

  EXPECT_TRUE(tools::approx_equal(q_d_nom, q_d_test, 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(2 * (q1.derivate() + q1.derivate()),
                                  2 * (dmat + dmat) * q1, 1.0e-9));
}

TEST(GLLSpline, Transpose_Left_Multiplication) {

  GLLSpline q1 = GaussLobattoLagrangeSpline::approximate(
      optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)), nglp,
      n_inter);

  GLLSpline q2 = GaussLobattoLagrangeSpline::approximate(
      optimization::minimum_jerk_path(Eigen::MatrixXd::Random(wpn, dim)), nglp,
      n_inter);

  GLLSpline q_nom =
      GaussLobattoLagrangeSpline::approximate(q1.dot(q2), nglp, n_inter);

  TransposeLeftMultiplication q1_t(q1);

  GLLSpline q_test = q1_t * q2;

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

// Test the evaluation of the continuity functions
// Confirm that the ouput of the interpolaiton is contiuous up the the first
// deriviate
TEST(Collocation, ContinuityError) {

  gsplines::basis::BasisLagrangeGaussLobatto bgl(6);

  Eigen::MatrixXd waypoints(Eigen::MatrixXd::Random(n_inter + 1, dim));

  Eigen::VectorXd tauv((Eigen::VectorXd::Random(n_inter).array() + 3.0) / 2.0);

  gsplines::GSpline curve_1 = gsplines::interpolate(tauv, waypoints, bgl);

  ContinuityError cerr(curve_1, 3);

  Eigen::VectorXd res = cerr(curve_1);

  EXPECT_TRUE(tools::approx_zero(cerr(curve_1), 1.0e-7))
      << cerr(curve_1).array().abs().maxCoeff() << std::endl;
}

TEST(Collocation, IfOpt) {

  n_inter = 2;
  dim = 1;
  nglp = 4;
  GLLSpline _in = GaussLobattoLagrangeSpline::approximate(
      optimization::minimum_jerk_path(
          Eigen::MatrixXd::Random(n_inter + 1, dim)),
      nglp, n_inter);

  Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
      n_inter + 1, _in.get_domain().first, _in.get_domain().second);

  Eigen::VectorXd tauv =
      Eigen::VectorXd::Ones(n_inter) * _in.get_domain_length() / n_inter;

  GLLSpline first_guess =
      interpolate(tauv, _in(time_spam), basis::BasisLagrangeGaussLobatto(nglp));

  std::shared_ptr<GLLSplineVariable> variable =
      std::make_shared<GLLSplineVariable>(first_guess);
  // 1.2 Constraints
  std::shared_ptr<ConstraintWrapper<ContinuityError>> continuity =
      std::make_shared<ConstraintWrapper<ContinuityError>>(0.0, 0.0,
                                                           first_guess, 1);

  std::shared_ptr<CostWrapper<SobolevDistance>> cost =
      std::make_shared<CostWrapper<SobolevDistance>>(first_guess, nglp, n_inter,
                                                     1);

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
  ipopt.SetOption("print_level", 0);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
}

GLLSpline vector = GaussLobattoLagrangeSpline::approximate(
    optimization::minimum_jerk_path(Eigen::MatrixXd::Random(n_inter + 1, dim)),
    nglp, n_inter);

GLLSpline scalar = GaussLobattoLagrangeSpline::approximate(
    optimization::minimum_jerk_path(Eigen::MatrixXd::Random(n_inter + 1, 1)),
    nglp, n_inter);

Eigen::VectorXd time_spam =
    legendre_gauss_lobatto_points({0, 1}, nglp, n_inter);

TEST(Collocation, MultiplicationConstConst) {

  GLLSpline res = vector * scalar;
  GLLSpline res1 = scalar * vector;
  EXPECT_TRUE(tools::approx_equal(res, res1, 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(res(time_spam), res1(time_spam), 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(vector(time_spam).array().colwise() *
                                      scalar(time_spam).col(0).array(),
                                  res1(time_spam), 1.0e-4))
      << "Fail value\n"
      << ((vector(time_spam).array().colwise() *
           scalar(time_spam).col(0).array())
              .matrix() -
          res1(time_spam))
             .array()
             .abs()
             .maxCoeff()
      << "  -- \n";
}
TEST(Collocation, MultiplicationConstNonConst) {

  GLLSpline res = vector * (scalar * scalar);
  GLLSpline res1 = (vector * scalar) * scalar;
  EXPECT_TRUE(tools::approx_equal(res, res1, 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(res(time_spam), res1(time_spam), 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(
      vector(time_spam).array().colwise() *
          (scalar(time_spam).col(0).array() * scalar(time_spam).col(0).array()),
      res1(time_spam), 1.0e-7))
      << "Fail Value";
}
TEST(Collocation, MultiplicationNonConstNonConst) {

  GLLSpline res1 = ((vector + vector) * scalar) * scalar;
  GLLSpline res2 = (vector * scalar + vector * scalar) * scalar;
  GLLSpline res3 = vector * scalar * scalar + vector * scalar * scalar;
  GLLSpline res4 = (vector + vector) * (scalar * scalar);
  EXPECT_TRUE(tools::approx_equal(res1, res2, 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(res2, res3, 1.0e-9));
  EXPECT_TRUE(tools::approx_equal(res3, res4, 1.0e-9));
}
*/
TEST(Collocation, Approximation) {

  std::size_t codom_dim = 8;
  std::size_t intervals = 20;
  std::size_t nc = 20;

  Eigen::VectorXd tau = Eigen::VectorXd::Ones(intervals);

  Eigen::MatrixXd wp = Eigen::MatrixXd::Random(intervals + 1, codom_dim);

  GaussLobattoLagrangeSpline path =
      interpolate(tau, wp, basis::BasisLagrangeGaussLobatto(nc));
  GaussLobattoLagrangeSpline diffeo =
      GaussLobattoLagrangeSpline::identity({0, intervals}, nc, intervals);
  gsplines::functions::FunctionExpression trj = path.compose(diffeo);
  GLLSpline res = GLLSpline::approximate(trj, nc, intervals);
  GLLSpline res2({0, intervals}, codom_dim, intervals, nc);

  res2 = GLLSpline::approximate(trj, nc, intervals);

  Eigen::VectorXd domain = res2.get_domain_discretization();
  EXPECT_TRUE(tools::approx_equal(res2(domain), trj(domain), 1.0e-9))
      << "error " << tools::last_error;
  res2 = trj;
  EXPECT_TRUE(tools::approx_equal(res2(domain), trj(domain), 1.0e-9));
}

TEST(Collocation, Operations) {

  Eigen::VectorXd tau = Eigen::VectorXd::Ones(n_inter);

  GLLSpline q1 = GaussLobattoLagrangeSpline::approximate(
      interpolate(tau, Eigen::MatrixXd::Random(n_inter + 1, dim),
                  basis::BasisLagrangeGaussLobatto(10)),
      nglp, n_inter);

  GLLSpline q2 = GaussLobattoLagrangeSpline::approximate(
      interpolate(tau, Eigen::MatrixXd::Random(n_inter + 1, dim),
                  basis::BasisLagrangeGaussLobatto(10)),
      nglp, n_inter);

  Derivative der(q1);
  Integral integral(q1);
  TransposeLeftMultiplication tr1_(q1 - q2);
  TransposeLeftMultiplication tr2_(q1 - der * q2);

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

  tr1_ = q1 - der * q2;
  EXPECT_TRUE(tools::approx_equal(tr1_ * q2, tr2_ * q2, 1.0e-9))
      << "error " << tools::last_error << std::endl
      << "relative" << tools::last_relative_error << std::endl;
}
TEST(Collocation, Value_At) {

  Eigen::VectorXd tau = Eigen::VectorXd::Ones(n_inter);

  GLLSpline q1 = GaussLobattoLagrangeSpline::approximate(
      interpolate(tau, Eigen::MatrixXd::Random(n_inter + 1, dim),
                  basis::BasisLagrangeGaussLobatto(10)),
      nglp, n_inter);
  Eigen::VectorXd time_spam =
      legendre_gauss_lobatto_points(q1.get_domain(), nglp, n_inter);

  Eigen::MatrixXd res = q1(time_spam);

  for (long i = 0; i < time_spam.size(); i++)
    EXPECT_TRUE(
        tools::approx_equal(q1.value_at(i), res.row(i).transpose(), 1.0e-9))
        << tools::last_error;
}

int main(int argc, char **argv) {

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
