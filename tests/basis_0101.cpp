

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_interface.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
/* Test that the Legendre polynomials respect the recursive relation */
using namespace gsplines;

/* Test that the derivtaive of the Legendre polynomials respect the recursive
 * relation */
TEST(Basis0101, Derivative) {
  gsplines::basis::Basis0101 basis(0.5);
  Eigen::VectorXd buff(6);
  Eigen::VectorXd buff_d1(6);
  Eigen::VectorXd buff_d2(6);
  double t = -1;
  basis.eval_derivative_on_window(t, 2, 1, buff);
  basis.eval_derivative_on_window(t, 2, 2, buff_d1);
  basis.eval_derivative_on_window(t, 2, 3, buff_d2);
}

TEST(Basis0101, InnderProductMatrixDerivatibeWRTTau) {
  Eigen::MatrixXd matrix_ground(6, 6);
  Eigen::MatrixXd matrix_test(6, 6);
  gsplines::basis::Basis0101 basis(0.5);
  std::array<double, 9> coeff = {1.0 / 280.0, -4.0 / 105.0, 1.0 / 5.0,
                                 -4.0 / 5.0,  0.0,          4.0 / 5.0,
                                 -1.0 / 5.0,  4.0 / 105.0,  -1.0 / 280.0};

  for (const auto& d_case : std::vector<int>{1, 3}) {
    const double delta = 1.0e-5;

    const double tau0 = 2.0;

    Eigen::MatrixXd buff(6, 6);
    matrix_ground.setZero();
    matrix_test.setZero();
    double tau = tau0 - 4 * delta;
    for (int i = 0; i < 9; i++) {
      buff.setZero();
      basis.add_derivative_matrix(tau, d_case, buff);
      matrix_ground += buff * coeff[i];
      tau += delta;
    }
    matrix_ground /= delta;

    basis.add_derivative_matrix_deriv_wrt_tau(tau0, d_case, matrix_test);
    // std::cout << matrix_test << "\n\n ---- \n";
    // std::cout << matrix_ground << "\n\n ---- \n";
    // std::cout << (matrix_ground - matrix_test).array().abs() << "\n\n ----
    // \n";

    EXPECT_TRUE(tools::approx_equal(matrix_ground, matrix_test, 1.0e-5));
  }
}
TEST(Basis0101, IpoptDerivativeTest) {
  Eigen::MatrixXd waypoints = Eigen::MatrixXd::Random(3, 2);
  ifopt::Problem nlp;
  double exec_time = 10;
  std::size_t num_intervals = waypoints.rows() - 1;
  std::size_t codom_dim = waypoints.cols();
  auto basis = basis::Basis0101(0.5);
  gsplines::Interpolator inter(codom_dim, num_intervals, basis);

  std::vector<std::pair<std::size_t, double>> weights = {{1, 0.5}, {3, 0.5}};
  // 1. Build problem objects
  // 1.1 Variables
  std::shared_ptr<optimization::TimeSegmentLenghtsVar> variable(
      new optimization::TimeSegmentLenghtsVar(num_intervals, exec_time));
  // 1.2 Constraints
  std::shared_ptr<optimization::ExecTimeConstraint> constraints(
      new optimization::ExecTimeConstraint(num_intervals, exec_time));
  // 1.3 Cost Function
  std::shared_ptr<optimization::SobolevNorm> cost_function(
      new optimization::SobolevNorm("Name", waypoints, basis, weights));

  // 2. Use the problem objects to build the problem
  nlp.AddVariableSet(variable);
  nlp.AddConstraintSet(constraints);
  nlp.AddCostSet(cost_function);
  // nlp.PrintCurrent();

  // 3. Instantiate ipopt solver
  ifopt::IpoptSolver ipopt;
  // 3.1 Customize the solver
  ipopt.SetOption("linear_solver", "mumps");
  ipopt.SetOption("jacobian_approximation", "exact");
  ipopt.SetOption("fast_step_computation", "yes");
  ipopt.SetOption("derivative_test", "first-order");
  ipopt.SetOption("hessian_approximation", "limited-memory");
  ipopt.SetOption("jac_c_constant", "yes");
  ipopt.SetOption("print_level", 5);

  ipopt.Solve(nlp);
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
