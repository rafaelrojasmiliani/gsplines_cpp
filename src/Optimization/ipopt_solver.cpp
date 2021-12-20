#include <iostream>

#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_interface.hpp>
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>

namespace gsplines {
namespace optimization {

::gsplines::GSpline
optimal_sobolev_norm(const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
                     const gsplines::basis::Basis &_basis,
                     std::vector<std::pair<std::size_t, double>> _weights,
                     double _exec_time) {

  ifopt::Problem nlp;
  std::size_t num_intervals = _waypoints.rows() - 1;
  std::size_t codom_dim = _waypoints.cols();
  gsplines::Interpolator inter(codom_dim, num_intervals, _basis);

  if (num_intervals == 1) {
    Eigen::VectorXd tauv_aux(1);
    tauv_aux(0) = _exec_time;
    return inter.interpolate(tauv_aux, _waypoints);
  }

  // 1. Build problem objects
  // 1.1 Variables
  std::shared_ptr<TimeSegmentLenghtsVar> variable(
      new TimeSegmentLenghtsVar(num_intervals, _exec_time));
  // 1.2 Constraints
  std::shared_ptr<ExecTimeConstraint> constraints(
      new ExecTimeConstraint(num_intervals, _exec_time));
  // 1.3 Cost Function
  std::shared_ptr<SobolevNorm> cost_function(
      new SobolevNorm("Name", _waypoints, _basis, _weights));

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
  // ipopt.SetOption("derivative_test", "first-order");
  ipopt.SetOption("hessian_approximation", "limited-memory");
  ipopt.SetOption("jac_c_constant", "yes");
  ipopt.SetOption("print_level", 0);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  // 5. Retrive problem solution
  Eigen::VectorXd tauv = nlp.GetOptVariables()->GetValues();
  // std::cout << tauv << '\n';
  // 6. Return solution
  return inter.interpolate(tauv, _waypoints);
}

::gsplines::GSpline
broken_lines_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(2),
                              {{1, 1.0}}, _waypoints.rows() - 1)
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline
minimum_acceleration_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(4),
                              {{2, 1.0}}, _waypoints.rows() - 1)
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline
minimum_jerk_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(6),
                              {{3, 1.0}}, _waypoints.rows() - 1)
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline
minimum_snap_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(8),
                              {{4, 1.0}}, _waypoints.rows() - 1)
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline
minimum_crackle_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(10),
                              {{5, 1.0}}, _waypoints.rows() - 1)
      .linear_scaling_new_execution_time(1.0);
}
} // namespace optimization
} // namespace gsplines
