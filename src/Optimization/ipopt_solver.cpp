#include "gsplines/Basis/BasisLagrange.hpp"
#include "gsplines/Collocation/GaussLobattoLagrange.hpp"
#include "gsplines/Collocation/GaussLobattoPointsWeights.hpp"
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Collocation/GaussLobattoLagrangeFunctionals.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_interface.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>
#include <ifopt/ipopt_solver.h>
#include <ifopt/problem.h>
#include <iostream>
#include <memory>

namespace gsplines {
namespace optimization {

std::optional<IpoptSolverOptions> IpoptSolverOptions::instance_ = std::nullopt;

IpoptSolverOptions::IpoptSolverOptions() = default;

IpoptSolverOptions& IpoptSolverOptions::instance() {
  if (!instance_.has_value()) {
    instance_ = IpoptSolverOptions();
  }
  return instance_.value();
}

void IpoptSolverOptions::set_option(const std::string& _option_name,
                                    const std::string& _option_value) {
  auto iter = std::find_if(
      instance().string_options_.begin(), instance().string_options_.end(),
      [_option_name](const auto& in) { return in.first == _option_name; });

  if (iter == instance().string_options_.end()) {
    instance().string_options_.emplace_back(_option_name, _option_value);
  } else {
    iter->second = _option_value;
  }
}

void IpoptSolverOptions::set_option(const std::string& _option_name,
                                    int _option_value) {
  auto iter = std::find_if(
      instance().int_options_.begin(), instance().int_options_.end(),
      [_option_name](const auto& in) { return in.first == _option_name; });

  if (iter == instance().int_options_.end()) {
    instance().int_options_.emplace_back(_option_name, _option_value);
  } else {
    iter->second = _option_value;
  }
}

void IpoptSolverOptions::set_option(const std::string& _option_name,
                                    double _option_value) {
  auto iter = std::find_if(
      instance().int_options_.begin(), instance().int_options_.end(),
      [_option_name](const auto& in) { return in.first == _option_name; });

  if (iter == instance().int_options_.end()) {
    instance().int_options_.emplace_back(_option_name, _option_value);
  } else {
    iter->second = _option_value;
  }
}
void IpoptSolverOptions::set_options_on_interface(ifopt::IpoptSolver& solver) {
  for (const auto& p : instance().string_options_) {
    solver.SetOption(p.first, p.second);
  }

  for (const auto& p : instance().int_options_) {
    solver.SetOption(p.first, p.second);
  }

  for (const auto& p : instance().double_options_) {
    solver.SetOption(p.first, p.second);
  }
}

::gsplines::GSpline optimal_sobolev_norm(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints,
    const gsplines::basis::Basis& _basis,
    const std::vector<std::pair<std::size_t, double>>& _weights,
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
  IpoptSolverOptions::set_options_on_interface(ipopt);

  // 4. Ask the solver to solve the problem
  ipopt.Solve(nlp);
  // 5. Retrive problem solution
  Eigen::VectorXd tauv = nlp.GetOptVariables()->GetValues();
  // std::cout << tauv << '\n';
  // 6. Return solution
  return inter.interpolate(tauv, _waypoints);
}

::gsplines::GSpline broken_lines_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints) {
  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(2),
                              {{1, 1.0}},
                              static_cast<double>(_waypoints.rows() - 1))
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline minimum_acceleration_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints) {
  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(4),
                              {{2, 1.0}},
                              static_cast<double>(_waypoints.rows() - 1))
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline minimum_jerk_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints) {
  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(6),
                              {{3, 1.0}},
                              static_cast<double>(_waypoints.rows() - 1))
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline minimum_snap_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints) {
  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(8),
                              {{4, 1.0}},
                              static_cast<double>(_waypoints.rows() - 1))
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline minimum_crackle_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints) {
  return optimal_sobolev_norm(_waypoints, gsplines::basis::BasisLegendre(10),
                              {{5, 1.0}},
                              static_cast<double>(_waypoints.rows() - 1))
      .linear_scaling_new_execution_time(1.0);
}

::gsplines::GSpline rojas_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints, double _k) {
  auto ni = static_cast<double>(_waypoints.rows() - 1);

  double execution_time = 4.0 * ni / std::sqrt(2.0);
  double aux = 4 * ni * _k / std::sqrt(2);
  double alpha = 1.0 / (1.0 + std::pow(execution_time / aux, 4));

  return optimal_sobolev_norm(_waypoints, gsplines::basis::Basis0101(alpha),
                              {{1, alpha}, {3, 1 - alpha}}, execution_time)
      .linear_scaling_new_execution_time(1.0);
}

}  // namespace optimization
}  // namespace gsplines
