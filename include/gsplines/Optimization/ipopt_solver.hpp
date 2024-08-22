#ifndef IPOPT_SOLVER_H
#define IPOPT_SOLVER_H

#include "gsplines/Functions/FunctionBase.hpp"
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_interface.hpp>
#include <cstddef>
#include <optional>

namespace ifopt {
class IpoptSolver;
}

namespace gsplines {

namespace optimization {

class IpoptSolverOptions {
 private:
  static std::optional<IpoptSolverOptions> instance_;
  IpoptSolverOptions();

  std::vector<std::pair<std::string, std::string>> string_options_ = {
      {"linear_solver", "mumps"},
      {"jacobian_approximation", "exact"},
      {"fast_step_computation", "yes"},
      {"derivative_test", "none"},
      {"hessian_approximation", "limited-memory"},
      {"jac_c_constant", "yes"},
      {"linear_solver", "mumps"},
      {"print_timing_statistics", "no"},
      {"dependency_detector", "mumps"},
      {"dependency_detection_with_rhs", "no"}};

  std::vector<std::pair<std::string, int>> int_options_ = {{"print_level", 0}};

  std::vector<std::pair<std::string, double>> double_options_ = {
      {"tol", 1.0e-3}};

 public:
  static IpoptSolverOptions& instance();

  static void set_option(const std::string& _option_name,
                         const std::string& _option_value);

  static void set_option(const std::string& _option_name, int _option_value);

  static void set_option(const std::string& _option_name, double _option_value);

  static void set_options_on_interface(ifopt::IpoptSolver& solver);
};

std::optional<gsplines::GSpline> optimal_sobolev_norm(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints,
    const gsplines::basis::Basis& _basis,
    const std::vector<std::pair<std::size_t, double>>& _weights,
    double _exec_time);

std::optional<gsplines::GSpline> broken_lines_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints);

std::optional<gsplines::GSpline> minimum_acceleration_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints);

std::optional<gsplines::GSpline> minimum_jerk_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints);

std::optional<gsplines::GSpline> minimum_snap_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints);

std::optional<gsplines::GSpline> minimum_crackle_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints);

std::optional<gsplines::GSpline> best_approximation(
    const functions::FunctionBase& _in, std::size_t _n_glp,
    std::size_t _n_inter);

std::optional<gsplines::GSpline> rojas_path(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints, double _k);
}  // namespace optimization
}  // namespace gsplines
#endif /* IPOPT_SOLVER_H */
