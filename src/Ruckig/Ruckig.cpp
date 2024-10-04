#include <gsplines/Ruckig/Ruckig.hpp>
#include <ruckig/result.hpp>
#include <ruckig/ruckig.hpp>
#include <iostream>

namespace gsplines::ruckig {

inline auto std_vector_to_eigen(std::vector<double>& _int) {
  return Eigen::Map<Eigen::VectorXd>(_int.data(),
                                     static_cast<long>(_int.size()));
}

inline auto std_vector_to_eigen(const std::vector<double>& _int) {
  return Eigen::Map<const Eigen::VectorXd>(_int.data(),
                                           static_cast<long>(_int.size()));
}

std::optional<RuckigCurve> interpolator(
    Eigen::Ref<const Eigen::MatrixXd> _waypoints,
    const std::vector<double>& _max_abs_vel,
    const std::vector<double>& _max_abs_acc,
    const std::vector<double>& _max_abs_jerk) {
  //
  auto dof = static_cast<std::size_t>(_waypoints.cols());
  ::ruckig::InputParameter<::ruckig::DynamicDOFs> ruckig_input(dof);

  auto number_of_waypoints = static_cast<std::size_t>(_waypoints.rows());

  if (number_of_waypoints < 2) {
    return std::nullopt;
  }
  if (_waypoints.cols() != dof) {
    return std::nullopt;
  }
  if (_max_abs_vel.size() != dof) {
    return std::nullopt;
  }
  if (_max_abs_acc.size() != dof) {
    return std::nullopt;
  }
  if (_max_abs_jerk.size() != dof) {
    return std::nullopt;
  }

  ruckig_input.max_velocity = _max_abs_vel;
  ruckig_input.max_acceleration = _max_abs_acc;
  ruckig_input.max_jerk = _max_abs_jerk;
  ruckig_input.current_velocity.assign(dof, 0.0);
  ruckig_input.current_acceleration.assign(dof, 0.0);

  for (int i = 0; i < dof; ++i) {
    ruckig_input.current_position[i] = _waypoints.topRows(1)(0, i);
    ruckig_input.target_position[i] = _waypoints.bottomRows(1)(0, i);
  }

  std::vector<double> vec(dof, 0.0);
  ruckig_input.intermediate_positions.reserve(number_of_waypoints - 2);
  for (int i = 1; i < number_of_waypoints - 1; ++i) {
    std_vector_to_eigen(vec) = _waypoints.row(i).transpose();
    ruckig_input.intermediate_positions.emplace_back(vec);
  }
  ::ruckig::Ruckig<::ruckig::DynamicDOFs> otg(dof, 0.001);
  ::ruckig::Trajectory<::ruckig::DynamicDOFs> trajectory(dof);
  ::ruckig::Result error = otg.calculate(ruckig_input, trajectory);

  std::map<::ruckig::Result, std::string> error_map = {
      {::ruckig::Result::Working, "Working"},
      {::ruckig::Result::Finished, "Finished"},
      {::ruckig::Result::Error, "Error"},
      {::ruckig::Result::ErrorInvalidInput, "ErrorInvalidInput"},
      {::ruckig::Result::ErrorTrajectoryDuration, "ErrorTrajectoryDuration"},
      {::ruckig::Result::ErrorPositionalLimits, "ErrorPositionalLimits"},
      {::ruckig::Result::ErrorZeroLimits, "ErrorZeroLimits"},
      {::ruckig::Result::ErrorExecutionTimeCalculation,
       "ErrorExecutionTimeCalculation"},
      {::ruckig::Result::ErrorSynchronizationCalculation,
       "ErrorSynchronizationCalculation"}};

  if (error == ::ruckig::Result::Working ||
      error == ::ruckig::Result::Finished) {
    return RuckigCurve(trajectory, RuckigCurve::type::POS);
  }
  std::cout << "GSplines Ruckig Error: " << error_map[error] << "\n\n";
  return std::nullopt;
}

RuckigCurve::RuckigCurve(::ruckig::Trajectory<::ruckig::DynamicDOFs> _trj,
                         type _type)
    : FunctionInheritanceHelper({0.0, _trj.get_duration()},
                                _trj.degrees_of_freedom, "Ruckig"),
      ruckig_(std::move(_trj)),
      type_(_type) {}

void RuckigCurve::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {
  std::vector<double> pos(get_codom_dim(), 0.0);
  std::vector<double> vel(get_codom_dim(), 0.0);
  std::vector<double> acc(get_codom_dim(), 0.0);
  std::vector<double> jerk(get_codom_dim(), 0.0);

  std::size_t sec = 0;

  for (int i = 0; i < _domain_points.size(); ++i) {
    ruckig_.at_time(_domain_points(i), pos, vel, acc, jerk, sec);
    switch (type_) {
      case type::POS:
        _result.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            pos.data(), static_cast<long>(pos.size()));
        break;
      case type::VEL:
        _result.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            vel.data(), static_cast<long>(vel.size()));

        break;
      case type::ACC:
        _result.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            acc.data(), static_cast<long>(acc.size()));

        break;
      case type::JERK:
        _result.row(i) = Eigen::Map<const Eigen::RowVectorXd>(
            jerk.data(), static_cast<long>(jerk.size()));

        break;
      case type::ZERO:
        _result.row(i).setZero();
        break;
    }
  }
}

RuckigCurve* RuckigCurve::deriv_impl(std::size_t _deg) const {
  int nex_deg = static_cast<int>(type_) + static_cast<int>(_deg);

  switch (nex_deg) {
    case 0:
      return new RuckigCurve(ruckig_, RuckigCurve::type::POS);  // NOLINT
    case 1:
      return new RuckigCurve(ruckig_, RuckigCurve::type::VEL);
    case 2:
      return new RuckigCurve(ruckig_, RuckigCurve::type::ACC);
    case 3:
      return new RuckigCurve(ruckig_, RuckigCurve::type::JERK);
    default:
      return new RuckigCurve(ruckig_, RuckigCurve::type::ZERO);
  }
}

}  // namespace gsplines::ruckig
