#include <gsplines/Ruckig/Ruckig.hpp>
#include <ruckig/ruckig.hpp>

namespace gsplines::ruckig {

std::optional<RuckigCurve> interpolator(
    Eigen::Ref<const Eigen::MatrixXd> _waypoints,
    const std::vector<double>& _max_abs_vel,
    const std::vector<double>& _max_abs_acc,
    const std::vector<double>& _max_abs_jerk) {
  //
  std::size_t dof = static_cast<std::size_t>(_waypoints.cols());
  ::ruckig::InputParameter<::ruckig::DynamicDOFs> result(dof);

  std::size_t number_of_waypoints = static_cast<std::size_t>(_waypoints.rows());

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

  result.max_velocity = _max_abs_vel;
  result.min_acceleration = _max_abs_acc;
  result.max_jerk = _max_abs_jerk;
  result.current_velocity.assign(dof, 0.0);
  result.current_acceleration.assign(dof, 0.0);

  for (int i = 0; i < dof; ++i) {
    result.current_position[i] = _waypoints.topRows(1)(0, i);
    result.target_position[i] = _waypoints.bottomRows(1)(0, i);
  }

  std::vector<double> vec(dof, 0.0);
  result.intermediate_positions.reserve(number_of_waypoints - 2);
  for (int i = 1; i < number_of_waypoints - 1; ++i) {
    for (int j = 0; j < dof; ++j) {
      vec[j] = _waypoints(i, j);
    }
    result.intermediate_positions.emplace_back(vec);
  }
  ::ruckig::Ruckig<::ruckig::DynamicDOFs> otg(dof, 0.001);
  ::ruckig::Trajectory<::ruckig::DynamicDOFs> trajectory(dof);
  otg.calculate(result, trajectory);
  return RuckigCurve(trajectory, RuckigCurve::type::POS);
}

RuckigCurve::RuckigCurve(::ruckig::Trajectory<::ruckig::DynamicDOFs> _trj,
                         type _type)
    : FunctionInheritanceHelper({0.0, _trj.get_duration()},
                                _trj.degrees_of_freedom, "Ruckig"),
      ruckig_(_trj) {}

void RuckigCurve::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {
  std::vector<double> pos(get_codom_dim(), 0.0);
  std::vector<double> vel(get_codom_dim(), 0.0);
  std::vector<double> acc(get_codom_dim(), 0.0);
  std::vector<double> jerk(get_codom_dim(), 0.0);

  std::size_t sec;

  for (int i = 0; i < _domain_points.size(); ++i) {
    ruckig_.at_time(_domain_points(i), pos, vel, acc, jerk, sec);
    switch (type_) {
      case type::POS:
        _result.row(i) =
            Eigen::Map<const Eigen::RowVectorXd>(pos.data(), pos.size());

        break;
      case type::VEL:
        _result.row(i) =
            Eigen::Map<const Eigen::RowVectorXd>(vel.data(), vel.size());

        break;
      case type::ACC:
        _result.row(i) =
            Eigen::Map<const Eigen::RowVectorXd>(acc.data(), acc.size());

        break;
      case type::JERK:
        _result.row(i) =
            Eigen::Map<const Eigen::RowVectorXd>(acc.data(), acc.size());

        break;
      case type::ZERO:
        break;
    }
  }
}

RuckigCurve* RuckigCurve::deriv_impl(std::size_t _deg) const {
  int nex_deg = static_cast<int>(type_) + _deg;

  switch (nex_deg) {
    case 0:
      return new RuckigCurve(ruckig_, RuckigCurve::type::POS);
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
