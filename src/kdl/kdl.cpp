
#include <gsplines/kdl/kdl.hpp>
#include <kdl/velocityprofile_trap.hpp>
#include <ruckig/result.hpp>
#include <iostream>

namespace gsplines::kdl {

std::pair<std::size_t, double> get_interval(
    double _domain_point, const Eigen::VectorXd& domain_interval_lengths_) {
  double left_breakpoint = 0.0;
  double right_breakpoint = 0.0;
  if (_domain_point <= left_breakpoint) {
    return {0, 0.0};
  }

  for (std::size_t i = 0; i < domain_interval_lengths_.size(); i++) {
    right_breakpoint =
        left_breakpoint + domain_interval_lengths_(static_cast<long>(i));
    if (left_breakpoint < _domain_point and _domain_point <= right_breakpoint) {
      double t = _domain_point - left_breakpoint;
      return {i, t};
    }
    left_breakpoint = right_breakpoint;
  }
  return {domain_interval_lengths_.size() - 1, domain_interval_lengths_.sum()};
}

struct KdlTrap::Impl {
  std::vector<std::vector<::KDL::VelocityProfile_Trap>> velocity_profiles;
  Eigen::MatrixXd waypoints;
  Eigen::VectorXd time_intervals;
  std::pair<double, double> domain;
};
KdlTrap::~KdlTrap() = default;
KdlTrap::KdlTrap(const KdlTrap& other)
    : FunctionInheritanceHelper(other.impl_->domain, other.get_codom_dim(),
                                "kdl"),
      impl_(std::make_unique<KdlTrap::Impl>(*other.impl_)),
      type_(other.type_)

{}

inline auto std_vector_to_eigen(std::vector<double>& _int) {
  return Eigen::Map<Eigen::VectorXd>(_int.data(),
                                     static_cast<long>(_int.size()));
}

inline auto std_vector_to_eigen(const std::vector<double>& _int) {
  return Eigen::Map<const Eigen::VectorXd>(_int.data(),
                                           static_cast<long>(_int.size()));
}

std::optional<KdlTrap> interpolator(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints,
    const std::vector<double>& _max_abs_vel,
    const std::vector<double>& _max_abs_acc) {
  //
  auto dof = static_cast<std::size_t>(_waypoints.cols());

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
  auto number_of_intervals = static_cast<std::size_t>(_waypoints.rows() - 1);

  std::vector<std::vector<::KDL::VelocityProfile_Trap>> velocity_profiles;

  for (int i = 0; i < number_of_intervals; ++i) {
    std::cout << "runig i = " << i << "\n";
    velocity_profiles.emplace_back();
    double max_time = 0;
    for (int j = 0; j < dof; ++j) {
      std::cout << "runig j = " << j << "\n";
      velocity_profiles.back().emplace_back(_max_abs_vel[j], _max_abs_acc[j]);
      velocity_profiles.back().back().SetProfile(_waypoints(i, j),
                                                 _waypoints(i + 1, j));
      double time = velocity_profiles.back().back().Duration();
      if (time > max_time) {
        max_time = time;
      }
      /// Sync joints, make all joints arrive to their goal at the same time.
      /// This is done by setting all joints to last the same amount of time
      /// during its trajectory.
    }
    std::cout << " starting with set duration\n";
    for (int j = 0; j < dof; ++j) {
      velocity_profiles.back()[j].SetProfileDuration(
          _waypoints(i, j), _waypoints(i + 1, j), max_time);
    }
  }

  KdlTrap::Impl impl;
  impl.time_intervals.setZero(static_cast<long>(number_of_intervals));
  for (int i = 0; i < number_of_intervals; ++i) {
    impl.time_intervals(i) = velocity_profiles[i].front().Duration();
  }
  impl.waypoints = _waypoints;
  impl.velocity_profiles = std::move(velocity_profiles);
  impl.domain = std::make_pair<double, double>(0.0, impl.time_intervals.sum());
  return KdlTrap(impl, KdlTrap::type::POS);
}

KdlTrap::KdlTrap(const Impl& _impl, type _type)
    : FunctionInheritanceHelper(_impl.domain,
                                _impl.velocity_profiles.front().size(), "kdl"),
      impl_(std::make_unique<KdlTrap::Impl>(_impl)),
      type_(_type) {}

void KdlTrap::value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                         Eigen::Ref<Eigen::MatrixXd> _result) const {
  std::vector<double> pos(get_codom_dim(), 0.0);
  std::vector<double> vel(get_codom_dim(), 0.0);
  std::vector<double> acc(get_codom_dim(), 0.0);

  std::size_t sec = 0;

  const std::size_t result_cols(_domain_points.size());
  std::size_t current_interval = 0;
  double s = 0.0;
  double tau = 0.0;

  std::size_t i = 0;
  std::size_t j = 0;

  for (i = 0; i < result_cols; i++) {
    auto [current_interval, t] = get_interval(
        _domain_points(static_cast<long>(i)), impl_->time_intervals);
    const auto& vel_profile_array = impl_->velocity_profiles[current_interval];
    for (j = 0; j < get_codom_dim(); j++) {
      switch (type_) {
        case type::POS:
          _result(static_cast<long>(i), static_cast<long>(j)) =
              impl_->velocity_profiles[current_interval][j].Pos(t);
          break;
        case type::VEL:
          _result(static_cast<long>(i), static_cast<long>(j)) =
              impl_->velocity_profiles[current_interval][j].Vel(t);

          break;
        case type::ACC:
          _result(static_cast<long>(i), static_cast<long>(j)) =
              impl_->velocity_profiles[current_interval][j].Acc(t);
          break;
        case type::ZERO:
          _result(static_cast<long>(i), static_cast<long>(j)) = 0.0;
          break;
      }
    }
  }
}

KdlTrap* KdlTrap::deriv_impl(std::size_t _deg) const {
  int nex_deg = static_cast<int>(type_) + static_cast<int>(_deg);

  switch (nex_deg) {
    case 0:
      return new KdlTrap(*impl_, KdlTrap::type::POS);  // NOLINT
    case 1:
      return new KdlTrap(*impl_, KdlTrap::type::VEL);
    case 2:
      return new KdlTrap(*impl_, KdlTrap::type::ACC);
    default:
      return new KdlTrap(*impl_, KdlTrap::type::ZERO);
  }
}

}  // namespace gsplines::kdl
