#pragma once
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <ruckig/ruckig.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace gsplines::ruckig {

class RuckigCurve : public functions::FunctionInheritanceHelper<
                        RuckigCurve, functions::Function, RuckigCurve> {
 public:
  friend std::optional<RuckigCurve> interpolator(
      Eigen::Ref<const Eigen::MatrixXd> _waypoints,
      const std::vector<double>& _max_abs_vel,
      const std::vector<double>& _max_abs_acc,
      const std::vector<double>& _max_abs_jerk);
  ~RuckigCurve() override = default;

 private:
  enum class type : int {
    POS = 0,
    VEL,
    ACC,
    JERK,
    ZERO

  };
  ::ruckig::Trajectory<::ruckig::DynamicDOFs> ruckig_;
  type type_;
  RuckigCurve(::ruckig::Trajectory<::ruckig::DynamicDOFs> _trj, type _type);
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;
  RuckigCurve* deriv_impl(std::size_t _deg) const override;
};

std::optional<RuckigCurve> interpolator(
    Eigen::Ref<const Eigen::MatrixXd> _waypoints,
    const std::vector<double>& _max_abs_vel,
    const std::vector<double>& _max_abs_acc,
    const std::vector<double>& _max_abs_jerk);
}  // namespace gsplines::ruckig
