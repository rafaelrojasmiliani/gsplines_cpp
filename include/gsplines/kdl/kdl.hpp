#pragma once
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

namespace gsplines::kdl {

class KdlTrap
    : public functions::FunctionInheritanceHelper<KdlTrap, functions::Function,
                                                  KdlTrap> {
 public:
  friend std::optional<KdlTrap> interpolator(
      const Eigen::Ref<const Eigen::MatrixXd>& _waypoints,
      const std::vector<double>& _max_abs_vel,
      const std::vector<double>& _max_abs_acc);
  ~KdlTrap() override;

  KdlTrap(const KdlTrap& other);

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
  enum class type : int {

    POS = 0,
    VEL,
    ACC,
    ZERO

  };
  type type_;
  KdlTrap(const Impl&, type _type);
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;
  [[nodiscard]] KdlTrap* deriv_impl(std::size_t _deg) const override;
};

std::optional<KdlTrap> interpolator(
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints,
    const std::vector<double>& _max_abs_vel,
    const std::vector<double>& _max_abs_acc);
}  // namespace gsplines::kdl
