#include <gsplines/Optimization/ipopt_interface.hpp>

namespace gsplines {
namespace optimization {

TimeSegmentLenghtsVar::TimeSegmentLenghtsVar(std::size_t _num_intervals,
                                             double _exec_time)
    : VariableSet(static_cast<int>(_num_intervals), "TimeSegmentLenghtsVar"),
      exec_time_(_exec_time),
      values_(_num_intervals),
      bounds_(GetRows()) {
  values_.setOnes();
  values_ *= exec_time_ / GetRows();
  ifopt::Bounds default_bound(0.0, ifopt::inf);
  bounds_ = ifopt::Component::VecBound(GetRows(), default_bound);
}

void TimeSegmentLenghtsVar::SetVariables(const Eigen::VectorXd& _vec) {
  values_ = _vec;
}

Eigen::VectorXd TimeSegmentLenghtsVar::GetValues() const { return values_; }

ifopt::Component::VecBound TimeSegmentLenghtsVar::GetBounds() const {
  return bounds_;
}

ExecTimeConstraint::ExecTimeConstraint(std::size_t _num_intervals,
                                       double _exec_time)
    : ConstraintSet(1, "ExecTimeConstraint"), values_(1) {
  ifopt::Bounds default_bound(_exec_time, _exec_time);
  bounds_ = ifopt::Component::VecBound(1, default_bound);
}

Eigen::VectorXd ExecTimeConstraint::GetValues() const {
  Eigen::VectorXd result(1);

  result(0) =
      GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues().sum();
  return result;
}
ifopt::Component::VecBound ExecTimeConstraint::GetBounds() const {
  return bounds_;
}

void ExecTimeConstraint::FillJacobianBlock(std::string _set_name,
                                           Jacobian& _jac_block) const {
  (void)_set_name;
  for (unsigned int i = 0;
       i < GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows();
       i++) {
    _jac_block.coeffRef(0, i) = 1.0;
  }
}

SobolevNorm::SobolevNorm(
    const std::string& _name,
    const Eigen::Ref<const Eigen::MatrixXd>& _waypoints,
    const gsplines::basis::Basis& _basis,
    const std::vector<std::pair<std::size_t, double>>& _weights)
    : CostTerm(_name),
      basis_(_basis.clone()),
      weights_(_weights),
      waypoints_(_waypoints),
      sobol_norm_(_waypoints, _basis, _weights),
      buff_1_(_waypoints.rows() - 1),
      buff_2_(_waypoints.rows() - 1) {}

double SobolevNorm::GetCost() const {
  buff_1_.noalias() =
      GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues();
  return sobol_norm_(buff_1_);
}
void SobolevNorm::FillJacobianBlock(std::string _var_set,
                                    Jacobian& _jac) const {
  (void)_var_set;
  buff_1_.noalias() =
      GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues();

  sobol_norm_.deriv_wrt_interval_len(buff_1_, buff_2_);

  for (unsigned int i = 0;
       i < GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows();
       i++) {
    _jac.coeffRef(0, i) = buff_2_(i);
  }
}
}  // namespace optimization
}  // namespace gsplines
