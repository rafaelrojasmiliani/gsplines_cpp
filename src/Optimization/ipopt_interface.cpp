#include <gsplines/Optimization/ipopt_interface.hpp>

namespace gsplines {
namespace optimization {

TimeSegmentLenghtsVar::TimeSegmentLenghtsVar(std::size_t _num_intervals,
                                             double _exec_time)
    : VariableSet(_num_intervals, "TimeSegmentLenghtsVar"),
      exec_time_(_exec_time), values_(_num_intervals), bounds_(GetRows()) {
  values_.setOnes();
  values_ *= exec_time_ / GetRows();
  ifopt::Bounds default_bound(0.0, ifopt::inf);
  bounds_ = ifopt::Component::VecBound(GetRows(), default_bound);
}
TimeSegmentLenghtsVar::~TimeSegmentLenghtsVar() {}

void TimeSegmentLenghtsVar::SetVariables(const Eigen::VectorXd &_vec) {
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

  static Eigen::VectorXd tauv(GetRows());
  static Eigen::VectorXd result(1);

  tauv = GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues();
  result(0) = tauv.sum();
  return result;
}
ifopt::Component::VecBound ExecTimeConstraint::GetBounds() const {
  return bounds_;
}

void ExecTimeConstraint::FillJacobianBlock(std::string _set_name,
                                           Jacobian &_jac_block) const {
  for (unsigned int i = 0;
       i < GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows();
       i++)
    _jac_block.coeffRef(0, i) = 1.0;
}

SobolevNorm::SobolevNorm(std::string _name,
                         const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
                         const gsplines::basis::Basis &_basis,
                         std::vector<std::pair<std::size_t, double>> _weights)
    : CostTerm(_name), basis_(_basis.clone()), weights_(_weights),
      waypoints_(_waypoints) {}

double SobolevNorm::GetCost() const {

  static Eigen::VectorXd tauv(GetVariables()->GetRows());
  static ::gsplines::functional_analysis::SobolevNorm sobol_norm(
      waypoints_, *basis_, weights_);
  tauv = GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues();
  return sobol_norm(tauv);
}
void SobolevNorm::FillJacobianBlock(std::string _var_set,
                                    Jacobian &_jac) const {
  static Eigen::VectorXd tauv(
      GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows());
  static ::gsplines::functional_analysis::SobolevNorm sobol_norm(
      waypoints_, *basis_, weights_);
  static Eigen::VectorXd result(
      GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows());

  tauv = GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetValues();

  sobol_norm.deriv_wrt_interval_len(tauv, result);

  for (unsigned int i = 0;
       i < GetVariables()->GetComponent("TimeSegmentLenghtsVar")->GetRows();
       i++)
    _jac.coeffRef(0, i) = result(i);
}
} // namespace optimization
} // namespace gsplines
