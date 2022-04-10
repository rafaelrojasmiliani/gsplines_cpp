#include <gsplines/Collocation/GaussLobattoLagrangeFunctionals.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
namespace gsplines {

namespace collocation {

SobolevError::SobolevError(const gsplines::functions::FunctionBase &_fun,
                           std::size_t _nglp, std::size_t _n_inter)
    : Functional(1),
      approx_(GaussLobattoLagrangeSpline::approximate(_fun, _nglp, _n_inter)),
      approx_d_(GaussLobattoLagrangeSpline::approximate(_fun.derivate(), _nglp,
                                                        _n_inter)),
      der_(1, _fun.get_codom_dim() * _nglp * _n_inter),
      int_(_fun.get_domain(), _nglp, _n_inter) {}

Eigen::VectorXd
SobolevError::operator()(const GaussLobattoLagrangeSpline &_in) {
  Derivative der_(_in);
  TransposeLeftMultiplication tr1_(approx_ - _in);
  TransposeLeftMultiplication tr2_(approx_d_ - der_ * _in);
  return int_(tr1_ * (approx_ - _in) + tr2_ * (approx_d_ - der_ * _in));
}

const Eigen::MatrixXd &
SobolevError::derivative(const GaussLobattoLagrangeSpline &_in) const {
  Derivative der(_in);
  TransposeLeftMultiplication tr1_(approx_ - _in);
  TransposeLeftMultiplication tr2_(approx_d_ - der * _in);
  der_ = int_.derivative(_in) * (tr1_ + tr2_ * der).to_matrix();
  return der_;
}
Integral::Integral(std::tuple<double, double> _domain, std::size_t _nglp,
                   std::size_t _n_intervals)
    : Functional(1),
      glw_(legendre_gauss_lobatto_weights(_domain, _nglp, _n_intervals)
               .transpose()) {}

Eigen::VectorXd
Integral::operator()(const GaussLobattoLagrangeSpline &_in) const {
  return glw_ * _in.get_coefficients();
}
const Eigen::MatrixXd &
Integral::derivative(const GaussLobattoLagrangeSpline & /*_in*/) const {
  return glw_;
}

GLLSplineVariable::GLLSplineVariable(const GaussLobattoLagrangeSpline &_in)
    : ifopt::VariableSet(_in.get_coefficients().size(), "GLL"),
      GaussLobattoLagrangeSpline(_in) {

  ifopt::Bounds default_bound(-ifopt::inf, ifopt::inf);
  bounds_ = ifopt::Component::VecBound(GetRows(), default_bound);
}
void GLLSplineVariable::SetVariables(const Eigen::VectorXd &_vec) {
  coefficients_ = _vec;
}
Eigen::VectorXd GLLSplineVariable::GetValues() const {
  return get_coefficients();
}

ConstraintWrapper::ConstraintWrapper(Functional<T> _fun, double _min,
                                     double _max);
ConstraintWrapper::ConstraintWrapper(Functional<T> _fun,
                                     const std::vector<double> &_min,
                                     double _max);
ConstraintWrapper::ConstraintWrapper(Functional<T> _fun, double _min,
                                     const std::vector<double> &_max);
ConstraintWrapper::ConstraintWrapper(Functional<T> _fun,
                                     const std::vector<double> &_min,
                                     const std::vector<double> &_max);

Eigen::VectorXd ConstraintWrapper::GetValues() const override;
ifopt::Component::VecBound ConstraintWrapper::GetBounds() const override;

void ConstraintWrapper::FillJacobianBlock(std::string _set_name,
                                          Jacobian &_jac_block) const override;

} // namespace collocation
} // namespace gsplines
