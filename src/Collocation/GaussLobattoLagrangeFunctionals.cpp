#include <gsplines/Collocation/GaussLobattoLagrangeFunctionals.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
namespace gsplines {

namespace collocation {
/*
LinearFunctionalSparse::operator LinearFunctionalDense() const {
  return LinearFunctionalDense(mat_.toDense());
}
LinearFunctionalDense::operator LinearFunctionalSparse() const {
  return LinearFunctionalSparse(mat_.sparseView());
}
*/

Integral::Integral(std::tuple<double, double> _domain, std::size_t _nglp,
                   std::size_t _n_intervals)
    : LinearFunctionalDense(
          legendre_gauss_lobatto_weights(_domain, _nglp, _n_intervals)
              .transpose()) {}

Integral::Integral(const GaussLobattoLagrangeSpline &_in)
    : Integral(_in.get_domain(), _in.get_basis().get_dim(),
               _in.get_number_of_intervals()) {}

SobolevDistance::SobolevDistance(const gsplines::functions::FunctionBase &_fun,
                                 std::size_t _nglp, std::size_t _n_inter,
                                 std::size_t /*_deg*/)
    : Functional(1),
      approx_(GaussLobattoLagrangeSpline::approximate(_fun, _nglp, _n_inter)),
      approx_d_(GaussLobattoLagrangeSpline::approximate(_fun.derivate(), _nglp,
                                                        _n_inter)),
      int_(approx_), der_(approx_),
      diff_(std::make_shared<LinearFunctionalDense>(
          1u, _nglp * _n_inter * _fun.get_codom_dim())) {}

Eigen::VectorXd
SobolevDistance::operator()(const GaussLobattoLagrangeSpline &_in) const {
  Derivative der_(_in);
  TransposeLeftMultiplication tr1_(approx_ - _in);
  TransposeLeftMultiplication tr2_(approx_d_ - der_ * _in);
  return int_(tr1_ * (approx_ - _in) + tr2_ * (approx_d_ - der_ * _in));
}

std::shared_ptr<LinearFunctionalBase>
SobolevDistance::differential(const GaussLobattoLagrangeSpline &_in) const {
  TransposeLeftMultiplication tr1_(approx_ - _in);
  TransposeLeftMultiplication tr2_(approx_d_ - _in.derivate());
  diff_->set(int_ * (tr1_ + tr2_ * der_) * der_);
  return diff_;
}

/*
const Eigen::MatrixXd &
Integral::derivative(const GaussLobattoLagrangeSpline & _in) const {
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
*/

} // namespace collocation
} // namespace gsplines
