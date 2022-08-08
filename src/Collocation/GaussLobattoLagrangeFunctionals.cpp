#include <algorithm>
#include <gsplines/Collocation/GaussLobattoLagrangeFunctionals.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <iostream>
namespace gsplines {

namespace collocation {

// ----------------
// Derivative
// ----------------

Derivative::Derivative(std::pair<double, double> _domain,
                       std::size_t _codom_dim, std::size_t _n_glp,
                       std::size_t _n_intervals, std::size_t _deg)
    : Derivative(_codom_dim, _n_glp,
                 Eigen::VectorXd::Ones(_n_intervals) *
                     (_domain.second - _domain.first) / _n_intervals,
                 _deg) {}

Derivative::Derivative(std::size_t _codom_dim, std::size_t _n_glp,
                       const Eigen::VectorXd &_interval_lengths,
                       std::size_t _deg)
    : LinearFunctionalSparse(_n_glp * _interval_lengths.size() * _codom_dim,
                             _n_glp * _interval_lengths.size() * _codom_dim),
      basis_(_n_glp), n_intervals_(_interval_lengths.size()),
      codom_dim_(_codom_dim), deg_(_deg) {

  mat_ = basis_.gspline_derivative_matrix(_interval_lengths.size(), codom_dim_,
                                          deg_, _interval_lengths);
  // ---
  mat_.makeCompressed();
}

Derivative::Derivative(const GaussLobattoLagrangeSpline &_that)
    : Derivative(_that.get_codom_dim(), _that.get_basis().get_dim(),
                 _that.get_interval_lengths()) {}

void Derivative::update(const Eigen::VectorXd &_interval_lengths) {

  mat_ = basis_.gspline_derivative_matrix(_interval_lengths.size(), codom_dim_,
                                          deg_, _interval_lengths);

  mat_.makeCompressed();
}

void Derivative::update(double _interval_length) {

  mat_ = basis_.gspline_derivative_matrix(
      n_intervals_, codom_dim_, deg_,
      Eigen::VectorXd::Constant(n_intervals_, _interval_length));
  mat_.makeCompressed();
}

// ----------------
// Integral
// ----------------
Integral::Integral(std::tuple<double, double> _domain, std::size_t _nglp,
                   std::size_t _n_intervals)
    : LinearFunctionalDense(
          legendre_gauss_lobatto_weights(_domain, _nglp, _n_intervals)
              .transpose()),
      nglp_(_nglp) {}

Integral::Integral(const GaussLobattoLagrangeSpline &_in)
    : Integral(_in.get_domain(), _in.get_basis().get_dim(),
               _in.get_number_of_intervals()) {}

void Integral::update(std::tuple<double, double> _domain, std::size_t _nglp,
                      std::size_t _n_intervals) {
  mat_.noalias() =
      legendre_gauss_lobatto_weights(_domain, _nglp, _n_intervals).transpose();
}

void Integral::update(const Eigen::VectorXd &_intervals) {
  if (static_cast<long>(_intervals.size() * nglp_) != mat_.rows()) {
    throw std::invalid_argument(" Invaled size of intervals vector");
  }
  mat_.noalias() =
      legendre_gauss_lobatto_weights(nglp_, _intervals).transpose();
}

// ----------------
// Sobolev Distance
// ----------------
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
  TransposeLeftMultiplication tr1_(approx_ - _in);
  TransposeLeftMultiplication tr2_(approx_d_ - der_ * _in);

  return int_(tr1_ * (approx_ - _in) + tr2_ * (approx_d_ - der_ * _in));
}

std::shared_ptr<LinearFunctionalBase>
SobolevDistance::differential(const GaussLobattoLagrangeSpline &_in) const {
  TransposeLeftMultiplication tr1_(approx_ - _in);
  TransposeLeftMultiplication tr2_(approx_d_ - der_ * _in);

  diff_->set(-2 * int_ * (tr1_ + tr2_ * der_));

  return diff_;
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

} // namespace collocation
} // namespace gsplines
