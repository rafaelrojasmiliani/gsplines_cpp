#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Tools.hpp>
#include <iostream>
#include <random>

namespace gsplines {

GSplineBase::GSplineBase(const GSplineBase &that)
    : FunctionInheritanceHelper(that), coefficients_(that.coefficients_),
      domain_interval_lengths_(that.domain_interval_lengths_),
      basis_(that.basis_->clone()), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(get_number_of_intervals() * basis_->get_dim() * get_codom_dim())) {
    throw std::invalid_argument(
        "GSplineBase instantation Error: The number of coefficients is "
        "incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(that.get_codom_dim()) +
        " req num inter: " + std::to_string(that.get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

GSplineBase::GSplineBase(GSplineBase &&that)
    : FunctionInheritanceHelper(std::move(that)),
      coefficients_(std::move(that.coefficients_)),
      domain_interval_lengths_(std::move(that.domain_interval_lengths_)),
      basis_(that.basis_->move_clone()), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(get_number_of_intervals() * basis_->get_dim() * get_codom_dim())) {
    throw std::invalid_argument(
        "GSplineBase instantation Error: The number of coefficients is "
        "incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(that.get_codom_dim()) +
        " req num inter: " + std::to_string(that.get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

GSplineBase::GSplineBase(std::pair<double, double> _domain,
                         std::size_t _codom_dim, std::size_t _n_intervals,
                         const basis::Basis &_basis,
                         const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                         const Eigen::Ref<const Eigen::VectorXd> _tauv,
                         const std::string &_name)
    : FunctionInheritanceHelper(_domain, _codom_dim, _name),
      coefficients_(_coefficents), domain_interval_lengths_(_tauv),
      basis_(_basis.clone()), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(_n_intervals * basis_->get_dim() * _codom_dim)) {
    throw std::invalid_argument(
        "GSplineBase instantation Error: The number of coefficients is "
        "incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(get_codom_dim()) +
        " req num inter: " + std::to_string(get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

GSplineBase::GSplineBase(std::pair<double, double> _domain,
                         std::size_t _codom_dim, std::size_t _n_intervals,
                         const basis::Basis &_basis,
                         Eigen::VectorXd &&_coefficents,
                         Eigen::VectorXd &&_tauv, const std::string &_name)
    : FunctionInheritanceHelper(_domain, _codom_dim, _name),
      coefficients_(std::move(_coefficents)),
      domain_interval_lengths_(std::move(_tauv)), basis_(_basis.clone()),
      basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(_n_intervals * basis_->get_dim() * _codom_dim)) {
    throw std::invalid_argument(
        "GSplineBase instantation Error: The number of coefficients is "
        "incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(get_codom_dim()) +
        " req num inter: " + std::to_string(get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

void GSplineBase::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {

  std::size_t result_cols(_domain_points.size());
  std::size_t current_interval = 0;
  double s, tau;

  std::size_t i, j;

  for (i = 0; i < result_cols; i++) {
    current_interval = get_interval(_domain_points(i));
    s = interval_to_window(_domain_points(i), current_interval);
    tau = domain_interval_lengths_(current_interval);
    basis_->eval_on_window(s, tau, basis_buffer_);
    for (j = 0; j < get_codom_dim(); j++) {
      _result(i, j) =
          coefficient_segment(current_interval, j).adjoint() * basis_buffer_;
    }
  }
}

std::size_t GSplineBase::get_interval(double _domain_point) const {
  double left_breakpoint = get_domain().first;
  double right_breakpoint;
  if (_domain_point <= left_breakpoint)
    return 0;

  for (std::size_t i = 0; i < get_number_of_intervals(); i++) {
    right_breakpoint = left_breakpoint + domain_interval_lengths_(i);
    if (left_breakpoint < _domain_point and _domain_point <= right_breakpoint) {
      return i;
    }
    left_breakpoint = right_breakpoint;
  }
  return get_number_of_intervals() - 1;
}

Eigen::Ref<const Eigen::VectorXd>
GSplineBase::coefficient_segment(std::size_t _interval,
                                 std::size_t _component) const {
  int i0 = _interval * basis_->get_dim() * get_codom_dim() +
           basis_->get_dim() * _component;
  return coefficients_.segment(i0, basis_->get_dim());
}

Eigen::Ref<Eigen::VectorXd>
GSplineBase::coefficient_segment(std::size_t _interval,
                                 std::size_t _component) {
  int i0 = _interval * basis_->get_dim() * get_codom_dim() +
           basis_->get_dim() * _component;
  return coefficients_.segment(i0, basis_->get_dim());
}

double GSplineBase::interval_to_window(double _domain_point,
                                       std::size_t _interval) const {
  double left_breakpoint =
      get_domain().first +
      domain_interval_lengths_.segment(0, _interval).array().sum();
  return 2.0 * (_domain_point - left_breakpoint) /
             domain_interval_lengths_[_interval] -
         1.0;
}

bool GSplineBase::operator==(const GSplineBase &_that) const {
  return tools::approx_equal(*this, _that, 1.0e-7);
}
bool GSplineBase::operator!=(const GSplineBase &_that) const {
  return not(*this == _that);
}

const Eigen::Ref<const Eigen::VectorXd>
get_coefficient_segment(const Eigen::Ref<const Eigen::VectorXd> _coefficients,
                        basis::Basis &_basis, std::size_t /*_num_interval*/,
                        std::size_t _codom_dim, std::size_t _interval,
                        std::size_t _component) {

  int i0 =
      _interval * _basis.get_dim() * _codom_dim + _basis.get_dim() * _component;
  return _coefficients.segment(i0, _basis.get_dim());
}
Eigen::VectorXd GSplineBase::get_domain_breakpoints() const {

  double time_instant = get_domain().first;
  Eigen::VectorXd result(get_number_of_intervals() + 1);
  result(0) = time_instant;
  for (std::size_t i = 1; i < get_number_of_intervals() + 1; i++) {
    time_instant += domain_interval_lengths_(i - 1);
    result(i) = time_instant;
  }
  return result;
}

Eigen::MatrixXd GSplineBase::get_waypoints() const {
  Eigen::MatrixXd result(get_number_of_intervals() + 1, get_codom_dim());
  GSplineBase::value(get_domain_breakpoints(), result);
  return result;
}

bool GSplineBase::same_vector_space(const GSplineBase &_that) const {
  return get_basis() == _that.get_basis() and
         get_codom_dim() == _that.get_codom_dim() and
         tools::approx_equal(get_interval_lengths(),
                             _that.get_interval_lengths(), 1.0e-8);
}

GSpline operator*(double _a, const GSpline &_that) {
  GSpline result(_that);
  result.coefficients_ *= _a;
  return result;
}
GSpline operator*(double _a, GSpline &&_that) {
  _that.coefficients_ *= _a;
  return std::move(_that);
}
GSpline operator*(const GSpline &_that, double _a) { return _a * _that; }
GSpline operator*(GSpline &&_that, double _a) { return _a * std::move(_that); }
GSpline operator-(const GSpline &_that) {

  GSpline result(_that);
  result.coefficients_ *= -1.0;
  return result;
}
GSpline operator-(GSpline &&_that) {

  _that.coefficients_ *= -1.0;
  return std::move(_that);
}

GSpline random_gspline(std::pair<double, double> _domain,
                       std::size_t _codom_dim) {

  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);
  std::size_t number_of_intervals = uint_dist(mt);

  Eigen::VectorXd tau =
      1.0 + (Eigen::VectorXd::Random(number_of_intervals).array()) / 2.0;

  Eigen::MatrixXd wp =
      Eigen::MatrixXd::Random(number_of_intervals + 1, _codom_dim);

  GSpline result = interpolate(tau, wp, *basis::BasisLegendre::get(6));
  return result.linear_scaling_new_execution_time(_domain.second);
}

} // namespace gsplines
