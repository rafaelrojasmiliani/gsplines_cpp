#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Tools.hpp>
#include <iostream>

namespace gsplines {

GSpline *GSpline::deriv_impl(std::size_t _deg) const {

  Eigen::VectorXd result_coeff(coefficients_);
  std::size_t interval_coor;
  std::size_t codom_coor;
  std::size_t i0;

  if (_deg > 0) {
    for (interval_coor = 0; interval_coor < get_number_of_intervals();
         interval_coor++) {
      for (codom_coor = 0; codom_coor < get_codom_dim(); codom_coor++) {
        i0 = interval_coor * basis_->get_dim() * get_codom_dim() +
             basis_->get_dim() * codom_coor;
        result_coeff.segment(i0, basis_->get_dim()) =
            basis_->get_derivative_matrix_block(_deg) *
            result_coeff.segment(i0, basis_->get_dim()) *
            std::pow(2 / domain_interval_lengths_(interval_coor), _deg);
      }
    }
  }

  return new GSpline(get_domain(), get_codom_dim(), get_number_of_intervals(),
                     *basis_, result_coeff, domain_interval_lengths_,
                     get_name());
}

GSpline::GSpline(const GSpline &that)
    : FunctionInheritanceHelper(that), coefficients_(that.coefficients_),
      domain_interval_lengths_(that.domain_interval_lengths_),
      basis_(that.basis_->clone()), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(get_number_of_intervals() * basis_->get_dim() * get_codom_dim())) {
    throw std::invalid_argument(
        "GSpline instantation Error: The number of coefficients is incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(that.get_codom_dim()) +
        " req num inter: " + std::to_string(that.get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

GSpline::GSpline(GSpline &&that)
    : FunctionInheritanceHelper(std::move(that)),
      coefficients_(std::move(that.coefficients_)),
      domain_interval_lengths_(std::move(that.domain_interval_lengths_)),
      basis_(that.basis_->move_clone()), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(get_number_of_intervals() * basis_->get_dim() * get_codom_dim())) {
    throw std::invalid_argument(
        "GSpline instantation Error: The number of coefficients is incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(that.get_codom_dim()) +
        " req num inter: " + std::to_string(that.get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

GSpline::GSpline(std::pair<double, double> _domain, std::size_t _codom_dim,
                 std::size_t _n_intervals, const basis::Basis &_basis,
                 const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                 const Eigen::Ref<const Eigen::VectorXd> _tauv,
                 const std::string &_name)
    : FunctionInheritanceHelper(_domain, _codom_dim, _name),
      coefficients_(_coefficents), domain_interval_lengths_(_tauv),
      basis_(_basis.clone()), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      (long)(_n_intervals * basis_->get_dim() * _codom_dim)) {
    throw std::invalid_argument(
        "GSpline instantation Error: The number of coefficients is incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(get_codom_dim()) +
        " req num inter: " + std::to_string(get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
}

void GSpline::value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
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

GSpline
GSpline::linear_scaling_new_execution_time(double _new_exec_time) const {
  assert(_new_exec_time > 0);
  Eigen::VectorXd new_domain_interva_lengths =
      domain_interval_lengths_ * _new_exec_time / get_domain_length();

  gsplines::Interpolator inter(get_codom_dim(), get_number_of_intervals(),
                               *basis_);

  return inter.interpolate(new_domain_interva_lengths, get_waypoints());
}

std::size_t GSpline::get_interval(double _domain_point) const {
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
GSpline::coefficient_segment(std::size_t _interval,
                             std::size_t _component) const {
  int i0 = _interval * basis_->get_dim() * get_codom_dim() +
           basis_->get_dim() * _component;
  return coefficients_.segment(i0, basis_->get_dim());
}

Eigen::Ref<Eigen::VectorXd>
GSpline::coefficient_segment(std::size_t _interval, std::size_t _component) {
  int i0 = _interval * basis_->get_dim() * get_codom_dim() +
           basis_->get_dim() * _component;
  return coefficients_.segment(i0, basis_->get_dim());
}

double GSpline::interval_to_window(double _domain_point,
                                   std::size_t _interval) const {
  double left_breakpoint =
      get_domain().first +
      domain_interval_lengths_.segment(0, _interval).array().sum();
  return 2.0 * (_domain_point - left_breakpoint) /
             domain_interval_lengths_[_interval] -
         1.0;
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
Eigen::VectorXd GSpline::get_domain_breakpoints() const {

  double time_instant = get_domain().first;
  Eigen::VectorXd result(get_number_of_intervals() + 1);
  result(0) = time_instant;
  for (std::size_t i = 1; i < get_number_of_intervals() + 1; i++) {
    time_instant += domain_interval_lengths_(i - 1);
    result(i) = time_instant;
  }
  return result;
}

Eigen::MatrixXd GSpline::get_waypoints() const {
  Eigen::MatrixXd result(get_number_of_intervals() + 1, get_codom_dim());
  GSpline::value(get_domain_breakpoints(), result);
  return result;
}

bool GSpline::compatible(const GSpline &_that) const {
  return get_basis() == _that.get_basis() and
         tools::approx_equal(get_interval_lengths(),
                             _that.get_interval_lengths(), 1.0e-8);
}

GSpline GSpline::operator+(const GSpline &that) const & {
  if (not compatible(that))
    throw std::invalid_argument("Cannot sum Incompatible Gspline");
  GSpline result = GSpline(*this);
  result.coefficients_ += that.coefficients_;
  return result;
}
GSpline GSpline::operator+(GSpline &&that) const & {
  if (not compatible(that))
    throw std::invalid_argument("Cannot sum Incompatible Gspline");
  that.coefficients_ += coefficients_;
  return std::move(that);
}
GSpline GSpline::operator+(const GSpline &that) && {

  if (not compatible(that))
    throw std::invalid_argument("Cannot sum Incompatible Gspline");
  coefficients_ += that.coefficients_;
  return std::move(*this);
}

GSpline GSpline::operator-(const GSpline &that) const & {

  if (not compatible(that))
    throw std::invalid_argument("Cannot sum Incompatible Gspline");
  GSpline result = GSpline(*this);
  result.coefficients_ -= that.coefficients_;
  return result;
}
GSpline GSpline::operator-(GSpline &&that) const & {

  if (not compatible(that))
    throw std::invalid_argument("Cannot sum Incompatible Gspline");
  that.coefficients_ -= coefficients_;
  return std::move(that);
}
GSpline GSpline::operator-(const GSpline &that) && {

  if (not compatible(that))
    throw std::invalid_argument("Cannot sum Incompatible Gspline");
  coefficients_ -= that.coefficients_;
  return std::move(*this);
}
} // namespace gsplines
