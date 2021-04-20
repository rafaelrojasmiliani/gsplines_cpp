#include <gsplines++/PiecewiseFunction.hpp>
#include <iostream>

namespace gsplines {

PiecewiseFunction PiecewiseFunction::deriv(std::size_t _deg) {
  Eigen::VectorXd result_coeff(coefficients_);
  for (int d = 1; d <= _deg; d++) {
    for (int interval = 0; interval < number_of_intervals_; interval++) {
      for (int comp = 0; comp < codom_dim_; comp++) {
        int i0 = interval * comp + comp * basis_->get_dim();
        result_coeff.segment(i0, basis_->get_dim()) =
            basis_->get_derivative_matrix().transpose() *
            result_coeff.segment(i0, basis_->get_dim()) * 2 /
            domain_interval_lengths_(interval);
      }
    }
  }

  return PiecewiseFunction(codom_dim_, number_of_intervals_, *basis_,
                           result_coeff, domain_interval_lengths_);
}

PiecewiseFunction::PiecewiseFunction(const PiecewiseFunction &that)
    : coefficients_(that.coefficients_), codom_dim_(that.codom_dim_),
      number_of_intervals_(that.number_of_intervals_),
      domain_break_points_(that.domain_break_points_),
      domain_interval_lengths_(that.domain_interval_lengths_),
      basis_(that.basis_), basis_buffer_(that.basis_->get_dim()) {
  if (coefficients_.size() !=
      number_of_intervals_ * basis_->get_dim() * codom_dim_) {
    printf("Error: The number of coefficients is incorrect\n");
    fflush(stdout);
  }
}

PiecewiseFunction::PiecewiseFunction(
    std::size_t _codom_dim, std::size_t _n_intervals, basis::Basis &_basis,
    const Eigen::Ref<const Eigen::VectorXd> _coefficents,
    const Eigen::Ref<const Eigen::VectorXd> _tauv)
    : coefficients_(_coefficents), codom_dim_(_codom_dim),
      number_of_intervals_(_n_intervals),
      domain_break_points_(_n_intervals + 1), domain_interval_lengths_(_tauv),
      basis_(&_basis), basis_buffer_(_basis.get_dim()) {
  if (coefficients_.size() != _n_intervals * basis_->get_dim() * _codom_dim) {
    printf("Error: The number of coefficients is incorrect\n");
    fflush(stdout);
  }
  double time_instant = 0.0;
  domain_break_points_(0) = time_instant;
  for (int i = 1; i < number_of_intervals_ + 1; i++) {
    time_instant += domain_interval_lengths_(i - 1);
    domain_break_points_(i) = time_instant;
  }
}

PiecewiseFunction::~PiecewiseFunction() {}

Eigen::MatrixXd PiecewiseFunction::operator()(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(), codom_dim_);
  std::size_t result_cols(_domain_points.size());
  std::size_t current_interval = 0;
  double s, tau;

  for (int i = 0; i < result_cols; i++) {
    current_interval = get_interval(_domain_points(i));
    s = interval_to_window(_domain_points(i), current_interval);
    tau = domain_interval_lengths_(current_interval);
    basis_->eval_on_window(s, tau, basis_buffer_);
    for (int j = 0; j < codom_dim_; j++) {
      result(i, j) =
          coefficient_segment(current_interval, j).adjoint() * basis_buffer_;
    }
  }
  return result;
}

std::size_t PiecewiseFunction::get_interval(double _domain_point) const {
  if (_domain_point <= domain_break_points_(0))
    return 0.0;
  for (int i = 0; i < number_of_intervals_; i++) {
    if (domain_break_points_(i) < _domain_point and
        _domain_point <= domain_break_points_(i + 1)) {
      return i;
    }
  }
  if (domain_break_points_.tail(1)(0) < _domain_point)
    return number_of_intervals_ - 1;
}

Eigen::Ref<Eigen::VectorXd>
PiecewiseFunction::coefficient_segment(std::size_t _interval,
                                       std::size_t _component) {
  int i0 = _interval * _component + _component * basis_->get_dim();
  return coefficients_.segment(i0, basis_->get_dim());
}

double PiecewiseFunction::interval_to_window(double _domain_point,
                                             std::size_t _interval) const {
  return 2.0 * (_domain_point - domain_break_points_[_interval]) /
             domain_interval_lengths_[_interval] -
         1.0;
}
} // namespace gsplines
