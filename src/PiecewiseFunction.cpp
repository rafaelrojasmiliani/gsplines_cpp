#include <gsplines++/PiecewiseFunction.hpp>
#include <iostream>

namespace gsplines {

PiecewiseFunction PiecewiseFunction::deriv(std::size_t _deg) {
  Eigen::VectorXd result_coeff(coefficients_);
  int der_coor;
  int interval_coor;
  int codom_coor;
  int i0;

  for (der_coor = 1; der_coor <= _deg; der_coor++) {
    for (interval_coor = 0; interval_coor < number_of_intervals_;
         interval_coor++) {
      for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
        i0 = interval_coor * basis_->get_dim() * codom_dim_ +
             basis_->get_dim() * codom_coor;
        result_coeff.segment(i0, basis_->get_dim()) =
            basis_->get_derivative_matrix().transpose() *
            result_coeff.segment(i0, basis_->get_dim()) * 2 /
            domain_interval_lengths_(interval_coor);
      }
    }
  }

  return PiecewiseFunction(codom_dim_, number_of_intervals_, *basis_,
                           result_coeff, domain_interval_lengths_);
}

PiecewiseFunction::PiecewiseFunction(const PiecewiseFunction &that)
    : codom_dim_(that.codom_dim_),
      number_of_intervals_(that.number_of_intervals_), basis_(that.basis_),
      coefficients_(that.coefficients_),
      domain_break_points_(that.domain_break_points_),
      domain_interval_lengths_(that.domain_interval_lengths_),
      basis_buffer_(that.basis_->get_dim()) {

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
    : codom_dim_(_codom_dim), number_of_intervals_(_n_intervals),
      basis_(&_basis), coefficients_(_coefficents),
      domain_break_points_(_n_intervals + 1), domain_interval_lengths_(_tauv),
      basis_buffer_(_basis.get_dim()) {

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

Eigen::MatrixXd PiecewiseFunction::
operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(), codom_dim_);
  std::size_t result_cols(_domain_points.size());
  std::size_t current_interval = 0;
  double s, tau;

  int i, j;

  for (i = 0; i < result_cols; i++) {
    current_interval = get_interval(_domain_points(i));
    s = interval_to_window(_domain_points(i), current_interval);
    tau = domain_interval_lengths_(current_interval);
    basis_->eval_on_window(s, tau, basis_buffer_);
    for (j = 0; j < codom_dim_; j++) {
      result(i, j) =
          coefficient_segment(current_interval, j).adjoint() * basis_buffer_;
    }
  }
  return result;
}

std::size_t PiecewiseFunction::get_interval(double _domain_point) const {
  if (_domain_point <= domain_break_points_(0))
    return 0;
  for (int i = 0; i < number_of_intervals_; i++) {
    if (domain_break_points_(i) < _domain_point and
        _domain_point <= domain_break_points_(i + 1)) {
      return i;
    }
  }
  return number_of_intervals_ - 1;
}

Eigen::Ref<Eigen::VectorXd>
PiecewiseFunction::coefficient_segment(std::size_t _interval,
                                       std::size_t _component) {
  int i0 = _interval * basis_->get_dim() * codom_dim_ +
           basis_->get_dim() * _component;
  return coefficients_.segment(i0, basis_->get_dim());
}

double PiecewiseFunction::interval_to_window(double _domain_point,
                                             std::size_t _interval) const {
  return 2.0 * (_domain_point - domain_break_points_[_interval]) /
             domain_interval_lengths_[_interval] -
         1.0;
}

Eigen::VectorXd PiecewiseFunction::get_coeff() { return coefficients_; }

const Eigen::Ref<const Eigen::VectorXd>
get_coefficient_segment(const Eigen::Ref<const Eigen::VectorXd> _coefficients,
                        basis::Basis &_basis, std::size_t _num_interval,
                        std::size_t _codom_dim, std::size_t _interval,
                        std::size_t _component) {

  int i0 =
      _interval * _basis.get_dim() * _codom_dim + _basis.get_dim() * _component;
  return _coefficients.segment(i0, _basis.get_dim());
}
} // namespace gsplines
