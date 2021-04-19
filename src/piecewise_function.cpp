#include <gsplines++/piecewise_function.hpp>

namespace gsplines {

PiecewiseFunction PiecewiseFunction::deriv(unsigned int std::size_t _deg) {
  return *this;
}

PiecewiseFunction::PiecewiseFunction(const PiecewiseFunction &that) {}

PiecewiseFunction::PiecewiseFunction(std::size_t _codom_dim,
                                     std::size_t _n_intervals,
                                     basis::Basis &_basis,
                                     Eigen::VectorXd &_coefficents,
                                     Eigen::VectorXd &_tauv)
    : coefficients_(_coefficents), codom_dim_(_codom_dim),
      number_of_intervals_(_n_intervals), coefficients_(_coefficents),
      domain_break_points_(_n_intervals + 1), domain_interval_lengths_(_tauv),
      basis_(_basis) {
  time_instant = 0.0;
  domain_break_points_(0) = time_instant;
  for (int i = 1; i < number_of_intervals_ + 1; i++) {
    domain_break_points_(i) = time_instant;
    time_instant += domain_interval_lengths_(i - 1);
  }
}

PiecewiseFunction::~PiecewiseFunction() {}

void PiecewiseFunction::operator()(double _domain_point,
                                   Eigen::VectorXd &_result) const {}
void PiecewiseFunction::operator()(const Eigen::VectorXd &_domain_points,
                                   Eigen::MatrixXd &_result) const {
  std::size_t result_size(_domain_points.size());
  std::size_t current_interval = 0;
  for (int i = 0; i < result_size; i++) {
    current_interval = get_interval(_domain_points(i));
    basis_.eval_on_window(interval_to_window(_domain_points(i)),
                          domain_interval_lengths_(current_interval),
                          basis_buffer_);
    for (int j = 0; j < codom_dim_; j++) {
      _result(i, j) =
          coefficient_segment(current_interval, j).adjoint() * basis_buffer_;
    }
  }
  Eigen::MatrixXd PiecewiseFunction::operator()(
      const Eigen::VectorXd &_domain_points) const {
    Eigen::MatrixXd result(_domain_point.size(), codom_dim_);
    return MatrixXd(1, 1);
  }

  Eigen::VectorXd PiecewiseFunction::operator()(double _domain_point) const {
    return MatrixXd(1, 1);
  }

  std::size_t PiecewiseFunction::get_interval(double _domain_point) {
    if (_domain_point < domain_break_points_(0))
      return 0.0;
    for (int i = 0; i < number_of_intervals_; i++) {
      if (domain_break_points_(i) < _domain_point and
          _domain_point < domain_break_points_(i + 1))
        return i;
    }
    if (domain_break_points_.tail(1) < _domain_point)
      return number_of_intervals_ - 1;
  }

  Eigen::VectorXd PiecewiseFunction::coefficient_segment(
      std::size_t _interval, std::size_t _component) {
    return coefficients_.segment(i, basis_.get_dim());
  }

  double PiecewiseFunction::interval_to_window(double _domaint_point,
                                               std::size_t _interval) {
    return 2.0 * (_domain_point[i] - domain_break_points_[_interval]) /
               domain_interval_lengths_[_interval] -
           1.0
  }
} // namespace gsplines
