#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <iostream>

namespace gsplines {

GSpline *GSpline::deriv_impl(std::size_t _deg) const {

  Eigen::VectorXd result_coeff(coefficients_);
  int der_coor;
  int interval_coor;
  int codom_coor;
  int i0;

  if (_deg > 0) {
    for (interval_coor = 0; interval_coor < number_of_intervals_;
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

  return new GSpline(get_domain(), get_codom_dim(), number_of_intervals_,
                     *basis_, result_coeff, domain_interval_lengths_,
                     get_name());
}

GSpline::GSpline(const GSpline &that)
    : FunctionInheritanceHelper(that),
      number_of_intervals_(that.number_of_intervals_),
      basis_(that.basis_->clone()), coefficients_(that.coefficients_),
      domain_break_points_(that.domain_break_points_),
      domain_interval_lengths_(that.domain_interval_lengths_),
      waypoints_(that.waypoints_), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      number_of_intervals_ * basis_->get_dim() * get_codom_dim()) {
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
      number_of_intervals_(that.number_of_intervals_),
      basis_(that.basis_->move_clone()),
      coefficients_(std::move(that.coefficients_)),
      domain_break_points_(std::move(that.domain_break_points_)),
      domain_interval_lengths_(std::move(that.domain_interval_lengths_)),
      waypoints_(std::move(that.waypoints_)), basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() !=
      number_of_intervals_ * basis_->get_dim() * get_codom_dim()) {
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
      number_of_intervals_(_n_intervals), basis_(_basis.clone()),
      coefficients_(_coefficents), domain_break_points_(_n_intervals + 1),
      domain_interval_lengths_(_tauv), waypoints_(_n_intervals + 1, _codom_dim),
      basis_buffer_(basis_->get_dim()) {

  if (coefficients_.size() != _n_intervals * basis_->get_dim() * _codom_dim) {
    throw std::invalid_argument(
        "GSpline instantation Error: The number of coefficients is incorrect "
        "req base dim: " +
        std::to_string(basis_->get_dim()) +
        " req codom dim: " + std::to_string(get_codom_dim()) +
        " req num inter: " + std::to_string(get_number_of_intervals()) +
        ". However, the number of coeff was " +
        std::to_string(coefficients_.size()));
  }
  double time_instant = 0.0;
  domain_break_points_(0) = time_instant;
  for (int i = 1; i < number_of_intervals_ + 1; i++) {
    time_instant += domain_interval_lengths_(i - 1);
    domain_break_points_(i) = time_instant;
  }

  GSpline::value(domain_break_points_, waypoints_);
}

void GSpline::value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                         Eigen::Ref<Eigen::MatrixXd> _result) const {

  std::size_t result_cols(_domain_points.size());
  std::size_t current_interval = 0;
  double s, tau;

  int i, j;

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

  return inter.interpolate(new_domain_interva_lengths, waypoints_);
}

std::size_t GSpline::get_interval(double _domain_point) const {
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
  return 2.0 * (_domain_point - domain_break_points_[_interval]) /
             domain_interval_lengths_[_interval] -
         1.0;
}

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
