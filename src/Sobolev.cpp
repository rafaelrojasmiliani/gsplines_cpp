
#include <gsplines++/PiecewiseFunction.hpp>
#include <gsplines++/Sobolev.hpp>

namespace gsplines {
SobolevNorm::SobolevNorm(const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
                         basis::Basis &_basis,
                         std::vector<std::pair<std::size_t, double>> _weights)
    : basis_(&_basis), num_intervals_(_waypoints.rows() - 1),
      codom_dim_(_waypoints.cols()),
      interpolator_(codom_dim_, num_intervals_, _basis), weights_(_weights),
      waypoints_(_waypoints), matrix_(_basis.get_dim(), _basis.get_dim()),
      matrix_2_(_basis.get_dim(), _basis.get_dim()) {}

double SobolevNorm::
operator()(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths) {
  const Eigen::Ref<const Eigen::VectorXd> coeff =
      interpolator_.solve_interpolation(_interval_lengths, waypoints_);

  return inner_prod(_interval_lengths, coeff, coeff);
}

void SobolevNorm::deriv_wrt_interval_len(
    Eigen::Ref<Eigen::VectorXd> _interval_lengths,
    Eigen::Ref<Eigen::VectorXd> _buff) {
  /* y = coefficients
   * tau = _interval_lengths
   *
   * returns [] = y^T [ dQdtau_i y + 2 Q dy_dtau_i]*/

  unsigned int interval_coor;
  unsigned int codom_coor;
  double tau;
  const Eigen::Ref<const Eigen::VectorXd> coeff =
      interpolator_.solve_interpolation(_interval_lengths, waypoints_);
  double result = 0.0;

  for (interval_coor = 0; interval_coor < num_intervals_; interval_coor++) {
    // get the derivaties of y wrt tau_i
    const Eigen::VectorXd dy_dtau_i =
        interpolator_.get_coeff_derivative_wrt_tau(coeff, _interval_lengths,
                                                   interval_coor);
    // get the value of the current interval length
    tau = _interval_lengths(interval_coor);
    // initialize the matrix buffers and compute the de matrix values
    matrix_ = Eigen::MatrixXd::Zero(basis_->get_dim(), basis_->get_dim());
    for (std::pair<std::size_t, double> w : weights_) {
      basis_->add_derivative_matrix_deriv_wrt_tau(tau, w.first, matrix_);
      matrix_ *= w.second;
    }
    // compute y^T dQdtau_i y
    for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
      Eigen::Ref<Eigen::VectorXd> v1 =
          get_coefficient_segment(coeff, *basis_, num_intervals_, codom_dim_,
                                  interval_coor, codom_coor);
      result += v1.transpose() * matrix_ * v1;
    }
    result += inner_prod(_interval_lengths, coeff, dy_dtau_i);
    _buff[interval_coor] = result;
  }
}

double SobolevNorm::inner_prod(
    const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
    const Eigen::Ref<const Eigen::VectorXd> _v1,
    const Eigen::Ref<const Eigen::VectorXd> _v2) {
  if (_v1 != _v2) {
    return 0.0;
  }

  unsigned int interval_coor;
  unsigned int codom_coor;
  double tau;
  double result = 0.0;
  for (interval_coor = 0; interval_coor < num_intervals_; interval_coor++) {
    matrix_ = Eigen::MatrixXd::Zero(basis_->get_dim(), basis_->get_dim());
    tau = _interval_lengths(interval_coor);
    for (std::pair<std::size_t, double> w : weights_) {
      basis_->add_derivative_matrix(tau, w.first, matrix_);
      matrix_ *= w.second;
    }
    for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
      Eigen::Ref<Eigen::VectorXd> v1 = get_coefficient_segment(
          _v1, *basis_, num_intervals_, codom_dim_, interval_coor, codom_coor);
      Eigen::Ref<Eigen::VectorXd> v2 = get_coefficient_segment(
          _v2, *basis_, num_intervals_, codom_dim_, interval_coor, codom_coor);
      result += v1.transpose() * matrix_ * v2;
    }
  }
  return result;
}
} // namespace gsplines
