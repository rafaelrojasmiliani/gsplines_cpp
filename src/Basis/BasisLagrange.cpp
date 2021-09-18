#include <gsplines/Basis/BasisLagrange.hpp>
#include <iostream>
#include <math.h>
namespace gsplines {

namespace basis {

bool almost_equal(double _a, double _b, double _epsilon) {

  if (std::fabs(_a) < _epsilon or std::fabs(_b) < _epsilon) {
    if (std::fabs(_a - _b) < 2 * _epsilon)
      return true;
    return false;
  }
  if (std::fabs(_a - _b) < _epsilon * std::fabs(_a) and
      std::fabs(_a - _b) < _epsilon * std::fabs(_b))
    return true;
  return false;
}

void gsplines_lagrange_dmat(size_t _dim, Eigen::MatrixXd &_dmat);

BasisLagrange::BasisLagrange(Eigen::Ref<const Eigen::VectorXd> _domain_points)
    : Basis(_domain_points.size(),
            "lagrange_" + std::to_string(_domain_points.size())),
      domain_points_(_domain_points),
      barycentric_weights_(barycentric_weights(domain_points_)) {

  derivative_matrix_ = derivative_matrix(domain_points_);
  /*
    derivative_matrix(get_dim(), derivative_matrix_);

    Eigen::MatrixXd dmat(derivative_matrix_);

    Eigen::MatrixXd l2normbase_matrix(get_dim(), get_dim());
    l2normbase_matrix.setZero();
    for (int i = 0; i < _dim; i++) {
      l2normbase_matrix(i, i) = 2.0 / (2.0 * i + 1.0);
    }

    derivative_matrices_buffer_.push_back(l2normbase_matrix);
    for (int i = 1; i < _dim + 1; i++) {
      derivative_matrices_buffer_.push_back(dmat * l2normbase_matrix *
                                            dmat.transpose());
      dmat *= derivative_matrix_;
    }*/
}

BasisLagrange::BasisLagrange(const BasisLagrange &that)
    : Basis(that), domain_points_(that.domain_points_),
      barycentric_weights_(that.barycentric_weights_) {}

BasisLagrange::BasisLagrange(BasisLagrange &&that)
    : Basis(std::move(that)), domain_points_(std::move(that.domain_points_)),
      barycentric_weights_(std::move(that.barycentric_weights_)) {}

void BasisLagrange::eval_derivative_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd> _buff) const {

  eval_on_window(_s, _tau, _buff);
  for (std::size_t i = 1; i <= _deg; i++)
    _buff = derivative_matrix_.transpose() * _buff * (2.0 / _tau);
}

void BasisLagrange::eval_derivative_wrt_tau_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd> _buff) const {

  eval_derivative_on_window(_s, _tau, _deg, _buff);
  _buff *= -0.5 * _deg * (2.0 / _tau);
}

void BasisLagrange::eval_on_window(double _s, double _tau,
                                   Eigen::Ref<Eigen::VectorXd> _buff) const {

  /*
   *  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 34: LagrangeInterpolatingPolynomials
   * */
  _buff.setConstant(0.0);
  for (std::size_t j = 0; j < get_dim(); j++) {
    if (almost_equal(_s, domain_points_(j), 1.0e-9)) {
      _buff(j) = 1.0;
      return;
    }
  }
  double t_book, s_book;
  s_book = 0.0;
  for (std::size_t j = 0; j < get_dim(); j++) {
    t_book = barycentric_weights_(j) / (_s - domain_points_(j));
    _buff(j) = t_book;
    s_book = s_book - t_book;
  }
  _buff /= -s_book;
}

void BasisLagrange::add_derivative_matrix(double tau, std::size_t _deg,
                                          Eigen::Ref<Eigen::MatrixXd> _mat) {}

void BasisLagrange::add_derivative_matrix_deriv_wrt_tau(
    double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) {}

Eigen::VectorXd
BasisLagrange::barycentric_weights(Eigen::Ref<const Eigen::VectorXd> _points) {
  /*
   *  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 30: BarycentricWeights: Weights for Lagrange Interpolation
   * */
  Eigen::VectorXd result(Eigen::VectorXd::Ones(_points.size()));

  for (std::size_t uicj = 1; uicj < _points.size(); uicj++) {
    for (std::size_t uick = 0; uick < uicj; uick++) {
      result(uick) = result(uick) * (_points(uick) - _points(uicj));
      result(uicj) = result(uicj) * (_points(uicj) - _points(uick));
    }
  }

  result = (1.0 / result.array()).matrix();

  return std::move(result);
}

Eigen::MatrixXd
BasisLagrange::derivative_matrix(Eigen::Ref<const Eigen::VectorXd> _points) {
  /*
   *  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 37: PolynomialDerivativeMatrix: First Derivative Approximation
   * */

  Eigen::MatrixXd result(Eigen::MatrixXd::Zero(_points.size(), _points.size()));

  Eigen::VectorXd bw = barycentric_weights(_points);

  for (std::size_t uici = 0; uici < _points.size(); uici++) {
    result(uici, uici) = 0;
    for (std::size_t uicj = 0; uicj < _points.size(); uicj++)
      if (uici != uicj) {
        result(uici, uicj) =
            bw(uicj) / bw(uici) * 1.0 / (_points(uici) - _points(uicj));
        result(uici, uici) += -result(uici, uicj);
      }
  }

  return std::move(result);
}

Eigen::MatrixXd
BasisLagrange::derivative_matrix(Eigen::Ref<const Eigen::VectorXd> _points,
                                 std::size_t _deg) {
  /*
   *  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 38: mthOrderPolynomialDerivativeMatrix
   * */

  if (_deg == 0)
    return Eigen::MatrixXd::Identity(_points.size(), _points.size());
  if (_deg == 1)
    return derivative_matrix(_points);

  Eigen::MatrixXd result(derivative_matrix(_points));

  Eigen::VectorXd bw = barycentric_weights(_points);

  for (std::size_t uick = 2; uick <= _deg; uick++) {
    for (std::size_t uici = 0; uici < _points.size(); uici++) {
      result(uick, uick) = 0;
      for (std::size_t uicj = 0; uicj < _points.size(); uicj++) {
        if (uici != uicj) {
          result(uici, uicj) =
              static_cast<double>(uick) / (_points(uici) - _points(uicj)) *
              (bw(uicj) / bw(uici) * result(uici, uici) - result(uici, uicj));
          result(uici, uici) += -result(uici, uicj);
        }
      }
    }
  }

  return std::move(result);
}
} // namespace basis
} // namespace gsplines
