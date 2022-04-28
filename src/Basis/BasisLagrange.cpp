#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <iostream>
#include <math.h>
#include <memory>
namespace gsplines {

namespace basis {

std::shared_ptr<BasisLagrange>
BasisLagrange::get(const Eigen::Ref<const Eigen::VectorXd> &_domain_points) {
  return std::make_shared<BasisLagrange>(_domain_points);
}

std::shared_ptr<BasisLagrange>
BasisLagrange::get(const std::vector<double> &_domain_points) {
  return std::make_shared<BasisLagrange>(_domain_points);
}

bool almost_equal(double _a, double _b, double _epsilon) {

  if (std::fabs(_a) < _epsilon or std::fabs(_b) < _epsilon) {
    if (std::fabs(_a - _b) < _epsilon)
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
    : Basis(_domain_points.size(), "lagrange", _domain_points),
      domain_points_(_domain_points),
      barycentric_weights_(barycentric_weights(domain_points_)),
      deriv_buff_1_(derivative_matrix(domain_points_)),
      deriv_buff_2_(derivative_matrix(domain_points_)) {
  derivative_matrix_ = derivative_matrix(domain_points_);

  Eigen::MatrixXd dmat(derivative_matrix_);

  Eigen::MatrixXd points_to_glp_matrix(change_interpolation_points(
      _domain_points, collocation::legendre_gauss_lobatto_points(get_dim())));

  Eigen::MatrixXd l2normbase_matrix(get_dim(), get_dim());

  l2normbase_matrix.diagonal() =
      collocation::legendre_gauss_lobatto_weights(get_dim());

  derivative_matrices_buffer_.push_back(l2normbase_matrix);

  for (std::size_t i = 1; i < get_dim() + 1; i++) {
    derivative_matrices_buffer_.push_back(dmat * l2normbase_matrix *
                                          dmat.transpose());
    dmat = derivative_matrix(domain_points_, i + 1);
  }
  // get_derivative_matrix(get_dim());
}

BasisLagrange::BasisLagrange(const std::vector<double> &_domain_points)
    : BasisLagrange(Eigen::Map<const Eigen::VectorXd>(_domain_points.data(),
                                                      _domain_points.size())) {}

BasisLagrange::BasisLagrange(const BasisLagrange &that)
    : Basis(that), domain_points_(that.domain_points_),
      barycentric_weights_(that.barycentric_weights_),
      derivative_matrices_buffer_(that.derivative_matrices_buffer_),
      deriv_buff_1_(that.deriv_buff_1_), deriv_buff_2_(that.deriv_buff_2_) {}

BasisLagrange::BasisLagrange(BasisLagrange &&that)
    : Basis(std::move(that)), domain_points_(std::move(that.domain_points_)),
      barycentric_weights_(std::move(that.barycentric_weights_)),
      derivative_matrices_buffer_(std::move(that.derivative_matrices_buffer_)),
      deriv_buff_1_(std::move(that.deriv_buff_1_)),
      deriv_buff_2_(std::move(that.deriv_buff_2_)) {}

void BasisLagrange::eval_derivative_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const {
  eval_on_window(_s, _tau, _buff);
  if (_deg == 0)
    return;
  double term = std::pow(2.0 / _tau, _deg);
  /*

    const Eigen::VectorXd &x = domain_points_;

    for (std::size_t uick = 2; uick <= _deg; uick++) {
      term *= term;
      for (std::size_t uici = 0; uici < get_dim(); uici++) {
        deriv_buff_1_(uici, uici) = 0.0;
        for (std::size_t uicj = 0; uicj < get_dim(); uicj++) {
          if (uici != uicj) {
            double wj = barycentric_weights_(uicj);
            double wi = barycentric_weights_(uici);

            deriv_buff_1_(uici, uicj) = (((double)uick) / (x(uici) - x(uicj))) *
                                        ((wj / wi * deriv_buff_2_(uici, uici)) -
                                         deriv_buff_2_(uici, uicj));

            deriv_buff_1_(uici, uici) =
                deriv_buff_1_(uici, uici) - deriv_buff_1_(uici, uicj);
          }
        }
      }
      deriv_buff_2_.noalias() = deriv_buff_1_;
    }*/
  _buff = get_derivative_matrix_block(_deg).transpose() * _buff * term;
}

void BasisLagrange::eval_derivative_wrt_tau_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const {
  eval_derivative_on_window(_s, _tau, _deg, _buff);
  _buff *= -0.5 * _deg * (2.0 / _tau);
}

void BasisLagrange::eval_on_window(
    double _s, double /*_tau*/,
    Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const {
  /*  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 34: LagrangeInterpolatingPolynomials */
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
    s_book += t_book;
  }
  _buff /= s_book;
}

void BasisLagrange::add_derivative_matrix(double tau, std::size_t _deg,
                                          Eigen::Ref<Eigen::MatrixXd> _mat) {
  double scale = _deg > 0 ? pow(2.0 / tau, 2 * _deg - 1) : tau / 2.0;

  if (_deg < get_dim() + 1)
    _mat.noalias() += derivative_matrices_buffer_[_deg] * scale;
}

void BasisLagrange::add_derivative_matrix_deriv_wrt_tau(
    double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) {
  double scale = _deg > 0 ? -0.5 * (2.0 * _deg - 1.0) * pow(2.0 / tau, 2 * _deg)
                          : 1.0 / 2.0;
  if (_deg < get_dim() + 1)
    _mat.noalias() += derivative_matrices_buffer_[_deg] * scale;
}

Eigen::VectorXd
BasisLagrange::barycentric_weights(Eigen::Ref<const Eigen::VectorXd> _points) {
  /*  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 30: BarycentricWeights: Weights for Lagrange Interpolation*/
  Eigen::VectorXd result(Eigen::VectorXd::Ones(_points.size()));

  for (long uicj = 1; uicj < _points.size(); uicj++) {
    for (long uick = 0; uick < uicj; uick++) {
      result(uick) = result(uick) * (_points(uick) - _points(uicj));
      result(uicj) = result(uicj) * (_points(uicj) - _points(uick));
    }
  }

  result = (1.0 / result.array()).matrix();

  return result;
}

Eigen::MatrixXd
BasisLagrange::derivative_matrix(Eigen::Ref<const Eigen::VectorXd> _points) {
  /*  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 37: PolynomialDerivativeMatrix: First Derivative Approximation*/

  Eigen::MatrixXd result(Eigen::MatrixXd::Zero(_points.size(), _points.size()));

  Eigen::VectorXd bw = barycentric_weights(_points);

  for (long uici = 0; uici < _points.size(); uici++) {
    result(uici, uici) = 0;
    for (long uicj = 0; uicj < _points.size(); uicj++)
      if (uici != uicj) {
        result(uici, uicj) =
            bw(uicj) / bw(uici) * 1.0 / (_points(uici) - _points(uicj));
        result(uici, uici) += -result(uici, uicj);
      }
  }

  return result;
}

Eigen::MatrixXd
BasisLagrange::derivative_matrix(Eigen::Ref<const Eigen::VectorXd> _points,
                                 std::size_t _deg) {
  /*  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 38: mthOrderPolynomialDerivativeMatrix*/

  if (_deg == 0)
    return Eigen::MatrixXd::Identity(_points.size(), _points.size());
  if (_deg == 1)
    return derivative_matrix(_points);

  Eigen::MatrixXd result(derivative_matrix(_points));
  Eigen::MatrixXd buff(derivative_matrix(_points));

  const Eigen::VectorXd bw = barycentric_weights(_points);

  for (std::size_t uick = 2; uick <= _deg; uick++) {
    for (long uici = 0; uici < _points.size(); uici++) {
      result(uici, uici) = 0;
      for (long uicj = 0; uicj < _points.size(); uicj++) {
        if (uici != uicj) {
          result(uici, uicj) =
              static_cast<double>(uick) / (_points(uici) - _points(uicj)) *
              (bw(uicj) / bw(uici) * buff(uici, uici) - buff(uici, uicj));
          result(uici, uici) += -result(uici, uicj);
        }
      }
    }
    buff.noalias() = result;
  }

  return result;
}

Eigen::MatrixXd BasisLagrange::change_interpolation_points(
    Eigen::Ref<const Eigen::VectorXd> _old_points,
    Eigen::Ref<const Eigen::VectorXd> _new_points) {
  /*  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 32: Matrix for Interpolation Between Two Sets of Points*/

  std::size_t book_m = _new_points.size();
  std::size_t book_n = _old_points.size();
  double t_book, s_book;
  bool book_row_has_match = false;

  Eigen::MatrixXd result(book_m, book_n);
  Eigen::VectorXd bw = barycentric_weights(_old_points);

  for (std::size_t k = 0; k < book_m; k++) {
    for (std::size_t j = 0; j < book_n; j++) {
      result(k, j) = 0.0;
      if (almost_equal(_new_points(k), _old_points(j), 1.0e-9)) {
        book_row_has_match = true;
        result(k, j) = 1.0;
      }
    }
    if (not book_row_has_match) {
      s_book = 0.0;
      for (std::size_t j = 0; j < book_n; j++) {
        t_book = bw(j) / (_new_points(k) - _old_points(j));
        result(k, j) = t_book;
        s_book += t_book;
      }
      result.row(k) /= s_book;
    }
  }

  return result;
}

Eigen::MatrixXd BasisLagrange::derivative_matrix_impl(std::size_t _deg) const {
  if (_deg == 0) {
    return Eigen::MatrixXd::Identity(get_dim(), get_dim());
  }
  return derivative_matrix(domain_points_, _deg);
}

BasisLagrangeGaussLobatto::BasisLagrangeGaussLobatto(std::size_t _dim)
    : BasisLagrange(
          gsplines::collocation::legendre_gauss_lobatto_points(_dim)) {}
} // namespace basis
} // namespace gsplines
