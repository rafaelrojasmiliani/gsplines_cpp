#include <gsplines++/BasisLagrange.hpp>
#include <iostream>
#include <math.h>
namespace gsplines {

namespace basis {

void gsplines_lagrange_dmat(size_t _dim, Eigen::MatrixXd &_dmat);

BasisLagrange::BasisLagrange(Eigen::Ref<const Eigen::VectorXd> _domain_points)
    : Basis(_domain_points.size()), domain_points_(_domain_points) {
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
    : Basis(that), domain_points_(that.domain_points_) {}

BasisLagrange::BasisLagrange(BasisLagrange &&that)
    : Basis(std::move(that)), domain_points_(std::move(that.domain_points_)) {}

BasisLagrange::~BasisLagrange() {}

void BasisLagrange::eval_derivative_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd> _buff) const {}

void BasisLagrange::eval_derivative_wrt_tau_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd> _buff) const {}

void BasisLagrange::eval_on_window(double _s, double _tau,
                                   Eigen::Ref<Eigen::VectorXd> _buff) const {}

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
    for (std::size_t uick = 0; uick < _points.size(); uick++) {
      result(uick) = result(uick) * (_points(uick) - _points(uicj));
      result(uicj) = result(uicj) * (_points(uicj) - _points(uick));
    }
  }

  result = (1.0 / result.array()).matrix();

  return std::move(result);
}

Eigen::MatrixXd
BasisLagrange::derivative_matrix(Eigen::Ref<const Eigen::VectorXd> _points) {

  Eigen::MatrixXd result(Eigen::MatrixXd::Zero(_points.size(), _points.size()));

  Eigen::VectorXd bw = barycentric_weights(_points);

  for (std::size_t uici = 0; uici < _points.size(); uici++) {
    result(uici, uici) = 0;
    for (std::size_t uicj = 0; uicj < _points.size(); uicj++)
      if (uici != uicj) {
        result(uici, uicj) =
            bw(uicj) / bw(uici) * 1.0 / (_points(uici) - _points(uicj));
        result(uici, uici) += -result(uicj, uici);
      }
  }

  return std::move(result);
}
} // namespace basis
} // namespace gsplines
