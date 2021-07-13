#include <gsplines++/BasisLegendre.hpp>
#include <iostream>
#include <math.h>
namespace gsplines {

namespace basis {

void gsplines_legendre_dmat(size_t _dim, Eigen::MatrixXd &_dmat);
BasisLegendre::BasisLegendre(std::size_t _dim) : Basis(_dim) {
  gsplines_legendre_dmat(_dim, derivative_matrix_);

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
  }
}

BasisLegendre::BasisLegendre(const BasisLegendre &that)
    : Basis(that),
      derivative_matrices_buffer_(that.derivative_matrices_buffer_) {}
BasisLegendre::~BasisLegendre() {}

void BasisLegendre::eval_derivative_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd> _buff) {
  Eigen::VectorXd buff_next(get_dim());
  double term = 0;
  double aux = 0;
  double mutiplier = 1.0;
  _buff(0) = 1.0;
  _buff(1) = _s;
  for (unsigned int i = 1; i < get_dim() - 1; i++) {
    _buff(i + 1) =
        1.0 / ((double)i + 1.0) *
        ((2.0 * (double)i + 1.0) * _s * _buff(i) - (double)i * _buff(i - 1));
  }
  aux = 1.0;
  for (int d = 1; d <= _deg; d++) {
    buff_next(d - 1) = 0.0;
    buff_next(d) = aux;

    for (int i = d; i < get_dim() - 1; i++) {
      term =
          (2.0 * (double)i + 1.0) * ((double)d * _buff(i) + _s * buff_next(i));
      buff_next(i + 1) =
          1.0 / ((double)i + 1.0) * (term - (double)i * buff_next(i - 1));
    }
    aux = (2.0 * (double)d + 1.0) / ((double)d + 1.0) * ((double)d + 1.0) * aux;
    _buff = buff_next;
    mutiplier *= (2.0 / _tau);
  }
  _buff *= mutiplier;
}

void BasisLegendre::eval_derivative_wrt_tau_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd> _buff) {

  eval_derivative_on_window(_s, _tau, _deg, _buff);
  _buff *= -0.5 * _deg * (2.0 / _tau);
}

void BasisLegendre::eval_on_window(double _s, double _tau,
                                   Eigen::Ref<Eigen::VectorXd> _buff) const {
  _buff(0) = 1.0;
  _buff(1) = _s;
  for (int i = 1; i < get_dim() - 1; i++) {
    _buff(i + 1) =
        1.0 / ((double)i + 1.0) *
        ((2.0 * (double)i + 1.0) * _s * _buff(i) - (double)i * _buff(i - 1));
  }
}
double alpha(int i) { return ((double)(i + 1.0)) / ((double)(2.0 * i + 1.0)); }

double gamma(int i) { return ((double)i) / ((double)(2.0 * i + 1.0)); }

void gsplines_legendre_dmat(size_t _dim, Eigen::MatrixXd &_dmat) {

  _dmat = Eigen::MatrixXd::Zero(_dim, _dim);
  double firstTerm, secondTerm, thirdTerm, fourthTerm;

  for (int i = 1; i < _dim; i++) {       // for on i
    for (int j = 0; j < _dim - 1; j++) { // for on j
      if (i == j + 1) {
        _dmat(i, j) = ((double)i) / alpha(i - 1);
      } else if (i > j + 1) {
        firstTerm = 0.0;
        secondTerm = (j >= 1) ? alpha(j - 1) * _dmat(i - 1, j - 1) : 0.0;
        thirdTerm = (i >= 1) ? gamma(j + 1) * _dmat(i - 1, j + 1) : 0.0;
        fourthTerm = (i >= 2) ? gamma(i - 1) * _dmat(i - 2, j) : 0.0;

        _dmat(i, j) = 1.0 / alpha(i - 1) *
                      (firstTerm + secondTerm + thirdTerm - fourthTerm);
      } else {
        _dmat(i, j) = 0.0;
      }
    } // for on j
  }   // for on i
}

void BasisLegendre::add_derivative_matrix(double tau, std::size_t _deg,
                                          Eigen::Ref<Eigen::MatrixXd> _mat) {

  unsigned int i, j;
  double scale = _deg > 0 ? pow(2.0 / tau, 2 * _deg - 1) : tau / 2.0;
  if (_deg < get_dim() + 1)
    _mat.noalias() += derivative_matrices_buffer_[_deg] * scale;
}

void BasisLegendre::add_derivative_matrix_deriv_wrt_tau(
    double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) {
  unsigned int i, j;
  double scale = _deg > 0 ? -0.5 * (2.0 * _deg - 1.0) * pow(2.0 / tau, 2 * _deg)
                          : 1.0 / 2.0;
  if (_deg < get_dim() + 1)
    _mat.noalias() += derivative_matrices_buffer_[_deg] * scale;
}
} // namespace basis
} // namespace gsplines
