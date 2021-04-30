#include <gsplines++/BasisLegendre.hpp>
#include <iostream>
#include <math.h>
namespace gsplines {

namespace basis {

void gsplines_legendre_dmat(size_t _dim, Eigen::MatrixXd &_dmat);
BasisLegendre::BasisLegendre(std::size_t _dim) : Basis(_dim) {
  gsplines_legendre_dmat(_dim, derivative_matrix_);
}

BasisLegendre::BasisLegendre(const BasisLegendre &that) : Basis(that) {}
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
                                   Eigen::Ref<Eigen::VectorXd> _buff) {
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
} // namespace basis
} // namespace gsplines
