#include <gsplines++/BasisLegendre.hpp>
#include <math.h>
namespace gsplines {

namespace basis {

void gsplines_legendre_dmat(size_t _dim, Eigen::MatrixXd &_dmat);
BasisLegendre::BasisLegendre(std::size_t _dim)
    : Basis(_dim), dmat_(_dim, _dim) {

  gsplines_legendre_dmat(_dim, dmat_);
}

BasisLegendre::~BasisLegendre() {}

void BasisLegendre::eval_derivative_on_window(double _s, double _tau,
                                              unsigned int _deg,
                                              double _buff[]) {
  double buff_next[20];
  double term = 0;
  double aux = 0;
  _buff[0] = 1.0;
  _buff[1] = _s;
  for (int i = 1; i < get_dim() - 1; i++) {
    _buff[i + 1] =
        1.0 / ((double)i + 1.0) *
        ((2.0 * (double)i + 1.0) * _s * _buff[i] - (double)i * _buff[i - 1]);
  }

  buff_next[0] = 0;
  buff_next[1] = 1.0;
  for (int i = 1; i < get_dim() - 1; i++) {
    term = _buff[i] + _s * buff_next[i];
    buff_next[i + 1] =
        1.0 / ((double)i + 1.0) *
        ((2.0 * (double)i + 1.0) * term - (double)i * buff_next[i - 1]);
  }
  memcpy(_buff, buff_next, get_dim() * sizeof(double));
  buff_next[1] = 0.0;
  buff_next[2] = 3.0;
  for (int d = 2; d <= _deg; d++) {
    for (int i = d; i < get_dim() - 1; i++) {
      term = (1.0 + _s) * buff_next[i] + (d - 1.0) * _buff[i];
      buff_next[i + 1] =
          1.0 / ((double)i + 1.0) *
          ((2.0 * (double)i + 1.0) * term - (double)i * buff_next[i - 1]);
    }
    term = (2.0) * buff_next[d + 1] + (d - 1) * _buff[d];
    memcpy(_buff, buff_next, get_dim() * sizeof(double));
    buff_next[d] = 0.0;
    buff_next[d + 1] =
        1.0 / ((double)d + 1.0) * ((2.0 * (double)d + 1.0) * term);
  }
}

void BasisLegendre::eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                                      unsigned int _deg,
                                                      double _buff[]) {}

void BasisLegendre::eval_on_window(double _s, double _tau, double _buff[]) {

  _buff[0] = 1.0;
  _buff[1] = _s;
  for (int i = 1; i < get_dim() - 1; i++) {
    _buff[i + 1] =
        1.0 / ((double)i + 1.0) *
        ((2.0 * (double)i + 1.0) * _s * _buff[i] - (double)i * _buff[i - 1]);
  }
}
double alpha(int i) { return (double)(i + 1.0) / (double)(2.0 * i + 1.0); }

double gamma(int i) { return (double)i / (double)(2.0 * i + 1.0); }

void gsplines_legendre_dmat(size_t _dim, Eigen::MatrixXd &_dmat) {

  double firstTerm, secondTerm, thirdTerm, fourthTerm;

  for (int i = 0; i < _dim; i++) {
    for (int j = 0; j < _dim - 1; j++) {
      if (i == j + 1) {
        _dmat(i, j) = (i) / alpha(i - 1);
      } else if (i > j + 1) {
        firstTerm = 0.0;
        secondTerm = (j >= 1) ? alpha(j - 1) * _dmat(i - 1, j - 1) : 0.0;
        thirdTerm = (i >= 1) ? gamma(j + 1) * _dmat(i - 1, j + 1) : 0.0;
        fourthTerm = (i >= 2) ? gamma(i - 1) * _dmat(i - 2, j) : 0.0;

        _dmat(i, j) = 1.0f / alpha(i - 1) *
                      (firstTerm + secondTerm + thirdTerm - fourthTerm);
      }
    }
  }
}
} // namespace basis
} // namespace gsplines
