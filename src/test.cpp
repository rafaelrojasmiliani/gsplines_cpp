#include "basis.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_bessel.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;
class BasisLagrange : public Basis {
private:
  double *gen_points_;
  double **dd_coeff_;
  double *deriv_mat_;

public:
  BasisLagrange(size_t _dim, double _gen_points[]) : Basis(_dim) {
    gen_points_ = (double *)calloc(_dim, sizeof(double));
    memcpy((void *)gen_points_, (const void *)_gen_points,
           sizeof(double) * _dim);
    double *codom_buff = (double *)calloc(_dim, sizeof(double));
    dd_coeff_ = (double **)calloc(_dim, sizeof(double *));
    for (int i = 0; i < _dim; i++) {
      *(codom_buff + i) = 1.0;
      *(dd_coeff_ + i) = (double *)calloc(_dim, sizeof(double));
      gsl_poly_dd_init(*(dd_coeff_ + i), _gen_points, codom_buff, _dim);
      *(codom_buff + i) = 0.0;
    }
    free(codom_buff);
  }
  ~BasisLagrange() {

    for (int i = 0; i < get_dim(); i++)
      free(*(dd_coeff_ + i));

    free(dd_coeff_);
    free(gen_points_);
  }
  void eval_on_window(double _s, double _tau, double _buff[]) {
    for (int i = 0; i < get_dim(); i++) {
      _buff[i] = gsl_poly_dd_eval(*(dd_coeff_ + i), gen_points_, get_dim(), _s);
    }
  }
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 double _buff[]) {}

  void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                         unsigned int _deg, double _buff[]) {}
};
int main(void) {
  int i = 0;
  int a;
  double k;
  gsl_integration_fixed_workspace *ws = gsl_integration_fixed_alloc(
      gsl_integration_fixed_legendre, 5, -1.0, 1.0, 0, 0);
  a = i + 100;
  printf("----\n");
  BasisLagrange test(5, ws->x);
  printf("----\n");
  gsl_integration_fixed_free(ws);
  return 0;
}
