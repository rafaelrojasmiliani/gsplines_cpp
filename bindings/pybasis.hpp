#ifndef PYBASIS_H
#define PYBASIS_H
#include <gsplines++/basis.hpp>

class PyBasis : public Basis {
public:
  /* Inherit the constructors */
  using Basis::Basis;
  void eval_on_window(double _s, double _tau, double _buff[]) {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_on_window, _s, _tau, _buff)
  }
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 double _buff[]) {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_on_window, _s, _tau,
                           _deg, _buff)
  }

  void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                         unsigned int _deg, double _buff[]) {

    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_wrt_tau_on_window, _s,
                           _tau, _deg, _buff)
  }
};

#endif /* PYBASIS_H */
