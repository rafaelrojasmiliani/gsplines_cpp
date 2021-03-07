#include <gsplines++/basis0010.hpp>
namespace gsplines {
namespace basis {
Basis0010::Basis0010() : Basis(6) {}

Basis0010::~Basis0010() {}

void Basis0010::eval_derivative_on_window(double _s, double _tau,
                                          unsigned int _deg, double _buff[]) {}

void Basis0010::eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                                  unsigned int _deg,
                                                  double _buff[]) {}

void Basis0010::eval_on_window(double _s, double _tau, double _buff[]) {}
} // namespace basis
} // namespace gsplines
