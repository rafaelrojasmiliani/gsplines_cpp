#include <gsplines++/Basis.hpp>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace gsplines {
namespace basis {

class PyBasis : public Basis {
public:
  /* Inherit the constructors */
  using Basis::Basis;
  void eval_on_window(double _s, double _tau, double _buff[]) override {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_on_window, _s, _tau, _buff);
  }
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 double _buff[]) override {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_on_window, _s, _tau,
                           _deg, _buff);
  }

  void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                         unsigned int _deg,
                                         double _buff[]) override {

    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_wrt_tau_on_window, _s,
                           _tau, _deg, _buff);
  }
};

namespace py = pybind11;
PYBIND11_MODULE(pygsplines, m) {
  py::class_<Basis, PyBasis>(m, "Basis")
      .def(py::init<std::size_t>())
      .def("get_dim", &Basis::get_dim)
      .def("eval_on_window", &Basis::eval_on_window)
      .def("eval_derivative_on_window", &Basis::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &Basis::eval_derivative_wrt_tau_on_window);
}

} // namespace basis
} // namespace gsplines
