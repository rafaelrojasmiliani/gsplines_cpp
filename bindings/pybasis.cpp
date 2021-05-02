#include <eigen3/Eigen/Core>
#include <gsplines++/Basis.hpp>
#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/Interpolator.hpp>
#include <gsplines++/PiecewiseFunction.hpp>
#include <gsplines++/Sobolev.hpp>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace gsplines {
namespace basis {

class PyBasis : public Basis {
public:
  /* Inherit the constructors */
  using Basis::Basis;
  void eval_on_window(double _s, double _tau,
                      Eigen::Ref<Eigen::VectorXd> _buff) override {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_on_window, _s, _tau, _buff);
  }
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 Eigen::Ref<Eigen::VectorXd> _buff) override {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_on_window, _s, _tau,
                           _deg, _buff);
  }

  void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd> _buff) override {

    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_wrt_tau_on_window, _s,
                           _tau, _deg, _buff);
  }
  void add_derivative_matrix(double _tau, std::size_t _deg,
                             Eigen::Ref<Eigen::MatrixXd> _mat) override {

    PYBIND11_OVERRIDE_PURE(void, Basis, add_derivative_matrix, _tau, _deg,
                           _mat);
  }
  void add_derivative_matrix_deriv_wrt_tau(
      double _tau, std::size_t _deg,
      Eigen::Ref<Eigen::MatrixXd> _mat) override {

    PYBIND11_OVERRIDE_PURE(void, Basis, add_derivative_matrix_deriv_wrt_tau,
                           _tau, _deg, _mat);
  }
};

class PyInterpolator : public Interpolator {
public:
  using Interpolator::Interpolator;

  PiecewiseFunction
  py_interpolate(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
                 const py::EigenDRef<const Eigen::MatrixXd> _waypoints) {
    return interpolate(_interval_lengths, _waypoints);
  }
  const Eigen::Ref<const Eigen::VectorXd> py_solve_interpolation(
      const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
      const py::EigenDRef<const Eigen::MatrixXd> _waypoints) {
    return solve_interpolation(_interval_lengths, _waypoints);
  }
};

PYBIND11_MODULE(pygsplines, m) {
  py::class_<Basis, PyBasis>(m, "Basis")
      .def(py::init<std::size_t>())
      .def("get_dim", &Basis::get_dim)
      .def("eval_on_window", &Basis::eval_on_window)
      .def("eval_derivative_on_window", &Basis::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &Basis::eval_derivative_wrt_tau_on_window);

  py::class_<BasisLegendre, Basis>(m, "BasisLegendre")
      .def(py::init<std::size_t>())
      .def("eval_on_window", &BasisLegendre::eval_on_window)
      .def("eval_derivative_on_window",
           &BasisLegendre::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &BasisLegendre::eval_derivative_wrt_tau_on_window);

  py::class_<PiecewiseFunction>(m, "PiecewiseFunction")
      .def(py::init<std::size_t, std::size_t, Basis &,
                    const Eigen::Ref<const Eigen::VectorXd>,
                    const Eigen::Ref<const Eigen::VectorXd>>())
      .def("__call__", &PiecewiseFunction::operator())
      .def("get_exec_time", &PiecewiseFunction::get_exec_time)
      .def("get_codom_dim", &PiecewiseFunction::get_codom_dim)
      .def("get_domain_breakpoints", &PiecewiseFunction::get_domain_breakpoints)
      .def("get_coeff", &PiecewiseFunction::get_coeff)
      .def("deriv", &PiecewiseFunction::deriv);

  py::class_<Interpolator>(m, "Interpolator")
      .def(py::init<std::size_t, std::size_t, Basis &>());

  py::class_<PyInterpolator, Interpolator>(m, "PyInterpolator")
      .def(py::init<std::size_t, std::size_t, Basis &>())
      .def("interpolate", &PyInterpolator::py_interpolate)
      .def("solve_interpolation", &PyInterpolator::py_solve_interpolation)
      .def("print_interpolating_matrix",
           &PyInterpolator::print_interpolating_matrix)
      .def("get_coeff_derivative_wrt_tau",
           &PyInterpolator::get_coeff_derivative_wrt_tau)
      .def("print_interpolating_vector",
           &PyInterpolator::print_interpolating_vector);

  py::class_<SobolevNorm>(m, "SobolevNorm")
      .def(py::init<const Eigen::Ref<const Eigen::MatrixXd>, basis::Basis &,
                    std::vector<std::pair<std::size_t, double>>>())
      .def("__call__", &SobolevNorm::operator())
      .def("deriv_wrt_interval_len", &SobolevNorm::deriv_wrt_interval_len);
}

} // namespace basis
} // namespace gsplines
