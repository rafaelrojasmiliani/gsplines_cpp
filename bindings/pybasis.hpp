
#ifndef PYBASIS
#define PYBASIS
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/FunctionalAnalysis/Sobolev.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
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
                      Eigen::Ref<Eigen::VectorXd> _buff) const override {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_on_window, _s, _tau, _buff);
  }
  void
  eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                            Eigen::Ref<Eigen::VectorXd> _buff) const override {
    PYBIND11_OVERRIDE_PURE(void, Basis, eval_derivative_on_window, _s, _tau,
                           _deg, _buff);
  }

  void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd> _buff) const override {

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
  Eigen::MatrixXd derivative_matrix_impl(std::size_t _deg) const override {

    PYBIND11_OVERRIDE_PURE(Eigen::MatrixXd, Basis, derivative_matrix_impl,
                           _deg);
  }
  std::unique_ptr<Basis> clone() const override { return nullptr; }
  std::unique_ptr<Basis> move_clone() override { return nullptr; }
};

} // namespace basis

class PyInterpolator : public Interpolator {
public:
  using Interpolator::Interpolator;

  GSpline
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

class PySobolevNorm : public functional_analysis::SobolevNorm {
public:
  PySobolevNorm(const py::EigenDRef<const Eigen::MatrixXd> _waypoints,
                basis::Basis &_basis,
                std::vector<std::pair<std::size_t, double>> _weights)
      : SobolevNorm(_waypoints, _basis, _weights) {}
};
} // namespace gsplines
#endif /* ifndef PYBASIS */
