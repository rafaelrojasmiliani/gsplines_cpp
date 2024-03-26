
#ifndef PYFUNCTION
#define PYFUNCTION
#include <eigen3/Eigen/Core>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionBase.hpp>
#include <gsplines/GSpline.hpp>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <cstddef>

namespace gsplines {
namespace functions {

class PyFunctionBase : public FunctionBase {
 public:
  using FunctionBase::FunctionBase;
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override {
    PYBIND11_OVERRIDE_PURE(void, FunctionBase, value_impl, _domain_points,
                           _result);
  }

  [[nodiscard]] FunctionBase* clone_impl() const override {
    PYBIND11_OVERRIDE_PURE(FunctionBase*, FunctionBase, clone_impl);
  }

  FunctionBase* move_clone_impl() override {
    PYBIND11_OVERRIDE_PURE(FunctionBase*, FunctionBase, move_clone_impl);
  };

  [[nodiscard]] FunctionBase* deriv_impl(std::size_t /*_deg*/) const override {
    PYBIND11_OVERRIDE_PURE(FunctionBase*, FunctionBase, deriv_impl);
  };

  ~PyFunctionBase() override = default;
};
class PyFunction : public Function {
 public:
  using Function::Function;
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override {
    PYBIND11_OVERRIDE_PURE(void, Function, value_impl, _domain_points, _result);
  }

  [[nodiscard]] FunctionBase* clone_impl() const override {
    PYBIND11_OVERRIDE_PURE(FunctionBase*, Function, clone_impl);
  }

  FunctionBase* move_clone_impl() override {
    PYBIND11_OVERRIDE_PURE(FunctionBase*, Function, move_clone_impl);
  };

  [[nodiscard]] FunctionBase* deriv_impl(std::size_t /*_deg*/) const override {
    PYBIND11_OVERRIDE_PURE(FunctionBase*, Function, deriv_impl);
  };

  ~PyFunction() override = default;
};

class PyGSplineBase : public GSplineBase {
 public:
  using GSplineBase::GSplineBase;
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override {
    PYBIND11_OVERRIDE_PURE(void, GSplineBase, value_impl, _domain_points,
                           _result);
  }

  GSplineBase* deriv_impl(std::size_t /*_deg*/) const override {
    PYBIND11_OVERRIDE_PURE(GSplineBase*, GSplineBase, deriv_impl);
  };

  ~PyGSplineBase() override = default;
};

}  // namespace functions
}  // namespace gsplines
#endif /* ifndef PYFUNCTIOND */
