
#ifndef PYFUNCTION
#define PYFUNCTION
#include <eigen3/Eigen/Core>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <gsplines++/Functions/Functions.hpp>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace gsplines {
namespace functions {

class PyFunction : public Function {
public:
  /* Inherit the constructors */
  using Function::Function;
  virtual std::unique_ptr<Function> deriv(int _deg = 1) { return nullptr; }
  virtual std::unique_ptr<Function> clone() const { return nullptr; }
  // https://github.com/pybind/pybind11/issues/673#issuecomment-280755589
  Function *deriv_wrapper(int _deg) {
    PYBIND11_OVERRIDE_PURE(Function *, Function, deriv, _deg);
  }

  Function *clone_wrapper(int _deg) {
    PYBIND11_OVERRIDE_PURE(Function *, Function, clone, _deg);
  }

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

    PYBIND11_OVERRIDE_PURE(Eigen::MatrixXd, Function, operator(),
                           _domain_points);
  }
};

} // namespace functions
} // namespace gsplines
#endif /* ifndef PYFUNCTIOND */
