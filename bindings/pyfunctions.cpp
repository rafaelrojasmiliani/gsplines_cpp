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
}; // namespace functions
class PyFunctionExpression : public FunctionExpression {
public:
  /* Inherit the constructors */
  using FunctionExpression::FunctionExpression;

  PyFunctionExpression sum(const FunctionExpression &that) {
    return *this + that;
  }
}; // namespace functions
PYBIND11_MODULE(pygsplines, m) {
  // Function
  py::class_<PyFunction>(m, "Function")
      .def(py::init<std::pair<double, double>, std::size_t>())
      .def("get_codom_dim", &Function::get_codom_dim)
      .def("is_point_in_domain", &Function::is_point_in_domain)
      .def("get_domain", &Function::get_domain)
      .def("deriv", &PyFunction::deriv_wrapper)
      .def("__call__", &PyFunction::operator());

  // Function Expression
  py::class_<FunctionExpression>(m, "FunctionExpression")
      .def(py::init<std::pair<double, double>, std::size_t,
                    FunctionExpression::Type>())
      .def("get_codom_dim", &FunctionExpression::get_codom_dim)
      .def("is_point_in_domain", &FunctionExpression::is_point_in_domain)
      .def("get_domain", &FunctionExpression::get_domain)
      .def("deriv", &FunctionExpression::derivate)
      .def("__call__", &FunctionExpression::operator());
  // Operations
  // m.def("optimal_sobolev_norm", &gsplines_opt::optimal_sobolev_norm);
}
} // namespace functions
} // namespace gsplines
