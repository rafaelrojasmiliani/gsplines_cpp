#include <eigen3/Eigen/Core>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace gsplines {
namespace functions {

PYBIND11_MODULE(pygsplines, m) {
  // Function
  py::class_<FunctionExpression>(m, "FunctionExpression")
      .def(py::init<std::pair<double, double>, std::size_t>())
      .def("get_codom_dim", &FunctionExpression::get_codom_dim)
      .def("is_point_in_domain", &FunctionExpression::is_point_in_domain)
      .def("get_domain", &FunctionExpression::get_domain)
      .def("deriv", &FunctionExpression::derivate)
      .def("__call__", &FunctionExpression::operator());
}
} // namespace functions
} // namespace gsplines
