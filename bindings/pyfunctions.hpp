
#ifndef PYFUNCTION
#define PYFUNCTION
#include <eigen3/Eigen/Core>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <gsplines++/Functions/Function.hpp>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace gsplines {
namespace functions {} // namespace functions
} // namespace gsplines
#endif /* ifndef PYFUNCTIOND */
