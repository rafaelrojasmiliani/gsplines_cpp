
#include "pybasis.hpp"
#include "pyfunctions.hpp"
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/FunctionalAnalysis/Sobolev.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Optimization/ipopt_solver.hpp>

PYBIND11_MODULE(gsplines, gsplines_module) {
  gsplines_module.doc() = "Generalized Splines Library with Optimization tools";

  py::module basis_submodule = gsplines_module.def_submodule("basis");
  py::module functions_submodule = gsplines_module.def_submodule("functions");
  py::module optimization_submodule =
      gsplines_module.def_submodule("optimization");
  py::module collocation_submodule =
      gsplines_module.def_submodule("collocation");
  py::module functional_analysis_submodule =
      gsplines_module.def_submodule("functional_analysis");

  // ----------------
  // Basis Submodule
  // ----------------
  py::class_<gsplines::basis::Basis, gsplines::basis::PyBasis>(basis_submodule,
                                                               "Basis")
      .def(py::init<std::size_t, const std::string &>())
      .def("get_dim", &gsplines::basis::Basis::get_dim)
      .def("get_name", &gsplines::basis::Basis::get_name)
      .def("eval_on_window", &gsplines::basis::Basis::eval_on_window)
      .def("eval_derivative_on_window",
           &gsplines::basis::Basis::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &gsplines::basis::Basis::eval_derivative_wrt_tau_on_window);

  py::class_<gsplines::basis::BasisLegendre, gsplines::basis::Basis>(
      basis_submodule, "BasisLegendre")
      .def(py::init<std::size_t>())
      .def("eval_on_window", &gsplines::basis::BasisLegendre::eval_on_window)
      .def("eval_derivative_on_window",
           &gsplines::basis::BasisLegendre::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &gsplines::basis::BasisLegendre::eval_derivative_wrt_tau_on_window);

  basis_submodule.def("string_to_basis", gsplines::basis::string_to_basis);

  // --------------------
  // Functions Submodule
  // --------------------
  py::class_<gsplines::functions::FunctionBase,
             gsplines::functions::PyFunctionBase>(functions_submodule,
                                                  "FunctionBase")
      .def(py::init<std::pair<double, double>, std::size_t,
                    const std::string &>())
      .def("get_codom_dim", &gsplines::functions::FunctionBase::get_codom_dim)
      .def("__call__", &gsplines::functions::FunctionBase::operator())
      .def("is_point_in_domain",
           &gsplines::functions::FunctionBase::is_point_in_domain)
      .def("get_domain", &gsplines::functions::FunctionBase::get_domain)
      .def("print", &gsplines::functions::FunctionBase::print,
           py::arg("_indent") = 0);

  py::class_<gsplines::functions::FunctionExpression,
             gsplines::functions::FunctionBase>(functions_submodule,
                                                "FunctionExpression")
      .def(py::init<std::pair<double, double>, std::size_t>())
      .def("deriv", &gsplines::functions::FunctionExpression::derivate,
           py::arg("_deg") = 1)
      .def("concat",
           [](const gsplines::functions::FunctionExpression &_self,
              const gsplines::functions::FunctionExpression &_that) {
             return _self.concat(_that);
           })
      .def("compose",
           [](const gsplines::functions::FunctionExpression &_self,
              const gsplines::functions::FunctionExpression &_that) {
             return _self.compose(_that);
           })
      .def(
          "__add__",
          [](const gsplines::functions::FunctionExpression &_lhs,
             const gsplines::functions::FunctionExpression &_rhs) {
            return _lhs + _rhs;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const gsplines::functions::FunctionExpression &_lhs,
             const gsplines::functions::FunctionExpression &_rhs) {
            return _lhs - _rhs;
          },
          py::is_operator())
      .def(
          "__mul__",
          [](const gsplines::functions::FunctionExpression &_lhs,
             const gsplines::functions::FunctionExpression &_rhs) {
            return _lhs * _rhs;
          },
          py::is_operator())
      .def(
          "__neg__",
          [](const gsplines::functions::FunctionExpression &_this) {
            return -_this;
          },
          py::is_operator())
      .def("concat",
           [](const gsplines::functions::FunctionExpression &_self,
              const gsplines::functions::Function &_that) {
             return _self.concat(_that);
           })
      .def("compose",
           [](const gsplines::functions::FunctionExpression &_self,
              const gsplines::functions::Function &_that) {
             return _self.compose(_that);
           })
      .def(
          "__add__",
          [](const gsplines::functions::FunctionExpression &_lhs,
             const gsplines::functions::Function &_rhs) { return _lhs + _rhs; },
          py::is_operator())
      .def(
          "__sub__",
          [](const gsplines::functions::FunctionExpression &_lhs,
             const gsplines::functions::Function &_rhs) { return _lhs - _rhs; },
          py::is_operator())
      .def(
          "__mul__",
          [](const gsplines::functions::FunctionExpression &_lhs,
             const gsplines::functions::Function &_rhs) { return _lhs * _rhs; },
          py::is_operator());

  py::class_<gsplines::functions::Function, gsplines::functions::PyFunction>(
      gsplines_module, "Function")
      .def(py::init<std::pair<double, double>, std::size_t,
                    const std::string &>())
      .def("__call__", &gsplines::functions::FunctionBase::operator())
      .def("get_codom_dim", &gsplines::functions::FunctionBase::get_codom_dim)
      .def("is_point_in_domain",
           &gsplines::functions::FunctionBase::is_point_in_domain)
      .def("get_domain", &gsplines::functions::FunctionBase::get_domain)
      .def("print", &gsplines::functions::FunctionBase::print,
           py::arg("_indent") = 0)
      .def("concat",
           [](const gsplines::functions::Function &_self,
              const gsplines::functions::FunctionExpression &_that) {
             return _self.concat(_that);
           })
      .def("compose",
           [](const gsplines::functions::Function &_self,
              const gsplines::functions::FunctionExpression &_that) {
             return _self.compose(_that);
           })
      .def(
          "__add__",
          [](const gsplines::functions::Function &_lhs,
             const gsplines::functions::FunctionExpression &_rhs) {
            return _lhs + _rhs;
          },
          py::is_operator())
      .def(
          "__sub__",
          [](const gsplines::functions::Function &_lhs,
             const gsplines::functions::FunctionExpression &_rhs) {
            return _lhs - _rhs;
          },
          py::is_operator())
      .def(
          "__mul__",
          [](const gsplines::functions::Function &_lhs,
             const gsplines::functions::FunctionExpression &_rhs) {
            return _lhs * _rhs;
          },
          py::is_operator())
      .def(
          "__neg__",
          [](const gsplines::functions::Function &_this) { return -_this; },
          py::is_operator())
      .def("concat",
           [](const gsplines::functions::Function &_self,
              const gsplines::functions::Function &_that) {
             return _self.concat(_that);
           })
      .def("compose",
           [](const gsplines::functions::Function &_self,
              const gsplines::functions::Function &_that) {
             return _self.compose(_that);
           })
      .def(
          "__add__",
          [](const gsplines::functions::Function &_lhs,
             const gsplines::functions::Function &_rhs) { return _lhs + _rhs; },
          py::is_operator())
      .def(
          "__sub__",
          [](const gsplines::functions::Function &_lhs,
             const gsplines::functions::Function &_rhs) { return _lhs - _rhs; },
          py::is_operator())
      .def(
          "__mul__",
          [](const gsplines::functions::Function &_lhs,
             const gsplines::functions::Function &_rhs) { return _lhs * _rhs; },
          py::is_operator())
      .def("deriv",
           [](const gsplines::functions::Function &_self, std::size_t _deg) {
             return gsplines::functions::FunctionExpression(_self).deriv(_deg);
           });

  py::class_<gsplines::functions::Exponential, gsplines::functions::Function>(
      functions_submodule, "Exponential")
      .def(py::init<std::pair<double, double>>());

  py::class_<gsplines::functions::ConstFunction, gsplines::functions::Function>(
      functions_submodule, "ConstFunction")
      .def(py::init<std::pair<double, double>, std::size_t, double>());

  py::class_<gsplines::functions::Cos, gsplines::functions::Function>(
      functions_submodule, "Cos")
      .def(py::init<std::pair<double, double>>());

  py::class_<gsplines::functions::DomainLinearDilation,
             gsplines::functions::Function>(functions_submodule,
                                            "DomainLinearDilation")
      .def(py::init<std::pair<double, double>, double>());

  py::class_<gsplines::functions::Identity,
             gsplines::functions::DomainLinearDilation>(functions_submodule,
                                                        "Identity")
      .def(py::init<std::pair<double, double>>());

  py::class_<gsplines::functions::Sin, gsplines::functions::Function>(
      functions_submodule, "Sin")
      .def(py::init<std::pair<double, double>>());

  py::class_<gsplines::functions::CanonicPolynomial,
             gsplines::functions::Function>(functions_submodule,
                                            "CanonicPolynomial")
      .def(py::init<std::pair<double, double>,
                    const Eigen::Ref<const Eigen::VectorXd>>());

  // ----------------------
  // Collocation Submodule
  // ----------------------
  collocation_submodule.def(
      "legendre_gauss_lobatto_points_weights",
      &gsplines::collocation::legendre_gauss_lobatto_points_and_weights);
  collocation_submodule.def("q_and_evaluation",
                            &gsplines::collocation::q_and_evaluation);

  // -----------------------
  // Optimization Submodule
  // -----------------------
  optimization_submodule.def("broken_lines_path",
                             &gsplines::optimization::broken_lines_path);
  optimization_submodule.def(
      "minimum_acceleration_path",
      &gsplines::optimization::minimum_acceleration_path);
  optimization_submodule.def("minimum_jerk_path",
                             &gsplines::optimization::minimum_jerk_path);
  optimization_submodule.def("minimum_snap_path",
                             &gsplines::optimization::minimum_snap_path);
  optimization_submodule.def("minimum_crackle_path",
                             &gsplines::optimization::minimum_crackle_path);

  optimization_submodule.def("optimal_sobolev_norm",
                             &gsplines::optimization::optimal_sobolev_norm);

  // ---------------
  // Gsplines module
  // ---------------
  py::class_<gsplines::Interpolator>(gsplines_module, "InterpolatorBase")
      .def(py::init<std::size_t, std::size_t, gsplines::basis::Basis &>());

  py::class_<gsplines::PyInterpolator, gsplines::Interpolator>(gsplines_module,
                                                               "Interpolator")
      .def(py::init<std::size_t, std::size_t, gsplines::basis::Basis &>())
      .def("interpolate", &gsplines::PyInterpolator::py_interpolate)
      .def("solve_interpolation",
           &gsplines::PyInterpolator::py_solve_interpolation)
      .def("print_interpolating_matrix",
           &gsplines::PyInterpolator::print_interpolating_matrix)
      .def("get_coeff_derivative_wrt_tau",
           &gsplines::PyInterpolator::get_coeff_derivative_wrt_tau)
      .def("print_interpolating_vector",
           &gsplines::PyInterpolator::print_interpolating_vector);

  py::class_<gsplines::GSpline, gsplines::functions::Function>(gsplines_module,
                                                               "GSpline")
      .def(py::init<
           std::pair<double, double>, std::size_t, std::size_t,
           gsplines::basis::Basis &, const Eigen::Ref<const Eigen::VectorXd>,
           const Eigen::Ref<const Eigen::VectorXd>, const std::string &>())
      .def("get_exec_time", &gsplines::GSpline::get_exec_time)
      .def("get_domain_breakpoints", &gsplines::GSpline::get_domain_breakpoints)
      .def("get_number_of_intervals",
           &gsplines::GSpline::get_number_of_intervals)
      .def("get_interval_lengths", &gsplines::GSpline::get_interval_lengths)
      .def("get_waypoints", &gsplines::GSpline::get_waypoints)
      .def("get_basis_name", &gsplines::GSpline::get_basis_name)
      .def("linear_scaling_new_execution_time",
           &gsplines::GSpline::linear_scaling_new_execution_time)
      .def("get_coefficients", &gsplines::GSpline::get_coefficients);

  py::class_<gsplines::functional_analysis::SobolevNorm>(
      functional_analysis_submodule, "SobolevNormBase")
      .def(py::init<const Eigen::Ref<const Eigen::MatrixXd>,
                    gsplines::basis::Basis &,
                    std::vector<std::pair<std::size_t, double>>>())
      .def("__call__", &gsplines::functional_analysis::SobolevNorm::operator())
      .def("deriv_wrt_interval_len",
           &gsplines::functional_analysis::SobolevNorm::deriv_wrt_interval_len);

  py::class_<gsplines::PySobolevNorm>(functional_analysis_submodule,
                                      "SobolevNorm")
      .def(py::init<const py::EigenDRef<const Eigen::MatrixXd>,
                    gsplines::basis::Basis &,
                    std::vector<std::pair<std::size_t, double>>>())
      .def("__call__", &gsplines::PySobolevNorm::operator())
      .def("deriv_wrt_interval_len",
           &gsplines::PySobolevNorm::deriv_wrt_interval_len);
}
