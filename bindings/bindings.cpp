
#include "pybasis.hpp"
#include "pyfunctions.hpp"
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>

PYBIND11_MODULE(gsplines, m) {
  py::class_<gsplines::basis::Basis, gsplines::basis::PyBasis>(m, "Basis")
      .def(py::init<std::size_t, const std::string &>())
      .def("get_dim", &gsplines::basis::Basis::get_dim)
      .def("get_name", &gsplines::basis::Basis::get_name)
      .def("eval_on_window", &gsplines::basis::Basis::eval_on_window)
      .def("eval_derivative_on_window",
           &gsplines::basis::Basis::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &gsplines::basis::Basis::eval_derivative_wrt_tau_on_window);

  py::class_<gsplines::basis::BasisLegendre, gsplines::basis::Basis>(
      m, "BasisLegendre")
      .def(py::init<std::size_t>())
      .def("eval_on_window", &gsplines::basis::BasisLegendre::eval_on_window)
      .def("eval_derivative_on_window",
           &gsplines::basis::BasisLegendre::eval_derivative_on_window)
      .def("eval_derivative_wrt_tau_on_window",
           &gsplines::basis::BasisLegendre::eval_derivative_wrt_tau_on_window);

  // Function Expression
  py::class_<gsplines::functions::FunctionExpression>(m, "FunctionExpression")
      .def(py::init<std::pair<double, double>, std::size_t>())
      .def("get_codom_dim",
           &gsplines::functions::FunctionExpression::get_codom_dim)
      .def("is_point_in_domain",
           &gsplines::functions::FunctionExpression::is_point_in_domain)
      .def("get_domain", &gsplines::functions::FunctionExpression::get_domain)
      .def("deriv", &gsplines::functions::FunctionExpression::derivate,
           py::arg("_deg") = 1)
      .def("print", &gsplines::functions::FunctionExpression::print,
           py::arg("_indent") = 0)
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
      .def("__call__", &gsplines::functions::FunctionExpression::operator());

  py::class_<gsplines::functions::Exponential,
             gsplines::functions::FunctionExpression>(m, "Exponential")
      .def(py::init<std::pair<double, double>>())
      .def("__call__", &gsplines::functions::FunctionExpression::operator());
  py::class_<gsplines::functions::ConstFunction,
             gsplines::functions::FunctionExpression>(m, "ConstFunction")
      .def(py::init<std::pair<double, double>, std::size_t, double>())
      .def("__call__", &gsplines::functions::FunctionExpression::operator());

  py::class_<gsplines::functions::Cos, gsplines::functions::FunctionExpression>(
      m, "Cos")
      .def(py::init<std::pair<double, double>>())
      .def("__call__", &gsplines::functions::FunctionExpression::operator());

  py::class_<gsplines::functions::DomainLinearDilation,
             gsplines::functions::FunctionExpression>(m, "DomainLinearDilation")
      .def(py::init<std::pair<double, double>, double>())
      .def("__call__", &gsplines::functions::FunctionExpression::operator());

  py::class_<gsplines::functions::Identity,
             gsplines::functions::DomainLinearDilation>(m, "Identity")
      .def(py::init<std::pair<double, double>>())
      .def("__call__", &gsplines::functions::FunctionExpression::operator());

  py::class_<gsplines::functions::Sin, gsplines::functions::FunctionExpression>(
      m, "Sin")
      .def(py::init<std::pair<double, double>>())
      .def("__call__", &gsplines::functions::FunctionExpression::operator());

  py::class_<gsplines::GSpline, gsplines::functions::FunctionExpression>(
      m, "GSpline")
      .def(py::init<std::pair<double, double>, std::size_t, std::size_t,
                    gsplines::basis::Basis &,
                    const Eigen::Ref<const Eigen::VectorXd>,
                    const Eigen::Ref<const Eigen::VectorXd>>())
      .def("__call__", &gsplines::GSpline::operator())
      .def("get_exec_time", &gsplines::GSpline::get_exec_time)
      .def("get_domain_breakpoints", &gsplines::GSpline::get_domain_breakpoints)
      .def("get_number_of_intervals", &gsplines::GSpline::get_number_of_intervals)
      .def("get_interval_lengths", &gsplines::GSpline::get_interval_lengths)
      .def("get_waypoints", &gsplines::GSpline::get_waypoints)
      .def("get_basis_name", &gsplines::GSpline::get_basis_name)
      .def("linear_scaling_new_execution_time",
           &gsplines::GSpline::linear_scaling_new_execution_time)
      .def("get_coefficients", &gsplines::GSpline::get_coefficients);

  py::class_<gsplines::functions::CanonicPolynomial,
             gsplines::functions::FunctionExpression>(m, "CanonicPolynomial")
      .def(py::init<std::pair<double, double>,
                    const Eigen::Ref<const Eigen::VectorXd>>());
  /* Interpolator
   * */

  py::class_<gsplines::Interpolator>(m, "Interpolator")
      .def(py::init<std::size_t, std::size_t, gsplines::basis::Basis &>());

  py::class_<gsplines::PyInterpolator, gsplines::Interpolator>(m,
                                                               "PyInterpolator")
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

  py::class_<gsplines::SobolevNorm>(m, "SobolevNorm")
      .def(py::init<const Eigen::Ref<const Eigen::MatrixXd>,
                    gsplines::basis::Basis &,
                    std::vector<std::pair<std::size_t, double>>>())
      .def("__call__", &gsplines::SobolevNorm::operator())
      .def("deriv_wrt_interval_len",
           &gsplines::SobolevNorm::deriv_wrt_interval_len);

  py::class_<gsplines::PySobolevNorm>(m, "PySobolevNorm")
      .def(py::init<const py::EigenDRef<const Eigen::MatrixXd>,
                    gsplines::basis::Basis &,
                    std::vector<std::pair<std::size_t, double>>>())
      .def("__call__", &gsplines::SobolevNorm::operator())
      .def("deriv_wrt_interval_len",
           &gsplines::SobolevNorm::deriv_wrt_interval_len);

  m.def("broken_lines_path", &gsplines_opt::broken_lines_path);
  m.def("minimum_acceleration_path", &gsplines_opt::minimum_acceleration_path);
  m.def("minimum_jerk_path", &gsplines_opt::minimum_jerk_path);
  m.def("minimum_snap_path", &gsplines_opt::minimum_snap_path);
  m.def("minimum_crackle_path", &gsplines_opt::minimum_crackle_path);

  m.def("legendre_gauss_lobatto_points_weights",
        &gsplines::collocation::legendre_gauss_lobatto_points_and_weights);
  m.def("q_and_evaluation", &gsplines::collocation::q_and_evaluation);
  // Operations
  m.def("optimal_sobolev_norm", &gsplines_opt::optimal_sobolev_norm);
  m.def("string_to_basis", gsplines::basis::string_to_basis);
}
