
#include "pybasis.hpp"
#include "pyfunctions.hpp"

PYBIND11_MODULE(pygsplines, m) {
  py::class_<gsplines::basis::Basis, gsplines::basis::PyBasis>(m, "Basis")
      .def(py::init<std::size_t>())
      .def("get_dim", &gsplines::basis::Basis::get_dim)
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

  py::class_<gsplines::PiecewiseFunction>(m, "PiecewiseFunction")
      .def(py::init<std::size_t, std::size_t, gsplines::basis::Basis &,
                    const Eigen::Ref<const Eigen::VectorXd>,
                    const Eigen::Ref<const Eigen::VectorXd>>())
      .def("__call__", &gsplines::PiecewiseFunction::operator())
      .def("get_exec_time", &gsplines::PiecewiseFunction::get_exec_time)
      .def("get_codom_dim", &gsplines::PiecewiseFunction::get_codom_dim)
      .def("get_domain_breakpoints",
           &gsplines::PiecewiseFunction::get_domain_breakpoints)
      .def("get_coeff", &gsplines::PiecewiseFunction::get_coeff)
      .def("deriv", &gsplines::PiecewiseFunction::deriv);

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

  m.def("optimal_sobolev_norm", &gsplines_opt::optimal_sobolev_norm);

  // Function Expression
  /*
  py::class_<gsplines::functions::FunctionExpression>(m, "FunctionExpression")
      .def(py::init<
           std::pair<double, double>, std::size_t,
           gsplines::functions::FunctionExpression::Type,
           std::list<std::unique_ptr<gsplines::functions::FunctionExpression>>
               &,
           const std::string &>());
*/
  /*
        .def(py::init(
            [](std::pair<double, double> a, std::size_t b,
               gsplines::functions::FunctionExpression::Type c,
               std::list<std::unique_ptr<gsplines::functions::FunctionExpression>>
                   &d,
               const std::string &e) {
              return std::make_unique<gsplines::functions::FunctionExpression>(
                  a, b, c, d, e);
            }));
        .def("get_codom_dim",
             &gsplines::functions::FunctionExpression::get_codom_dim)
        .def("is_point_in_domain",
             &gsplines::functions::FunctionExpression::is_point_in_domain)
        .def("get_domain", &gsplines::functions::FunctionExpression::get_domain)
        .def("deriv", &gsplines::functions::FunctionExpression::derivate)
        .def(
            "__sum__",
            [](const gsplines::functions::FunctionExpression &_lhs,
               const gsplines::functions::FunctionExpression &_rhs) {
              return _lhs + _rhs;
            },
            py::is_operator())
        .def(
            "__mul__",
            [](const gsplines::functions::FunctionExpression &_lhs,
               const gsplines::functions::FunctionExpression &_rhs) {
              return _lhs * _rhs;
            },
            py::is_operator())
        .def("__call__",
     &gsplines::functions::FunctionExpression::operator());*/
  // Operations
  // m.def("optimal_sobolev_norm", &gsplines_opt::optimal_sobolev_norm);
}
