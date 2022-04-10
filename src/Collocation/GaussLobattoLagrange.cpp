
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/FunctionalAnalysis/Integral.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Tools.hpp>

namespace gsplines {

namespace collocation {
GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    std::pair<double, double> _domain, std::size_t _codom_dim,
    std::size_t _n_intervals, std::size_t _n_glp,
    const Eigen::Ref<const Eigen::VectorXd> _coefficents,
    const Eigen::Ref<const Eigen::VectorXd> _tauv)
    : FunctionInheritanceHelper(
          _domain, _codom_dim, _n_intervals,
          ::gsplines::basis::BasisLagrangeGaussLobatto(_n_glp), _coefficents,
          _tauv) {

  Eigen::MatrixXd mat = Eigen::Map<const Eigen::MatrixXd>(
      get_coefficients().data(), _n_glp * _n_intervals, _codom_dim);
}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    const GaussLobattoLagrangeSpline &_that)
    : FunctionInheritanceHelper(_that) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    GaussLobattoLagrangeSpline &&_that)
    : FunctionInheritanceHelper(std::move(_that)) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(const GSpline &_that)
    : FunctionInheritanceHelper(_that) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(GSpline &&_that)
    : FunctionInheritanceHelper(std::move(_that)) {}

GaussLobattoLagrangeSpline *
GaussLobattoLagrangeSpline::deriv_impl(std::size_t _deg) const {
  GSpline *aux = GSpline::deriv_impl(_deg);

  GLLSpline *result = new GaussLobattoLagrangeSpline(
      get_domain(), get_codom_dim(), get_number_of_intervals(),
      get_basis().get_dim(), aux->get_coefficients(), get_interval_lengths());

  delete aux;
  return result;
}

GaussLobattoLagrangeSpline GaussLobattoLagrangeSpline::approximate(
    const ::gsplines::functions::FunctionBase &_in, std::size_t _n_glp,
    std::size_t _n_intervals) {
  Eigen::VectorXd result_coeff(_in.get_codom_dim() * _n_glp * _n_intervals);
  Eigen::MatrixXd local_value(_n_glp, _in.get_codom_dim());

  Eigen::VectorXd glp = legendre_gauss_lobatto_points(_n_glp);
  double left_bound;
  double right_bound;
  std::tie(left_bound, right_bound) = _in.get_domain();

  std::size_t codom_dim = _in.get_codom_dim();

  double local_interval_length = (right_bound - left_bound) / _n_intervals;

  for (std::size_t interval = 0; interval < _n_intervals; interval++) {

    double local_left_bound = left_bound + interval * local_interval_length;
    // double local_right_bound =
    //    left_bound + (interval + 1) * local_interval_length;

    Eigen::VectorXd local_points =
        (glp.array() + 1.0) * local_interval_length / 2.0 + local_left_bound;

    _in.value(local_points, local_value);

    std::size_t i0 = interval * _n_glp * codom_dim;
    // std::size_t i1 = (interval + 1) * _n_glp * codom_dim;

    result_coeff.segment(i0, _n_glp * codom_dim) =
        Eigen::Map<Eigen::VectorXd>(local_value.data(), local_value.size());
  }
  return GaussLobattoLagrangeSpline(
      _in.get_domain(), codom_dim, _n_intervals, _n_glp, result_coeff,
      Eigen::VectorXd::Ones(_n_intervals) * local_interval_length);
  /*
  return new GaussLobattoLagrangeSpline(
      get_domain(), get_codom_dim(), get_number_of_intervals(),
  get_basis().get_dim(), result_coeff, get_interval_lengths());*/
}

GaussLobattoLagrangeSpline
GaussLobattoLagrangeSpline::identity(std::pair<double, double> _domain,
                                     std::size_t _n_glp,
                                     std::size_t _n_intervals) {
  ::gsplines::functions::Identity id(_domain);
  return approximate(id, _n_glp, _n_intervals);
}

double integral(::gsplines::functions::FunctionBase &_diffeo,
                ::gsplines::functions::FunctionBase &_path) {
  return ::gsplines::functional_analysis::integral(
      _path.dot(_path).compose(_diffeo.to_expression()));
}
} // namespace collocation
} // namespace gsplines
