
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/FunctionalAnalysis/Integral.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>

namespace gsplines {

namespace collocation {
GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    std::pair<double, double> _domain, std::size_t _codom_dim,
    const Eigen::Ref<const Eigen::VectorXd> _coefficents, std::size_t _n_glp,
    std::size_t _n_intervals)
    : FunctionInheritanceHelper(
          _domain, _codom_dim, _n_intervals,
          ::gsplines::basis::BasisLagrange(
              legendre_gauss_lobatto_points(_n_glp)),
          _coefficents,
          (_domain.second - _domain.first) *
              Eigen::VectorXd::Ones(_n_glp * _n_intervals * _codom_dim)),
      value_at_nodes_((_n_glp - 1) * (_n_intervals - 1) + _n_glp, _codom_dim) {

  Eigen::MatrixXd mat = Eigen::Map<const Eigen::MatrixXd>(
      get_coefficients().data(), _n_glp * _n_intervals, _codom_dim);

  for (std::size_t interval = 0; interval < _n_intervals - 1; interval++) {
    value_at_nodes_.middleRows(interval * (_n_glp - 1), _n_glp - 1) =
        mat.middleRows(interval * _n_glp, _n_glp - 1);
  }
  value_at_nodes_.bottomRows(_n_glp) = mat.bottomRows(_n_glp);
}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    const GaussLobattoLagrangeSpline &_that)
    : FunctionInheritanceHelper(_that), value_at_nodes_(_that.value_at_nodes_) {
}
GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    GaussLobattoLagrangeSpline &&_that)
    : FunctionInheritanceHelper(std::move(_that)),
      value_at_nodes_(std::move(_that.value_at_nodes_)) {}

GaussLobattoLagrangeSpline *
GaussLobattoLagrangeSpline::deriv_impl(std::size_t _deg) const {

  Eigen::VectorXd result_coeff(get_coefficients());
  int i0;

  for (std::size_t der_coor = 1; der_coor <= _deg; der_coor++) {
    for (std::size_t interval_coor = 0; interval_coor < get_intervals_num();
         interval_coor++) {
      for (std::size_t codom_coor = 0; codom_coor < get_codom_dim();
           codom_coor++) {
        i0 = interval_coor * get_basis_dim() * get_codom_dim() +
             get_basis_dim() * codom_coor;
        result_coeff.segment(i0, get_basis_dim()) =
            get_basis().get_derivative_matrix_block(_deg) *
            result_coeff.segment(i0, get_basis_dim());
      }
    }
  }

  return new GaussLobattoLagrangeSpline(get_domain(), get_codom_dim(),
                                        result_coeff, get_basis_dim(),
                                        get_intervals_num());
}

GaussLobattoLagrangeSpline
gauss_lobatto_lagrange_approximtion(::gsplines::functions::FunctionBase &_in,
                                    std::size_t _n_glp,
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
  return GaussLobattoLagrangeSpline(_in.get_domain(), codom_dim, result_coeff,
                                    _n_glp, _n_intervals);
}

double integral(::gsplines::functions::FunctionBase &_diffeo,
                ::gsplines::functions::FunctionBase &_path) {

  return ::gsplines::functional_analysis::integral(
      _path.dot(_path).compose(_diffeo.to_expression()));
}
} // namespace collocation
} // namespace gsplines
