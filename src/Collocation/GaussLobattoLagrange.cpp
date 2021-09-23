
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>

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
          Eigen::VectorXd::Zero(_n_glp * _n_intervals * _codom_dim),
          (_domain.second - _domain.first) *
              Eigen::VectorXd::Ones(_n_glp * _n_intervals * _codom_dim)) {}
GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    const GaussLobattoLagrangeSpline &_that)
    : FunctionInheritanceHelper(_that) {}
GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    GaussLobattoLagrangeSpline &&_that)
    : FunctionInheritanceHelper(_that) {}

GaussLobattoLagrangeSpline *
GaussLobattoLagrangeSpline::deriv_impl(std::size_t _deg) const {

  Eigen::VectorXd result_coeff(get_coefficients());
  int der_coor;
  int interval_coor;
  int codom_coor;
  int i0;

  for (der_coor = 1; der_coor <= _deg; der_coor++) {
    for (interval_coor = 0; interval_coor < get_intervals_num();
         interval_coor++) {
      for (codom_coor = 0; codom_coor < get_codom_dim(); codom_coor++) {
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
    double local_right_bound =
        left_bound + (interval + 1) * local_interval_length;

    Eigen::VectorXd local_points =
        (glp.array() + 1.0) * local_interval_length / 2.0 + local_left_bound;

    _in.value(local_points, local_value);

    std::size_t i0 = interval * _n_glp * codom_dim;
    std::size_t i1 = (interval + 1) * _n_glp * codom_dim;

    result_coeff.segment(i0, _n_glp * codom_dim) =
        Eigen::Map<Eigen::VectorXd>(local_value.data(), local_value.size());
  }
  return GaussLobattoLagrangeSpline(_in.get_domain(), codom_dim, result_coeff,
                                    _n_glp, _n_intervals);
}
} // namespace collocation
} // namespace gsplines
