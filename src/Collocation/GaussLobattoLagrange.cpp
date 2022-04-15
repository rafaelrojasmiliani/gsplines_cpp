
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/FunctionalAnalysis/Integral.hpp>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Tools.hpp>
#include <stdexcept>

namespace gsplines {

namespace collocation {

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    std::pair<double, double> _domain, std::size_t _codom_dim,
    std::size_t _n_intervals, std::size_t _n_glp,
    const Eigen::Ref<const Eigen::VectorXd> _coefficents,
    const Eigen::Ref<const Eigen::VectorXd> _tauv)
    : GSplineInheritanceHelper(
          _domain, _codom_dim, _n_intervals,
          ::gsplines::basis::BasisLagrangeGaussLobatto(_n_glp), _coefficents,
          _tauv) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    std::pair<double, double> _domain, std::size_t _codom_dim,
    std::size_t _n_intervals, std::size_t _n_glp,
    Eigen::VectorXd &&_coefficents, Eigen::VectorXd &&_tauv)
    : GSplineInheritanceHelper(
          _domain, _codom_dim, _n_intervals,
          ::gsplines::basis::BasisLagrangeGaussLobatto(_n_glp),
          std::move(_coefficents), std::move(_tauv)) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    const GaussLobattoLagrangeSpline &_that)
    : GSplineInheritanceHelper(_that) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(
    GaussLobattoLagrangeSpline &&_that)
    : GSplineInheritanceHelper(std::move(_that)) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(const GSpline &_that)
    : GSplineInheritanceHelper(_that) {}

GaussLobattoLagrangeSpline::GaussLobattoLagrangeSpline(GSpline &&_that)
    : GSplineInheritanceHelper(std::move(_that)) {}

GaussLobattoLagrangeSpline &GaussLobattoLagrangeSpline::operator=(
    const GaussLobattoLagrangeSpline &_that) & {

  if (not same_discretization(_that))
    throw std::invalid_argument("Discretization must be the same");
  if (not same_domain(_that))
    throw std::invalid_argument("GSplines must have same domain");
  if (not(same_codomain(_that)))
    throw std::invalid_argument("Assigmention requires same codomain");

  coefficients_ = _that.coefficients_;
  return *this;
}
GaussLobattoLagrangeSpline &
GaussLobattoLagrangeSpline::operator=(GaussLobattoLagrangeSpline &&_that) & {
  if (this == &_that) {
    throw std::invalid_argument(" Cannot move assign to itself");
  }
  if (not same_discretization(_that))
    throw std::invalid_argument("Discretization must be the same");
  if (not same_domain(_that))
    throw std::invalid_argument("GSplines must have same domain");
  if (not(same_codomain(_that)))
    throw std::invalid_argument("Assigmention requires same codomain");

  coefficients_ = std::move(_that.coefficients_);
  return *this;
}

Eigen::Map<Eigen::MatrixXd>
GaussLobattoLagrangeSpline::get_value_block(std::size_t _interval) {
  std::size_t nglp = get_nglp();
  std::size_t dim = get_codom_dim();
  std::size_t chunk_size = nglp * dim;
  if (_interval >= get_number_of_intervals()) {
    throw std::invalid_argument("Cannot accest to that interval");
  }
  double *data =
      coefficients_.segment(_interval * chunk_size, chunk_size).data();
  return Eigen::Map<Eigen::MatrixXd>(data, nglp, get_codom_dim());
}

Eigen::Map<const Eigen::MatrixXd>
GaussLobattoLagrangeSpline::get_value_block(std::size_t _interval) const {
  std::size_t nglp = get_nglp();
  std::size_t dim = get_codom_dim();
  std::size_t chunk_size = nglp * dim;
  if (_interval >= get_number_of_intervals()) {
    throw std::invalid_argument("Cannot accest to that interval");
  }
  const double *data =
      coefficients_.segment(_interval * chunk_size, chunk_size).data();
  return Eigen::Map<const Eigen::MatrixXd>(data, nglp, get_codom_dim());
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
bool GaussLobattoLagrangeSpline::same_discretization(
    const GaussLobattoLagrangeSpline &_that) const {
  return get_nglp() == _that.get_nglp() and
         get_number_of_intervals() == _that.get_number_of_intervals();
}

// Operations -----------------------------
// -----------------------------------------

GaussLobattoLagrangeSpline GaussLobattoLagrangeSpline::operator*(
    const GaussLobattoLagrangeSpline &_that) const & {

  std::cout << "----------------------\n";
  fflush(stdout);
  if (not same_discretization(_that))
    throw std::invalid_argument("Discretization must be the same");
  if (not same_domain(_that))
    throw std::invalid_argument("GSplines must have same domain");
  if (not(get_codom_dim() == 1 or _that.get_codom_dim() == 1))
    throw std::invalid_argument("Multiplication can be only done by scalars");

  const GaussLobattoLagrangeSpline &scalar =
      (_that.get_codom_dim() == 1) ? _that : *this;
  const GaussLobattoLagrangeSpline &vector =
      (_that.get_codom_dim() != 1) ? _that : *this;
  GaussLobattoLagrangeSpline result(vector);

  std::size_t n_inter = get_number_of_intervals();

  for (std::size_t i = 0; i < n_inter; i++) {

    Eigen::Map<Eigen::MatrixXd> matrix = result.get_value_block(i);

    matrix.array().colwise() *= scalar.get_value_block(i).col(0).array();
  }

  return result;
}

GaussLobattoLagrangeSpline GaussLobattoLagrangeSpline::operator*(
    const GaussLobattoLagrangeSpline &_that) && {

  if (not same_discretization(_that))
    throw std::invalid_argument("Discretization must be the same");
  if (not same_domain(_that))
    throw std::invalid_argument("GSplines must have same domain");
  if (not(get_codom_dim() == 1 or _that.get_codom_dim() == 1))
    throw std::invalid_argument("Multiplication can be only done by scalars");

  const GaussLobattoLagrangeSpline &scalar =
      (_that.get_codom_dim() == 1) ? _that : *this;
  const GaussLobattoLagrangeSpline &vector =
      (_that.get_codom_dim() != 1) ? _that : *this;

  if (&vector == this) {
    std::size_t n_inter = get_number_of_intervals();

    for (std::size_t i = 0; i < n_inter; i++) {

      Eigen::Map<Eigen::MatrixXd> matrix = get_value_block(i);

      matrix.array().colwise() *= scalar.get_value_block(i).col(0).array();
    }
    return std::move(*this);
  }
  GaussLobattoLagrangeSpline result(vector);

  std::size_t n_inter = get_number_of_intervals();

  for (std::size_t i = 0; i < n_inter; i++) {

    Eigen::Map<Eigen::MatrixXd> matrix = result.get_value_block(i);

    matrix.array().colwise() *= scalar.get_value_block(i).col(0).array();
  }

  return result;
}

GaussLobattoLagrangeSpline GaussLobattoLagrangeSpline::operator*(
    GaussLobattoLagrangeSpline &&_that) const & {

  if (not same_discretization(_that))
    throw std::invalid_argument("Discretization must be the same");
  if (not same_domain(_that))
    throw std::invalid_argument("GSplines must have same domain");
  if (not(get_codom_dim() == 1 or _that.get_codom_dim() == 1))
    throw std::invalid_argument("Multiplication can be only done by scalars");

  const GaussLobattoLagrangeSpline &scalar =
      (_that.get_codom_dim() == 1) ? _that : *this;
  const GaussLobattoLagrangeSpline &vector =
      (_that.get_codom_dim() != 1) ? _that : *this;

  if (&vector == &_that) {
    std::size_t n_inter = get_number_of_intervals();

    for (std::size_t i = 0; i < n_inter; i++) {

      Eigen::Map<Eigen::MatrixXd> matrix = _that.get_value_block(i);

      matrix.array().colwise() *= scalar.get_value_block(i).col(0).array();
    }
    return std::move(_that);
  }
  GaussLobattoLagrangeSpline result(vector);

  std::size_t n_inter = get_number_of_intervals();

  for (std::size_t i = 0; i < n_inter; i++) {

    Eigen::Map<Eigen::MatrixXd> matrix = result.get_value_block(i);

    matrix.array().colwise() *= scalar.get_value_block(i).col(0).array();
  }

  return result;
}

GaussLobattoLagrangeSpline
GaussLobattoLagrangeSpline::operator*(GaussLobattoLagrangeSpline &&_that) && {

  if (not same_discretization(_that))
    throw std::invalid_argument("Discretization must be the same");
  if (not same_domain(_that))
    throw std::invalid_argument("GSplines must have same domain");
  if (not(get_codom_dim() == 1 or _that.get_codom_dim() == 1))
    throw std::invalid_argument("Multiplication can be only done by scalars");

  GaussLobattoLagrangeSpline &scalar =
      (_that.get_codom_dim() == 1) ? _that : *this;
  GaussLobattoLagrangeSpline &vector =
      (_that.get_codom_dim() != 1) ? _that : *this;

  GaussLobattoLagrangeSpline &result = vector;

  std::size_t n_inter = get_number_of_intervals();

  for (std::size_t i = 0; i < n_inter; i++) {

    Eigen::Map<Eigen::MatrixXd> matrix = result.get_value_block(i);

    matrix.array().colwise() *= scalar.get_value_block(i).col(0).array();
  }

  return std::move(result);
}
GaussLobattoLagrangeSpline operator*(double _a,
                                     const GaussLobattoLagrangeSpline &_that) {
  GaussLobattoLagrangeSpline result(_that);
  result.coefficients_ *= _a;
  return result;
}
GaussLobattoLagrangeSpline operator*(double _a,
                                     GaussLobattoLagrangeSpline &&_that) {
  _that.coefficients_ *= _a;
  return std::move(_that);
}
GaussLobattoLagrangeSpline operator*(const GaussLobattoLagrangeSpline &_that,
                                     double _a) {
  return _a * _that;
}
GaussLobattoLagrangeSpline operator*(GaussLobattoLagrangeSpline &&_that,
                                     double _a) {
  return _a * std::move(_that);
}
GaussLobattoLagrangeSpline operator-(const GaussLobattoLagrangeSpline &_that) {

  GaussLobattoLagrangeSpline result(_that);
  result.coefficients_ *= -1.0;
  return result;
}
GaussLobattoLagrangeSpline operator-(GaussLobattoLagrangeSpline &&_that) {

  _that.coefficients_ *= -1.0;
  return std::move(_that);
}

} // namespace collocation
} // namespace gsplines
