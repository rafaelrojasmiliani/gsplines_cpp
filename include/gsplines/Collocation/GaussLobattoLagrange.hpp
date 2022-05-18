#ifndef GAUSSLOBATTOLAGRANGE_H
#define GAUSSLOBATTOLAGRANGE_H
#include <Eigen/SparseCore>
#include <cstddef>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <gsplines/GSpline.hpp>

namespace gsplines {

namespace collocation {
class GaussLobattoLagrangeSpline
    : public GSplineInheritanceHelper<GaussLobattoLagrangeSpline> {

  friend GSpline;
  friend GaussLobattoLagrangeSpline
  operator*(double _a, const GaussLobattoLagrangeSpline &_that);
  friend GaussLobattoLagrangeSpline
  operator*(double _a, GaussLobattoLagrangeSpline &&_that);
  friend GaussLobattoLagrangeSpline
  operator*(const GaussLobattoLagrangeSpline &_that, double _a);
  friend GaussLobattoLagrangeSpline
  operator*(GaussLobattoLagrangeSpline &&_that, double _a);
  friend GaussLobattoLagrangeSpline
  operator-(const GaussLobattoLagrangeSpline &_that);
  friend GaussLobattoLagrangeSpline
  operator-(GaussLobattoLagrangeSpline &&_that);

private:
  Eigen::Map<Eigen::MatrixXd> get_value_block(std::size_t _interval);
  Eigen::Map<const Eigen::MatrixXd>
  get_value_block(std::size_t _interval) const;

public:
  using GSplineInheritanceHelper<
      GaussLobattoLagrangeSpline>::GSplineInheritanceHelper;
  GaussLobattoLagrangeSpline(std::pair<double, double> _domain,
                             std::size_t _codom_dim, std::size_t _n_intervals,
                             std::size_t _n_glp,
                             const Eigen::VectorXd &_coefficents,
                             const Eigen::VectorXd &_tauv);

  // create moving data
  GaussLobattoLagrangeSpline(std::pair<double, double> _domain,
                             std::size_t _codom_dim, std::size_t _n_intervals,
                             std::size_t _n_glp, Eigen::VectorXd &&_coefficents,
                             Eigen::VectorXd &&_tauv);

  // create a zero polynomial with equal intervals
  GaussLobattoLagrangeSpline(std::pair<double, double> _domain,
                             std::size_t _codom_dim, std::size_t _n_intervals,
                             std::size_t _n_glp);

  GaussLobattoLagrangeSpline(const GaussLobattoLagrangeSpline &_that);
  GaussLobattoLagrangeSpline(GaussLobattoLagrangeSpline &&_that);
  GaussLobattoLagrangeSpline(const GSpline &_that);
  GaussLobattoLagrangeSpline(GSpline &&_that);
  virtual ~GaussLobattoLagrangeSpline() = default;

  static GaussLobattoLagrangeSpline identity(std::pair<double, double> _domain,
                                             std::size_t _n_glp,
                                             std::size_t _n_intervals);

  static GaussLobattoLagrangeSpline
  approximate(const ::gsplines::functions::FunctionBase &_in,
              std::size_t _n_glp, std::size_t _n_intervals);

  std::size_t get_nglp() const { return get_basis_dim(); }
  Eigen::VectorXd get_domain_discretization() const {
    return legendre_gauss_lobatto_points(get_domain(), get_nglp(),
                                         get_number_of_intervals());
  }

  GaussLobattoLagrangeSpline
  operator*(const GaussLobattoLagrangeSpline &_that) const &;
  GaussLobattoLagrangeSpline
  operator*(const GaussLobattoLagrangeSpline &_that) &&;
  GaussLobattoLagrangeSpline
  operator*(GaussLobattoLagrangeSpline &&_that) const &;
  GaussLobattoLagrangeSpline operator*(GaussLobattoLagrangeSpline &&_that) &&;

  GaussLobattoLagrangeSpline &
  operator=(const GaussLobattoLagrangeSpline &_that) &;
  GaussLobattoLagrangeSpline &
  operator=(const ::gsplines::functions::FunctionBase &_that) &;
  GaussLobattoLagrangeSpline &operator=(GaussLobattoLagrangeSpline &&_that) &;

  GaussLobattoLagrangeSpline &norm(const GaussLobattoLagrangeSpline &_that) &;

  bool same_discretization(const GaussLobattoLagrangeSpline &_that) const;
};

double integral(::gsplines::functions::FunctionBase &_diffeo,
                ::gsplines::functions::FunctionBase &_path);

using GLLSpline = GaussLobattoLagrangeSpline;
} // namespace collocation
} // namespace gsplines

#endif /* GAUSSLOBATTOLAGRANGE_H */
