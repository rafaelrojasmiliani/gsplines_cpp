#ifndef GAUSSLOBATTOLAGRANGE_H
#define GAUSSLOBATTOLAGRANGE_H
#include <Eigen/SparseCore>
#include <gsplines/Basis/BasisLagrange.hpp>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <gsplines/GSpline.hpp>

namespace gsplines {

namespace collocation {
class GaussLobattoLagrangeSpline
    : public ::gsplines::functions::FunctionInheritanceHelper<
          GaussLobattoLagrangeSpline, ::gsplines::GSpline,
          GaussLobattoLagrangeSpline> {
  friend GSpline;

public:
  GaussLobattoLagrangeSpline(
      std::pair<double, double> _domain, std::size_t _codom_dim,
      std::size_t _n_intervals, std::size_t _n_glp,
      const Eigen::Ref<const Eigen::VectorXd> _coefficents,
      const Eigen::Ref<const Eigen::VectorXd> _tauv);
  GaussLobattoLagrangeSpline(const GaussLobattoLagrangeSpline &_that);
  GaussLobattoLagrangeSpline(GaussLobattoLagrangeSpline &&_that);
  GaussLobattoLagrangeSpline(const GSpline &_that);
  GaussLobattoLagrangeSpline(GSpline &&_that);
  GaussLobattoLagrangeSpline &
  operator=(const GaussLobattoLagrangeSpline &) = delete;
  virtual ~GaussLobattoLagrangeSpline() = default;

  static GaussLobattoLagrangeSpline identity(std::pair<double, double> _domain,
                                             std::size_t _n_glp,
                                             std::size_t _n_intervals);

  static GaussLobattoLagrangeSpline
  approximate(const ::gsplines::functions::FunctionBase &_in,
              std::size_t _n_glp, std::size_t _n_intervals);

  std::size_t get_nglp() const { return get_basis_dim(); }

protected:
  GaussLobattoLagrangeSpline *deriv_impl(std::size_t _deg = 1) const override;
};

double integral(::gsplines::functions::FunctionBase &_diffeo,
                ::gsplines::functions::FunctionBase &_path);

using GLLSpline = GaussLobattoLagrangeSpline;
} // namespace collocation
} // namespace gsplines

#endif /* GAUSSLOBATTOLAGRANGE_H */
