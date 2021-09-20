#ifndef GAUSSLOBATTOLAGRANGE_H
#define GAUSSLOBATTOLAGRANGE_H
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <gsplines/GSpline.hpp>

namespace gsplines {

namespace collocation {
class GaussLobattoLagrangeSpline
    : public ::gsplines::functions::FunctionInheritanceHelper<
          GaussLobattoLagrangeSpline, ::gsplines::GSpline,
          GaussLobattoLagrangeSpline> {
private:
public:
  GaussLobattoLagrangeSpline(
      std::pair<double, double> _domain, std::size_t _codom_dim,
      const Eigen::Ref<const Eigen::VectorXd> _coefficents, std::size_t _n_glp,
      std::size_t _n_intervals);
  GaussLobattoLagrangeSpline(const GaussLobattoLagrangeSpline &_that);
  GaussLobattoLagrangeSpline(GaussLobattoLagrangeSpline &&_that);
  GaussLobattoLagrangeSpline &
  operator=(const GaussLobattoLagrangeSpline &) = delete;
  virtual ~GaussLobattoLagrangeSpline() = default;

protected:
  GaussLobattoLagrangeSpline *deriv_impl(std::size_t _deg = 1) const override;
};
} // namespace collocation
} // namespace gsplines

#endif /* GAUSSLOBATTOLAGRANGE_H */
