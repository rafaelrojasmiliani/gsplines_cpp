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
  Eigen::MatrixXd value_at_nodes_;

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
  const Eigen::MatrixXd &value_at_nodes() const { return value_at_nodes_; };

  static GaussLobattoLagrangeSpline identity(std::pair<double, double> _domain,
                                             std::size_t _n_glp,
                                             std::size_t _n_intervals);

  static GaussLobattoLagrangeSpline
  approximate(::gsplines::functions::FunctionBase &_in, std::size_t _n_glp,
              std::size_t _n_intervals);

protected:
  GaussLobattoLagrangeSpline *deriv_impl(std::size_t _deg = 1) const override;
};

double integral(::gsplines::functions::FunctionBase &_diffeo,
                ::gsplines::functions::FunctionBase &_path);

} // namespace collocation
} // namespace gsplines

#endif /* GAUSSLOBATTOLAGRANGE_H */
