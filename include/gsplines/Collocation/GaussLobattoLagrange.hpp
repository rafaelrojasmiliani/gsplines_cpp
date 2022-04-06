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
  approximate(const ::gsplines::functions::FunctionBase &_in,
              std::size_t _n_glp, std::size_t _n_intervals);

  std::size_t get_nglp() const { return get_basis_dim(); }

  template <typename T> class LinearOperator {
  protected:
    T mat_;

  public:
    template <typename... Ts> LinearOperator(Ts &&... _args) : mat_(_args...) {}

    const T &to_matrix() const { return mat_; }

    template <typename M> inline auto operator*(const LinearOperator<M> &_rhs) {
      return LinearOperator(to_matrix() * _rhs.to_matrix());
    }

    GaussLobattoLagrangeSpline
    operator*(const GaussLobattoLagrangeSpline &_in) {
      return GaussLobattoLagrangeSpline(_in.get_domain(), _in.get_codom_dim(),
                                        to_matrix() * _in.get_coefficients(),
                                        _in.get_basis().get_dim(),
                                        _in.get_intervals_num());
    }
  };

  class Derivative : public LinearOperator<Eigen::SparseMatrix<double>> {

  private:
    gsplines::basis::BasisLagrangeGaussLobatto basis_;

  public:
    Derivative(std::size_t _codom_dim, std::size_t _n_glp,
               std::size_t _n_intervals, const Eigen::VectorXd &_int_lengs)
        : LinearOperator(_n_glp * _n_intervals * _codom_dim,
                         _n_glp * _n_intervals * _codom_dim),
          basis_(_n_glp) {

      const Eigen::MatrixXd &d_mat = basis_.get_derivative_matrix_block();
      for (std::size_t uic_inter = 0; uic_inter < _n_intervals; uic_inter++) {
        for (std::size_t uic_dim = 0; uic_dim < _codom_dim; uic_dim++) {
          for (std::size_t i = 0; i < _n_glp; i++) {
            for (std::size_t j = 0; j < _n_glp; j++) {
              mat_.block(uic_inter * uic_dim * _n_glp,
                         uic_inter * uic_dim * _n_glp, _n_glp, _n_glp)
                  .coeffRef(i, j) = d_mat(i, j) / _int_lengs(uic_inter);
            }
          }
        }
      }
    }
    Derivative(const GaussLobattoLagrangeSpline &_that)
        : Derivative(_that.get_codom_dim(), _that.get_basis().get_dim(),
                     _that.get_number_of_intervals(),
                     _that.get_interval_lengths()) {}
  };

  class TransposeLeftMultiplication
      : public LinearOperator<Eigen::SparseMatrix<double>> {

  public:
    TransposeLeftMultiplication(const GaussLobattoLagrangeSpline &_that)
        : LinearOperator(
              _that.get_basis().get_dim() * _that.get_intervals_num(),
              _that.get_basis().get_dim() * _that.get_intervals_num() *
                  _that.get_codom_dim()) {
      const Eigen::VectorXd &vec = _that.get_coefficients();
      std::size_t nglp = _that.get_basis().get_dim();
      for (std::size_t i = 0;
           i < _that.get_number_of_intervals() * _that.get_basis().get_dim();
           i++) {
        for (std::size_t j = 0; j < _that.get_codom_dim(); j++) {
          mat_.coeffRef(i, j * nglp) = vec(j * nglp);
        }
      }
    }
  };

protected:
  GaussLobattoLagrangeSpline *deriv_impl(std::size_t _deg = 1) const override;
};

double integral(::gsplines::functions::FunctionBase &_diffeo,
                ::gsplines::functions::FunctionBase &_path);

using GLLSpline = GaussLobattoLagrangeSpline;
} // namespace collocation
} // namespace gsplines

#endif /* GAUSSLOBATTOLAGRANGE_H */
