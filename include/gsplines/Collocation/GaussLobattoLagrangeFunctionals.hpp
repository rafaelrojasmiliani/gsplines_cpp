#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/FunctionalAnalysis/Sobolev.hpp>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/variable_set.h>
#ifndef GAUSSLOBATTOLAGRANGEFUNCTIONALS_H
#define GAUSSLOBATTOLAGRANGEFUNCTIONALS_H

namespace gsplines {
namespace collocation {

template <typename T> class Functional {
private:
  std::size_t codom_dim_;

public:
  Functional(std::size_t codom_dim);
  virtual Eigen::VectorXd
  operator()(const GaussLobattoLagrangeSpline &_in) const = 0;
  virtual const T &derivative(const GaussLobattoLagrangeSpline &_in) const = 0;
  std::size_t get_codom_dim() { return codom_dim_; }
};

template <typename T> class LinearOperator : public Functional<T> {
protected:
  T mat_;

public:
  LinearOperator(const T &_args) : Functional<T>(_args.rows()), mat_(_args) {}
  LinearOperator(long _rows, long _cols)
      : Functional<T>(_rows), mat_(_rows, _cols) {}

  const T &to_matrix() const { return mat_; }

  template <typename M>
  inline auto operator*(const LinearOperator<M> &_rhs) const {
    return LinearOperator(to_matrix() * _rhs.to_matrix());
  }
  template <typename M>
  inline auto operator+(const LinearOperator<M> &_rhs) const {
    return LinearOperator(to_matrix() + _rhs.to_matrix());
  }

  template <typename M>
  inline auto operator*(const Eigen::MatrixBase<M> &_rhs) const {
    return to_matrix() * _rhs;
  }
  Eigen::VectorXd operator()(const GaussLobattoLagrangeSpline &_in) const {
    return ((*this) * _in).get_coefficients();
  }
  T derivative(const GaussLobattoLagrangeSpline & /*_in*/) { return mat_; }
  const T &derivative(const GaussLobattoLagrangeSpline & /*_in*/) const {
    return mat_;
  }

  GaussLobattoLagrangeSpline
  operator*(const GaussLobattoLagrangeSpline &_in) const {
    Eigen::VectorXd coeff = to_matrix() * _in.get_coefficients();
    std::size_t codom_dim = to_matrix().rows() / _in.get_number_of_intervals() /
                            _in.get_basis().get_dim();
    return GaussLobattoLagrangeSpline(
        _in.get_domain(), codom_dim, _in.get_number_of_intervals(),
        _in.get_basis().get_dim(), to_matrix() * _in.get_coefficients(),
        _in.get_interval_lengths());
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
    std::size_t total_size = _n_glp * _n_intervals * _codom_dim;
    // ---
    for (std::size_t uici = 0; uici < total_size; uici += _n_glp) {
      for (std::size_t i = 0; i < _n_glp; i++) {
        for (std::size_t j = 0; j < _n_glp; j++) {
          mat_.block(uici, uici, _n_glp, _n_glp).coeffRef(i, j) =
              2 * d_mat(i, j) / _int_lengs(0);
        }
      }
    }
    // ---
    mat_.makeCompressed();
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
            _that.get_basis().get_dim() * _that.get_number_of_intervals(),
            _that.get_basis().get_dim() * _that.get_number_of_intervals() *
                _that.get_codom_dim()) {

    std::size_t n_glp = _that.get_basis().get_dim();
    std::size_t n_inter = _that.get_number_of_intervals();
    std::size_t codom_dim = _that.get_codom_dim();

    const Eigen::VectorXd &vec = _that.get_coefficients();
    // ---
    for (std::size_t uic_interval = 0; uic_interval < n_inter; uic_interval++) {
      for (std::size_t uic_coor = 0; uic_coor < codom_dim; uic_coor++) {
        std::size_t i0 = uic_interval * n_glp;
        std::size_t j0 = uic_coor * n_glp + uic_interval * n_glp * codom_dim;

        for (std::size_t i = 0; i < n_glp; i++) {

          mat_.coeffRef(i0 + i, j0 + i) =
              vec(uic_interval * n_glp * codom_dim + uic_coor * n_glp + i);
        }
      }
    }
    mat_.makeCompressed();
  }
};

class ContinuityError
    : public LinearOperator<Eigen::SparseMatrix<double, Eigen::RowMajor>> {

public:
  ContinuityError(const GaussLobattoLagrangeSpline &_that, std::size_t _deg)
      : LinearOperator(_that.get_basis().continuity_matrix(
            _that.get_number_of_intervals(), _that.get_codom_dim(), _deg,
            _that.get_interval_lengths())) {}
};

class Integral : public Functional<Eigen::MatrixXd> {
  Eigen::MatrixXd glw_;

public:
  Integral(std::tuple<double, double> _domain, std::size_t _nglp,
           std::size_t _n_intervals);
  Eigen::VectorXd operator()(const GaussLobattoLagrangeSpline &_in) const;
  const Eigen::MatrixXd &
  derivative(const GaussLobattoLagrangeSpline &_in) const;
};
class SobolevError : public Functional<Eigen::MatrixXd> {
private:
  GaussLobattoLagrangeSpline approx_;
  GaussLobattoLagrangeSpline approx_d_;
  mutable Eigen::MatrixXd der_;
  Integral int_;

public:
  SobolevError(const gsplines::functions::FunctionBase &_fun, std::size_t _nglp,
               std::size_t _n_inter);
  Eigen::VectorXd operator()(const GaussLobattoLagrangeSpline &_in);
  const Eigen::MatrixXd &
  derivative(const GaussLobattoLagrangeSpline &_in) const;
};

class GLLSplineVariable : public ifopt::VariableSet,
                          public GaussLobattoLagrangeSpline {
private:
  ifopt::Component::Component::VecBound bounds_;

public:
  GLLSplineVariable(const GaussLobattoLagrangeSpline &_in);
  virtual ~GLLSplineVariable() = default;
  void SetVariables(const Eigen::VectorXd &_vec) override;
  Eigen::VectorXd GetValues() const override;
  ifopt::Component::VecBound GetBounds() const override { return bounds_; };
};

template <typename T> class CostWrapper : public ifopt::CostTerm {
private:
  std::unique_ptr<Functional<T>> fun_;

public:
  CostWrapper(Functional<T> _fun);
  double GetCost() const override;
  void FillJacobianBlock(std::string var_set, Jacobian &jac) const override;
};

template <typename T> class ConstraintWrapper : public ifopt::ConstraintSet {
private:
  std::unique_ptr<Functional<T>> fun_;
  ifopt::Component::Component::VecBound bounds_;

public:
  ConstraintWrapper(Functional<T> _fun, double _min, double _max);
  ConstraintWrapper(Functional<T> _fun, const std::vector<double> &_min,
                    double _max);
  ConstraintWrapper(Functional<T> _fun, double _min,
                    const std::vector<double> &_max);
  ConstraintWrapper(Functional<T> _fun, const std::vector<double> &_min,
                    const std::vector<double> &_max);

  ~ConstraintWrapper() = default;

  Eigen::VectorXd GetValues() const override;
  ifopt::Component::VecBound GetBounds() const override;

  void FillJacobianBlock(std::string _set_name,
                         Jacobian &_jac_block) const override;
};
} // namespace collocation
} // namespace gsplines
#endif /* GAUSSLOBATTOLAGRANGEFUNCTIONALS_H */
