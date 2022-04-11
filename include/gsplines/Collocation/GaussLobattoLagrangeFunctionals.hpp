#include <gsplines/Collocation/GaussLobattoLagrange.hpp>
#include <gsplines/FunctionalAnalysis/Sobolev.hpp>
#include <ifopt/constraint_set.h>
#include <ifopt/cost_term.h>
#include <ifopt/variable_set.h>
#ifndef GAUSSLOBATTOLAGRANGEFUNCTIONALS_H
#define GAUSSLOBATTOLAGRANGEFUNCTIONALS_H

namespace gsplines {
namespace collocation {

class LinearFunctionalBase;
class Functional {
private:
  std::size_t codom_dim_;

public:
  Functional(std::size_t _codom_dim) : codom_dim_(_codom_dim){};

  virtual Eigen::VectorXd
  operator()(const GaussLobattoLagrangeSpline &_in) const = 0;

  virtual std::shared_ptr<LinearFunctionalBase>
  differential(const GaussLobattoLagrangeSpline &_in) const = 0;

  std::size_t get_codom_dim() const { return codom_dim_; }

  virtual ~Functional() = default;
};

class LinearFunctionalBase : public Functional {
public:
  using Functional::Functional;

  LinearFunctionalBase(const LinearFunctionalBase &_that) = default;

  virtual Eigen::MatrixXd to_dense_matrix() const = 0;

  virtual Eigen::SparseMatrix<double> to_sparse_matrix() const = 0;

  ~LinearFunctionalBase() = default;
};

template <typename T> class LinearFunctional : public LinearFunctionalBase {
  template <typename M>
  friend LinearFunctional<M> operator*(double, const LinearFunctional<M> &);
  template <typename M>
  friend LinearFunctional<M> &&operator*(double, LinearFunctional<M> &&);

protected:
  T mat_;

public:
  LinearFunctional(const T &_args)
      : LinearFunctionalBase(_args.rows()), mat_(_args) {}

  LinearFunctional(T &&_args)
      : LinearFunctionalBase(_args.rows()), mat_(std::move(_args)) {}

  LinearFunctional(const LinearFunctional &_that)
      : LinearFunctionalBase(_that.get_codom_dim()), mat_(_that.mat_) {}

  LinearFunctional(LinearFunctional &&_that)
      : LinearFunctionalBase(_that.get_codom_dim()),
        mat_(std::move(_that.mat_)) {}

  LinearFunctional(long _rows, long _cols)
      : LinearFunctionalBase((std::size_t)_rows), mat_(_rows, _cols) {}

  const T &to_matrix() const { return mat_; }

  LinearFunctional &operator=(const T &_mat) {
    mat_ = _mat;
    return *this;
  }

  template <typename M>
  inline auto operator*(const LinearFunctional<M> &_rhs) const {
    return LinearFunctional(to_matrix() * _rhs.to_matrix());
  }

  template <typename M>
  inline auto operator+(const LinearFunctional<M> &_rhs) const {
    return LinearFunctional(to_matrix() + _rhs.to_matrix());
  }

  template <typename M>
  inline auto operator*(const Eigen::MatrixBase<M> &_rhs) const {
    return to_matrix() * _rhs;
  }
  Eigen::VectorXd
  operator()(const GaussLobattoLagrangeSpline &_in) const override {
    return mat_ * _in.get_coefficients();
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

  Eigen::MatrixXd to_dense_matrix() const override {
    if constexpr (std::is_same_v<Eigen::SparseMatrix<double>, T>) {
      return mat_.toDense();
    } else if constexpr (std::is_same_v<Eigen::MatrixXd, T>) {
      return mat_;
    }
    return Eigen::MatrixXd();
  }
  Eigen::SparseMatrix<double> to_sparse_matrix() const override {

    if constexpr (std::is_same_v<Eigen::MatrixXd, T>) {
      return mat_.sparseView();
    } else if constexpr (std::is_same_v<Eigen::MatrixXd, T>) {
      return mat_;
    }
    return Eigen::SparseMatrix<double>();
  }

  std::shared_ptr<LinearFunctionalBase>
  differential(const GaussLobattoLagrangeSpline & /*_in*/) const override {

    return std::make_shared<LinearFunctional<T>>(*this);
  }

  ~LinearFunctional() = default;
};

template <typename T>
LinearFunctional<T> operator*(double _num, const LinearFunctional<T> &_lm) {
  LinearFunctional result(_lm);
  result.mat_ *= _num;
  return result;
}

template <typename T>
LinearFunctional<T> &&operator*(double _num, LinearFunctional<T> &&_lm) {
  _lm.mat_ *= _num;
  return std::move(_lm);
}

class LinearFunctionalSparse
    : public LinearFunctional<Eigen::SparseMatrix<double>> {
public:
  using LinearFunctional<Eigen::SparseMatrix<double>>::LinearFunctional;
  ~LinearFunctionalSparse() = default;
};

class LinearFunctionalDense : public LinearFunctional<Eigen::MatrixXd> {
public:
  using LinearFunctional<Eigen::MatrixXd>::LinearFunctional;
  ~LinearFunctionalDense() = default;
};

class Derivative : public LinearFunctionalSparse {

private:
  gsplines::basis::BasisLagrangeGaussLobatto basis_;

public:
  Derivative(std::size_t _codom_dim, std::size_t _n_glp,
             std::size_t _n_intervals, const Eigen::VectorXd &_int_lengs)
      : LinearFunctionalSparse(_n_glp * _n_intervals * _codom_dim,
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

class TransposeLeftMultiplication : public LinearFunctionalSparse {

public:
  TransposeLeftMultiplication(const GaussLobattoLagrangeSpline &_that)
      : LinearFunctionalSparse(
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
    : public LinearFunctional<Eigen::SparseMatrix<double, Eigen::RowMajor>> {

public:
  ContinuityError(const GaussLobattoLagrangeSpline &_that, std::size_t _deg)
      : LinearFunctional(_that.get_basis().continuity_matrix(
            _that.get_number_of_intervals(), _that.get_codom_dim(), _deg,
            _that.get_interval_lengths())) {}
};

class Integral : public LinearFunctionalDense {
  Eigen::MatrixXd glw_;

public:
  Integral(std::tuple<double, double> _domain, std::size_t _nglp,
           std::size_t _n_intervals);
  Integral(const GaussLobattoLagrangeSpline &_in);
};

class SobolevDistance : public Functional {
private:
  GaussLobattoLagrangeSpline approx_;
  GaussLobattoLagrangeSpline approx_d_;
  Integral int_;
  Derivative der_;

  LinearFunctionalDense diff_;

public:
  SobolevDistance(const gsplines::functions::FunctionBase &_fun,
                  std::size_t _nglp, std::size_t _n_inter, std::size_t _deg);
  Eigen::VectorXd
  operator()(const GaussLobattoLagrangeSpline &_in) const override;

  std::shared_ptr<LinearFunctionalBase>
  differential(const GaussLobattoLagrangeSpline & /*_in*/) const override;
};

/*

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

class __Counter {
protected:
  static int counter__;
  __Counter() { counter__++; }
};

template <typename T>
class CostWrapper : public __Counter, public T, public ifopt::CostTerm {

public:
  template <typename... Ts>
  CostWrapper(Ts &&... args)
      : __Counter(), T(args...), CostTerm("cost_" + std::to_string(counter__)) {
  }
  double GetCost() const override {
    std::shared_ptr<GaussLobattoLagrangeSpline> spline =
        std::dynamic_pointer_cast<GaussLobattoLagrangeSpline>(
            GetVariables()->GetComponent("GLL"));
    if (not spline)
      throw std::logic_error("could not cast GLL into GLLSpline");
    Eigen::VectorXd result = (*this)(*spline);
    return result(0);
  }
  void FillJacobianBlock(std::string var_set,
                         Jacobian &jac) const override {
    std::shared_ptr<GaussLobattoLagrangeSpline> spline =
        std::dynamic_pointer_cast<GaussLobattoLagrangeSpline>(
            GetVariables()->GetComponent("GLL"));
    if (not spline)
      throw std::logic_error("could not cast GLL into GLLSpline");
    jac = this->derivate(*spline);
  }
};

template <typename T>
class ConstraintWrapper : public __Counter,
                          public T,
                          public ifopt::ConstraintSet {
private:
  ifopt::Component::Component::VecBound bounds_;
  static int counter__;

public:
  template <typename... Ts>
  ConstraintWrapper(double _min, double _max, Ts &&... args)
      : __Counter(), T(args...),
        ConstraintSet(T::get_codom_dim(),
                      "Wrapper_" + std::to_string(counter__)) {

    counter__++;
    ifopt::Bounds default_bound(_min, _max);
    bounds_ = ifopt::Component::VecBound(GetRows(), default_bound);
  }

  Eigen::VectorXd GetValues() const override {

    std::shared_ptr<GaussLobattoLagrangeSpline> spline =
        std::dynamic_pointer_cast<GaussLobattoLagrangeSpline>(
            GetVariables()->GetComponent("GLL"));
    if (not spline)
      throw std::logic_error("could not cast GLL into GLLSpline");

    return T::operator()(*spline);
  }

  ifopt::Component::VecBound GetBounds() const override { return bounds_; }

  void FillJacobianBlock(std::string _set_name,
                         Jacobian &_jac_block) const override {
    std::shared_ptr<GaussLobattoLagrangeSpline> spline =
        std::dynamic_pointer_cast<GaussLobattoLagrangeSpline>(
            GetVariables()->GetComponent("GLL"));
    if (_set_name == "GLL")
      _jac_block = T::differential(*spline);
  }
};*/
} // namespace collocation
} // namespace gsplines
#endif /* GAUSSLOBATTOLAGRANGEFUNCTIONALS_H */
