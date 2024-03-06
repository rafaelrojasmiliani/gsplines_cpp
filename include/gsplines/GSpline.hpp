#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include "gsplines/Functions/FunctionExpression.hpp"
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <stdexcept>

namespace gsplines {
class SobolevNorm;

class GSplineBase : public functions::FunctionInheritanceHelper<
                        GSplineBase, functions::Function, GSplineBase> {
  friend SobolevNorm;
  template <typename Current, typename Base>
  friend class GSplineInheritanceHelper;

protected:
  Eigen::VectorXd coefficients_;
  Eigen::VectorXd domain_interval_lengths_;

  std::unique_ptr<basis::Basis> basis_;
  mutable Eigen::VectorXd basis_buffer_;
  double interval_to_window(double _t, std::size_t _interval) const;

  Eigen::Ref<const Eigen::VectorXd>
  coefficient_segment(std::size_t _interval, std::size_t _component) const;

  Eigen::Ref<Eigen::VectorXd> coefficient_segment(std::size_t _interval,
                                                  std::size_t _component);

  std::size_t get_interval(double _domain_point) const;

public:
  GSplineBase(std::pair<double, double> _domain, std::size_t _codom_dim,
              std::size_t _n_intervals, const basis::Basis &_basis,
              const Eigen::Ref<const Eigen::VectorXd> _coefficents,
              const Eigen::Ref<const Eigen::VectorXd> _tauv,
              const std::string &_name = "GSplineBase");

  GSplineBase(std::pair<double, double> _domain, std::size_t _codom_dim,
              std::size_t _n_intervals, const basis::Basis &_basis,
              Eigen::VectorXd &&_coefficents, Eigen::VectorXd &&_tauv,
              const std::string &_name = "GSplineBase");

  GSplineBase(const GSplineBase &that);
  GSplineBase(GSplineBase &&that);

public:
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

  Eigen::VectorXd get_domain_breakpoints() const;

  Eigen::MatrixXd get_waypoints() const;

  virtual ~GSplineBase() = default;
  const Eigen::VectorXd &get_coefficients() const { return coefficients_; }

  const Eigen::VectorXd &get_interval_lengths() const {
    return domain_interval_lengths_;
  }

  const std::string &get_basis_name() const { return basis_->get_name(); }

  std::size_t get_number_of_intervals() const {
    return domain_interval_lengths_.size();
  }

  std::size_t get_basis_dim() const { return basis_->get_dim(); }

  const basis::Basis &get_basis() const { return *basis_; }

  bool same_vector_space(const GSplineBase &_that) const;

  bool operator==(const GSplineBase &_that) const;

  bool operator!=(const GSplineBase &_that) const;

protected:
  GSplineBase *deriv_impl(std::size_t) const override {
    throw std::runtime_error("Base class cannot be used");
    return nullptr;
  }
};

const Eigen::Ref<const Eigen::VectorXd>
get_coefficient_segment(const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                        basis::Basis &_basis, std::size_t _num_interval,
                        std::size_t _codom_dim, std::size_t _interval,
                        std::size_t _component);

template <typename Current, typename Base = GSplineBase>
class GSplineInheritanceHelper
    : public functions::FunctionInheritanceHelper<Current, Base, Current> {

protected:
  Base *deriv_impl(std::size_t _deg) const {

    Eigen::VectorXd result_coeff(Base::coefficients_);
    std::size_t interval_coor;
    std::size_t codom_coor;
    std::size_t i0;

    if (_deg > 0) {
      for (interval_coor = 0; interval_coor < Base::get_number_of_intervals();
           interval_coor++) {
        for (codom_coor = 0; codom_coor < Base::get_codom_dim(); codom_coor++) {
          i0 = interval_coor * Base::basis_->get_dim() * Base::get_codom_dim() +
               Base::basis_->get_dim() * codom_coor;

          double scale = 1.0;
          /// workaround to handle Basis0101, which derivative matrix depends
          /// on tau
          if (dynamic_cast<const basis::Basis0101 *>(this->basis_.get())) {
            scale = (_deg == 0) ? 0.0 : 1.0;
          } else {
            scale = std::pow(
                2.0 / Base::domain_interval_lengths_(interval_coor), _deg);
          }

          result_coeff.segment(i0, Base::basis_->get_dim()) =
              Base::basis_->get_derivative_matrix_block(_deg) *
              result_coeff.segment(i0, Base::basis_->get_dim()) * scale;
        }
      }
    }

    return new Current(Base::get_domain(), Base::get_codom_dim(),
                       Base::get_number_of_intervals(), *Base::basis_,
                       result_coeff, Base::domain_interval_lengths_,
                       Base::get_name());
  }

public:
  // Import constructor of Basis
  using functions::FunctionInheritanceHelper<
      Current, Base, Current>::FunctionInheritanceHelper;

  Current operator+(const Current &that) const & {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    Current result = Current(*this);
    result.coefficients_ += that.coefficients_;
    return result;
  }

  Current operator+(Current &&that) const & {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    that.coefficients_ += Base::coefficients_;
    return std::move(that);
  }

  Current operator+(const Current &that) && {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    Base::coefficients_ += that.coefficients_;
    return std::move(*this);
  }

  Current operator+(Current &&that) && {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    Base::coefficients_ += that.coefficients_;
    return std::move(*this);
  }

  Current operator-(const Current &that) const & {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    GSplineBase result = GSplineBase(*this);
    result.coefficients_ -= that.coefficients_;
    return result;
  }

  Current operator-(Current &&that) const & {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    that.coefficients_ = Base::coefficients_ - that.coefficients_;
    return std::move(that);
  }

  Current operator-(const Current &that) && {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    Base::coefficients_ -= that.coefficients_;
    return std::move(*this);
  }

  Current operator-(Current &&that) && {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    Base::coefficients_ -= that.coefficients_;
    return std::move(*this);
  }

  Current &operator+=(const Current &that) {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    Base::coefficients_ += that.coefficients_;
    return *this;
  }

  Current &operator-=(const Current &that) {
    if (not Base::same_vector_space(that))
      throw std::invalid_argument("Cannot sum Incompatible Gspline");

    Base::coefficients_ -= that.coefficients_;
    return *this;
  }

  Current linear_scaling_new_execution_time(double _new_exec_time) const {

    if (_new_exec_time < 1.0e-8) {
      throw std::invalid_argument("New time cannot be negative");
    }

    Eigen::VectorXd new_domain_interva_lengths =
        Base::domain_interval_lengths_ * _new_exec_time /
        Base::get_domain_length();

    std::pair<double, double> new_domain = Base::get_domain();
    new_domain.second = new_domain.first + _new_exec_time;

    return Current(new_domain, Base::get_codom_dim(),
                   Base::get_number_of_intervals(), *Base::basis_,
                   Base::coefficients_, new_domain_interva_lengths,
                   Base::get_name());
  }
  virtual ~GSplineInheritanceHelper() = default;
};

class GSpline : public GSplineInheritanceHelper<GSpline> {
  friend GSpline operator*(double _a, const GSpline &_that);
  friend GSpline operator*(double _a, GSpline &&_that);
  friend GSpline operator*(const GSpline &_that, double _a);
  friend GSpline operator*(GSpline &&_that, double _a);
  friend GSpline operator-(const GSpline &_that);
  friend GSpline operator-(GSpline &&_that);

public:
  using GSplineInheritanceHelper<GSpline>::GSplineInheritanceHelper;
  GSpline(const GSpline &_that) : GSplineInheritanceHelper(_that) {}
  GSpline(GSpline &&_that) : GSplineInheritanceHelper(std::move(_that)) {}
};

GSpline operator*(double _a, const GSpline &_that);
GSpline operator*(double _a, GSpline &&_that);
GSpline operator*(const GSpline &_that, double _a);
GSpline operator*(GSpline &&_that, double _a);
GSpline operator-(const GSpline &_that);
GSpline operator-(GSpline &&_that);

GSpline random_gspline(std::pair<double, double> _domain,
                       std::size_t _codom_dim);

} // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
