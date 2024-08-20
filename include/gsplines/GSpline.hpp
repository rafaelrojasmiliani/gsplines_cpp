#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include "gsplines/Functions/FunctionExpression.hpp"
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

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
  double interval_to_window(double _domain_point, std::size_t _interval) const;

  Eigen::Ref<const Eigen::VectorXd> coefficient_segment(
      std::size_t _interval, std::size_t _component) const;

  Eigen::Ref<Eigen::VectorXd> coefficient_segment(std::size_t _interval,
                                                  std::size_t _component);

  std::size_t get_interval(double _domain_point) const;

 public:
  GSplineBase(std::pair<double, double> _domain, std::size_t _codom_dim,
              std::size_t _n_intervals, const basis::Basis& _basis,
              Eigen::Ref<const Eigen::VectorXd> _coefficents,
              Eigen::Ref<const Eigen::VectorXd> _tauv,
              const std::string& _name = "GSplineBase");

  GSplineBase(std::pair<double, double> _domain, std::size_t _codom_dim,
              std::size_t _n_intervals, const basis::Basis& _basis,
              Eigen::VectorXd&& _coefficents, Eigen::VectorXd&& _tauv,
              const std::string& _name = "GSplineBase");

  GSplineBase(const GSplineBase& that);
  GSplineBase(GSplineBase&& that) noexcept;

 public:
  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

  Eigen::VectorXd get_domain_breakpoints() const;

  Eigen::MatrixXd get_waypoints() const;

  virtual ~GSplineBase() = default;
  const Eigen::VectorXd& get_coefficients() const { return coefficients_; }

  const Eigen::VectorXd& get_interval_lengths() const {
    return domain_interval_lengths_;
  }

  const std::string& get_basis_name() const { return basis_->get_name(); }

  std::size_t get_number_of_intervals() const {
    return domain_interval_lengths_.size();
  }

  std::size_t get_basis_dim() const { return basis_->get_dim(); }

  const basis::Basis& get_basis() const { return *basis_; }

  bool same_vector_space(const GSplineBase& _that) const;

  bool operator==(const GSplineBase& _that) const;

  bool operator!=(const GSplineBase& _that) const;

 protected:
  GSplineBase* deriv_impl(std::size_t unused_degree) const override {
    (void)unused_degree;
    throw std::runtime_error(
        "GSplineBase class has not a default derivative calculator.");
    return nullptr;
  }
};

Eigen::Ref<const Eigen::VectorXd> get_coefficient_segment(
    Eigen::Ref<const Eigen::VectorXd> _coefficients, basis::Basis& _basis,
    std::size_t _num_interval, std::size_t _codom_dim, std::size_t _interval,
    std::size_t _component);

template <typename Current, typename Base = GSplineBase>
class GSplineInheritanceHelper
    : public functions::FunctionInheritanceHelper<Current, Base, Current> {
 protected:
  Base* deriv_impl(std::size_t _deg) const override {
    Eigen::VectorXd result_coeff(Base::coefficients_);
    std::size_t interval_coor = 0;
    std::size_t codom_coor = 0;
    std::size_t i0 = 0;

    if (_deg > 0) {
      for (interval_coor = 0; interval_coor < Base::get_number_of_intervals();
           interval_coor++) {
        for (codom_coor = 0; codom_coor < Base::get_codom_dim(); codom_coor++) {
          i0 = interval_coor * Base::basis_->get_dim() * Base::get_codom_dim() +
               Base::basis_->get_dim() * codom_coor;

          double scale = 1.0;
          /// workaround to handle Basis0101, which derivative matrix depends
          /// on tau
          if (dynamic_cast<const basis::Basis0101*>(this->basis_.get()) !=
              nullptr) {
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

    return new Current(Base::get_domain(), Base::get_codom_dim(),  // NOLINT
                       Base::get_number_of_intervals(), *Base::basis_,
                       result_coeff, Base::domain_interval_lengths_,
                       Base::get_name());
  }

 public:
  // Import constructor of Basis
  using functions::FunctionInheritanceHelper<
      Current, Base, Current>::FunctionInheritanceHelper;

  Current operator+(const Current& that) const& {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    Current result = Current(*this);
    result.coefficients_ += that.coefficients_;
    return result;
  }

  Current operator+(Current&& that) const& {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    that.coefficients_ += Base::coefficients_;
    return std::move(that);
  }

  Current operator+(const Current& that) && {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    Base::coefficients_ += that.coefficients_;
    return std::move(*this);
  }

  Current operator+(Current&& that) && {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    Base::coefficients_ += that.coefficients_;
    return std::move(*this);
  }

  Current operator-(const Current& that) const& {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    GSplineBase result = GSplineBase(*this);
    result.coefficients_ -= that.coefficients_;
    return result;
  }

  Current operator-(Current&& that) const& {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    that.coefficients_ = Base::coefficients_ - that.coefficients_;
    return std::move(that);
  }

  Current operator-(const Current& that) && {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    Base::coefficients_ -= that.coefficients_;
    return std::move(*this);
  }

  Current operator-(Current&& that) && {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    Base::coefficients_ -= that.coefficients_;
    return std::move(*this);
  }

  Current& operator+=(const Current& that) {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }
    Base::coefficients_ += that.coefficients_;
    return *this;
  }

  Current& operator-=(const Current& that) {
    if (not Base::same_vector_space(that)) {
      throw std::invalid_argument("Cannot sum Incompatible Gspline");
    }

    Base::coefficients_ -= that.coefficients_;
    return *this;
  }

  Current linear_scaling_new_execution_time(double _new_exec_time) const {
    if (_new_exec_time < 1.0e-8) {
      throw std::invalid_argument("New time cannot be negative");
    }

    double time_scale_factor = _new_exec_time / Base::get_domain_length();

    Eigen::VectorXd new_domain_interva_lengths =
        Base::domain_interval_lengths_ * time_scale_factor;

    std::pair<double, double> new_domain = Base::get_domain();
    new_domain.second = new_domain.first + _new_exec_time;

    auto* b = dynamic_cast<basis::Basis0101*>(Base::basis_.get());
    if (b != nullptr) {
      const double alpha_0 = b->get_alpha();
      const double k_0 = std::sqrt(2.0) / 4.0 * std::pow(alpha_0, 0.25) /
                         std::pow((1.0 - alpha_0), 0.25);
      const double k = k_0 / time_scale_factor;
      const double k4 = std::pow(k, 4);

      const double alpha = k4 / (1.0 + k4);
      return Current(new_domain, Base::get_codom_dim(),
                     Base::get_number_of_intervals(), basis::Basis0101(alpha),
                     Base::coefficients_, new_domain_interva_lengths,
                     Base::get_name());
    }
    return Current(new_domain, Base::get_codom_dim(),
                   Base::get_number_of_intervals(), *Base::basis_,
                   Base::coefficients_, new_domain_interva_lengths,
                   Base::get_name());
  }

  Current linear_scaling_new_execution_time_max_velocity_max_acceleration(
      std::optional<double> _velocity_bound,
      std::optional<double> _acceleration_bound, double _dt = 0.01) const {
    std::optional<std::vector<double>> velocity_bound;
    std::optional<std::vector<double>> acceleration_bound;
    if (_velocity_bound.has_value()) {
      velocity_bound =
          std::vector<double>(this->get_codom_dim(), _velocity_bound.value());
    }
    if (_acceleration_bound.has_value()) {
      acceleration_bound = std::vector<double>(this->get_codom_dim(),
                                               _acceleration_bound.value());
    }
    return this
        ->linear_scaling_new_execution_time_max_velocity_max_acceleration(
            std::move(velocity_bound), std::move(acceleration_bound), _dt);
  }

  Current linear_scaling_new_execution_time_max_velocity_max_acceleration(
      std::optional<std::vector<double>> _velocity_bound,
      std::optional<std::vector<double>> _acceleration_bound,
      double _dt = 0.01) const {
    if (!_velocity_bound.has_value() && !_acceleration_bound.has_value()) {
      throw std::invalid_argument(
          "This method need a may value of velocity or acceleration");
    }
    if (_dt < 1.0e-8) {
      throw std::invalid_argument("dt cannot be negative or too small");
    }
    if (_velocity_bound.value().size() != this->get_codom_dim() ||
        _acceleration_bound.value().size() != this->get_codom_dim()) {
      throw std::invalid_argument(
          "Bounds for velocity and acceleration must correspont to codom dim");
    }

    std::size_t number_of_segments = this->get_domain_length() / _dt;
    const double t0 = this->get_domain().first;
    const double t1 = this->get_domain().second;
    const Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
        static_cast<long>(number_of_segments) + 1, t0, t1);

    const gsplines::functions::FunctionExpression gspline_diff_1 =
        this->derivate();
    const gsplines::functions::FunctionExpression gspline_diff_2 =
        gspline_diff_1.derivate();

    Eigen::MatrixXd gspline_diff_1_evaluated = gspline_diff_1(time_spam);
    Eigen::MatrixXd gspline_diff_2_evaluated = gspline_diff_2(time_spam);

    const Eigen::Map<const Eigen::VectorXd> velocity_bound(
        _velocity_bound.value().data(),
        static_cast<long>(_velocity_bound.value().size()));
    const Eigen::Map<const Eigen::VectorXd> acceleration_bound(
        _acceleration_bound.value().data(),
        static_cast<long>(_acceleration_bound.value().size()));

    // d / d t {f (\alpha t)} = df/dt \alpha
    // the best value of alpha is
    //
    // {df/dt}_max \alpha = v_bound
    //
    // alpha =  v_bound / {df/dt}_max
    //
    const double max_velocity_ratio =
        (velocity_bound.transpose().array() /
         gspline_diff_1_evaluated.array().abs().colwise().maxCoeff())
            .maxCoeff();

    const double max_acceleration_ratio =
        Eigen::sqrt(acceleration_bound.transpose().array() /
                    gspline_diff_2_evaluated.array().abs().colwise().maxCoeff())
            .maxCoeff();

    const double time_scale_factor =
        std::max(max_velocity_ratio, max_acceleration_ratio);

    Eigen::VectorXd new_domain_interva_lengths =
        Base::domain_interval_lengths_ * time_scale_factor;

    std::pair<double, double> new_domain = Base::get_domain();
    new_domain.second =
        new_domain.first + new_domain_interva_lengths.array().sum();

    /// work arround, as the basis0101 depens on tau, if we scale,
    /// we change the basis.
    auto* b = dynamic_cast<basis::Basis0101*>(Base::basis_.get());
    if (b != nullptr) {
      const double alpha_0 = b->get_alpha();
      const double k_0 = std::sqrt(2.0) / 4.0 * std::pow(alpha_0, 0.25) /
                         std::pow((1.0 - alpha_0), 0.25);
      const double k = k_0 / time_scale_factor;
      const double k4 = std::pow(k, 4);

      const double alpha = k4 / (1.0 + k4);
      return Current(new_domain, Base::get_codom_dim(),
                     Base::get_number_of_intervals(), basis::Basis0101(alpha),
                     Base::coefficients_, new_domain_interva_lengths,
                     Base::get_name());
    }
    return Current(new_domain, Base::get_codom_dim(),
                   Base::get_number_of_intervals(), *Base::basis_,
                   Base::coefficients_, new_domain_interva_lengths,
                   Base::get_name());
  }
  ~GSplineInheritanceHelper() override = default;
};

class GSpline : public GSplineInheritanceHelper<GSpline> {
  friend GSpline operator*(double _a, const GSpline& _that);
  friend GSpline operator*(double _a, GSpline&& _that);
  friend GSpline operator*(const GSpline& _that, double _a);
  friend GSpline operator*(GSpline&& _that, double _a);
  friend GSpline operator-(const GSpline& _that);
  friend GSpline operator-(GSpline&& _that);

 public:
  using GSplineInheritanceHelper<GSpline>::GSplineInheritanceHelper;
  GSpline(const GSpline& _that) : GSplineInheritanceHelper(_that) {}
  GSpline(GSpline&& _that) : GSplineInheritanceHelper(std::move(_that)) {}
};

GSpline operator*(double _a, const GSpline& _that);
GSpline operator*(double _a, GSpline&& _that);
GSpline operator*(const GSpline& _that, double _a);
GSpline operator*(GSpline&& _that, double _a);
GSpline operator-(const GSpline& _that);
GSpline operator-(GSpline&& _that);

GSpline random_gspline(std::pair<double, double> _domain,
                       std::size_t _codom_dim);

GSpline random_gspline(std::pair<double, double> _domain,
                       std::size_t _codom_dim, const basis::Basis& _basis);

}  // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
