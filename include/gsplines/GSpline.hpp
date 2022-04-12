#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>

namespace gsplines {
class SobolevNorm;
class GSpline
    : public functions::FunctionInheritanceHelper<GSpline, functions::Function,
                                                  GSpline> {
  friend SobolevNorm;
  friend GSpline operator*(double _a, const GSpline &_in);
  friend GSpline &operator*(double _a, GSpline &&_in);

protected:
  Eigen::VectorXd coefficients_;
  Eigen::VectorXd domain_interval_lengths_;

private:
  std::unique_ptr<basis::Basis> basis_;
  mutable Eigen::VectorXd basis_buffer_;
  double interval_to_window(double _t, std::size_t _interval) const;

  Eigen::Ref<const Eigen::VectorXd>
  coefficient_segment(std::size_t _interval, std::size_t _component) const;

  Eigen::Ref<Eigen::VectorXd> coefficient_segment(std::size_t _interval,
                                                  std::size_t _component);

  std::size_t get_interval(double _domain_point) const;

public:
  GSpline(std::pair<double, double> _domain, std::size_t _codom_dim,
          std::size_t _n_intervals, const basis::Basis &_basis,
          const Eigen::Ref<const Eigen::VectorXd> _coefficents,
          const Eigen::Ref<const Eigen::VectorXd> _tauv,
          const std::string &_name = "GSpline");

  GSpline(const GSpline &that);
  GSpline(GSpline &&that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

  Eigen::VectorXd get_domain_breakpoints() const;

  Eigen::MatrixXd get_waypoints() const;

  virtual ~GSpline() = default;
  const Eigen::VectorXd &get_coefficients() const { return coefficients_; }

  const Eigen::VectorXd &get_interval_lengths() const {
    return domain_interval_lengths_;
  }

  const std::string &get_basis_name() const { return basis_->get_name(); }

  std::size_t get_number_of_intervals() const {
    return domain_interval_lengths_.size();
  }

  std::size_t get_basis_dim() const { return basis_->get_dim(); }

  GSpline linear_scaling_new_execution_time(double _new_exec_time) const;

  const basis::Basis &get_basis() const { return *basis_; }

  bool same_vector_space(const GSpline &_that) const;

  GSpline operator+(const GSpline &that) const &;
  GSpline operator+(GSpline &&that) const &;
  GSpline operator+(const GSpline &that) &&;
  GSpline operator+(GSpline &&that) &&;

  GSpline operator-(const GSpline &that) const &;
  GSpline operator-(GSpline &&that) const &;
  GSpline operator-(const GSpline &that) &&;
  GSpline operator-(GSpline &&that) &&;

  GSpline &operator+=(const GSpline &that);
  GSpline &operator-=(const GSpline &that);
  // GSpline operator-() const &;
  // GSpline operator-() &&;

protected:
  virtual GSpline *deriv_impl(std::size_t _deg = 1) const override;
};

const Eigen::Ref<const Eigen::VectorXd>
get_coefficient_segment(const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                        basis::Basis &_basis, std::size_t _num_interval,
                        std::size_t _codom_dim, std::size_t _interval,
                        std::size_t _component);

GSpline operator*(double _a, const GSpline &_in);
GSpline &operator*(double _a, GSpline &&_in);
} // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
