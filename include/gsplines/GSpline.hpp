#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include <eigen3/Eigen/Core>
#include <gsplines/Basis.hpp>
#include <gsplines/Functions/Function.hpp>

namespace gsplines {
class SobolevNorm;
class GSpline : public functions::Function {
  friend SobolevNorm;

private:
  GSpline &operator=(const GSpline &);
  const std::size_t number_of_intervals_;
  std::unique_ptr<basis::Basis> basis_;
  Eigen::VectorXd coefficients_;
  Eigen::VectorXd domain_break_points_;
  Eigen::VectorXd domain_interval_lengths_;
  Eigen::MatrixXd waypoints_;
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
          const std::string &_name = "PieceWiseFunction");

  GSpline(const GSpline &that);
  GSpline(GSpline &&that);

  std::unique_ptr<FunctionExpression> deriv(int _deg = 1) const override;

  std::unique_ptr<FunctionExpression> clone() const override;

  std::unique_ptr<FunctionExpression> move_clone() override;

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;

  std::size_t get_intervals_num() { return number_of_intervals_; }
  double get_exec_time() {
    return domain_break_points_.tail(1)(0) - domain_break_points_(0);
  }
  const Eigen::Ref<const Eigen::VectorXd> get_domain_breakpoints() {
    return domain_break_points_;
  }
  const Eigen::Ref<const Eigen::MatrixXd> get_waypoints() const {
    return waypoints_;
  }

  virtual ~GSpline();
  const Eigen::VectorXd &get_coefficients() const { return coefficients_; }

  const Eigen::VectorXd &get_interval_lengths() const {
    return domain_interval_lengths_;
  }

  const std::string &get_basis_name() const { return basis_->get_name(); }

  std::size_t get_number_of_intervals() const { return number_of_intervals_; }

  GSpline linear_scaling_new_execution_time(double _new_exec_time) const;
};

const Eigen::Ref<const Eigen::VectorXd>
get_coefficient_segment(const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                        basis::Basis &_basis, std::size_t _num_interval,
                        std::size_t _codom_dim, std::size_t _interval,
                        std::size_t _component);

} // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
