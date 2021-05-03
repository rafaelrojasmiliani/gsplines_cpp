#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include <eigen3/Eigen/Core>
#include <gsplines++/Basis.hpp>

namespace gsplines {
class SobolevNorm;
class PiecewiseFunction {
  friend SobolevNorm;

private:
  PiecewiseFunction &operator=(const PiecewiseFunction &);
  const std::size_t codom_dim_;
  const std::size_t number_of_intervals_;
  std::unique_ptr<basis::Basis> basis_;
  Eigen::VectorXd coefficients_;
  Eigen::VectorXd domain_break_points_;
  Eigen::VectorXd domain_interval_lengths_;
  Eigen::VectorXd basis_buffer_;
  std::size_t get_interval(double _domain_point) const;
  double interval_to_window(double _t, std::size_t _interval) const;

  Eigen::Ref<Eigen::VectorXd> coefficient_segment(std::size_t _interval,
                                                  std::size_t _component);

public:
  PiecewiseFunction(std::size_t _codom_dim, std::size_t _n_intervals,
                    const basis::Basis &_basis,
                    const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                    const Eigen::Ref<const Eigen::VectorXd> _tauv);
  PiecewiseFunction(const PiecewiseFunction &that);
  PiecewiseFunction deriv(std::size_t _deg);
  virtual ~PiecewiseFunction();
  // Eigen::VectorXd operator()(double _domain_point);
  //  void operator()(double _domain_point, Eigen::VectorXd &_result);
  //  void operator()(const Eigen::VectorXd &_domain_points,
  //                  Eigen::MatrixXd &_result) const;
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points);

  std::size_t get_codom_dim() { return codom_dim_; }
  std::size_t get_intervals_num() { return number_of_intervals_; }
  double get_exec_time() {
    return domain_break_points_.tail(1)(0) - domain_break_points_(0);
  }
  const Eigen::Ref<const Eigen::VectorXd> get_domain_breakpoints() {
    return domain_break_points_;
  }
  Eigen::VectorXd get_coeff();
};

const Eigen::Ref<const Eigen::VectorXd>
get_coefficient_segment(const Eigen::Ref<const Eigen::VectorXd> _coefficents,
                        basis::Basis &_basis, std::size_t _num_interval,
                        std::size_t _codom_dim, std::size_t _interval,
                        std::size_t _component);

} // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
