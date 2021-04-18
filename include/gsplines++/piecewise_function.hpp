#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include <eigen3/Eigen/Core>
#include <gsplines++/Basis.hpp>

namespace gsplines {

class PiecewiseFunction {
private:
  PiecewiseFunction &operator=(const PiecewiseFunction &);
  Eigen::VectorXd coefficients_;
  Eigen::VectorXd domain_break_points_;
  Eigen::VectorXd domain_interval_lengths_;
  std::size_t codom_dim_;
  std::size_t number_of_intervals_;
  basis::Basis basis_;
  std::size_t get_interval(double _domain_point);
  double interval_to_window(double _t, std::size_t _interval);

public:
  PiecewiseFunction(std::size_t _codom_dim, std::size_t _n_intervals,
                    basis::Basis _basis);
  PiecewiseFunction(const PiecewiseFunction &that);
  PiecewiseFunction deriv(unsigned int std::size_t _deg);
  virtual ~PiecewiseFunction();
  void operator()(double _domain_point, Eigen::VectorXd &_result) const;
  Eigen::VectorXd operator()(double _domain_point) const;
  void operator()(const Eigen::VectorXd &_domain_points,
                  Eigen::MatrixXd &_result) const;
  Eigen::MatrixXd operator()(const Eigen::VectorXd &_domain_points) const;
  Eigen::VectorXd coefficient_segment(std::size_t _interval,
                                      std::size_t _component);
};
} // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
