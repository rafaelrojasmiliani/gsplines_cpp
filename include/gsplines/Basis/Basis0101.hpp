#pragma once

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>

namespace gsplines::basis {
/** Represent a set of function used to build a GSpline.
 * Let I be an interval of R. This class represent a set of
 * functions f_i: I -> R and contains the tools to compute them.*/
class Basis0101 : public Basis {
 private:
  mutable Eigen::MatrixXd qBuffer_ = Eigen::MatrixXd(6, 6);

 public:
  static std::shared_ptr<Basis0101> get(double _alpha);
  Basis0101(double _alpha);
  void eval_on_window(double _s, double _tau,
                      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>>
                          _buff) const override;
  void eval_derivative_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff)
      const override;

  void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff)
      const override;

  void add_derivative_matrix(double tau, std::size_t _deg,
                             Eigen::Ref<Eigen::MatrixXd> _mat) override;
  void add_derivative_matrix_deriv_wrt_tau(
      double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) override;

  double get_alpha() const;

  std::unique_ptr<Basis> clone() const override;
  std::unique_ptr<Basis> move_clone() override;
  Eigen::MatrixXd derivative_matrix_impl(std::size_t _deg) const override;
};

}  // namespace gsplines::basis
