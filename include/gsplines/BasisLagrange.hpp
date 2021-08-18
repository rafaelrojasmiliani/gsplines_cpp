#ifndef BASISLAGRANGE_H
#define BASISLAGRANGE_H

#include <gsplines/Basis.hpp>
namespace gsplines {

namespace basis {

class BasisLagrange : public Basis {

private:
  BasisLagrange &operator=(const BasisLagrange &);
  Eigen::VectorXd domain_points_;
  Eigen::VectorXd barycentric_weights_;

public:
  BasisLagrange(Eigen::Ref<const Eigen::VectorXd> _domain_points);
  BasisLagrange(const BasisLagrange &that);
  BasisLagrange(BasisLagrange &&that);

  virtual ~BasisLagrange();
  void eval_on_window(double _s, double _tau,
                      Eigen::Ref<Eigen::VectorXd> _buff) const override;

  void
  eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                            Eigen::Ref<Eigen::VectorXd> _buff) const override;

  void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd> _buff) const override;

  const Eigen::MatrixXd &get_derivative_matrix(double _tau, std::size_t _deg);

  void add_derivative_matrix(double tau, std::size_t _deg,
                             Eigen::Ref<Eigen::MatrixXd> _mat) override;

  void add_derivative_matrix_deriv_wrt_tau(
      double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) override;

  virtual std::unique_ptr<Basis> clone() const override {
    return std::make_unique<BasisLagrange>(*this);
  };

  virtual std::unique_ptr<Basis> move_clone() override {
    return std::make_unique<BasisLagrange>(std::move(*this));
  };

  static Eigen::VectorXd
  barycentric_weights(Eigen::Ref<const Eigen::VectorXd> _points);
  static Eigen::MatrixXd
  derivative_matrix(Eigen::Ref<const Eigen::VectorXd> _points);
};

} // namespace basis
} // namespace gsplines
#endif /* BASISLAGRANGE_H */
