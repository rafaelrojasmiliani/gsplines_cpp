#ifndef BASISLAGRANGE_H
#define BASISLAGRANGE_H

#include <gsplines++/Basis.hpp>
namespace gsplines {

namespace basis {

class BasisLagrange : public Basis {

private:
  BasisLagrange &operator=(const BasisLagrange &);
  std::vector<Eigen::MatrixXd> derivative_matrices_buffer_;
  Eigen::VectorXd domain_points_;

public:
  BasisLagrange(const Eigen::Ref<const Eigen::VectorXd> _domain_points);
  BasisLagrange(const BasisLagrange &that);
  virtual ~BasisLagrange();
  void eval_on_window(double _s, double _tau,
                      Eigen::Ref<Eigen::VectorXd> _buff);
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 Eigen::Ref<Eigen::VectorXd> _buff);

  void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                         unsigned int _deg,
                                         Eigen::Ref<Eigen::VectorXd> _buff);
  const Eigen::MatrixXd &get_derivative_matrix(double _tau, std::size_t _deg);
  void add_derivative_matrix(double tau, std::size_t _deg,
                             Eigen::Ref<Eigen::MatrixXd> _mat);
  void add_derivative_matrix_deriv_wrt_tau(double tau, std::size_t _deg,
                                           Eigen::Ref<Eigen::MatrixXd> _mat);

  virtual std::unique_ptr<Basis> clone() const override {
    return std::make_unique<BasisLagrange>(*this);
  };
};

} // namespace basis
} // namespace gsplines
#endif /* BASISLAGRANGE_H */
