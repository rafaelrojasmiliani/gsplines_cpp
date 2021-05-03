#ifndef BASISLEGENDRE_H
#define BASISLEGENDRE_H

#include <gsplines++/Basis.hpp>
namespace gsplines {

namespace basis {

class BasisLegendre : public Basis {

private:
  BasisLegendre &operator=(const BasisLegendre &);
  std::vector<Eigen::MatrixXd> derivative_matrices_buffer_;

public:
  BasisLegendre(std::size_t _dim);
  BasisLegendre(const BasisLegendre &that);
  virtual ~BasisLegendre();
  void eval_on_window(double _s, double _tau,
                      Eigen::Ref<Eigen::VectorXd> _buff) override;
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 Eigen::Ref<Eigen::VectorXd> _buff) override;

  void
  eval_derivative_wrt_tau_on_window(double _s, double _tau, unsigned int _deg,
                                    Eigen::Ref<Eigen::VectorXd> _buff) override;
  void add_derivative_matrix(double tau, std::size_t _deg,
                             Eigen::Ref<Eigen::MatrixXd> _mat) override;
  void add_derivative_matrix_deriv_wrt_tau(
      double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) override;

  virtual std::unique_ptr<Basis> clone() const override {
    return std::make_unique<BasisLegendre>(*this);
  };
};

} // namespace basis
} // namespace gsplines
#endif /* BASISLEGENDRE_H */
