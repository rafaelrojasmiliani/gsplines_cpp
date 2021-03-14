#ifndef BASISLEGENDRE_H
#define BASISLEGENDRE_H

#include <eigen3/Eigen/Core>
#include <gsplines++/Basis.hpp>
namespace gsplines {

namespace basis {

class BasisLegendre : public Basis {
  Eigen::MatrixXd dmat_;

private:
  BasisLegendre(const BasisLegendre &that);
  BasisLegendre &operator=(const BasisLegendre &);

public:
  BasisLegendre(std::size_t _dim);
  virtual ~BasisLegendre();
  void eval_on_window(double _s, double _tau,
                      Eigen::Ref<Eigen::VectorXd> _buff);
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 Eigen::Ref<Eigen::VectorXd> _buff);

  void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                         unsigned int _deg,
                                         Eigen::Ref<Eigen::VectorXd> _buff);
};

} // namespace basis
} // namespace gsplines
#endif /* BASISLEGENDRE_H */
