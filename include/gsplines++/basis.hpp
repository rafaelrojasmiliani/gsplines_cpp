#ifndef BASIS_H
#define BASIS_H
#include <cstddef>

namespace gsplines {
namespace basis {
class Basis {
private:
  size_t dim_;
  Basis(const Basis &);
  Basis &operator=(const Basis &);

public:
  Basis(std::size_t _dim) : dim_(_dim){};
  virtual ~Basis() {}
  size_t get_dim() const { return dim_; }
  virtual void eval_on_window(double _s, double _tau,
                              Eigen::Ref<Eigen::VectorXd> _buff) = 0;
  virtual void eval_derivative_on_window(double _s, double _tau,
                                         unsigned int _deg,
                                         Eigen::Ref<Eigen::VectorXd> _buff) = 0;

  virtual void
  eval_derivative_wrt_tau_on_window(double _s, double _tau, unsigned int _deg,
                                    Eigen::Ref<Eigen::VectorXd> _buff) = 0;
};
} // namespace basis
} // namespace gsplines

#endif /* BASIS_H */
