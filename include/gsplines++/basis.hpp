#ifndef BASIS_H
#define BASIS_H
#include <cstddef>
#include <eigen3/Eigen/Core>

namespace gsplines {
namespace basis {
class Basis {
private:
  size_t dim_;
  Basis &operator=(const Basis &);

protected:
  Eigen::MatrixXd derivative_matrix_;

public:
  Basis(std::size_t _dim) : dim_(_dim), derivative_matrix_(dim_, dim_) {}
  Basis(const Basis &that)
      : dim_(that.get_dim()), derivative_matrix_(that.get_derivative_matrix()) {
  }
  virtual ~Basis() {}
  std::size_t get_dim() const { return dim_; }
  virtual void eval_on_window(double _s, double _tau,
                              Eigen::Ref<Eigen::VectorXd> _buff) = 0;
  virtual void eval_derivative_on_window(double _s, double _tau,
                                         unsigned int _deg,
                                         Eigen::Ref<Eigen::VectorXd> _buff) = 0;

  virtual void
  eval_derivative_wrt_tau_on_window(double _s, double _tau, unsigned int _deg,
                                    Eigen::Ref<Eigen::VectorXd> _buff) = 0;

  const Eigen::MatrixXd &get_derivative_matrix() const {
    return derivative_matrix_;
  }
};
} // namespace basis
} // namespace gsplines

#endif /* BASIS_H */
