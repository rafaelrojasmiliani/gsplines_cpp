#ifndef BASIS_H
#define BASIS_H
#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <string>
#include <vector>

namespace gsplines {
namespace basis {
class Basis {
private:
  size_t dim_;
  Basis &operator=(const Basis &);
  const std::string name_;

protected:
  Eigen::MatrixXd derivative_matrix_;

public:
  Basis(std::size_t _dim, const std::string &_name)
      : dim_(_dim), derivative_matrix_(dim_, dim_), name_(_name) {}
  Basis(const Basis &that)
      : dim_(that.get_dim()), derivative_matrix_(that.derivative_matrix_),
        name_(that.name_) {}

  Basis(Basis &&that)
      : dim_(that.get_dim()),
        derivative_matrix_(std::move(that.derivative_matrix_)),
        name_(that.name_) {}

  virtual ~Basis() {}
  std::size_t get_dim() const { return dim_; }

  virtual void eval_on_window(double _s, double _tau,
                              Eigen::Ref<Eigen::VectorXd> _buff) const = 0;

  virtual void
  eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                            Eigen::Ref<Eigen::VectorXd> _buff) const = 0;

  virtual void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd> _buff) const = 0;

  const Eigen::MatrixXd &get_derivative_matrix() { return derivative_matrix_; };

  virtual void add_derivative_matrix(double tau, std::size_t _deg,
                                     Eigen::Ref<Eigen::MatrixXd> _mat) = 0;
  virtual void
  add_derivative_matrix_deriv_wrt_tau(double tau, std::size_t _deg,
                                      Eigen::Ref<Eigen::MatrixXd> _mat) = 0;

  virtual std::unique_ptr<Basis> clone() const = 0;
  virtual std::unique_ptr<Basis> move_clone() = 0;

  const std::string &get_name() const { return name_; };
};

std::unique_ptr<Basis> string_to_basis(const std::string &_basis_name);
} // namespace basis
} // namespace gsplines

#endif /* BASIS_H */
