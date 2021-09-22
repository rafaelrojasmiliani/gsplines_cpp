#ifndef BASIS_H
#define BASIS_H
#include <cstddef>
#include <eigen3/Eigen/Core>
#include <iostream>
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
  Eigen::VectorXd parameters_;

  mutable std::vector<Eigen::MatrixXd> derivative_matrix_array_;

protected:
  Eigen::MatrixXd derivative_matrix_;

public:
  Basis(std::size_t _dim, const std::string &_name)
      : dim_(_dim), derivative_matrix_(dim_, dim_), name_(_name) {

    derivative_matrix_array_.push_back(
        Eigen::MatrixXd::Identity(get_dim(), get_dim()));
  }
  Basis(const Basis &that)
      : dim_(that.get_dim()), derivative_matrix_(that.derivative_matrix_),
        name_(that.name_),
        derivative_matrix_array_(that.derivative_matrix_array_) {}

  Basis(Basis &&that)
      : dim_(that.get_dim()),
        derivative_matrix_(std::move(that.derivative_matrix_)),
        name_(that.name_),
        derivative_matrix_array_(std::move(that.derivative_matrix_array_)) {}

  virtual ~Basis() = default;
  std::size_t get_dim() const { return dim_; }

  virtual void eval_on_window(
      double _s, double _tau,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const = 0;

  virtual void eval_derivative_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const = 0;

  virtual void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const = 0;

  const Eigen::MatrixXd &get_derivative_matrix(std::size_t _deg = 1) const {
    std::size_t current_deriv_calc = derivative_matrix_array_.size();
    if (current_deriv_calc <= _deg) {
      while (current_deriv_calc != _deg + 1) {
        derivative_matrix_array_.push_back(derivative_matrix_impl(_deg));
        current_deriv_calc++;
      }
    }
    return derivative_matrix_array_[_deg];
  }

  virtual void add_derivative_matrix(double tau, std::size_t _deg,
                                     Eigen::Ref<Eigen::MatrixXd> _mat) = 0;
  virtual void
  add_derivative_matrix_deriv_wrt_tau(double tau, std::size_t _deg,
                                      Eigen::Ref<Eigen::MatrixXd> _mat) = 0;

  virtual std::unique_ptr<Basis> clone() const = 0;
  virtual std::unique_ptr<Basis> move_clone() = 0;

  const std::string &get_name() const { return name_; };

  virtual Eigen::MatrixXd derivative_matrix_impl(std::size_t _deg) const = 0;
  /*
  const Eigen::MatrixXd& get_sobolev_matrix(std::size_t _deg);
  virtual Eigen::MatrixXd sobolev_matrix_impl(std::size_t _deg) const = 0;
  const Eigen::MatrixXd& left_continuity_block(std::size_t _deg);
  const Eigen::MatrixXd& right_continuity_block(std::size_t _deg);
  */
};

std::unique_ptr<Basis> string_to_basis(const std::string &_basis_name);
} // namespace basis
} // namespace gsplines

#endif /* BASIS_H */
