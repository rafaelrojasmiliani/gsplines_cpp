#ifndef BASIS_H
#define BASIS_H
#include <cstddef>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/SparseCore>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace gsplines {
namespace basis {
/** Represent a set of function used to build a GSpline */
class Basis {
private:
  /* Avoid copy basis **/
  Basis &operator=(const Basis &);

  size_t dim_;
  const std::string name_;
  Eigen::VectorXd parameters_float_;
  Eigen::VectorXi parameters_int_;

  mutable std::vector<Eigen::MatrixXd> derivative_matrix_array_;

  mutable std::map<
      std::size_t,
      std::map<
          std::size_t,
          std::map<std::size_t, Eigen::SparseMatrix<double, Eigen::RowMajor>>>>
      continuity_matrix_buff_;

  mutable std::map<
      std::size_t,
      std::map<
          std::size_t,
          std::map<std::size_t, Eigen::SparseMatrix<double, Eigen::RowMajor>>>>
      continuity_matrix_dynamic_buff_;

  mutable std::map<
      std::size_t,
      std::map<
          std::size_t,
          std::map<std::size_t, Eigen::SparseMatrix<double, Eigen::RowMajor>>>>
      derivative_matrix_buff_;

  mutable std::map<
      std::size_t,
      std::map<
          std::size_t,
          std::map<std::size_t, Eigen::SparseMatrix<double, Eigen::RowMajor>>>>
      derivative_matrix_dynamic_buff_;

protected:
  Eigen::MatrixXd derivative_matrix_;

public:
  /**
   * @brief Constructor
   *
   * @param _dim Dimension of the base (number of functions that constitute the
   * base)
   * @param _name name of the base: Any label
   */
  Basis(std::size_t _dim, const std::string &_name,
        const Eigen::VectorXd &_parameters_float = Eigen::VectorXd(),
        const Eigen::VectorXi &_parameters_int = Eigen::VectorXi())
      : dim_(_dim), name_(_name), parameters_float_(_parameters_float),
        parameters_int_(_parameters_int), derivative_matrix_(dim_, dim_) {

    derivative_matrix_array_.push_back(
        Eigen::MatrixXd::Identity(get_dim(), get_dim()));
  }

  Basis(const Basis &that)
      : dim_(that.get_dim()), name_(that.name_),
        parameters_float_(that.parameters_float_),
        parameters_int_(that.parameters_int_),
        derivative_matrix_array_(that.derivative_matrix_array_),
        derivative_matrix_(that.derivative_matrix_) {}

  Basis(Basis &&that)
      : dim_(that.get_dim()), name_(that.name_),
        parameters_float_(std::move(that.parameters_float_)),
        parameters_int_(std::move(that.parameters_int_)),
        derivative_matrix_array_(std::move(that.derivative_matrix_array_)),
        derivative_matrix_(std::move(that.derivative_matrix_)) {}

  virtual ~Basis() = default;
  /**
   * @brief Get the dimension of the base, i.e. the number of functions that
   * constitute the base
   */
  std::size_t get_dim() const { return dim_; }

  /**
   * @brief Evaluates the basis functions in the window (windows is the canonic
   * interval, [-1, 1])
   *
   * @param _s Value inside the window [-1, 1]
   * @param _tau scaling factor, actual length of the interval in the GSpline
   * [t_i, t_{i+2})
   * @param _buff Buffer where the output is stored.
   */
  virtual void eval_on_window(
      double _s, double _tau,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const = 0;

  /**
   * @brief Evaluates the derivative of the basis functions in the window
   * (windows is the canonic interval, [-1, 1])
   *
   * @param _s Value inside the window [-1, 1]
   * @param _tau scaling factor, actual length of the interval in the GSpline
   * [t_i, t_{i+2})
   * @param _deg degree of the derivative
   * @param _buff Buffer where the output is stored.
   */
  virtual void eval_derivative_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const = 0;

  /**
   * @brief Evaluate the derivative of the basis with respect to tau, the actual
   * interval length inside the gspline
   *
   * @param _s Value inside the window [-1, 1]
   * @param _tau scaling factor, actual length of the interval in the GSpline
   * [t_i, t_{i+2})
   * @param _deg degree of the derivative
   * @param _buff Buffer where the output is stored.
   */
  virtual void eval_derivative_wrt_tau_on_window(
      double _s, double _tau, unsigned int _deg,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const = 0;

  const Eigen::MatrixXd &
  get_derivative_matrix_block(std::size_t _deg = 1) const {

    std::size_t current_deriv_calc = derivative_matrix_array_.size();
    if (current_deriv_calc <= _deg) {
      while (current_deriv_calc != _deg + 1) {
        derivative_matrix_array_.push_back(
            derivative_matrix_impl(current_deriv_calc));
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

  const Eigen::VectorXd &get_parameters() const { return parameters_float_; }

  virtual Eigen::MatrixXd derivative_matrix_impl(std::size_t _deg) const = 0;

  /*
  const Eigen::MatrixXd& get_sobolev_matrix(std::size_t _deg);
  virtual Eigen::MatrixXd sobolev_matrix_impl(std::size_t _deg) const = 0;
  const Eigen::MatrixXd& left_continuity_block(std::size_t _deg);
  const Eigen::MatrixXd& right_continuity_block(std::size_t _deg);
  */
  const Eigen::SparseMatrix<double, Eigen::RowMajor> &
  continuity_matrix(std::size_t _number_of_intervals, std::size_t _codom_dim,
                    std::size_t _deriv_order,
                    Eigen::Ref<const Eigen::VectorXd> _interval_lengths) const;

  const Eigen::SparseMatrix<double, Eigen::RowMajor> &gspline_derivative_matrix(
      std::size_t _number_of_intervals, std::size_t _codom_dim,
      std::size_t _deriv_order,
      Eigen::Ref<const Eigen::VectorXd> _interval_lengths) const;

  bool operator==(const Basis &_that) const;
  bool operator!=(const Basis &_that) const { return not(*this == _that); }
};

std::shared_ptr<Basis> get_basis(const std::string &_basis_name,
                                 std::size_t _dim,
                                 const std::vector<double> &_params);

std::shared_ptr<Basis>
get_basis(const std::string &_basis_name, std::size_t _dim,
          const Eigen::Ref<const Eigen::VectorXd> _params);
} // namespace basis
} // namespace gsplines

#endif /* BASIS_H */
