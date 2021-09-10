#ifndef FUNCTION_BASE
#define FUNCTION_BASE

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class FunctionExpression;
class FunctionBase {
private:
  std::size_t codom_dim_;
  std::pair<double, double> window_;
  std::pair<double, double> domain_;
  const std::string name_;

public:
  static const double dom_tollerance_;
  FunctionBase(std::pair<double, double> _domain, std::size_t _codom_dim,
               const std::string &_name);
  FunctionBase(const FunctionBase &that);

  const std::pair<double, double> &get_domain() const { return domain_; };

  std::size_t get_codom_dim() const { return codom_dim_; };

  static bool same_domain(const FunctionBase &_f1, const FunctionBase &_f2);

  static bool same_codomain(const FunctionBase &_f1, const FunctionBase &_f2);

  virtual ~FunctionBase() = default;

  bool is_point_in_domain(double _domain_point) const;

  virtual void print(std::size_t _indent = 0) const;

  const std::string &get_name() const { return name_; }

  double get_domain_length() const { return domain_.second - domain_.first; }

  std::unique_ptr<FunctionBase> clone() const & {
    return std::unique_ptr<FunctionBase>(this->clone_impl());
  }

  std::unique_ptr<FunctionBase> clone() && {
    return std::unique_ptr<FunctionBase>(this->move_clone_impl());
  }

  std::unique_ptr<FunctionBase> move_clone() {
    return std::unique_ptr<FunctionBase>(this->move_clone_impl());
  }

  std::unique_ptr<FunctionBase> deriv(std::size_t _deg = 1) const & {
    return std::unique_ptr<FunctionBase>(this->deriv_impl(_deg));
  }
  virtual void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                     Eigen::Ref<Eigen::MatrixXd> _result) const = 0;

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) const {
    Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());
    value(_domain_points, result);
    return std::move(result);
  }

  Eigen::MatrixXd
  value(const Eigen::Ref<const Eigen::VectorXd> _domain_points) const {
    return operator()(_domain_points);
  }

  FunctionExpression derivate(std::size_t _deg = 1);

protected:
  void set_domain(double _t0, double _t1) {
    domain_.first = _t0;
    domain_.second = _t1;
  }

  virtual FunctionBase *clone_impl() const = 0;

  virtual FunctionBase *move_clone_impl() = 0;

  virtual FunctionBase *deriv_impl(std::size_t _deg) const = 0;
};

} // namespace functions
} // namespace gsplines
#endif
