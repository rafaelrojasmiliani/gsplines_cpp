#ifndef FUNCTION_BASE
#define FUNCTION_BASE

#include <eigen3/Eigen/Core>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class FunctionExpression;
class DotProduct;
class FunctionBase {
 private:
  std::size_t codom_dim_;
  std::pair<double, double> window_;
  std::pair<double, double> domain_;
  const std::string name_;

 public:
  static const double dom_tollerance_;
  FunctionBase(std::pair<double, double> _domain, std::size_t _codom_dim,
               const std::string& _name);
  FunctionBase(const FunctionBase& that);

  [[nodiscard]] const std::pair<double, double>& get_domain() const {
    return domain_;
  };

  [[nodiscard]] std::size_t get_codom_dim() const { return codom_dim_; };

  [[nodiscard]] bool same_domain(const FunctionBase& _f1) const;
  [[nodiscard]] bool same_codomain(const FunctionBase& _f1) const;

  static bool same_domain(const FunctionBase& _f1, const FunctionBase& _f2);

  static bool same_codomain(const FunctionBase& _f1, const FunctionBase& _f2);

  virtual ~FunctionBase() = default;

  [[nodiscard]] bool is_point_in_domain(double _domain_point) const;

  virtual void print(std::size_t _indent = 0) const;

  [[nodiscard]] const std::string& get_name() const { return name_; }

  [[nodiscard]] double get_domain_length() const {
    return domain_.second - domain_.first;
  }

  [[nodiscard]] std::unique_ptr<FunctionBase> clone() const& {
    return std::unique_ptr<FunctionBase>(this->clone_impl());
  }

  std::unique_ptr<FunctionBase> clone() && {
    return std::unique_ptr<FunctionBase>(this->move_clone_impl());
  }

  std::unique_ptr<FunctionBase> move_clone() {
    return std::unique_ptr<FunctionBase>(this->move_clone_impl());
  }

  [[nodiscard]] std::unique_ptr<FunctionBase> deriv(
      std::size_t _deg = 1) const& {
    return std::unique_ptr<FunctionBase>(this->deriv_impl(_deg));
  }
  virtual void value_impl(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points,
      Eigen::Ref<Eigen::MatrixXd> _result) const = 0;

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const {
    value_impl(_domain_points, _result);
  }

  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const {
    Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());
    value_impl(_domain_points, result);
    return result;
  }

  [[nodiscard]] Eigen::MatrixXd value(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const {
    return operator()(_domain_points);
  }

  DotProduct dot(const FunctionBase& _that) const&;
  DotProduct dot(FunctionBase&& _that) const&;
  DotProduct dot(const FunctionBase& _that) &&;
  DotProduct dot(FunctionBase&& _that) &&;

  [[nodiscard]] FunctionExpression derivate(std::size_t _deg = 1) const;

  [[nodiscard]] FunctionExpression to_expression() const&;
  FunctionExpression to_expression() &&;

  virtual FunctionExpression operator+(const FunctionExpression& that) const&;
  virtual FunctionExpression operator+(FunctionExpression&& that) const&;
  virtual FunctionExpression operator+(const FunctionExpression& that) &&;
  virtual FunctionExpression operator+(FunctionExpression&& that) &&;

  virtual FunctionExpression operator-(const FunctionExpression& that) const&;
  virtual FunctionExpression operator-(FunctionExpression&& that) const&;
  virtual FunctionExpression operator-(const FunctionExpression& that) &&;
  virtual FunctionExpression operator-(FunctionExpression&& that) &&;
  virtual FunctionExpression operator-() const&;
  virtual FunctionExpression operator-() &&;

  virtual FunctionExpression operator*(const FunctionExpression& that) const&;
  virtual FunctionExpression operator*(FunctionExpression&& that) const&;
  virtual FunctionExpression operator*(const FunctionExpression& that) &&;
  virtual FunctionExpression operator*(FunctionExpression&& that) &&;

  [[nodiscard]] virtual FunctionExpression compose(
      const FunctionExpression& that) const&;
  virtual FunctionExpression compose(FunctionExpression&& that) const&;

  virtual FunctionExpression compose(const FunctionExpression& that) &&;
  virtual FunctionExpression compose(FunctionExpression&& that) &&;

  [[nodiscard]] virtual FunctionExpression concat(
      const FunctionExpression& that) const&;
  virtual FunctionExpression concat(FunctionExpression&& that) const&;

  virtual FunctionExpression concat(const FunctionExpression& that) &&;
  virtual FunctionExpression concat(FunctionExpression&& that) &&;

 protected:
  void set_domain(double _t0, double _t1) {
    domain_.first = _t0;
    domain_.second = _t1;
  }
  void set_domain(const std::pair<double, double>& _new_domain) {
    domain_ = _new_domain;
  }

  virtual FunctionBase* clone_impl() const = 0;

  virtual FunctionBase* move_clone_impl() = 0;

  virtual FunctionBase* deriv_impl(std::size_t _deg) const = 0;
};

}  // namespace functions
}  // namespace gsplines
#endif
