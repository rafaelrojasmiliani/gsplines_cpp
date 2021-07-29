

#ifndef FUNCTION_EXPRESSION
#define FUNCTION_EXPRESSION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <functional>
#include <gsplines++/Functions/FunctionBase.hpp>
#include <list>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

/*
FunctionExpression compose(const FunctionExpression &_f1, const
FunctionExpression &_f2);

*/

class FunctionExpression : public FunctionBase {

public:
  enum Type {
    SUM = 0,
    MULTIPLICATION,
    COMPOSITION,
    CONCATENATION,
    SINGLE,
    UNIQUE
  };

  std::list<std::unique_ptr<FunctionExpression>> function_array_;

private:
  typedef void(Eval_Function_Type)(
      const std::list<std::unique_ptr<FunctionExpression>> &,
      const Eigen::Ref<const Eigen::VectorXd> &,
      Eigen::Ref<Eigen::MatrixXd> _result);
  typedef std::unique_ptr<FunctionExpression>(Deriv_Function_Type)(
      const std::list<std::unique_ptr<FunctionExpression>> &, std::size_t);

  std::function<Eval_Function_Type> eval_operation_;
  std::function<Deriv_Function_Type> deriv_operation_;

  Type type_;

  Eigen::VectorXd domain_break_points_;

  static std::size_t num_call_constructor_;
  static std::size_t num_call_copy_constructor_;
  static std::size_t num_call_simple_constructor_;
  static std::size_t num_call_move_constructor_;

public:
  FunctionExpression(
      std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
      const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
      const std::string &_name = "FunctionExpression");

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim);

  FunctionExpression(
      std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
      std::list<std::unique_ptr<FunctionExpression>> &&_function_array,
      const std::string &_name = "FunctionExpression");

  FunctionExpression(const FunctionExpression &that);

  FunctionExpression(FunctionExpression &&that);

  Eigen::MatrixXd
  value(const Eigen::Ref<const Eigen::VectorXd> _domain_points) const {
    return operator()(_domain_points);
  }

  virtual void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                     Eigen::Ref<Eigen::MatrixXd> _result) const {

    eval_operation_(function_array_, _domain_points, _result);
  };

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) const {
    Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());
    value(_domain_points, result);
    return std::move(result);
  }

  virtual std::unique_ptr<FunctionExpression> deriv(int _deg = 1) const {
    return deriv_operation_(function_array_, _deg);
  }

  virtual FunctionExpression derivate(int _deg = 1) const {
    std::unique_ptr<FunctionExpression> result =
        deriv_operation_(function_array_, _deg);
    return FunctionExpression(std::move(*result));
  }

  virtual std::unique_ptr<FunctionExpression> clone() const {
    return std::make_unique<FunctionExpression>(*this);
  }

  virtual std::unique_ptr<FunctionExpression> move_clone() {
    return std::make_unique<FunctionExpression>(std::move(*this));
  }

  const Type &get_type() const { return type_; }

  void print_performace();

  const Type &get_type() { return type_; }

  virtual FunctionExpression operator+(const FunctionExpression &that) const &;
  virtual FunctionExpression operator+(FunctionExpression &&that) const &;
  virtual FunctionExpression operator+(const FunctionExpression &that) &&;
  virtual FunctionExpression operator+(FunctionExpression &&that) &&;

  virtual FunctionExpression operator-(const FunctionExpression &that) const &;
  virtual FunctionExpression operator-(FunctionExpression &&that) const &;
  virtual FunctionExpression operator-(const FunctionExpression &that) &&;
  virtual FunctionExpression operator-(FunctionExpression &&that) &&;
  virtual FunctionExpression operator-() const &;
  virtual FunctionExpression operator-() &&;

  virtual FunctionExpression operator*(const FunctionExpression &that) const &;
  virtual FunctionExpression operator*(FunctionExpression &&that) const &;
  virtual FunctionExpression operator*(const FunctionExpression &that) &&;
  virtual FunctionExpression operator*(FunctionExpression &&that) &&;

  virtual FunctionExpression compose(const FunctionExpression &that) const &;
  virtual FunctionExpression compose(FunctionExpression &&that) const &;

  virtual FunctionExpression compose(const FunctionExpression &that) &&;
  virtual FunctionExpression compose(FunctionExpression &&that) &&;

  virtual FunctionExpression concat(const FunctionExpression &that) const &;
  virtual FunctionExpression concat(FunctionExpression &&that) const &;

  virtual FunctionExpression concat(const FunctionExpression &that) &&;
  virtual FunctionExpression concat(FunctionExpression &&that) &&;

  virtual void print(std::size_t _indent = 0) const;

  std::string type_to_str() const;

  std::vector<std::pair<double, double>> get_arg_domains() const;

  void initialize();

  virtual ~FunctionExpression() {}
};

FunctionExpression operator*(double, const FunctionExpression &);
FunctionExpression operator*(double, FunctionExpression &&);

void eval_unique_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_single_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_sum_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_mul_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_compose_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_concat_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

std::unique_ptr<FunctionExpression> deriv_single_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg);

std::unique_ptr<FunctionExpression> deriv_unique_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg);

std::unique_ptr<FunctionExpression> deriv_sum_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg);

std::unique_ptr<FunctionExpression> deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg);

std::unique_ptr<FunctionExpression> deriv_compose_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg);

std::unique_ptr<FunctionExpression> deriv_concat_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg);

} // namespace functions
} // namespace gsplines
#endif
