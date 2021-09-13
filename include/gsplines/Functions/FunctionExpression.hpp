

#ifndef FUNCTION_EXPRESSION
#define FUNCTION_EXPRESSION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <functional>
#include <gsplines/Functions/FunctionBase.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <list>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {
class Function;

class FunctionExpression
    : public FunctionInheritanceHelper<FunctionExpression, FunctionBase,
                                       FunctionExpression> {

public:
  enum Type {
    SUM = 0,
    MULTIPLICATION,
    COMPOSITION,
    CONCATENATION,
    UNIQUE,
    NEGATIVE,
    EMPTY
  };

  std::list<std::unique_ptr<FunctionBase>> function_array_;

private:
  typedef void(Eval_Function_Type)(
      const std::list<std::unique_ptr<FunctionBase>> &,
      const Eigen::Ref<const Eigen::VectorXd> &,
      Eigen::Ref<Eigen::MatrixXd> _result);

  typedef FunctionExpression *(Deriv_Function_Type)(
      const std::list<std::unique_ptr<FunctionBase>> &, std::size_t);

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
      const std::list<std::unique_ptr<FunctionBase>> &_function_array,
      const std::string &_name = "FunctionExpression");

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim);

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim,
                     Type _type,
                     std::list<std::unique_ptr<FunctionBase>> &&_function_array,
                     const std::string &_name = "FunctionExpression");

  FunctionExpression(const FunctionExpression &that);
  FunctionExpression(FunctionExpression &&that);

  FunctionExpression(const Function &that);
  FunctionExpression(Function &&that);

  virtual void
  value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override {

    eval_operation_(function_array_, _domain_points, _result);
  };

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

  virtual void print(std::size_t _indent = 0) const override;

  std::string type_to_str() const;

  std::vector<std::pair<double, double>> get_arg_domains() const;

  void initialize();

  virtual ~FunctionExpression() = default;

  std::unique_ptr<FunctionExpression> deriv(std::size_t _deg = 1) const {
    return std::unique_ptr<FunctionExpression>(this->deriv_impl(_deg));
  }

protected:
  //
  virtual FunctionExpression *deriv_impl(std::size_t _deg = 1) const override {
    FunctionExpression *result = deriv_operation_(function_array_, _deg);
    assert(result != nullptr);
    return result;
  }

public:
  static std::list<std::unique_ptr<FunctionBase>>
  const_const_operation_handler(const FunctionExpression &_first,
                                const FunctionExpression &_second,
                                FunctionExpression::Type _opt_type);

  static std::list<std::unique_ptr<FunctionBase>>
  const_nonconst_operation_handler(const FunctionExpression &_first,
                                   FunctionExpression &&_second,
                                   FunctionExpression::Type _opt_type);

  static std::list<std::unique_ptr<FunctionBase>>
  nonconst_const_operation_handler(FunctionExpression &&_first,
                                   const FunctionExpression &_second,
                                   FunctionExpression::Type _opt_type);

  static std::list<std::unique_ptr<FunctionBase>>
  nonconst_nonconst_operation_handler(FunctionExpression &&_first,
                                      FunctionExpression &&_second,
                                      FunctionExpression::Type _opt_type);
};

FunctionExpression operator*(double, const FunctionExpression &);
FunctionExpression operator*(double, FunctionExpression &&);

void eval_unique_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_sum_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_mul_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_compose_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_concat_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

void eval_negative_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

FunctionExpression *deriv_unique_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

FunctionExpression *deriv_sum_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

FunctionExpression *deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

FunctionExpression *deriv_compose_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

FunctionExpression *deriv_concat_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

FunctionExpression *deriv_negative_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

} // namespace functions
} // namespace gsplines
#endif
