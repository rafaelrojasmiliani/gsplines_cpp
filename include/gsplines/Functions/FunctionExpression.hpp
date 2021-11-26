

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

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
class FunctionExpression
    : public FunctionInheritanceHelper<FunctionExpression<Domain, DDIM, DIM>,
                                       FunctionBase<Domain, DDIM, DIM>,
                                       FunctionExpression<Domain, DDIM, DIM>> {

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

  typedef FunctionBase<Domain, DDIM, DIM> FunctionBaseType;
  typedef FunctionExpression<Domain, DDIM, DIM> FunctionExpressionType;

  std::list<std::unique_ptr<FunctionBaseType>> function_array_;

private:
  typedef void(Eval_Function_Type)(
      const std::list<std::unique_ptr<FunctionBaseType>> &,
      const Eigen::Ref<const Eigen::VectorXd> &,
      Eigen::Ref<Eigen::MatrixXd> _result);

  typedef FunctionExpression *(Deriv_Function_Type)(
      const std::list<std::unique_ptr<FunctionBaseType>> &, std::size_t);

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
      const std::list<std::unique_ptr<FunctionBaseType>> &_function_array,
      const std::string &_name = "FunctionExpression");

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim);

  FunctionExpression(
      std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
      std::list<std::unique_ptr<FunctionBaseType>> &&_function_array,
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

  virtual FunctionExpressionType
  operator+(const FunctionExpressionType &that) const &;
  virtual FunctionExpressionType
  operator+(FunctionExpressionType &&that) const &;
  virtual FunctionExpressionType
  operator+(const FunctionExpressionType &that) &&;
  virtual FunctionExpressionType operator+(FunctionExpressionType &&that) &&;

  virtual FunctionExpressionType
  operator-(const FunctionExpressionType &that) const &;
  virtual FunctionExpressionType
  operator-(FunctionExpressionType &&that) const &;
  virtual FunctionExpressionType
  operator-(const FunctionExpressionType &that) &&;
  virtual FunctionExpressionType operator-(FunctionExpressionType &&that) &&;
  virtual FunctionExpressionType operator-() const &;
  virtual FunctionExpressionType operator-() &&;

  virtual FunctionExpressionType
  operator*(const FunctionExpressionType &that) const &;
  virtual FunctionExpressionType
  operator*(FunctionExpressionType &&that) const &;
  virtual FunctionExpressionType
  operator*(const FunctionExpressionType &that) &&;
  virtual FunctionExpressionType operator*(FunctionExpressionType &&that) &&;

  virtual FunctionExpressionType
  compose(const FunctionExpressionType &that) const &;
  virtual FunctionExpressionType compose(FunctionExpressionType &&that) const &;

  virtual FunctionExpressionType compose(const FunctionExpressionType &that) &&;
  virtual FunctionExpressionType compose(FunctionExpressionType &&that) &&;

  virtual FunctionExpressionType
  concat(const FunctionExpressionType &that) const &;
  virtual FunctionExpressionType concat(FunctionExpressionType &&that) const &;

  virtual FunctionExpressionType concat(const FunctionExpressionType &that) &&;
  virtual FunctionExpressionType concat(FunctionExpressionType &&that) &&;

  void operator+=(FunctionExpressionType &&that);
  void operator+=(const FunctionExpressionType &that);

  void operator*=(FunctionExpressionType &&that);
  void operator*=(const FunctionExpressionType &that);

  virtual void print(std::size_t _indent = 0) const override;

  std::string type_to_str() const;

  std::vector<std::pair<double, double>> get_arg_domains() const;

  void initialize();

  virtual ~FunctionExpression() = default;

  std::unique_ptr<FunctionExpressionType> deriv(std::size_t _deg = 1) const {
    return std::unique_ptr<FunctionExpressionType>(this->deriv_impl(_deg));
  }

protected:
  //
  virtual FunctionExpressionType *
  deriv_impl(std::size_t _deg = 1) const override {
    FunctionExpressionType *result = deriv_operation_(function_array_, _deg);
    assert(result != nullptr);
    return result;
  }

public:
  static std::list<std::unique_ptr<FunctionBaseType>>
  const_const_operation_handler(const FunctionExpressionType &_first,
                                const FunctionExpressionType &_second,
                                FunctionExpressionType::Type _opt_type);

  static std::list<std::unique_ptr<FunctionBaseType>>
  const_nonconst_operation_handler(const FunctionExpressionType &_first,
                                   FunctionExpressionType &&_second,
                                   FunctionExpressionType::Type _opt_type);

  static std::list<std::unique_ptr<FunctionBaseType>>
  nonconst_const_operation_handler(FunctionExpressionType &&_first,
                                   const FunctionExpressionType &_second,
                                   FunctionExpressionType::Type _opt_type);

  static std::list<std::unique_ptr<FunctionBaseType>>
  nonconst_nonconst_operation_handler(FunctionExpressionType &&_first,
                                      FunctionExpressionType &&_second,
                                      FunctionExpressionType::Type _opt_type);
};

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression<Domain, DDIM, DIM>
operator*(double, const FunctionExpression<Domain, DDIM, DIM> &);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression<Domain, DDIM, DIM>
operator*(double, FunctionExpression<Domain, DDIM, DIM> &&);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
void eval_unique_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
void eval_sum_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
void eval_mul_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
void eval_compose_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
void eval_concat_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
void eval_negative_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression *deriv_unique_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression *deriv_sum_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression *deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression *deriv_compose_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression *deriv_concat_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

template <template <std::size_t> typename Domain, std::size_t DDIM,
          std::size_t DIM>
FunctionExpression *deriv_negative_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg);

} // namespace functions
} // namespace gsplines
#endif
