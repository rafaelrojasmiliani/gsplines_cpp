

#ifndef FUNCTION_EXPRESSION
#define FUNCTION_EXPRESSION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <functional>
#include <gsplines++/Functions/Functions.hpp>
#include <memory>
#include <utility>
#include <vector>
namespace gsplines {
namespace functions {

class FunctionExpression : public Function {
  friend FunctionExpression operator+(Function &_f1, Function &_f2);

public:
  enum Type { SUM = 0, MULTIPLICATION, COMPOSITION, CONCATENATION };

private:
  std::vector<std::unique_ptr<Function>> function_array_;
  typedef Eigen::MatrixXd(Eval_Function_Type)(
      std::vector<std::unique_ptr<Function>> &,
      const Eigen::Ref<const Eigen::VectorXd> &);
  typedef std::unique_ptr<Function>(Deriv_Function_Type)(
      std::vector<std::unique_ptr<Function>> &);

  std::function<Eval_Function_Type> eval_operation_;
  std::function<Deriv_Function_Type> deriv_operation_;

  Type type_;

public:
  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim,
                     Type _type,
                     std::vector<std::unique_ptr<Function>> &_function_array);

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim,
                     Type _type);

  FunctionExpression(const FunctionExpression &that);
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override {
    return eval_operation_(function_array_, _domain_points);
  }
  std::unique_ptr<Function> deriv(int _deg = 1) override {
    return deriv_operation_(function_array_);
  }

  FunctionExpression derivate(int _deg = 1) {
    std::unique_ptr<Function> result = deriv_operation_(function_array_);
    return FunctionExpression(*static_cast<FunctionExpression *>(result.get()));
  }

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<FunctionExpression>(*this);
  }

  const Type &get_type() const { return type_; }
};

Eigen::MatrixXd
eval_sum_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                   const Eigen::Ref<const Eigen::VectorXd> _domain_points);

Eigen::MatrixXd
eval_mul_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                   const Eigen::Ref<const Eigen::VectorXd> _domain_points);

Eigen::MatrixXd
eval_compose_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                       const Eigen::Ref<const Eigen::VectorXd> _domain_points);

Eigen::MatrixXd
eval_concat_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                      const Eigen::Ref<const Eigen::VectorXd> _domain_points);

std::unique_ptr<Function>
deriv_sum_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                    std::size_t _deg);

std::unique_ptr<Function>
deriv_mul_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                    std::size_t _deg);

std::unique_ptr<Function>
deriv_compose_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                        std::size_t _deg);

std::unique_ptr<Function>
deriv_concat_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                       std::size_t _deg);

FunctionExpression operator+(const Function &_f1, const Function &_f2);
FunctionExpression operator*(const Function &_f1, const Function &_f2);
FunctionExpression compose(const Function &_f1, const Function &_f2);
} // namespace functions
} // namespace gsplines
#endif
