

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

class FunctionExpression;

FunctionExpression operator+(const Function &_f1, const Function &_f2);
FunctionExpression operator+(FunctionExpression &&_f1, const Function &_f2);
// cause ambiguity FunctionExpression operator+(const Function &_f2,
// FunctionExpression &&_f1);
FunctionExpression operator+(const FunctionExpression &_f1,
                             const Function &_f2);

FunctionExpression operator*(const Function &_f1, const Function &_f2);
FunctionExpression operator*(FunctionExpression &&_f1, const Function &_f2);
// FunctionExpression operator*(const Function &_f1, FunctionExpression &&_f2);

FunctionExpression compose(const Function &_f1, const Function &_f2);

class FunctionExpression : public Function {

public:
  enum Type { SUM = 0, MULTIPLICATION, COMPOSITION, CONCATENATION };
  // Consider use a deque
  // https://stackoverflow.com/questions/18811948/vector-vs-deque-operator
  // What is more critical to this calss? scalar multiplication or
  // operations that require random access?
  //
  // Consider to use a list. NOt that we do not desire a random access. We
  // deseire an
  // - Fast ordered access 
  // - fast concatenation
  // - Fast insertion at the begining
  // - Fast  inserion at the end.
  std::vector<std::unique_ptr<Function>> function_array_;

private:
  typedef Eigen::MatrixXd(Eval_Function_Type)(
      std::vector<std::unique_ptr<Function>> &,
      const Eigen::Ref<const Eigen::VectorXd> &);
  typedef std::unique_ptr<Function>(Deriv_Function_Type)(
      std::vector<std::unique_ptr<Function>> &, std::size_t);

  std::function<Eval_Function_Type> eval_operation_;
  std::function<Deriv_Function_Type> deriv_operation_;

  Type type_;

  static std::size_t num_call_constructor_;
  static std::size_t num_call_copy_constructor_;
  static std::size_t num_call_simple_constructor_;
  static std::size_t num_call_move_constructor_;

public:
  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim,
                     Type _type,
                     std::vector<std::unique_ptr<Function>> &_function_array);

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim,
                     Type _type);

  FunctionExpression(std::pair<double, double> _domain, std::size_t _codom_dim,
                     Type _type,
                     std::vector<std::unique_ptr<Function>> &&_function_array);

  FunctionExpression(const FunctionExpression &that);

  FunctionExpression(FunctionExpression &&that);

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override {
    printf("we are in operator ()\n");
    fflush(stdout);
    return eval_operation_(function_array_, _domain_points);
  }
  std::unique_ptr<Function> deriv(int _deg = 1) override {
    return deriv_operation_(function_array_, _deg);
  }

  FunctionExpression derivate(int _deg = 1) {
    std::unique_ptr<Function> result = deriv_operation_(function_array_, _deg);
    return FunctionExpression(*static_cast<FunctionExpression *>(result.get()));
  }

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<FunctionExpression>(*this);
  }

  const Type &get_type() const { return type_; }

  FunctionExpression sum(const FunctionExpression &that) {
    return *this + that;
  }

  FunctionExpression mul(const FunctionExpression &that) {
    return (*this) * (that);
  }

  FunctionExpression copose(const FunctionExpression &that) {
    return compose(*this, that);
  }

  void print_performace();

  Type get_type() { return type_; }

  FunctionExpression &operator+=(const Function &that);
  FunctionExpression &operator+=(const FunctionExpression &that);
  FunctionExpression &operator+=(FunctionExpression &&that);

  FunctionExpression &operator*=(const Function &that);
  FunctionExpression &operator*=(const FunctionExpression &that);
  FunctionExpression &operator*=(FunctionExpression &&that);
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

} // namespace functions
} // namespace gsplines
#endif
