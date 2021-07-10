#include <algorithm>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <iostream>

namespace gsplines {
namespace functions {

/* -----
 *  Function Construction
 * -----*/

std::size_t FunctionExpression::num_call_constructor_ = 0;
std::size_t FunctionExpression::num_call_simple_constructor_ = 0;
std::size_t FunctionExpression::num_call_copy_constructor_ = 0;
std::size_t FunctionExpression::num_call_move_constructor_ = 0;

FunctionExpression::FunctionExpression(
    std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
    std::list<std::unique_ptr<Function>> &_function_array)
    : Function(_domain, _codom_dim, "FunctionExpression"), type_(_type),
      function_array_(), eval_operation_(nullptr), deriv_operation_(nullptr) {

  for (const std::unique_ptr<Function> &f : function_array_) {
    function_array_.push_back(f->clone());
  }

  switch (type_) {
  case SUM:
    eval_operation_ = eval_sum_functions;
    deriv_operation_ = deriv_sum_functions;
    break;

  case MULTIPLICATION:
    eval_operation_ = eval_mul_functions;
    deriv_operation_ = deriv_mul_functions;
    break;

  case COMPOSITION:
    eval_operation_ = eval_compose_functions;
    deriv_operation_ = deriv_compose_functions;
    break;

  case CONCATENATION:
    eval_operation_ = eval_concat_functions;
    deriv_operation_ = deriv_concat_functions;
    break;
  default:
    throw std::invalid_argument("Function Expression Type not defined");
  }
  num_call_constructor_++;
}

FunctionExpression::FunctionExpression(
    std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
    std::list<std::unique_ptr<Function>> &&_function_array)
    : Function(_domain, _codom_dim, "FunctionExpression"), type_(_type),
      eval_operation_(nullptr), deriv_operation_(nullptr),
      function_array_(std::move(_function_array)) {

  switch (type_) {
  case SUM:
    eval_operation_ = eval_sum_functions;
    deriv_operation_ = deriv_sum_functions;
    break;

  case MULTIPLICATION:
    eval_operation_ = eval_mul_functions;
    deriv_operation_ = deriv_mul_functions;
    break;

  case COMPOSITION:
    eval_operation_ = eval_compose_functions;
    deriv_operation_ = deriv_compose_functions;
    break;

  case CONCATENATION:
    eval_operation_ = eval_concat_functions;
    deriv_operation_ = deriv_concat_functions;
    break;
  default:
    throw std::invalid_argument("Function Expression Type not defined");
  }
  num_call_simple_constructor_++;
}

FunctionExpression::FunctionExpression(const FunctionExpression &that)
    : Function(that), type_(that.type_), eval_operation_(that.eval_operation_),
      deriv_operation_(that.deriv_operation_) {

  printf("copy  CONSTRUCTOR ----------\n");
  for (const std::unique_ptr<Function> &f : that.function_array_) {
    function_array_.push_back(f->clone());
  }
  num_call_copy_constructor_++;
}

FunctionExpression::FunctionExpression(FunctionExpression &&that)
    : Function(that), type_(that.type_), eval_operation_(that.eval_operation_),
      deriv_operation_(that.deriv_operation_),
      function_array_(std::move(that.function_array_)) {

  printf("MOOOOOOOOOOOOOOOOOVEEE  CONSTRUCTOR ----------\n");
  fflush(stdout);
  num_call_move_constructor_++;
}

void FunctionExpression::print_performace() {
  printf("Num call of the constructor %zu\n", num_call_constructor_);
  printf("Num call of the COPY constructor %zu\n", num_call_copy_constructor_);
  printf("Num call of the MOVE constructor %zu\n", num_call_move_constructor_);
  printf("Num call of the SIMPLE constructor %zu\n",
         num_call_simple_constructor_);
}

std::vector<std::pair<double, double>>
FunctionExpression::get_arg_domains() const {

  std::vector<std::pair<double, double>> result;

  std::transform(function_array_.begin(), function_array_.end(),
                 std::back_inserter(result),
                 [](const std::unique_ptr<Function> &element) {
                   return element->get_domain();
                 });

  return result;
}
} // namespace functions
} // namespace gsplines
