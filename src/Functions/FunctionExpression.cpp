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
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const std::string &_name)
    : FunctionBase(_domain, _codom_dim, _name), type_(_type), function_array_(),
      eval_operation_(nullptr), deriv_operation_(nullptr) {

  for (const std::unique_ptr<FunctionExpression> &f : _function_array) {
    function_array_.push_back(f->clone());
  }

  assert(not(get_type() == SINGLE and get_name() == ""));
  initialize();

  num_call_constructor_++;
}

FunctionExpression::FunctionExpression(std::pair<double, double> _domain,
                                       std::size_t _codom_dim)
    : FunctionBase(_domain, _codom_dim, "DUMMY"), type_(SINGLE),
      function_array_(), eval_operation_(nullptr), deriv_operation_(nullptr) {

  assert(not(get_type() == SINGLE and get_name() == ""));
  initialize();

  num_call_constructor_++;
}

FunctionExpression::FunctionExpression(
    std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
    std::list<std::unique_ptr<FunctionExpression>> &&_function_array,
    const std::string &_name)
    : FunctionBase(_domain, _codom_dim, _name), type_(_type),
      eval_operation_(nullptr), deriv_operation_(nullptr),
      function_array_(std::move(_function_array)) {

  assert(not(get_type() == SINGLE and get_name() == ""));
  initialize();

  num_call_simple_constructor_++;
}

FunctionExpression::FunctionExpression(const FunctionExpression &that)
    : FunctionBase(that), type_(that.type_),
      eval_operation_(that.eval_operation_),
      deriv_operation_(that.deriv_operation_) {

  assert(not(get_type() == SINGLE and get_name() == ""));
  printf("lllllllll\n");
  for (const std::unique_ptr<FunctionExpression> &f : that.function_array_) {
    function_array_.push_back(f->clone());
  }
  num_call_copy_constructor_++;
}

FunctionExpression::FunctionExpression(FunctionExpression &&that)
    : FunctionBase(that), type_(that.type_),
      eval_operation_(that.eval_operation_),
      deriv_operation_(that.deriv_operation_),
      function_array_(std::move(that.function_array_)) {

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
                 [](const std::unique_ptr<FunctionExpression> &element) {
                   return element->get_domain();
                 });

  return result;
}

void FunctionExpression::initialize() {

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

  case SINGLE:
    eval_operation_ = eval_single_functions;
    deriv_operation_ = deriv_single_functions;
    break;

  case UNIQUE:
    eval_operation_ = eval_unique_functions;
    deriv_operation_ = deriv_unique_functions;
    break;

  default:
    throw std::invalid_argument("Function Expression Type not defined");
  }
}

std::string FunctionExpression::type_to_str() const {

  switch (type_) {
  case SUM:
    return "SUM";

  case MULTIPLICATION:
    return "MULTIPLICATION";

  case COMPOSITION:
    return "COMPOSITION";

  case CONCATENATION:
    return "CONCATENATION";
  case SINGLE:
    return "SINGLE: " + get_name();
  case UNIQUE:
    return "UNIQUE: " + get_name();
  default:
    throw std::invalid_argument(
        "FunctionExpression Expression Type not defined");
  }
}

void FunctionExpression::print(std::size_t _indent) const {

  printf("%*s- %s  %s\n", 4 * (int)_indent, "", get_name().c_str(),
         type_to_str().c_str());
  for (const std::unique_ptr<FunctionExpression> &f : function_array_) {
    f->print(_indent + 1);
  }
}

std::unique_ptr<FunctionExpression> deriv_single_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg) {
  throw std::invalid_argument("Single functions can't implement this method");
}

Eigen::MatrixXd eval_single_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
  throw std::invalid_argument("Single functions can't implement this method");
}

std::unique_ptr<FunctionExpression> deriv_unique_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg) {
  return _function_array.front()->deriv(_deg);
}

Eigen::MatrixXd eval_unique_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
  printf("EVAL UNIQUE !°°° .............. \n\n");
  Eigen::MatrixXd result = _function_array.front()->value(_domain_points);
  std::cout << "in eval uniqune ..... \n" << result << "----\n";
  return result;
}

} // namespace functions
} // namespace gsplines
