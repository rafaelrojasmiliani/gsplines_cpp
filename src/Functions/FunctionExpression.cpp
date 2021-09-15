#include <algorithm>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
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
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const std::string &_name)
    : FunctionInheritanceHelper(_domain, _codom_dim, _name), type_(_type),
      function_array_(), eval_operation_(nullptr), deriv_operation_(nullptr) {

  for (const std::unique_ptr<FunctionBase> &f : _function_array) {
    function_array_.push_back(f->clone());
  }

  assert(not(get_type() == UNIQUE and get_name() == ""));
  // assert(not(get_type() == UNIQUE and function_array_.size() == 0));
  initialize();

  num_call_constructor_++;
}
FunctionExpression::FunctionExpression(std::pair<double, double> _domain,
                                       std::size_t _codom_dim)
    : FunctionExpression(_domain, _codom_dim, UNIQUE, {}, "ZERO") {

  function_array_.push_back(
      std::make_unique<ConstFunction>(_domain, _codom_dim, 0.0));
}

FunctionExpression::FunctionExpression(
    std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
    std::list<std::unique_ptr<FunctionBase>> &&_function_array,
    const std::string &_name)
    : FunctionInheritanceHelper(_domain, _codom_dim, _name), type_(_type),
      eval_operation_(nullptr), deriv_operation_(nullptr),
      function_array_(std::move(_function_array)) {

  assert(not(get_type() == UNIQUE and get_name() == ""));
  // assert(not(get_type() == UNIQUE and function_array_.size() == 0));
  initialize();

  num_call_simple_constructor_++;
}

FunctionExpression::FunctionExpression(const FunctionExpression &that)
    : FunctionInheritanceHelper(that), type_(that.type_),
      eval_operation_(that.eval_operation_),
      deriv_operation_(that.deriv_operation_) {

  assert(not(get_type() == UNIQUE and get_name() == ""));
  // printf("lllllllll\n");
  for (const std::unique_ptr<FunctionBase> &f : that.function_array_) {
    function_array_.push_back(f->clone());
  }
  num_call_copy_constructor_++;
}

FunctionExpression::FunctionExpression(FunctionExpression &&that)
    : FunctionInheritanceHelper(that), type_(that.type_),
      eval_operation_(that.eval_operation_),
      deriv_operation_(that.deriv_operation_),
      function_array_(std::move(that.function_array_)) {

  num_call_move_constructor_++;
}

FunctionExpression::FunctionExpression(const Function &_that)
    : FunctionExpression(_that.get_domain(), _that.get_codom_dim(), UNIQUE,
                         {}) {
  function_array_.push_back(_that.clone());
  num_call_copy_constructor_++;
}
FunctionExpression::FunctionExpression(Function &&_that)
    : FunctionExpression(_that.get_domain(), _that.get_codom_dim(), UNIQUE,
                         {}) {

  function_array_.push_back(_that.move_clone());
  num_call_copy_constructor_++;
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
                 [](const std::unique_ptr<FunctionBase> &element) {
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

  case UNIQUE:
    eval_operation_ = eval_unique_functions;
    deriv_operation_ = deriv_unique_functions;
    break;

  case NEGATIVE:
    eval_operation_ = eval_negative_functions;
    deriv_operation_ = deriv_negative_functions;
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
  case UNIQUE:
    return "UNIQUE: " + get_name();
  case NEGATIVE:
    return "NEGATIVE: " + get_name();
  default:
    throw std::invalid_argument(
        "FunctionExpression Expression Type not defined");
  }
}

void FunctionExpression::print(std::size_t _indent) const {

  FunctionBase::print(_indent);
  printf("%*s %s  %s\n", 4 * (int)_indent, "", "Expression type",
         type_to_str().c_str());
  for (const std::unique_ptr<FunctionBase> &f : function_array_) {
    f->print(_indent + 1);
  }
}

FunctionExpression *deriv_unique_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg) {

  assert(_function_array.size() == 1);
  std::pair<double, double> domain = _function_array.front()->get_domain();
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::list<std::unique_ptr<FunctionBase>> result_array;
  result_array.push_back(_function_array.front()->deriv(_deg));

  return new FunctionExpression(domain, codom_dim,
                                FunctionExpression::Type::UNIQUE,
                                std::move(result_array));
}

void eval_unique_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) {
  assert(_function_array.size() == 1);
  // printf("EVAL UNIQUE !°°° .............. \n\n");

  _function_array.front()->value(_domain_points, _result);
  // std::cout << "in eval uniqune ..... \n" << _result << "----\n";
}

FunctionExpression *deriv_negative_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg) {

  assert(_function_array.size() == 1);
  std::pair<double, double> domain = _function_array.front()->get_domain();
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::list<std::unique_ptr<FunctionBase>> result_array;
  result_array.push_back(_function_array.front()->deriv(_deg));

  return new FunctionExpression(domain, codom_dim,
                                FunctionExpression::Type::NEGATIVE,
                                std::move(result_array));
}

void eval_negative_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) {
  assert(_function_array.size() == 1);

  _function_array.front()->value(_domain_points, _result);
  _result *= -1.0;
}

} // namespace functions
} // namespace gsplines
