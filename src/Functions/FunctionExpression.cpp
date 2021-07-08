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
    : Function(_domain, _codom_dim), type_(_type), function_array_(),
      eval_operation_(nullptr), deriv_operation_(nullptr) {

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
    : Function(_domain, _codom_dim), type_(_type), eval_operation_(nullptr),
      deriv_operation_(nullptr), function_array_(std::move(_function_array)) {

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

FunctionExpression::FunctionExpression(std::pair<double, double> _domain,
                                       std::size_t _codom_dim, Type _type)
    : Function(_domain, _codom_dim), type_(_type) {

  printf("calling wrong constructor \n-----\n");
  fflush(stdout);
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

/* -----
 *  Function Evaluation
 * -----*/

Eigen::MatrixXd
eval_compose_functions(std::list<std::unique_ptr<Function>> &_function_array,
                       const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array.back()->get_codom_dim());
  Eigen::VectorXd domain_ponts = _domain_points;
  std::list<std::unique_ptr<Function>>::const_iterator it;
  std::list<std::unique_ptr<Function>>::const_iterator it_limit =
      std::next(_function_array.end(), -1);

  for (it = _function_array.begin(); it != it_limit; it++) {
    domain_ponts = (*it)->value(domain_ponts);
  }
  result = _function_array.back()->value(_domain_points);
  return result;
}

Eigen::MatrixXd
eval_concat_functions(std::list<std::unique_ptr<Function>> &_function_array,
                      const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array.front()->get_codom_dim());

  return result;
}

/* -----
 *  Function Derivation
 * -----*/
std::unique_ptr<Function> first_deriv_compose_functions(
    std::list<std::unique_ptr<Function>> &_function_array) {

  std::list<std::unique_ptr<Function>> result_array;

  result_array.push_back(_function_array.front()->deriv());

  std::list<std::unique_ptr<Function>>::const_iterator it;
  //
  // oritiginal compisition
  // +-----+-----+-----+-----+-----+
  // | f1  | f2  | f3  | f4  |  f5 |
  // +-----+-----+-----+-----+-----+
  //
  // This returs f(t) = f5 \circ f4 \circ f3 \circ \f2 \circ \f1(t)
  //
  // df dt =
  //
  // ( d f5 d t  \circ f4 \circ f3 \circ \f2 \circ \f1(t) )
  // ( d f4 ft \circ f3 \circ \f2 \circ \f1(t) )
  // ( d f3 ft \circ f2 \circ \f1(t) )
  // ( d f2 ft \circ \f1(t) )
  // ( d f1 dt (t) )
  //
  // MUL
  // +--------------------------------------------------+
  // | d f5 d t  \circ f4 \circ f3 \circ \f2 \circ \f1  |
  // +--------------------------------------------------+
  // | d f4 ft \circ f3 \circ \f2 \circ \f1             |
  // +--------------------------------------------------+
  // |         d f3 ft \circ f2 \circ \f1               |
  // +--------------------------------------------------+
  // |              d f2 ft \circ \f1                   |
  // +--------------------------------------------------+
  // |                d f1 dt (t)                       |
  // +--------------------------------------------------+

  for (it = std::next(_function_array.begin(), 1); it != _function_array.end();
       it++) {

    std::list<std::unique_ptr<Function>> elem_array;
    std::list<std::unique_ptr<Function>>::const_iterator it_elem;

    for (it_elem = _function_array.begin(); it_elem != it; it_elem++) {
      elem_array.push_back((*it_elem)->clone());
    }

    elem_array.push_back((*it)->deriv());

    std::pair<double, double> domain = (*it)->get_domain();
    std::size_t codom_dim = (*it)->get_codom_dim();

    result_array.push_front(std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::COMPOSITION,
        std::move(elem_array)));
  }
  std::size_t codom_dim = _function_array.back()->get_codom_dim();
  std::pair<double, double> domain = _function_array.back()->get_domain();
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
      std::move(result_array));
}

std::unique_ptr<Function>
deriv_compose_functions(std::list<std::unique_ptr<Function>> &_function_array,
                        std::size_t _deg) {

  std::size_t codom_dim = _function_array.back()->get_codom_dim();
  std::pair<double, double> domain = _function_array.back()->get_domain();
  if (_deg == 0) {
    return std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::COMPOSITION,
        _function_array);
  }

  std::unique_ptr<Function> result =
      first_deriv_compose_functions(_function_array);

  for (std::size_t k = 1; k <= _deg; k++)
    result = std::move(result->deriv());

  return result;
}

std::unique_ptr<Function>
deriv_concat_functions(std::list<std::unique_ptr<Function>> &_function_array,
                       std::size_t _deg) {

  std::list<std::unique_ptr<Function>> result_array;
  for (std::unique_ptr<Function> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::CONCATENATION, result_array);
}
} // namespace functions
} // namespace gsplines
