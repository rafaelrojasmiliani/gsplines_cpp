#include <gsplines++/Functions/FunctionExpression.hpp>
#include<iostream>
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
    std::vector<std::unique_ptr<Function>> &_function_array)
    : Function(_domain, _codom_dim), type_(_type), function_array_(), eval_operation_(nullptr), deriv_operation_(nullptr) {

  for (std::size_t i = 0; i < _function_array.size() ; i++) {
    function_array_.push_back(_function_array[i]->clone());
  }
  switch(type_){
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
    std::vector<std::unique_ptr<Function>> &&_function_array)
    : Function(_domain, _codom_dim), type_(_type), eval_operation_(nullptr), deriv_operation_(nullptr), function_array_(std::move(_function_array)) {

  switch(type_){
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

        printf("calling wrong constructor \n-----\n"); fflush(stdout);
    }

FunctionExpression::FunctionExpression(const FunctionExpression &that)
    : Function(that), type_(that.type_), eval_operation_(that.eval_operation_), deriv_operation_(that.deriv_operation_){

        printf("copy  CONSTRUCTOR ----------\n");
  for (std::size_t i = 0; i < that.function_array_.size(); i++) {
    function_array_.push_back(that.function_array_[i]->clone());
   }
  num_call_copy_constructor_++;

    }



FunctionExpression::FunctionExpression(FunctionExpression &&that)
    : Function(that), type_(that.type_), eval_operation_(that.eval_operation_), deriv_operation_(that.deriv_operation_), function_array_(std::move(that.function_array_)){

        printf("MOOOOOOOOOOOOOOOOOVEEE  CONSTRUCTOR ----------\n");
        fflush(stdout);
        num_call_move_constructor_++;
    }

void FunctionExpression::print_performace(){
    printf("Num call of the constructor %zu\n", num_call_constructor_);
    printf("Num call of the COPY constructor %zu\n", num_call_copy_constructor_);
    printf("Num call of the MOVE constructor %zu\n", num_call_move_constructor_);
    printf("Num call of the SIMPLE constructor %zu\n", num_call_simple_constructor_);

}

/* -----
 *  Function Evaluation
 * -----*/

Eigen::MatrixXd
eval_compose_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                       const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array.back()->get_codom_dim());
  Eigen::VectorXd domain_ponts = _domain_points;
  for (std::size_t i = 0; i < _function_array.size() - 1; i++) {
    domain_ponts = _function_array[i]->value(domain_ponts);
  }
  result = _function_array.back()->value(_domain_points);
  return result;
}

Eigen::MatrixXd
eval_concat_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                      const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array[0]->get_codom_dim());

  std::size_t j = 0;
  for (std::size_t i = 0; i < _domain_points.size(); i++) {
    if (_domain_points(i) >= _function_array[j]->get_domain().second)
      j++;
    result.row(i) = _function_array[j]->value(_domain_points.row(i));
  }
  return result;
}

/* -----
 *  Function Derivation
 * -----*/
std::unique_ptr<Function> first_deriv_compose_functions(
    std::vector<std::unique_ptr<Function>> &_function_array) {

  std::vector<std::unique_ptr<Function>> result_array;

  result_array.push_back(_function_array[0]->deriv());
  for (std::size_t i = 1; i < _function_array.size(); i++) {
    std::vector<std::unique_ptr<Function>> elem_array;
    for (std::size_t j = 0; j < i; j++) {
      elem_array.push_back(_function_array[j]->clone());
    }
    std::size_t codom_dim = _function_array[i]->get_codom_dim();
    std::pair<double, double> domain = _function_array[i]->get_domain();
    elem_array.push_back(_function_array[i]->deriv());
    result_array.push_back(std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::COMPOSITION, std::move(elem_array)));
  }
  std::size_t codom_dim = _function_array.back()->get_codom_dim();
  std::pair<double, double> domain = _function_array.back()->get_domain();
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
      std::move(result_array));
}

std::unique_ptr<Function>
deriv_compose_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                        std::size_t _deg) {

  std::size_t codom_dim = _function_array[0]->get_codom_dim();
  std::pair<double, double> domain = _function_array[0]->get_domain();
  if (_deg == 0) {
    return std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::COMPOSITION,
        _function_array);
  }
  std::unique_ptr<Function> result =
      first_deriv_compose_functions(_function_array);
  for (std::size_t k = 1; k <= _deg; k++)
    result.reset(result->deriv().get());

  return result;
}

std::unique_ptr<Function>
deriv_concat_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                    std::size_t _deg) {

  std::vector<std::unique_ptr<Function>> result_array;
  for (std::unique_ptr<Function> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }
  std::size_t codom_dim = _function_array[0]->get_codom_dim();
  std::pair<double, double> domain = _function_array[0]->get_domain();
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::CONCATENATION, result_array);
}
} // namespace functions
} // namespace gsplines
