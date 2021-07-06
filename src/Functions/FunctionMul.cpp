
#include <gsplines++/Functions/FunctionExpression.hpp>
#include<iostream>
namespace gsplines {
namespace functions {

FunctionExpression operator*(const Function &_f1, const Function &_f2) {

  if (not FunctionBase::same_domain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different domains");
  }

  if (_f1.get_codom_dim() > 1 and _f2.get_codom_dim() > 1) {
    throw std::invalid_argument(
        "At most one function can have vectorial value");
  }

  std::pair<double, double> domain = _f1.get_domain();
  std::size_t codom_dim;
  std::vector<std::unique_ptr<Function>> result_array;

  if (_f1.get_codom_dim() > 1) {
    codom_dim = _f1.get_codom_dim();
    result_array.push_back(_f1.clone());
    result_array.push_back(_f2.clone());
  } else {
    codom_dim = _f2.get_codom_dim();
    result_array.push_back(_f2.clone());
    result_array.push_back(_f1.clone());
  }

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}
  
FunctionExpression operator*(FunctionExpression &&_f1, const Function &_f2) {

  if (not FunctionBase::same_domain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different domains");
  }

  if (_f1.get_codom_dim() > 1 and _f2.get_codom_dim() > 1) {
    throw std::invalid_argument(
        "At most one function can have vectorial value");
  }

  std::pair<double, double> domain = _f1.get_domain();
  std::size_t codom_dim;
  std::vector<std::unique_ptr<Function>> result_array;

  if (_f1.get_codom_dim() > 1) {
    codom_dim = _f1.get_codom_dim();
    result_array.push_back(std::make_unique<FunctionExpression>(std::move(_f1)));
    result_array.push_back(_f2.clone());
  } else {
    codom_dim = _f2.get_codom_dim();
    result_array.push_back(_f2.clone());
    result_array.push_back(std::make_unique<FunctionExpression>(std::move(_f1)));
  }

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}
FunctionExpression operator*(const Function &_f1, FunctionExpression &&_f2){
    return std::move(_f2) * _f1;
}
/* -----
 *  Function Evaluation
 * -----*/

Eigen::MatrixXd
eval_mul_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                   const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
  // NOTE: the first element of _function_array has larger codomain dimension

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array[0]->get_codom_dim());
  result = _function_array[0]->value(_domain_points);
  for (std::size_t i = 1; i < _function_array.size(); i++) {
    Eigen::MatrixXd scalar_res = _function_array[i]->value(_domain_points);
    result = result.array().colwise() * scalar_res.col(0).array();
  }
  return result;
}

// https://scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1352&context=rhumj
std::unique_ptr<Function> first_deriv_mul_functions(
    std::vector<std::unique_ptr<Function>> &_function_array) {

  std::vector<std::unique_ptr<Function>> result_array;
  std::size_t codom_dim = _function_array[0]->get_codom_dim();
  std::pair<double, double> domain = _function_array[0]->get_domain();
  for (std::size_t i = 0; i < _function_array.size(); i++) {

    std::vector<std::unique_ptr<Function>> elem_array;
    elem_array.push_back(_function_array[i]->deriv());
    for (std::size_t k = 0; k < _function_array.size(); k++) {
      if (k != i)
        elem_array.push_back(_function_array[k]->clone());
    }
    result_array.push_back(std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
        result_array));
  }
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::SUM, std::move(result_array));
}


std::unique_ptr<Function>
deriv_mul_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                    std::size_t _deg) {

  std::size_t codom_dim = _function_array[0]->get_codom_dim();
  std::pair<double, double> domain = _function_array[0]->get_domain();
  if (_deg == 0) {
    return std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
        _function_array);
  }
  std::unique_ptr<Function> result = first_deriv_mul_functions(_function_array);
  for (std::size_t k = 1; k <= _deg; k++)
    result.reset(result->deriv().get());

  return result;
}
} // namespace functions
} // namespace gsplines
