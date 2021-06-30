#include <boost/math/special_functions/binomial.hpp>
#include <gsplines++/Functions/FunctionExpression.hpp>
namespace gsplines {
namespace functions {

/* -----
 *  Function Construction
 * -----*/

FunctionExpression::FunctionExpression(
    std::pair<double, double> _domain, std::size_t _codom_dim, Type _type,
    std::vector<std::unique_ptr<Function>> &_function_array)
    : Function(_domain, _codom_dim), type_(_type), function_array_() {

  for (std::size_t i = 0; i < _function_array.size() - 1; i++) {
    function_array_.push_back(_function_array[i]->clone());
  }
}

FunctionExpression operator+(const Function &_f1, const Function &_f2) {

  if (not FunctionBase::same_domain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different domains");
  }
  if (not FunctionBase::same_codomain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different codomains");
  }

  std::vector<std::unique_ptr<Function>> result_array;
  std::size_t codom_dim = _f1.get_codom_dim();
  std::pair<double, double> domain = _f1.get_domain();

  result_array.push_back(_f1.clone());
  result_array.push_back(_f2.clone());

  return FunctionExpression(domain, codom_dim, FunctionExpression::Type::SUM,
                            result_array);
}

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
                            result_array);
}
/* -----
 *  Function Evaluation
 * -----*/
Eigen::MatrixXd
eval_sum_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                   const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array[0]->get_codom_dim());
  result.setZero();
  for (std::unique_ptr<Function> &f : _function_array) {
    result += f->value(_domain_points);
  }
  return result;
}

Eigen::MatrixXd
eval_mul_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                   const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array[0]->get_codom_dim());
  result = _function_array[0]->value(_domain_points);
  for (std::size_t i = 1; i < _function_array.size(); i++) {
    Eigen::MatrixXd scalar_res = _function_array[i]->value(_domain_points);
    result = result.array().colwise() * scalar_res.col(0).array();
  }
  return result;
}

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
std::unique_ptr<Function>
deriv_sum_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                    std::size_t _deg) {

  std::vector<std::unique_ptr<Function>> result_array;
  for (std::unique_ptr<Function> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }
  std::size_t codom_dim = _function_array[0]->get_codom_dim();
  std::pair<double, double> domain = _function_array[0]->get_domain();
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::SUM, result_array);
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
        elem_array.push_back(_function_array[k]);
    }
    result_array.push_back(std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
        result_array));
  }
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::SUM, result_array);
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

std::unique_ptr<Function> first_deriv_compose_functions(
    std::vector<std::unique_ptr<Function>> &_function_array) {

  std::vector<std::unique_ptr<Function>> result_array;
  std::size_t codom_dim = _function_array.back()->get_codom_dim();
  std::pair<double, double> domain = _function_array.back()->get_domain();

  for (std::size_t i = _function_array.size() - 1; i >= 0; i--) {
    result_array.push_back(_function_array[i]->deriv());
  }
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
      result_array);
}
std::unique_ptr<Function>
deriv_compose_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                        std::size_t _deg) {}

std::unique_ptr<Function>
deriv_concat_functions(std::vector<std::unique_ptr<Function>> &_function_array,
                       std::size_t _deg);

/*
 */
} // namespace functions
} // namespace gsplines
