
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <iostream>
namespace gsplines {
namespace functions {

void sum_throw(const FunctionExpression &_f1, const FunctionExpression &_f2) {

  if (not FunctionBase::same_domain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different domains");
  }
  if (not FunctionBase::same_codomain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different codomains");
  }
}

FunctionExpression
FunctionExpression::operator+(const FunctionExpression &_that) const {

  sum_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;
  if (get_type() == SUM) {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(this->clone());
  }

  if (_that.get_type() == SUM) {
    std::transform(_that.function_array_.begin(), _that.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(_that.clone());
  }
  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator+(FunctionExpression &&_that) const {

  sum_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (get_type() == SUM) {

    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {

    result_array.push_back(this->clone());
  }

  if (_that.get_type() == SUM) {
    std::move(_that.function_array_.begin(), _that.function_array_.end(),
              std::back_inserter(result_array));
  } else {
    result_array.push_front(
        std::make_unique<FunctionExpression>(std::move(_that)));
  }
  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator-(const FunctionExpression &_that) const {

  sum_throw(*this, _that);
  std::list<std::unique_ptr<FunctionExpression>> result_array;
  if (get_type() != SUM) {
    result_array.push_back(this->clone());
  } else {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  }

  std::list<std::unique_ptr<FunctionExpression>> aux_array;
  aux_array.push_back(std::make_unique<FunctionExpression>(_that));
  aux_array.push_back(std::make_unique<ConstFunction>(get_domain(), 1, -1.0));

  result_array.push_back(std::make_unique<FunctionExpression>(
      get_domain(), get_codom_dim(), FunctionExpression::Type::MULTIPLICATION,
      std::move(aux_array)));

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator-(FunctionExpression &&_that) const {

  sum_throw(*this, _that);
  std::list<std::unique_ptr<FunctionExpression>> result_array;
  if (get_type() != SUM) {
    result_array.push_back(this->clone());
  } else {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  }

  std::list<std::unique_ptr<FunctionExpression>> aux_array;
  aux_array.push_back(std::make_unique<FunctionExpression>(std::move(_that)));
  aux_array.push_back(std::make_unique<ConstFunction>(get_domain(), 1, -1.0));

  result_array.push_back(std::make_unique<FunctionExpression>(
      get_domain(), get_codom_dim(), FunctionExpression::Type::MULTIPLICATION,
      std::move(aux_array)));

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

/* -----
 *  FunctionExpression Evaluation
 * -----*/
Eigen::MatrixXd eval_sum_functions(
    std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array.front()->get_codom_dim());
  result.setZero();
  for (const std::unique_ptr<FunctionExpression> &f : _function_array) {
    result += f->value(_domain_points);
  }
  return result;
}

/* -----
 *  FunctionExpression Derivation
 * -----*/
std::unique_ptr<FunctionExpression> deriv_sum_functions(
    std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg) {

  std::list<std::unique_ptr<FunctionExpression>> result_array;
  for (std::unique_ptr<FunctionExpression> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  return std::make_unique<FunctionExpression>(domain, codom_dim,
                                              FunctionExpression::Type::SUM,
                                              std::move(result_array));
}

} // namespace functions
} // namespace gsplines
