
#include <algorithm>
#include <gsplines++/Functions/FunctionExpression.hpp>
namespace gsplines {
namespace functions {

void concat_throw(const FunctionExpression &_f1,
                  const FunctionExpression &_f2) {

  if (not FunctionBase::same_codomain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different codomains");
  }

  if (abs(_f1.get_domain().second - _f2.get_domain().first) >
      FunctionBase::dom_tollerance_)
    throw std::invalid_argument("Functions with different codomains");
}

FunctionExpression
FunctionExpression::concat(const FunctionExpression &_that) const & {

  concat_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (get_type() == CONCATENATION) {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(this->clone());
  }

  if (_that.get_type() == CONCATENATION) {
    std::transform(_that.function_array_.begin(), _that.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(_that.clone());
  }
  return FunctionExpression(
      {get_domain().first, _that.get_domain().second}, get_codom_dim(),
      FunctionExpression::Type::CONCATENATION, std::move(result_array));
}

FunctionExpression
FunctionExpression::concat(FunctionExpression &&_that) const & {

  concat_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (get_type() == CONCATENATION) {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(this->clone());
  }

  if (_that.get_type() == CONCATENATION) {

    std::move(_that.function_array_.begin(), _that.function_array_.end(),
              std::back_inserter(result_array));
  } else {
    result_array.push_back(_that.move_clone());
  }
  return FunctionExpression(
      {get_domain().first, _that.get_domain().second}, get_codom_dim(),
      FunctionExpression::Type::CONCATENATION, std::move(result_array));
}

FunctionExpression
FunctionExpression::concat(const FunctionExpression &_that) && {

  concat_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  std::list<std::unique_ptr<FunctionExpression>> &target_array =
      (get_type() == CONCATENATION) ? function_array_ : result_array;

  if (_that.get_type() == CONCATENATION) {

    std::transform(_that.function_array_.begin(), _that.function_array_.end(),
                   std::back_inserter(target_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    target_array.push_back(_that.clone());
  }

  if (get_type() == CONCATENATION)
    return std::move(*this);

  result_array.push_front(this->move_clone());

  return FunctionExpression(
      {get_domain().first, _that.get_domain().second}, get_codom_dim(),
      FunctionExpression::Type::CONCATENATION, std::move(target_array));
}

FunctionExpression FunctionExpression::concat(FunctionExpression &&_that) && {

  concat_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  std::list<std::unique_ptr<FunctionExpression>> &target_array =
      (get_type() == CONCATENATION) ? function_array_ : result_array;

  if (_that.get_type() == CONCATENATION) {

    std::move(_that.function_array_.begin(), _that.function_array_.end(),
              std::back_inserter(target_array));
  } else {
    target_array.push_back(_that.move_clone());
  }

  if (get_type() == CONCATENATION)
    return std::move(*this);

  result_array.push_front(this->move_clone());
  return FunctionExpression(
      {get_domain().first, _that.get_domain().second}, get_codom_dim(),
      FunctionExpression::Type::CONCATENATION, std::move(target_array));
}

//  -----------------------------
//  FunctionExpression Evaluation
//  -----------------------------

std::list<std::unique_ptr<FunctionExpression>>::const_iterator
get_interval_function(
    double _domain_point,
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array) {

  std::list<std::unique_ptr<FunctionExpression>>::const_iterator result =
      std::find_if(
          _function_array.begin(), _function_array.end(),
          [&_domain_point](const std::unique_ptr<FunctionExpression> &element) {
            return element->is_point_in_domain(_domain_point);
          });

  return result;
}

std::unique_ptr<FunctionExpression> deriv_concat_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg) {

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  for (const std::unique_ptr<FunctionExpression> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }

  return std::make_unique<FunctionExpression>(
      _function_array.front()->get_domain(),
      _function_array.front()->get_codom_dim(),
      FunctionExpression::Type::CONCATENATION, result_array);
}

void eval_concat_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) {

  Eigen::VectorXd vector_one(1);

  Eigen::MatrixXd temp(1, _function_array.front()->get_codom_dim());

  for (std::size_t i = 0; i < _domain_points.size(); i++) {

    std::list<std::unique_ptr<FunctionExpression>>::const_iterator f =
        get_interval_function(_domain_points[i], _function_array);

    if (f != _function_array.end()) {
      (*f)->value(_domain_points.segment(i, 1), temp);
      _result.row(i) = temp.row(0);
    } else if (_domain_points[i] <
               _function_array.front()->get_domain().first) {
      _function_array.front()->value(Eigen::VectorXd::Constant(
          1, _function_array.front()->get_domain().first));
    } else {

      _function_array.front()->value(Eigen::VectorXd::Constant(
          1, _function_array.back()->get_domain().second));
    }
  }
}
} // namespace functions
} // namespace gsplines
