

#include <gsplines++/Functions/FunctionExpression.hpp>
namespace gsplines {
namespace functions {

void comp_throw(const FunctionExpression &_f1, const FunctionExpression &_f2) {

  if (_f2.get_codom_dim() > 1) {
    throw std::invalid_argument(
        "Cannont compose a function with a vector value function");
  }
  /*
  if (not FunctionBase::same_codomain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different codomains");
  }
  */
}

FunctionExpression
FunctionExpression::compose(const FunctionExpression &_that) const & {

  comp_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (get_type() == COMPOSITION) {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(this->clone());
  }

  if (_that.get_type() == COMPOSITION) {
    std::transform(_that.function_array_.rbegin(), _that.function_array_.rend(),
                   std::front_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_front(_that.clone());
  }

  return FunctionExpression(_that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::compose(FunctionExpression &&_that) const & {

  comp_throw(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (get_type() == COMPOSITION) {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(this->clone());
  }

  if (_that.get_type() == COMPOSITION) {
    std::move(_that.function_array_.rbegin(), _that.function_array_.rend(),
              std::front_inserter(result_array));
  } else {
    result_array.push_front(_that.move_clone());
  }

  return FunctionExpression(_that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::compose(const FunctionExpression &_that) && {

  comp_throw(*this, _that);

  if (get_type() == COMPOSITION) {
    if (_that.get_type() == COMPOSITION) {
      std::transform(_that.function_array_.rbegin(),
                     _that.function_array_.rend(),
                     std::front_inserter(function_array_),
                     [](const std::unique_ptr<FunctionExpression> &element) {
                       return element->clone();
                     });

    } else {
      function_array_.push_front(_that.clone());
    }
    return std::move(*this);
  }

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  result_array.push_back(this->move_clone());

  if (_that.get_type() == COMPOSITION) {

    std::transform(_that.function_array_.rbegin(), _that.function_array_.rend(),
                   std::front_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_front(_that.clone());
  }

  return FunctionExpression(_that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression FunctionExpression::compose(FunctionExpression &&_that) && {

  comp_throw(*this, _that);

  if (get_type() == COMPOSITION) {
    if (_that.get_type() == COMPOSITION) {
      std::move(_that.function_array_.rbegin(), _that.function_array_.rend(),
                std::front_inserter(function_array_));

    } else {
      function_array_.push_front(_that.move_clone());
    }
    return std::move(*this);
  }

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  result_array.push_back(this->move_clone());

  std::pair<double, double> domain = _that.get_domain();
  std::size_t codom_dim = get_codom_dim();

  if (_that.get_type() == COMPOSITION) {

    std::move(_that.function_array_.rbegin(), _that.function_array_.rend(),
              std::front_inserter(result_array));

  } else {
    result_array.push_front(_that.move_clone());
  }

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

/* -----
 *  Evaluation method
 * -----*/
Eigen::MatrixXd eval_compose_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array.back()->get_codom_dim());

  Eigen::VectorXd domain_points_copy = _domain_points;

  std::list<std::unique_ptr<FunctionExpression>>::const_iterator it;

  std::list<std::unique_ptr<FunctionExpression>>::const_iterator it_limit =
      std::next(_function_array.end(), -1);

  for (it = _function_array.begin(); it != it_limit; it++) {
    domain_points_copy = (*it)->value(domain_points_copy);
  }

  result = _function_array.back()->value(domain_points_copy);

  return result;
}

/* -----
 *  FunctionExpression Derivation
 * -----*/
std::unique_ptr<FunctionExpression> first_deriv_compose_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array) {

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  result_array.push_back(_function_array.front()->deriv());

  std::list<std::unique_ptr<FunctionExpression>>::const_iterator it;
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

    std::list<std::unique_ptr<FunctionExpression>> elem_array;
    std::list<std::unique_ptr<FunctionExpression>>::const_iterator it_elem;

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

std::unique_ptr<FunctionExpression> deriv_compose_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg) {

  std::size_t codom_dim = _function_array.back()->get_codom_dim();
  std::pair<double, double> domain = _function_array.back()->get_domain();
  if (_deg == 0) {
    return std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::COMPOSITION,
        _function_array);
  }

  std::unique_ptr<FunctionExpression> result =
      first_deriv_compose_functions(_function_array);

  for (std::size_t k = 1; k < _deg; k++)
    result = result->deriv()->move_clone();

  return result;
}
} // namespace functions
} // namespace gsplines
