

#include <gsplines++/Functions/FunctionExpression.hpp>
namespace gsplines {
namespace functions {

void comp_throw(const Function &_f1, const Function &_f2) {

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

FunctionExpression Function::compose(const Function &that) const {

  comp_throw(*this, that);

  std::list<std::unique_ptr<Function>> result_array;

  result_array.push_back(that.clone());
  result_array.push_back(this->clone());

  return FunctionExpression(that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression Function::compose(const FunctionExpression &that) const {

  comp_throw(*this, that);

  std::list<std::unique_ptr<Function>> result_array;

  if (that.get_type() == FunctionExpression::Type::COMPOSITION) {
    for (const std::unique_ptr<Function> &f : that.function_array_)
      result_array.push_back(f->clone());

  } else {
    result_array.push_back(that.clone());
  }

  result_array.push_back(this->clone());

  return FunctionExpression(that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression Function::compose(FunctionExpression &&that) const {

  comp_throw(*this, that);

  std::list<std::unique_ptr<Function>> result_array;

  if (that.get_type() == FunctionExpression::Type::COMPOSITION) {

    std::move(std::begin(that.function_array_), std::end(that.function_array_),
              std::back_inserter(result_array));

  } else {
    result_array.push_back(that.clone());
  }

  result_array.push_back(this->clone());

  return FunctionExpression(that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression FunctionExpression::compose(const Function &that) const {

  comp_throw(*this, that);

  if (get_type() == COMPOSITION) {
    FunctionExpression self_copy(*this);
    self_copy.function_array_.push_front(that.clone());
    return self_copy;
  }

  std::list<std::unique_ptr<Function>> result_array;
  result_array.push_back(that.clone());
  result_array.push_back(this->clone());

  return FunctionExpression(that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::compose(const FunctionExpression &that) const {

  comp_throw(*this, that);

  if (get_type() == COMPOSITION) {
    FunctionExpression self_copy(*this);
    if (that.get_type() == COMPOSITION) {
      std::list<std::unique_ptr<Function>>::const_reverse_iterator it;
      for (it = that.function_array_.rbegin();
           it != that.function_array_.rend(); ++it) {
        self_copy.function_array_.push_front((*it)->clone());
      }
    } else {
      self_copy.function_array_.push_front(that.clone());
    }
    return self_copy;
  }

  std::list<std::unique_ptr<Function>> result_array;

  if (that.get_type() == COMPOSITION) {
    for (const std::unique_ptr<Function> &f : that.function_array_) {
      result_array.push_back(f->clone());
    }
  } else {
    result_array.push_back(that.clone());
  }

  return FunctionExpression(that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::compose(FunctionExpression &&that) const {

  comp_throw(*this, that);

  if (get_type() == COMPOSITION) {
    FunctionExpression self_copy(*this);
    if (that.get_type() == COMPOSITION) {
      std::move(std::rbegin(that.function_array_),
                std::rend(that.function_array_),
                std::front_inserter(self_copy.function_array_));
    } else {
      self_copy.function_array_.push_front(
          std::make_unique<FunctionExpression>(std::move(that)));
    }
    return self_copy;
  }

  std::list<std::unique_ptr<Function>> result_array;

  if (that.get_type() == COMPOSITION) {
    std::move(std::begin(that.function_array_), std::end(that.function_array_),
              std::back_inserter(result_array));
  } else {
    result_array.push_front(
        std::make_unique<FunctionExpression>(std::move(that)));
  }

  return FunctionExpression(that.get_domain(), get_codom_dim(),
                            FunctionExpression::Type::COMPOSITION,
                            std::move(result_array));
}

/* -----
 *  Evaluation method
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
} // namespace functions
} // namespace gsplines
