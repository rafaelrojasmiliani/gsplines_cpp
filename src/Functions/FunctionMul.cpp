
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <iostream>
namespace gsplines {
namespace functions {
void compatibility_mul(const FunctionExpression &_f1,
                       const FunctionExpression &_f2) {

  if (not FunctionBase::same_domain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different domains");
  }

  if (_f1.get_codom_dim() > 1 and _f2.get_codom_dim() > 1) {
    throw std::invalid_argument(
        "At most one function can have vectorial value");
  }
}

const FunctionExpression &
return_first_or_max_codom_dim(const FunctionExpression &_f1,
                              const FunctionExpression &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f1;

  return _f2;
}

FunctionExpression &return_first_or_max_codom_dim(FunctionExpression &_f1,
                                                  FunctionExpression &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f1;

  return _f2;
}

const FunctionExpression &
return_second_or_mim_codom_dim(const FunctionExpression &_f1,
                               const FunctionExpression &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f2;

  return _f1;
}

FunctionExpression
FunctionExpression::operator*(const FunctionExpression &_that) const & {

  printf("HEEEEEEEEEEE ---- \n ");
  compatibility_mul(*this, _that);

  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);
  const FunctionExpression &f_scalar =
      return_second_or_mim_codom_dim(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (f_vector.get_type() == MULTIPLICATION) {
    std::transform(f_vector.function_array_.begin(),
                   f_vector.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(f_vector.clone());
  }

  if (f_scalar.get_type() == MULTIPLICATION) {
    std::transform(f_scalar.function_array_.begin(),
                   f_scalar.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionExpression> &element) {
                     return element->clone();
                   });
  } else {
    result_array.push_back(f_scalar.clone());
  }

  return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator*(FunctionExpression &&_that) const & {

  printf("KKKKKKKKKKKKKKKKKKKKKEEEEEEEEEEE ---- \n ");
  compatibility_mul(*this, _that);

  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);
  const FunctionExpression &f_scalar =
      return_second_or_mim_codom_dim(*this, _that);

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  if (_that.get_type() == FunctionExpression::Type::MULTIPLICATION) {

    std::move(std::begin(_that.function_array_),
              std::end(_that.function_array_),
              std::back_inserter(result_array));
  } else {
    result_array.push_back(_that.move_clone());
  }

  if (get_codom_dim() >= _that.get_codom_dim()) {
    if (get_type() == MULTIPLICATION) {
      std::transform(function_array_.begin(), function_array_.end(),
                     std::front_inserter(result_array),
                     [](const std::unique_ptr<FunctionExpression> &element) {
                       return element->clone();
                     });
    } else {
      result_array.push_front(this->clone());
    }
  } else {
    if (get_type() == MULTIPLICATION) {
      std::transform(function_array_.begin(), function_array_.end(),
                     std::back_inserter(result_array),
                     [](const std::unique_ptr<FunctionExpression> &element) {
                       return element->clone();
                     });
    } else {
      result_array.push_back(this->clone());
    }
  }

  return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator*(const FunctionExpression &_that) && {

  compatibility_mul(*this, _that);

  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);
  const FunctionExpression &f_scalar =
      return_second_or_mim_codom_dim(*this, _that);

  if (get_codom_dim() >= _that.get_codom_dim()) {
    if (get_type() == MULTIPLICATION) {
      if (_that.get_type() == MULTIPLICATION) {
        std::transform(_that.function_array_.begin(),
                       _that.function_array_.end(),
                       std::back_inserter(function_array_),
                       [](const std::unique_ptr<FunctionExpression> &element) {
                         return element->clone();
                       });
      } else {
        function_array_.push_back(_that.clone());
      }
      return std::move(*this);

    } else {
      std::list<std::unique_ptr<FunctionExpression>> result_array;

      std::pair<double, double> domain = get_domain();
      std::size_t codom_dim = get_codom_dim();

      result_array.push_back(this->move_clone());

      if (_that.get_type() == MULTIPLICATION) {
        std::transform(_that.function_array_.begin(),
                       _that.function_array_.end(),
                       std::back_inserter(result_array),
                       [](const std::unique_ptr<FunctionExpression> &element) {
                         return element->clone();
                       });
      } else {
        result_array.push_back(_that.clone());
      }
      return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                                FunctionExpression::Type::MULTIPLICATION,
                                std::move(result_array));
    }
  } else {
    if (get_type() == MULTIPLICATION) {
      if (_that.get_type() == MULTIPLICATION) {
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

    } else {
      std::list<std::unique_ptr<FunctionExpression>> result_array;
      result_array.push_back(this->move_clone());
      if (_that.get_type() == MULTIPLICATION) {
        std::transform(_that.function_array_.rbegin(),
                       _that.function_array_.rend(),
                       std::front_inserter(result_array),
                       [](const std::unique_ptr<FunctionExpression> &element) {
                         return element->clone();
                       });
      } else {
        result_array.push_front(_that.clone());
      }
      return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                                FunctionExpression::Type::MULTIPLICATION,
                                std::move(result_array));
    }
  }
  throw(std::invalid_argument(
      "This condition should not happen."
      " All Multiplication cases should be addressed FunctionMul.cpp"));
}

FunctionExpression
FunctionExpression::operator*(FunctionExpression &&_that) && {

  compatibility_mul(*this, _that);

  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);
  const FunctionExpression &f_scalar =
      return_second_or_mim_codom_dim(*this, _that);

  if (get_codom_dim() >= _that.get_codom_dim()) {
    if (get_type() == MULTIPLICATION) {
      if (_that.get_type() == MULTIPLICATION) {
        std::move(_that.function_array_.begin(), _that.function_array_.end(),
                  std::back_inserter(function_array_));
      } else {
        function_array_.push_back(_that.move_clone());
      }
      return std::move(*this);

    } else {
      std::list<std::unique_ptr<FunctionExpression>> result_array;
      result_array.push_back(this->move_clone());
      if (_that.get_type() == MULTIPLICATION) {
        std::move(_that.function_array_.begin(), _that.function_array_.end(),
                  std::back_inserter(result_array));
      } else {
        result_array.push_back(_that.move_clone());
      }
      return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                                FunctionExpression::Type::MULTIPLICATION,
                                std::move(result_array));
    }
  } else {
    if (get_type() == MULTIPLICATION) {
      if (_that.get_type() == MULTIPLICATION) {
        std::move(_that.function_array_.rbegin(), _that.function_array_.rend(),
                  std::front_inserter(function_array_));
      } else {
        function_array_.push_front(_that.move_clone());
      }
      return std::move(*this);

    } else {
      std::list<std::unique_ptr<FunctionExpression>> result_array;
      result_array.push_back(this->move_clone());

      if (_that.get_type() == MULTIPLICATION) {
        std::move(_that.function_array_.rbegin(), _that.function_array_.rend(),
                  std::front_inserter(result_array));
      } else {
        result_array.push_front(_that.move_clone());
      }
      return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                                FunctionExpression::Type::MULTIPLICATION,
                                std::move(result_array));
    }
  }
  throw(std::invalid_argument(
      "This condition should not happen."
      " All Multiplication cases should be addressed FunctionMul.cpp"));
}

FunctionExpression FunctionExpression::operator-() const & {

  return ConstFunction(get_domain(), 1, -1.0) * (*this);
}

FunctionExpression FunctionExpression::operator-() && {

  return ConstFunction(get_domain(), 1, -1.0) * std::move(*this);
}

FunctionExpression operator*(double _value, const FunctionExpression &_that) {

  return ConstFunction(_that.get_domain(), 1, _value) * _that;
}
FunctionExpression operator*(double _value, FunctionExpression &&_that) {

  return ConstFunction(_that.get_domain(), 1, _value) * std::move(_that);
}

/* -----
 *  FunctionExpression Evaluation
 * -----*/

void eval_mul_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) {
  // NOTE: the first element of _function_array has larger codomain dimension

  Eigen::MatrixXd temp(_domain_points.size(), 1);

  _function_array.front()->value(_domain_points, _result);

  std::list<std::unique_ptr<FunctionExpression>>::const_iterator it;
  for (it = std::next(_function_array.begin(), 1); it != _function_array.end();
       it++) {
    (*it)->value(_domain_points, temp);
    _result.array().colwise() *= temp.col(0).array();
  }
}

// http://
// scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1352&context=rhumj

std::unique_ptr<FunctionExpression> first_deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array) {

  std::list<std::unique_ptr<FunctionExpression>> result_array;
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  std::list<std::unique_ptr<FunctionExpression>>::const_iterator it_1 =
      _function_array.begin();
  std::list<std::unique_ptr<FunctionExpression>>::const_iterator it_2;

  std::list<std::unique_ptr<FunctionExpression>> elem_array_1;

  elem_array_1.push_back((*it_1)->deriv());

  for (it_2 = std::next(_function_array.begin(), 1);
       it_2 != _function_array.end(); it_2++) {
    elem_array_1.push_back((*it_2)->clone());
  }

  result_array.push_back(std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
      std::move(elem_array_1)));

  for (it_1 = std::next(_function_array.begin(), 1);
       it_1 != _function_array.end(); it_1++) {

    std::list<std::unique_ptr<FunctionExpression>> elem_array;

    for (it_2 = _function_array.begin(); it_2 != _function_array.end();
         it_2++) {
      if (it_1 != it_2)
        elem_array.push_back((*it_2)->clone());
    }

    elem_array.push_back((*it_1)->deriv());

    result_array.push_back(std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
        std::move(elem_array)));
  }
  return std::make_unique<FunctionExpression>(domain, codom_dim,
                                              FunctionExpression::Type::SUM,
                                              std::move(result_array));
}
std::unique_ptr<FunctionExpression> deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionExpression>> &_function_array,
    std::size_t _deg) {

  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  if (_deg == 0) {
    return std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
        _function_array);
  }

  std::unique_ptr<FunctionExpression> result =
      first_deriv_mul_functions(_function_array);

  for (std::size_t k = 1; k < _deg; k++) {
    result = result->deriv()->move_clone();
  }

  return result;
}
} // namespace functions
} // namespace gsplines
