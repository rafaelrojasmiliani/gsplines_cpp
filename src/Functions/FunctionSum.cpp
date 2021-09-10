
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
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
FunctionExpression::operator+(const FunctionExpression &_that) const & {

  sum_throw(*this, _that);

  // printf("AAAAAAAAAAAAAAAAA\n");

  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (get_type() == SUM) {
    // printf("this sum, other const& \n");
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  } else {
    // printf("this not sum, other const& \n");
    // printf("this name %s, that name %s\n", get_name().c_str(),
    //       _that.get_name().c_str());
    result_array.push_back(this->clone());
    // printf("cloned name %s \n", result_array.front()->get_name().c_str());
  }

  if (_that.get_type() == SUM) {
    std::transform(_that.function_array_.begin(), _that.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
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
FunctionExpression::operator+(FunctionExpression &&_that) const & {

  sum_throw(*this, _that);
  // printf("BBBBBBBBBBBBBBBBBB\n");

  std::list<std::unique_ptr<FunctionBase>> result_array;

  if (get_type() == SUM) {

    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  } else {

    result_array.push_back(this->clone());
  }

  if (_that.get_type() == SUM) {
    std::move(_that.function_array_.begin(), _that.function_array_.end(),
              std::back_inserter(result_array));
  } else {
    result_array.push_back(_that.move_clone());
  }
  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator+(const FunctionExpression &_that) && {

  sum_throw(*this, _that);
  // printf("FunctionExpression &&*this, const FunctionExpression &_that \n");
  if (get_type() == FunctionExpression::SUM) {
    if (_that.get_type() == FunctionExpression::SUM) {
      // printf("a-------\n");
      std::transform(_that.function_array_.begin(), _that.function_array_.end(),
                     std::back_inserter(function_array_),
                     [](const std::unique_ptr<FunctionBase> &element) {
                       return element->clone();
                     });
    } else {
      this->function_array_.push_back(_that.clone());
    }

    return std::move(*this);
  }
  std::list<std::unique_ptr<FunctionBase>> result_array;

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = get_codom_dim();

  result_array.push_back(this->move_clone());

  if (_that.get_type() == FunctionExpression::SUM) {
    std::transform(_that.function_array_.begin(), _that.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  } else
    result_array.push_back(_that.clone());

  return FunctionExpression(domain, codom_dim, FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator+(FunctionExpression &&_that) && {

  sum_throw(*this, _that);

  if (get_type() == FunctionExpression::SUM) {
    if (_that.get_type() == FunctionExpression::SUM) {
      std::move(_that.function_array_.begin(), _that.function_array_.end(),
                std::back_inserter(function_array_));
    } else {
      function_array_.push_back(this->move_clone());
    }

    return std::move(*this);
  }
  std::list<std::unique_ptr<FunctionBase>> result_array;

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = get_codom_dim();
  result_array.push_back(this->move_clone());

  // printf("WHERE ARE WE?? \n");
  if (_that.get_type() == FunctionExpression::SUM) {
    std::move(_that.function_array_.begin(), _that.function_array_.end(),
              std::back_inserter(result_array));
  } else
    result_array.push_back(_that.move_clone());

  return FunctionExpression(domain, codom_dim, FunctionExpression::Type::SUM,
                            std::move(result_array));
}
/**
 * SUBSTRACTION
 */

FunctionExpression
FunctionExpression::operator-(const FunctionExpression &_that) const & {

  sum_throw(*this, _that);
  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (get_type() != SUM) {
    result_array.push_back(this->clone());
  } else {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  }

  std::list<std::unique_ptr<FunctionBase>> aux_array;
  aux_array.push_back(_that.clone());
  aux_array.push_back(std::make_unique<ConstFunction>(get_domain(), 1, -1.0));

  result_array.push_back(std::make_unique<FunctionExpression>(
      get_domain(), get_codom_dim(), FunctionExpression::Type::MULTIPLICATION,
      std::move(aux_array)));

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator-(FunctionExpression &&_that) const & {

  sum_throw(*this, _that);
  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (get_type() != SUM) {
    result_array.push_back(this->clone());
  } else {
    std::transform(function_array_.begin(), function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  }

  std::list<std::unique_ptr<FunctionBase>> aux_array;
  aux_array.push_back(_that.move_clone());
  aux_array.push_back(std::make_unique<ConstFunction>(get_domain(), 1, -1.0));

  result_array.push_back(std::make_unique<FunctionExpression>(
      get_domain(), get_codom_dim(), FunctionExpression::Type::MULTIPLICATION,
      std::move(aux_array)));

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator-(const FunctionExpression &_that) && {

  sum_throw(*this, _that);

  if (get_type() == FunctionExpression::SUM) {
    function_array_.push_back(
        std::make_unique<FunctionExpression>(std::move(-_that)));

    return std::move(*this);
  }
  std::list<std::unique_ptr<FunctionBase>> result_array;

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = get_codom_dim();

  result_array.push_back(this->move_clone());

  function_array_.push_back(
      std::make_unique<FunctionExpression>(std::move(-_that)));

  return FunctionExpression(domain, codom_dim, FunctionExpression::Type::SUM,
                            std::move(result_array));
}
FunctionExpression
FunctionExpression::operator-(FunctionExpression &&_that) && {

  sum_throw(*this, _that);

  if (get_type() == FunctionExpression::SUM) {
    function_array_.push_back(
        std::make_unique<FunctionExpression>(std::move(-_that)));

    return std::move(*this);
  }
  std::list<std::unique_ptr<FunctionBase>> result_array;

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = get_codom_dim();
  result_array.push_back(this->move_clone());
  function_array_.push_back(
      std::make_unique<FunctionExpression>(std::move(-_that)));

  return FunctionExpression(domain, codom_dim, FunctionExpression::Type::SUM,
                            std::move(result_array));
}
/* -----
 *  FunctionExpression Evaluation
 * -----*/
void eval_sum_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) {
  Eigen::MatrixXd temp(_domain_points.size(),
                       _function_array.front()->get_codom_dim());
  _result.setZero();
  for (const std::unique_ptr<FunctionBase> &f : _function_array) {
    f->value(_domain_points, temp);
    _result += temp;
    // std::cout << "result \n " << _result << "\n ---\n";
    // printf("kkk\n");
  }
}

/* -----
 *  FunctionExpression Derivation
 * -----*/
FunctionExpression *deriv_sum_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg) {
  std::list<std::unique_ptr<FunctionBase>> result_array;
  for (const std::unique_ptr<FunctionBase> &f : _function_array) {
    // printf("func name = %s \n", f->get_name().c_str());
    result_array.push_back(f->deriv(_deg));
    // printf("func name = %s \n", result_array.back()->get_name().c_str());
  }
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  return new FunctionExpression(domain, codom_dim,
                                FunctionExpression::Type::SUM,
                                std::move(result_array));
}

} // namespace functions
} // namespace gsplines
