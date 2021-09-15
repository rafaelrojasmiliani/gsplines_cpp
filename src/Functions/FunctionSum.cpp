
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

  std::list<std::unique_ptr<FunctionBase>> result_array =
      const_const_operation_handler(*this, _that,
                                    FunctionExpression::Type::SUM);
  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator+(FunctionExpression &&_that) const & {

  sum_throw(*this, _that);

  std::list<std::unique_ptr<FunctionBase>> result_array =
      const_nonconst_operation_handler(*this, std::move(_that),
                                       FunctionExpression::Type::SUM);

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator+(const FunctionExpression &_that) && {

  sum_throw(*this, _that);

  std::list<std::unique_ptr<FunctionBase>> result_array =
      nonconst_const_operation_handler(std::move(*this), _that,
                                       FunctionExpression::Type::SUM);

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator+(FunctionExpression &&_that) && {

  sum_throw(*this, _that);

  std::list<std::unique_ptr<FunctionBase>> result_array =
      nonconst_nonconst_operation_handler(std::move(*this), std::move(_that),
                                          FunctionExpression::Type::SUM);

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::SUM,
                            std::move(result_array));
}
/**
 * SUBSTRACTION
 */

FunctionExpression
FunctionExpression::operator-(const FunctionExpression &_that) const & {

  sum_throw(*this, _that);
  return *this + (-_that);
}

FunctionExpression
FunctionExpression::operator-(FunctionExpression &&_that) const & {

  sum_throw(*this, _that);
  return *this + (-std::move(_that));
}

FunctionExpression
FunctionExpression::operator-(const FunctionExpression &_that) && {

  sum_throw(*this, _that);
  return std::move(*this) + (-_that);
}
FunctionExpression
FunctionExpression::operator-(FunctionExpression &&_that) && {

  sum_throw(*this, _that);
  return std::move(*this) + (-std::move(_that));
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

void FunctionExpression::operator+=(FunctionExpression &&that) {
  Type initial_type = get_type();
  if (initial_type == UNIQUE) {
    type_ = SUM;
    initialize();
  } else if (initial_type != SUM) {
    FunctionExpression aux = FunctionExpression(std::move(*this));
    function_array_.clear();
    type_ = SUM;
    initialize();
    function_array_.push_back(aux.move_clone());
  }
  function_array_.push_back(that.move_clone());
}

void FunctionExpression::operator+=(const FunctionExpression &that) {

  Type initial_type = get_type();
  if (initial_type == UNIQUE) {
    type_ = SUM;
    initialize();
  } else if (initial_type != SUM) {
    FunctionExpression aux = FunctionExpression(std::move(*this));
    function_array_.clear();
    type_ = SUM;
    initialize();
    function_array_.push_back(aux.move_clone());
  }
  function_array_.push_back(that.clone());
}
} // namespace functions
} // namespace gsplines
