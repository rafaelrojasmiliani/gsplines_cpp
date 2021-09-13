
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
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

  printf("FunctionExpression::operator*(const FunctionExpression &_that) const "
         "&\n");
  compatibility_mul(*this, _that);

  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);
  const FunctionExpression &f_scalar =
      return_second_or_mim_codom_dim(*this, _that);

  std::list<std::unique_ptr<FunctionBase>> result_array =
      const_const_operation_handler(f_vector, f_scalar,
                                    FunctionExpression::Type::MULTIPLICATION);

  return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator*(FunctionExpression &&_that) const & {

  compatibility_mul(*this, _that);

  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = f_vector.get_codom_dim();

  std::list<std::unique_ptr<FunctionBase>> result_array =
      (get_codom_dim() >= _that.get_codom_dim())
          ? const_nonconst_operation_handler(
                *this, std::move(_that),
                FunctionExpression::Type::MULTIPLICATION)
          : nonconst_const_operation_handler(
                std::move(_that), *this,
                FunctionExpression::Type::MULTIPLICATION);

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator*(const FunctionExpression &_that) && {

  compatibility_mul(*this, _that);
  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = f_vector.get_codom_dim();

  std::list<std::unique_ptr<FunctionBase>> result_array =
      (get_codom_dim() >= _that.get_codom_dim())
          ? nonconst_const_operation_handler(
                std::move(*this), _that,
                FunctionExpression::Type::MULTIPLICATION)
          : const_nonconst_operation_handler(
                _that, std::move(*this),
                FunctionExpression::Type::MULTIPLICATION);

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}

FunctionExpression
FunctionExpression::operator*(FunctionExpression &&_that) && {

  compatibility_mul(*this, _that);
  const FunctionExpression &f_vector =
      return_first_or_max_codom_dim(*this, _that);

  std::pair<double, double> domain = get_domain();
  std::size_t codom_dim = f_vector.get_codom_dim();

  std::list<std::unique_ptr<FunctionBase>> result_array =
      (get_codom_dim() >= _that.get_codom_dim())
          ? nonconst_nonconst_operation_handler(
                std::move(*this), std::move(_that),
                FunctionExpression::Type::MULTIPLICATION)
          : nonconst_nonconst_operation_handler(
                std::move(_that), std::move(*this),
                FunctionExpression::Type::MULTIPLICATION);

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}
FunctionExpression FunctionExpression::operator-() const & {

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::NEGATIVE,
                            function_array_);
}

FunctionExpression FunctionExpression::operator-() && {

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::NEGATIVE,
                            std::move(function_array_));
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
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) {
  // NOTE: the first element of _function_array has larger codomain dimension

  Eigen::MatrixXd temp(_domain_points.size(), 1);

  _function_array.front()->value(_domain_points, _result);

  std::list<std::unique_ptr<FunctionBase>>::const_iterator it;
  for (it = std::next(_function_array.begin(), 1); it != _function_array.end();
       it++) {
    (*it)->value(_domain_points, temp);
    _result.array().colwise() *= temp.col(0).array();
  }
}

// http://
// scholar.rose-hulman.edu/cgi/viewcontent.cgi?article=1352&context=rhumj

FunctionExpression *first_deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array) {

  std::list<std::unique_ptr<FunctionBase>> result_array;
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  std::list<std::unique_ptr<FunctionBase>>::const_iterator it_1 =
      _function_array.begin();
  std::list<std::unique_ptr<FunctionBase>>::const_iterator it_2;

  std::list<std::unique_ptr<FunctionBase>> elem_array_1;

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

    std::list<std::unique_ptr<FunctionBase>> elem_array;

    for (it_2 = _function_array.begin(); it_2 != _function_array.end();
         it_2++) {
      if (it_1 != it_2)
        elem_array.push_back((*it_2)->clone());
    }

    elem_array.push_back((*it_1)->deriv());

    // printf("multiplication derivative domain = [%+11.3lf %+11.3lf]\n",
    //       domain.first, domain.second);
    // fflush(stdout);
    result_array.push_back(std::make_unique<FunctionExpression>(
        domain, codom_dim, FunctionExpression::Type::MULTIPLICATION,
        std::move(elem_array)));
  }
  return new FunctionExpression(domain, codom_dim,
                                FunctionExpression::Type::SUM,
                                std::move(result_array));
}

FunctionExpression *deriv_mul_functions(
    const std::list<std::unique_ptr<FunctionBase>> &_function_array,
    std::size_t _deg) {

  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  if (_deg == 0) {
    return new FunctionExpression(domain, codom_dim,
                                  FunctionExpression::Type::MULTIPLICATION,
                                  _function_array);
  }

  FunctionExpression *result = first_deriv_mul_functions(_function_array);
  FunctionExpression *buffer_pointer = result;

  for (std::size_t k = 1; k < _deg; k++) {
    result = result->deriv()->move_clone().release();
    delete buffer_pointer;
    buffer_pointer = result;
  }
  buffer_pointer = nullptr;

  return result;
}
} // namespace functions
} // namespace gsplines
