
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <iostream>
namespace gsplines {
namespace functions {

void compatibility_mul(const Function &_f1, const Function &_f2) {

  if (not FunctionBase::same_domain(_f1, _f2)) {
    throw std::invalid_argument("Functions with different domains");
  }

  if (_f1.get_codom_dim() > 1 and _f2.get_codom_dim() > 1) {
    throw std::invalid_argument(
        "At most one function can have vectorial value");
  }
}

const Function &return_first_or_max_codom_dim(const Function &_f1,
                                              const Function &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f1;

  return _f2;
}

Function &return_first_or_max_codom_dim(Function &_f1, Function &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f1;

  return _f2;
}

const Function &return_second_or_mim_codom_dim(const Function &_f1,
                                               const Function &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f2;

  return _f1;
}

Function &return_second_or_mim_codom_dim(Function &_f1, Function &_f2) {

  if (_f1.get_codom_dim() >= _f2.get_codom_dim())
    return _f2;

  return _f1;
}

void push_back_to_array(std::unique_ptr<Function> &_self_copy, Function &_that,
                        std::vector<std::unique_ptr<Function>> _fa) {

  if (_self_copy->get_codom_dim() >= _that.get_codom_dim()) {
    _fa.push_back(std::move(_self_copy));
    _fa.push_back(_that.clone());
  } else {
    _fa.push_back(_that.clone());
    _fa.push_back(std::move(_self_copy));
  }
}

FunctionExpression &FunctionExpression::operator*=(const Function &that) {

  compatibility_mul(*this, that);

  if (type_ == MULTIPLICATION) {
    if (this->get_codom_dim() >= that.get_codom_dim()) {
      function_array_.push_back(that.clone());
    } else {
      function_array_.insert(function_array_.begin(), that.clone());
    }
  } else {

    type_ = MULTIPLICATION;

    std::unique_ptr<Function> self_copy = std::make_unique<FunctionExpression>(
        get_domain(), get_codom_dim(), type_, std::move(function_array_));

    function_array_.clear();

    if (self_copy->get_codom_dim() >= that.get_codom_dim()) {
      function_array_.push_back(std::move(self_copy));
      function_array_.push_back(that.clone());
    } else {
      function_array_.push_back(that.clone());
      function_array_.push_back(std::move(self_copy));
    }

    eval_operation_ = eval_mul_functions;
    deriv_operation_ = deriv_mul_functions;
  }

  return *this;
}

FunctionExpression &
FunctionExpression::operator*=(const FunctionExpression &that) {

  compatibility_mul(*this, that);

  if (type_ == MULTIPLICATION and that.get_type() == MULTIPLICATION) {
    if (this->get_codom_dim() >= that.get_codom_dim()) {
      for (const std::unique_ptr<Function> &f : that.function_array_) {
        function_array_.push_back(f->clone());
      }
    } else {
      codom_dim_ = that.get_codom_dim();
      // concatenate lists 1. copy the that list, then concatenate at front of
      // function_array_.
    }
  } else if (type_ == MULTIPLICATION and that.get_type() != MULTIPLICATION) {

    if (this->get_codom_dim() >= that.get_codom_dim()) {
      function_array_.push_back(that.clone());
    } else {
      codom_dim_ = that.get_codom_dim();
      // concatenate lists push_front
    }

  } else {
    type_ = MULTIPLICATION;

    std::unique_ptr<Function> self_copy = std::make_unique<FunctionExpression>(
        get_domain(), get_codom_dim(), type_, std::move(function_array_));

    function_array_.clear();
  }

  return *this;
}

FunctionExpression &FunctionExpression::operator*=(FunctionExpression &&that) {

  compatibility_mul(*this, that);

  if (type_ != MULTIPLICATION) {

    type_ = MULTIPLICATION;

    std::unique_ptr<Function> self_copy = std::make_unique<FunctionExpression>(
        get_domain(), get_codom_dim(), MULTIPLICATION,
        std::move(function_array_));

    function_array_.clear();

    if (that.get_type() != MULTIPLICATION) {
      if (self_copy->get_codom_dim() >= that.get_codom_dim()) {
        function_array_.push_back(std::move(self_copy));
        function_array_.push_back(
            std::make_unique<FunctionExpression>(std::move(that)));
      } else {
        function_array_.push_back(
            std::make_unique<FunctionExpression>(std::move(that)));
        function_array_.push_back(std::move(self_copy));
      }
    } else {
      if (self_copy->get_codom_dim() >= that.get_codom_dim()) {
        function_array_.push_back(std::move(self_copy));
        function_array_.reserve(function_array_.size() +
                                that.function_array_.size());
        std::move(std::begin(that.function_array_),
                  std::end(that.function_array_),
                  std::back_inserter(function_array_));
      } else {
        function_array_.reserve(function_array_.size() +
                                that.function_array_.size());
        std::move(std::begin(that.function_array_),
                  std::end(that.function_array_),
                  std::back_inserter(function_array_));
        function_array_.push_back(std::move(self_copy));
      }
    }
  } else {
    if (that.type_ == MULTIPLICATION) {
      if (this->get_codom_dim() >= that.get_codom_dim()) {
        function_array_.reserve(function_array_.size() +
                                that.function_array_.size());
        std::move(std::begin(that.function_array_),
                  std::end(that.function_array_),
                  std::back_inserter(function_array_));
      } else {
        function_array_.insert(function_array_.begin(),
                               std::move(that.function_array_[0]));
        function_array_.reserve(function_array_.size() +
                                that.function_array_.size() - 1);
        std::move(std::begin(that.function_array_) + 1,
                  std::end(that.function_array_),
                  std::back_inserter(function_array_));
      }
    } else {
      if (this->get_codom_dim() >= that.get_codom_dim()) {
        function_array_.push_back(
            std::make_unique<FunctionExpression>(std::move(that)));
      } else {
        function_array_.insert(function_array_.begin(),
                               std::move(that.function_array_[0]));
        function_array_.reserve(function_array_.size() +
                                that.function_array_.size() - 1);
        std::move(std::begin(that.function_array_) + 1,
                  std::end(that.function_array_),
                  std::back_inserter(function_array_));
      }
    }
  }
  fflush(stdout);

  return *this;
}

FunctionExpression operator*(const Function &_f1, const Function &_f2) {

  compatibility_mul(_f1, _f2);

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

  compatibility_mul(_f1, _f2);

  std::pair<double, double> domain = _f1.get_domain();
  std::size_t codom_dim;
  std::vector<std::unique_ptr<Function>> result_array;

  if (_f1.get_type() == FunctionExpression::Type::MULTIPLICATION) {

    if (_f1.get_codom_dim() > 1) {
      codom_dim = _f1.get_codom_dim();
      result_array.push_back(
          std::make_unique<FunctionExpression>(std::move(_f1)));
      result_array.push_back(_f2.clone());
    } else {
      codom_dim = _f2.get_codom_dim();
      result_array.push_back(_f2.clone());
      result_array.push_back(
          std::make_unique<FunctionExpression>(std::move(_f1)));
    }
  }
  if (_f1.get_codom_dim() > 1) {
    codom_dim = _f1.get_codom_dim();
    result_array.push_back(
        std::make_unique<FunctionExpression>(std::move(_f1)));
    result_array.push_back(_f2.clone());
  } else {
    codom_dim = _f2.get_codom_dim();
    result_array.push_back(_f2.clone());
    result_array.push_back(
        std::make_unique<FunctionExpression>(std::move(_f1)));
  }

  return FunctionExpression(domain, codom_dim,
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
}

FunctionExpression operator*(const FunctionExpression &_f1,
                             const Function &_f2) {

  compatibility_mul(_f1, _f2);

  std::vector<std::unique_ptr<Function>> result_array;

  const Function &f_vector = return_first_or_max_codom_dim(_f1, _f2);
  const Function &f_scalar = return_second_or_mim_codom_dim(_f1, _f2);

  if (_f1.get_type() == FunctionExpression::Type::MULTIPLICATION) {

    if (_f1.get_codom_dim() > _f2.get_codom_dim()) {

      for (const std::unique_ptr<Function> &f : _f1.function_array_) {
        result_array.push_back(f->clone());
      }
      result_array.push_back(_f2.clone());

    } else {
      result_array.push_back(_f2.clone());
      for (const std::unique_ptr<Function> &f : _f1.function_array_) {
        result_array.push_back(f->clone());
      }
    }
  } else {

    result_array.push_back(f_vector.clone());
    result_array.push_back(f_scalar.clone());
  }

  return FunctionExpression(f_vector.get_domain(), f_vector.get_codom_dim(),
                            FunctionExpression::Type::MULTIPLICATION,
                            std::move(result_array));
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
  return std::make_unique<FunctionExpression>(domain, codom_dim,
                                              FunctionExpression::Type::SUM,
                                              std::move(result_array));
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
