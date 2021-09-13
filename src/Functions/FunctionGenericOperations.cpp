
#include <gsplines/Functions/FunctionExpression.hpp>

namespace gsplines {
namespace functions {

std::list<std::unique_ptr<FunctionBase>>
FunctionExpression::const_const_operation_handler(
    const FunctionExpression &_first, const FunctionExpression &_second,
    FunctionExpression::Type _opt_type) {
  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (_first.get_type() == _opt_type or
      _first.get_type() == FunctionExpression::UNIQUE) {

    std::transform(_first.function_array_.begin(), _first.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  } else {

    result_array.push_back(_first.clone());
  }

  if (_second.get_type() == _opt_type or
      _second.get_type() == FunctionExpression::UNIQUE) {
    std::transform(_second.function_array_.begin(),
                   _second.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });

  } else {

    result_array.push_back(_second.clone());
  }

  return std::move(result_array);
}

std::list<std::unique_ptr<FunctionBase>>
FunctionExpression::const_nonconst_operation_handler(
    const FunctionExpression &_first, FunctionExpression &&_second,
    FunctionExpression::Type _opt_type) {

  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (_first.get_type() == _opt_type or
      _first.get_type() == FunctionExpression::UNIQUE) {

    std::transform(_first.function_array_.begin(), _first.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  } else {

    result_array.push_back(_first.clone());
  }

  if (_second.get_type() == _opt_type or
      _second.get_type() == FunctionExpression::UNIQUE) {

    std::move(_second.function_array_.begin(), _second.function_array_.end(),
              std::back_inserter(result_array));
  } else {

    result_array.push_back(_second.move_clone());
  }

  return std::move(result_array);
}

std::list<std::unique_ptr<FunctionBase>>
FunctionExpression::nonconst_const_operation_handler(
    FunctionExpression &&_first, const FunctionExpression &_second,
    FunctionExpression::Type _opt_type) {

  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (_first.get_type() == _opt_type or
      _first.get_type() == FunctionExpression::UNIQUE) {

    std::move(_first.function_array_.begin(), _first.function_array_.end(),
              std::back_inserter(result_array));
  } else {

    result_array.push_back(_first.move_clone());
  }
  if (_second.get_type() == _opt_type or
      _second.get_type() == FunctionExpression::UNIQUE) {

    std::transform(_second.function_array_.begin(),
                   _second.function_array_.end(),
                   std::back_inserter(result_array),
                   [](const std::unique_ptr<FunctionBase> &element) {
                     return element->clone();
                   });
  } else {

    result_array.push_back(_second.clone());
  }

  return std::move(result_array);
}

std::list<std::unique_ptr<FunctionBase>>
FunctionExpression::nonconst_nonconst_operation_handler(
    FunctionExpression &&_first, FunctionExpression &&_second,
    FunctionExpression::Type _opt_type) {

  std::list<std::unique_ptr<FunctionBase>> result_array;
  if (_first.get_type() == _opt_type or
      _first.get_type() == FunctionExpression::UNIQUE) {

    std::move(_first.function_array_.begin(), _first.function_array_.end(),
              std::back_inserter(result_array));
  } else {

    result_array.push_back(_first.move_clone());
  }
  if (_second.get_type() == _opt_type or
      _second.get_type() == FunctionExpression::UNIQUE) {

    std::move(_second.function_array_.begin(), _second.function_array_.end(),
              std::back_inserter(result_array));
  } else {

    result_array.push_back(_second.move_clone());
  }

  return std::move(result_array);
}

} // namespace functions
} // namespace gsplines
