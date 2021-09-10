#ifndef FUNCTION
#define FUNCTION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <gsplines/Functions/FunctionBase.hpp>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class FunctionExpression;
class Function : public FunctionBase {

public:
  using FunctionBase::FunctionBase;
  Function(const Function &_that) : FunctionBase(_that) {}
  Function(Function &&_that) : FunctionBase(std::move(_that)) {}
  ~Function() = default;
  FunctionExpression operator+(const FunctionExpression &that) const &;
  FunctionExpression operator+(FunctionExpression &&that) const &;
  FunctionExpression operator+(const FunctionExpression &that) &&;
  FunctionExpression operator+(FunctionExpression &&that) &&;

  FunctionExpression operator-(const FunctionExpression &that) const &;
  FunctionExpression operator-(FunctionExpression &&that) const &;
  FunctionExpression operator-(const FunctionExpression &that) &&;
  FunctionExpression operator-(FunctionExpression &&that) &&;
  FunctionExpression operator-() const &;
  FunctionExpression operator-() &&;

  FunctionExpression operator*(const FunctionExpression &that) const &;
  FunctionExpression operator*(FunctionExpression &&that) const &;
  FunctionExpression operator*(const FunctionExpression &that) &&;
  FunctionExpression operator*(FunctionExpression &&that) &&;

  FunctionExpression compose(const FunctionExpression &that) const &;
  FunctionExpression compose(FunctionExpression &&that) const &;

  FunctionExpression compose(const FunctionExpression &that) &&;
  FunctionExpression compose(FunctionExpression &&that) &&;

  FunctionExpression concat(const FunctionExpression &that) const &;
  FunctionExpression concat(FunctionExpression &&that) const &;

  FunctionExpression concat(const FunctionExpression &that) &&;
  FunctionExpression concat(FunctionExpression &&that) &&;
};

} // namespace functions
} // namespace gsplines
#endif
