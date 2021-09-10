
#include <cmath>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <stdexcept>

namespace gsplines {
namespace functions {

FunctionExpression Function::operator+(const FunctionExpression &that) const & {
  return FunctionExpression(*this) + that;
}
FunctionExpression Function::operator+(FunctionExpression &&that) const & {
  return FunctionExpression(*this) + std::move(that);
}
FunctionExpression Function::operator+(const FunctionExpression &that) && {
  return FunctionExpression(std::move(*this)) + that;
}
FunctionExpression Function::operator+(FunctionExpression &&that) && {
  return FunctionExpression(std::move(*this)) + std::move(that);
}

FunctionExpression Function::operator-(const FunctionExpression &that) const & {

  return FunctionExpression(*this) - that;
}
FunctionExpression Function::operator-(FunctionExpression &&that) const & {
  return FunctionExpression(*this) - std::move(that);
}
FunctionExpression Function::operator-(const FunctionExpression &that) && {
  return FunctionExpression(std::move(*this)) - that;
}
FunctionExpression Function::operator-(FunctionExpression &&that) && {
  return FunctionExpression(std::move(*this)) - std::move(that);
}
FunctionExpression Function::operator-() const & {
  return -FunctionExpression(*this);
}
FunctionExpression Function::operator-() && {

  return -FunctionExpression(std::move(*this));
}

FunctionExpression Function::operator*(const FunctionExpression &that) const & {
  return FunctionExpression(*this) * that;
}
FunctionExpression Function::operator*(FunctionExpression &&that) const & {
  return FunctionExpression(*this) + std::move(that);
}
FunctionExpression Function::operator*(const FunctionExpression &that) && {
  return FunctionExpression(std::move(*this)) * that;
}
FunctionExpression Function::operator*(FunctionExpression &&that) && {
  return FunctionExpression(std::move(*this)) + std::move(that);
}

FunctionExpression Function::compose(const FunctionExpression &that) const & {
  return FunctionExpression(*this).compose(that);
}
FunctionExpression Function::compose(FunctionExpression &&that) const & {
  return FunctionExpression(*this).compose(std::move(that));
}

FunctionExpression Function::compose(const FunctionExpression &that) && {
  return FunctionExpression(std::move(*this)).compose(that);
}
FunctionExpression Function::compose(FunctionExpression &&that) && {
  return FunctionExpression(std::move(*this)).compose(std::move(that));
}

FunctionExpression Function::concat(const FunctionExpression &that) const & {
  return FunctionExpression(*this).concat(that);
}
FunctionExpression Function::concat(FunctionExpression &&that) const & {
  return FunctionExpression(*this).concat(std::move(that));
}

FunctionExpression Function::concat(const FunctionExpression &that) && {
  return FunctionExpression(std::move(*this)).concat(that);
}
FunctionExpression Function::concat(FunctionExpression &&that) && {
  return FunctionExpression(std::move(*this)).concat(std::move(that));
}
} // namespace functions
} // namespace gsplines
