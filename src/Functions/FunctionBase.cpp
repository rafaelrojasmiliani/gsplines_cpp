#include <cmath>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionBase.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <stdexcept>

namespace gsplines {
namespace functions {
/**  Function Base */
const double FunctionBase::dom_tollerance_ = 1.0e-9;

FunctionBase::FunctionBase(std::pair<double, double> _domain,
                           std::size_t _codom_dim, const std::string &_name)
    : codom_dim_(_codom_dim), window_(_domain), domain_(_domain), name_(_name) {

  assert(_domain.first <= _domain.second);
}

FunctionBase::FunctionBase(const FunctionBase &that)
    : codom_dim_(that.codom_dim_), window_(that.window_), domain_(that.domain_),
      name_(that.name_) {

  assert(domain_.first <= domain_.second);
}

bool FunctionBase::same_domain(const FunctionBase &_f1,
                               const FunctionBase &_f2) {
  bool err1 = abs(_f1.domain_.first - _f2.domain_.first) <
              FunctionBase::dom_tollerance_;
  bool err2 = abs(_f1.domain_.second - _f2.domain_.second) <
              FunctionBase::dom_tollerance_;

  return err1 and err2;
}
bool FunctionBase::same_codomain(const FunctionBase &_f1,
                                 const FunctionBase &_f2) {
  return _f1.get_codom_dim() == _f2.get_codom_dim();
}

bool FunctionBase::is_point_in_domain(double _domain_point) const {
  return domain_.first <= _domain_point and _domain_point < domain_.second;
}

void FunctionBase::print(std::size_t _indent) const {
  printf("%*s- %s domain = [ %+11.3lf, %+11.3lf] codomain dim = %zu\n",
         4 * (int)_indent, "", name_.c_str(), domain_.first, domain_.second,
         codom_dim_);
}
DotProduct FunctionBase::dot(const FunctionBase &_that) const & {
  return DotProduct(*this, _that);
}
DotProduct FunctionBase::dot(FunctionBase &&_that) const & {

  return DotProduct(*this, std::move(_that));
}
DotProduct FunctionBase::dot(const FunctionBase &_that) && {

  return DotProduct(_that, std::move(*this));
}
DotProduct FunctionBase::dot(FunctionBase &&_that) && {

  return DotProduct(std::move(_that), std::move(*this));
}

FunctionExpression FunctionBase::to_expression() const & {
  std::list<std::unique_ptr<FunctionBase>> result_array;
  result_array.push_back(clone());

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::UNIQUE,
                            std::move(result_array));
}
FunctionExpression FunctionBase::to_expression() && {
  std::list<std::unique_ptr<FunctionBase>> result_array;
  result_array.push_back(move_clone());

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::UNIQUE,
                            std::move(result_array));
}

FunctionExpression FunctionBase::derivate(std::size_t _deg) const {
  std::list<std::unique_ptr<FunctionBase>> result_array;
  result_array.push_back(std::unique_ptr<FunctionBase>(deriv_impl(_deg)));

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::UNIQUE,
                            std::move(result_array));
}
} // namespace functions
} // namespace gsplines
