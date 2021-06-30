
#include <cmath>
#include <gsplines++/Functions/Functions.hpp>
#include <stdexcept>

namespace gsplines {
namespace functions {
/**  Function Base */
const double FunctionBase::dom_tollerance_ = 1.0e-5;
FunctionBase::FunctionBase(std::pair<double, double> _domain,
                           std::size_t _codom_dim)
    : codom_dim_(_codom_dim), window_(_domain), domain_(_domain) {}

FunctionBase::FunctionBase(const FunctionBase &that)
    : codom_dim_(that.codom_dim_), window_(that.window_),
      domain_(that.domain_) {}

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

bool FunctionBase::is_point_in_domain(double _domain_point) {
  return domain_.first <= _domain_point and _domain_point <= domain_.second;
}
/**  Function */
Function::Function(std::pair<double, double> _domain, std::size_t _codom_dim)
    : FunctionBase(_domain, _codom_dim) {}

Function::Function(const Function &that) : FunctionBase(that) {}

} // namespace functions
} // namespace gsplines
