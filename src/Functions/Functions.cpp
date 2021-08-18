
#include <cmath>
#include <gsplines/Functions/Function.hpp>
#include <stdexcept>

namespace gsplines {
namespace functions {

Function::Function(std::pair<double, double> _domain, std::size_t _codom_dim,
                   const std::string &_name)
    : FunctionExpression(_domain, _codom_dim, FunctionExpression::SINGLE,
                         std::list<std::unique_ptr<FunctionExpression>>(),
                         _name) {}

Function::Function(const Function &that) : FunctionExpression(that) {}

FunctionExpression Function::derivate(int _deg) const {

  std::list<std::unique_ptr<FunctionExpression>> result_array;

  result_array.push_back(this->deriv(_deg));
  const std::string name = get_name();

  return FunctionExpression(get_domain(), get_codom_dim(),
                            FunctionExpression::Type::UNIQUE,
                            std::move(result_array), name);
}
} // namespace functions
} // namespace gsplines
