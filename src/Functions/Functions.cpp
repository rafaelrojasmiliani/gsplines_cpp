
#include <cmath>
#include <gsplines++/Functions/Function.hpp>
#include <stdexcept>

namespace gsplines {
namespace functions {

Function::Function(std::pair<double, double> _domain, std::size_t _codom_dim,
                   const std::string &_name)
    : FunctionExpression(_domain, _codom_dim, FunctionExpression::SINGLE,
                         std::list<std::unique_ptr<FunctionExpression>>(),
                         _name) {}

Function::Function(const Function &that) : FunctionExpression(that) {}

} // namespace functions
} // namespace gsplines
