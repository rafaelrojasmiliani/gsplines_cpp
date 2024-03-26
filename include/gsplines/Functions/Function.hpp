#ifndef FUNCTION
#define FUNCTION

#include <eigen3/Eigen/Core>
#include <gsplines/Functions/FunctionBase.hpp>

namespace gsplines {
namespace functions {

class FunctionExpression;
class Function : public FunctionBase {
 public:
  using FunctionBase::FunctionBase;
  Function(const Function& _that) = default;
  Function(Function&& _that) = default;
  ~Function() override = default;
};

}  // namespace functions
}  // namespace gsplines
#endif
