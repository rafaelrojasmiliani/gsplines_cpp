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

class Function : public FunctionBase {

public:
  using FunctionBase::FunctionBase;
  Function(const Function &_that) : FunctionBase(_that) {}
  Function(Function &&_that) : FunctionBase(std::move(_that)) {}
  ~Function() = default;
};

} // namespace functions
} // namespace gsplines
#endif
