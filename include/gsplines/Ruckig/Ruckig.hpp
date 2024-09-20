#ifndef PIECEWISE_FUNCTION_H
#define PIECEWISE_FUNCTION_H

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace gsplines {

class RuckigCurve : public functions::FunctionInheritanceHelper<
                        RuckigCurve, functions::Function, RuckigCurve> {
 public:
 private:
};

}  // namespace gsplines

#endif /* PIECEWISE_FUNCTION_H */
