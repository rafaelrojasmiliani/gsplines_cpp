
#include <cmath>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

int main() {
  Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);
  Sin sin({-M_PI_2, M_PI_2});

  DomainLinearDilation linear_dilation({-1.0, 1.0}, 2.0);

  FunctionExpression f = sin.compose(linear_dilation);

  assert((f(time_span) - Eigen::sin(2 * time_span.array()).matrix()).norm() <
         1.0e-9);

  return 0;
}
