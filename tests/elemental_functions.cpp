
#include <cmath>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

int main() {
  const std::pair<double, double> base_domain{-1.0, 1.0};
  Exponential exp(base_domain);
  Sin sin(base_domain);
  Cos cos(base_domain);

  FunctionExpression f_nom = (exp + sin + cos) * cos;
  FunctionExpression f_d_nom =
      (exp + sin + cos) * (-sin) + (exp + cos - sin) * cos;

  for (int _ = 0; _ < 100; _++) {
    Eigen::VectorXd time_span = Eigen::VectorXd::Random(10000);

    assert(
        ((exp(time_span) + sin(time_span) + cos(time_span)) * cos(time_span) -
         f_nom(time_span))
            .norm() < 1.0e-9);

    FunctionExpression f_dot = f_nom.derivate();

    assert((f_dot(time_span) - exp(time_span) + cos(time_span) - sin(time_span))
               .norm() < 1.0e-9);
  }

  return 0;
}
