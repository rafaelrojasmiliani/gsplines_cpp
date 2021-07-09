
#include <cmath>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

int main() {
  const std::pair<double, double> base_domain{-1.0, 1.0};
  Exponential exp(base_domain);
  Sin sin(base_domain);
  Cos cos(base_domain);

  FunctionExpression f_nom = (exp + sin + cos);

  FunctionExpression f_d_nom = (exp + cos - sin);

  for (int _ = 0; _ < 100; _++) {
    Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

    assert((exp(time_span) + sin(time_span) + cos(time_span) - f_nom(time_span))
               .norm() < 1.0e-9);

    FunctionExpression f_dot = f_nom.derivate();

    f_dot.print();

    printf("function name %s\n", f_dot.get_name().c_str());

    std::cout << f_dot(time_span) << "\n";
    std::cout << exp(time_span) + cos(time_span) - sin(time_span) << "\n"
              << (f_dot(time_span) - exp(time_span) + cos(time_span) -
                  sin(time_span))
                     .norm()
              << "\n";

    assert((f_dot(time_span) - exp(time_span) - cos(time_span) + sin(time_span))
               .norm() < 1.0e-9);
  }

  return 0;
}
