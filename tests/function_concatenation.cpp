#include <cmath>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

int main() {
  Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);
  Sin sin({-M_PI_2, M_PI_2});
  Cos cos({-M_PI_2, M_PI_2});
  Identity identity({-1.0, 1.0});

  DomainLinearDilation linear_dilation({-1, 1}, 2);
  FunctionExpression f = cos.compose(linear_dilation.compose(linear_dilation));

  FunctionExpression f_dot_nom =
      -4 * sin.compose(DomainLinearDilation({-1, 1}, 4));

  assert((f(time_span) - Eigen::cos(4 * time_span.array()).matrix()).norm() <
         1.0e-9);

  assert((f.derivate().value(time_span) - f_dot_nom(time_span)).norm() <
         1.0e-9);

  /*
  FunctionExpression g =
      (sin + identity + sin.compose(sin)).compose(sin).compose(cos);

  FunctionExpression g_dot =
      (cos + ConstFunction({-1.0, 1.0}, 1, 0.0) + cos.compose(sin) * cos)
          .compose(sin)
          .compose(cos) *
      cos.compose(cos) * (-sin);

  printf("------------\n");

  assert((g.derivate().value(time_span) - g_dot(time_span)).norm() < 1.0e-9);*/
  return 0;
}
