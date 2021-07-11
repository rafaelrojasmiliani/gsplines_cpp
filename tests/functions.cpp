#include <cmath>
#include <gsplines++/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

int main() {
  Exponential s({-1.0, 1.0}), g({-1.0, 1.0});
  Eigen::VectorXd time_span = Eigen::VectorXd::Random(4);
  Eigen::MatrixXd exp_value = s(time_span);

  gsplines::functions::FunctionExpression f = s + g + s + g;

  printf(".................f = s + g + s +g\n");
  f.print_performace();
  printf(".................\n");
  gsplines::functions::FunctionExpression z = f;
  printf(".................z=f\n");
  f.print_performace();
  printf(".................\n");
  gsplines::functions::FunctionExpression m = f + g + s;
  printf(".................m = f + g + s;\n");
  f.print_performace();
  printf(".................\n");

  printf("|||||||||||||||||||||||||");
  assert((4 * exp_value - f(time_span)).norm() < 1.0e-9);
  printf("|||||||||||||||||||||||||");
  assert((m(time_span) - 6 * exp_value).norm() < 1.0e-9);
  printf("|||||||||||||||||||||||||");

  f.print_performace();

  return 0;
}
