#include <cmath>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <iostream>
using namespace gsplines::functions;

TEST(Function, Addition) {

  Exponential s({-1.0, 1.0}), g({-1.0, 1.0});
  Eigen::VectorXd time_spam = Eigen::VectorXd::Random(4);
  Eigen::MatrixXd exp_value = s(time_spam);
  gsplines::functions::FunctionExpression f = s + g + s + g;
  gsplines::functions::FunctionExpression m = f + g + s;

  EXPECT_TRUE(
      gsplines::tools::approx_equal(4 * exp_value, f(time_spam), 1.0e-9));
  EXPECT_TRUE(
      gsplines::tools::approx_equal(m(time_spam), 6 * exp_value, 1.0e-9));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
