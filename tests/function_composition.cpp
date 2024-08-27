#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <cmath>
using namespace gsplines::functions;

TEST(Function, Composition) {
  Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);
  Sin sin({-1.0, 1.0});
  Cos cos({-1.0, 1.0});
  Identity identity({-1.0, 1.0});

  DomainLinearDilation linear_dilation({-1, 1}, 2);
  FunctionExpression f = cos.compose(linear_dilation);

  FunctionExpression f_dot_nom =
      -2 * sin.compose(DomainLinearDilation({-1, 1}, 2));

  f.print();
  f.derivate().print();
  EXPECT_TRUE(
      (f(time_span) - Eigen::cos(2 * time_span.array()).matrix()).norm() <
      1.0e-9);
  (void)identity.value(time_span);

  EXPECT_TRUE((f.derivate().value(time_span) -
               -2 * Eigen::sin(2 * time_span.array()).matrix())
                  .norm() < 1.0e-9);

  FunctionExpression g =
      (sin + identity + sin.compose(sin)).compose(sin).compose(cos);

  FunctionExpression g_dot =
      (cos + ConstFunction({-1.0, 1.0}, 1, 1.0) + cos.compose(sin) * cos)
          .compose(sin)
          .compose(cos) *
      cos.compose(cos) * (-sin);

  EXPECT_TRUE((g.derivate().value(time_span) - g_dot(time_span)).norm() <
              1.0e-9);
}
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
