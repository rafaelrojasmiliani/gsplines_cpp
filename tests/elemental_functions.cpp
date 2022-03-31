
#include <cmath>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <iostream>
using namespace gsplines::functions;

TEST(Functions, Value) {
  const std::pair<double, double> base_domain{-1.0, 1.0};
  Exponential exp(base_domain);
  Sin sin(base_domain);
  Cos cos(base_domain);

  {
    FunctionExpression f_nom = (exp + sin + cos);

    FunctionExpression f_d_nom = (exp + cos - sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      EXPECT_TRUE(
          (exp(time_span) + sin(time_span) + cos(time_span) - f_nom(time_span))
              .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      f_dot.print();

      EXPECT_TRUE(
          (f_dot(time_span) - exp(time_span) - cos(time_span) + sin(time_span))
              .norm() < 1.0e-9);
    }
  }

  {
    FunctionExpression f_nom = (exp * sin);

    FunctionExpression f_d_nom = (exp * cos) + (exp * sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      f_nom.print();
      EXPECT_TRUE(((exp(time_span).array() * sin(time_span).array()).matrix() -
                   f_nom(time_span))
                      .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      EXPECT_TRUE((f_dot(time_span) - f_d_nom(time_span)).norm() < 1.0e-9);
    }
  }

  {
    FunctionExpression f_nom = (exp + sin + cos);

    FunctionExpression f_d_nom = (exp + cos - sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      EXPECT_TRUE(((exp(time_span) + sin(time_span) + cos(time_span)).matrix() -
                   f_nom(time_span))
                      .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      EXPECT_TRUE((f_dot(time_span) - f_d_nom(time_span)).norm() < 1.0e-9);
    }
  }

  {
    FunctionExpression f_nom = (exp + sin + cos) * cos;

    FunctionExpression f_d_nom =
        (exp + cos - sin) * cos + (exp + sin + cos) * (-sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      EXPECT_TRUE((((exp(time_span) + sin(time_span) + cos(time_span)).array() *
                    cos(time_span).array())
                       .matrix() -
                   f_nom(time_span))
                      .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      EXPECT_TRUE((f_dot(time_span) - f_d_nom(time_span)).norm() < 1.0e-9);
    }
  }
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
