
#include <cmath>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

int main() {
  const std::pair<double, double> base_domain{-1.0, 1.0};
  Exponential exp(base_domain);
  Sin sin(base_domain);
  Cos cos(base_domain);

  {
    FunctionExpression f_nom = (exp + sin + cos);

    FunctionExpression f_d_nom = (exp + cos - sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      assert(
          (exp(time_span) + sin(time_span) + cos(time_span) - f_nom(time_span))
              .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      f_dot.print();

      assert(
          (f_dot(time_span) - exp(time_span) - cos(time_span) + sin(time_span))
              .norm() < 1.0e-9);
    }
  }

  {
    FunctionExpression f_nom = (exp * sin);

    printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
    FunctionExpression f_d_nom = (exp * cos) + (exp * sin);
    printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      assert(((exp(time_span).array() * sin(time_span).array()).matrix() -
              f_nom(time_span))
                 .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      f_d_nom.print();
      f_dot.print();
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      printf("Evalueation of f_dot\n\n");
      std::cout << f_dot(time_span) << "<< \n";
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      printf("Evaluation of f_d_nom \n\n");
      std::cout << f_d_nom(time_span) << "<< \n";
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      std::cout << (f_dot(time_span) - f_d_nom(time_span)).norm() << "<< \n";
      assert((f_dot(time_span) - f_d_nom(time_span)).norm() < 1.0e-9);
      printf("------------ assertion paased\n\n");
      fflush(stdout);
    }
  }

  {
    FunctionExpression f_nom = (exp + sin + cos);

    FunctionExpression f_d_nom = (exp + cos - sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      assert(((exp(time_span) + sin(time_span) + cos(time_span)).matrix() -
              f_nom(time_span))
                 .norm() < 1.0e-9);

      FunctionExpression f_dot = f_nom.derivate();

      f_dot.print();

      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      f_d_nom.print();
      f_dot.print();
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      printf("Evalueation of f_dot\n\n");
      std::cout << f_dot(time_span) << "<< \n";
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      printf("Evaluation of f_d_nom \n\n");
      std::cout << f_d_nom(time_span) << "<< \n";
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      std::cout << (f_dot(time_span) - f_d_nom(time_span)).norm() << "<< \n";
      fflush(stdout);
      assert((f_dot(time_span) - f_d_nom(time_span)).norm() < 1.0e-9);
    }
  }

  {
    FunctionExpression f_nom = (exp + sin + cos) * cos;

    FunctionExpression f_d_nom =
        (exp + cos - sin) * cos + (exp + sin + cos) * (-sin);

    for (int _ = 0; _ < 100; _++) {
      Eigen::VectorXd time_span = Eigen::VectorXd::Random(5);

      assert((((exp(time_span) + sin(time_span) + cos(time_span)).array() *
               cos(time_span).array())
                  .matrix() -
              f_nom(time_span))
                 .norm() < 1.0e-9);
      printf("--------------------iiiiiiiiiiiiiiiiiiiii\n");

      FunctionExpression f_dot = f_nom.derivate();

      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      f_d_nom.print();
      f_dot.print();
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      printf("Evalueation of f_dot\n\n");
      std::cout << f_dot(time_span) << "<< \n";
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      printf("Evaluation of f_d_nom \n\n");
      std::cout << f_d_nom(time_span) << "<< \n";
      printf("+++++++++++++++++ ########## ------------ ++++++++++++ \n\n");
      std::cout << (f_dot(time_span) - f_d_nom(time_span)).norm() << "<< \n";
      fflush(stdout);
      assert((f_dot(time_span) - f_d_nom(time_span)).norm() < 1.0e-9);
    }
  }

  return 0;
}
