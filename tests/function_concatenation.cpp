#include <cmath>
#include <gsplines/Functions/ElementalFunctions.hpp>
#include <iostream>
using namespace gsplines::functions;

std::size_t number_of_wp = 3;
std::size_t codom_dim = 3;

double exec_time = (double)number_of_wp - 1.0;
double T = exec_time * 0.9;
double ti = 0.5 * exec_time;
double Ts = T - ti;
double sf = ti * (Ts / ti + 1);

void compare_assert(Eigen::MatrixXd &_m_nom, Eigen::MatrixXd &_m_test) {

  // std::cout << _m_nom << "\n----\n";
  // std::cout << _m_test << "\n----\n";
  // std::cout << (_m_nom - _m_test).rowwise().norm() << "\n----\n";
  if (_m_nom.array().abs().maxCoeff() < 1.0e-9) {
    assert((_m_nom - _m_test).rowwise().norm().maxCoeff() < 1.0e-9);
  } else {
    double err = (_m_nom - _m_test).rowwise().norm().maxCoeff() /
                 _m_nom.rowwise().norm().maxCoeff();

    assert(err < 1.0e-9);
  }
}

int main() {

  Eigen::VectorXd pol_coeff(6);
  pol_coeff << ti, Ts, 0.0, -6.0 * Ts + 10.0 * sf - 10.0 * ti,
      8.0 * Ts - 15.0 * sf + 15.0 * ti, -3.0 * Ts + 6.0 * sf - 6.0 * ti;

  FunctionExpression tau = (1.0 / Ts) * (Identity({ti, ti + Ts}) +
                                         ConstFunction({ti, ti + Ts}, 1, -ti));

  CanonicPolynomial pol({0, 1}, pol_coeff);

  FunctionExpression pol_t = pol.compose(tau);
  FunctionExpression diffeo =
      gsplines::functions::Identity({0, ti}).concat(pol_t);
  Identity id = Identity({0, ti});

  Eigen::VectorXd interval_1 = Eigen::VectorXd::LinSpaced(10, 0, ti - 1.0e-5);
  Eigen::VectorXd interval_2 =
      Eigen::VectorXd::LinSpaced(10, ti, ti + Ts - 1.0e-5);

  Eigen::MatrixXd tau_val_test = tau(interval_2);
  Eigen::MatrixXd tau_val_nom = ((interval_2.array() - ti).matrix() / Ts);

  compare_assert(tau_val_nom, tau_val_test);

  for (std::size_t i = 0; i < 5; i++) {
    FunctionExpression piecewise = diffeo.derivate(i);
    FunctionExpression piece_1 = id.derivate(i);
    FunctionExpression piece_2 = pol_t.derivate(i);

    Eigen::MatrixXd interval_1_p1 = piece_1(interval_1);
    printf("--------------------- calling piece 2 ------------------\n");
    Eigen::MatrixXd interval_2_p2 = piece_2(interval_2);

    Eigen::MatrixXd interval_1_pw = piecewise(interval_1);
    printf("--------------------- calling pw in 2 ------------------\n");
    Eigen::MatrixXd interval_2_pw = piecewise(interval_2);

    compare_assert(interval_1_pw, interval_1_p1);
    compare_assert(interval_2_pw, interval_2_p2);
  }
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
