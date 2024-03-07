

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

/* Test that the Legendre polynomials respect the recursive relation */
using namespace gsplines;

/* Test that the derivtaive of the Legendre polynomials respect the recursive
 * relation */
TEST(Basis0101, Derivative) {

  gsplines::basis::Basis0101 basis(0.5);
  Eigen::VectorXd buff(6);
  Eigen::VectorXd buff_d1(6);
  Eigen::VectorXd buff_d2(6);
  double t = -1;
  basis.eval_derivative_on_window(t, 2, 1, buff);
  basis.eval_derivative_on_window(t, 2, 2, buff_d1);
  basis.eval_derivative_on_window(t, 2, 3, buff_d2);
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
