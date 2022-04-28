

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

/* Test that the Legendre polynomials respect the recursive relation */
using namespace gsplines;
TEST(BasisLegendre, Value) {

  for (int i = 3; i < 17; i++) {

    basis::BasisLegendre basis(i);

    EXPECT_TRUE(basis::BasisLegendre(i) == basis::BasisLegendre(i));
    EXPECT_TRUE(basis::BasisLegendre(i) == *basis::BasisLegendre::get(i));
    EXPECT_FALSE(basis::BasisLegendre(i) != basis::BasisLegendre(i));
    EXPECT_FALSE(basis::BasisLegendre(i) != *basis::BasisLegendre::get(i));

    for (int j = 3; j < 17; j++) {
      if (i != j) {
        EXPECT_FALSE(basis::BasisLegendre(i) == basis::BasisLegendre(j));
        EXPECT_TRUE(basis::BasisLegendre(i) != basis::BasisLegendre(j));
      }
    }
    Eigen::VectorXd buff(i);
    for (const double &t :
         {-1.0, -0.9, -0.5, -0.33, 0.0, 0.33, 0.5, 0.9, 1.0}) {
      basis.eval_on_window(t, 2, buff);
      for (int n = 2; n < i - 1; n++) {
        double lhs = (n + 1.0) * buff(n + 1);
        double rhs = (2.0 * n + 1.0) * t * buff(n) - n * buff(n - 1);
        EXPECT_TRUE(gsplines::tools::approx_equal(lhs, rhs, 1.0e-9))
            << buff << "\n --- " << lhs << "   " << rhs;
      }
    }
  }
}

/* Test that the derivtaive of the Legendre polynomials respect the recursive
 * relation */
TEST(BasisLegendre, Derivative) {

  for (int i = 3; i < 17; i++) {

    gsplines::basis::BasisLegendre basis(i);
    Eigen::VectorXd buff(i);
    Eigen::VectorXd buff_d1(i);
    Eigen::VectorXd buff_d2(i);
    for (const double &t :
         {-1.0, -0.9, -0.5, -0.33, 0.0, 0.33, 0.5, 0.9, 1.0}) {
      basis.eval_derivative_on_window(t, 2, 0, buff);
      basis.eval_derivative_on_window(t, 2, 1, buff_d1);
      basis.eval_derivative_on_window(t, 2, 2, buff_d2);

      for (int n = 2; n < i - 1; n++) {
        double lhs = (1.0 - t * t) * buff_d2(n) - 2 * t * buff_d1(n) +
                     n * (n + 1) * buff(n);
        EXPECT_TRUE(gsplines::tools::approx_equal(lhs, 0.0, 1.0e-9))
            << buff << "\n --- " << lhs;
      }
    }
  }
}
int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
