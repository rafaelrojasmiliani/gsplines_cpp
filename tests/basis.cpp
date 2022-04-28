

#include "gsplines/Basis/Basis.hpp"
#include "gsplines/Basis/BasisLagrange.hpp"
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>

/* Test that we get the correct basis*/
using namespace gsplines;
TEST(Basis, Get_Basis) {

  for (int i = 3; i < 17; i++) {

    Eigen::VectorXd points = Eigen::VectorXd::Random(i);
    basis::BasisLegendre basis_legendre(i);
    basis::BasisLagrange basis_lagrange(points);

    EXPECT_TRUE(basis_legendre == *basis::get_basis("legendre", i, {}));
    for (int j = 3; j < 17; j++) {
      if (j != i)
        EXPECT_TRUE(basis_legendre != *basis::get_basis("legendre", j, {}));
    }

    EXPECT_TRUE(basis_lagrange == *basis::get_basis("lagrange", i, points));
    for (int j = 3; j < 17; j++) {
      EXPECT_TRUE(basis_lagrange !=
                  *basis::get_basis("lagrange", i, Eigen::VectorXd::Random(j)));
    }
  }
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
