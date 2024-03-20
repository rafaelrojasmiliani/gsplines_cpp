#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/GSpline.hpp>
#include <gtest/gtest.h>
#include <cstddef>
#include <memory>
#include <random>
#include <vector>
using namespace gsplines;
/** Test that the linear scaling works.
 **/
TEST(LinearScaling, Scale) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  std::vector<std::unique_ptr<basis::Basis>> basis_vec;
  basis_vec.push_back(basis::Basis0101(0.5).move_clone());
  basis_vec.push_back(basis::BasisLegendre(6).move_clone());

  for (const auto& basis : basis_vec) {
    for (int i = 0; i < 100; i++) {
      const std::size_t dim = uint_dist(mt);

      auto g1 = random_gspline({0.0, 10.0}, dim, *basis);

      auto g2 = g1.linear_scaling_new_execution_time(20.0);
      EXPECT_FLOAT_EQ(g2.get_domain_length(), 20.0);

      auto g3 =
          g1.linear_scaling_new_execution_time_max_velocity_max_acceleration(
              0.4, 0.6);

      const Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
          100, g3.get_domain().first, g3.get_domain().second);

      const Eigen::MatrixXd gspline_diff_1_evaluated = g3.derivate()(time_spam);
      const Eigen::MatrixXd gspline_diff_2_evaluated =
          g3.derivate(2)(time_spam);

      EXPECT_LE(gspline_diff_1_evaluated.array().abs().maxCoeff(), 0.4);
      EXPECT_LE(gspline_diff_2_evaluated.array().abs().maxCoeff(), 0.6);
    }
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
