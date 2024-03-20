#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/Basis0101.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <cstddef>
#include <memory>
#include <random>
#include <vector>
using namespace gsplines;
/** Test that the linear scaling works.
 **/
TEST(LinearScaling, BasisLegendre) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  std::vector<std::unique_ptr<basis::Basis>> basis_vec;

  for (int i = 0; i < 2; i++) {
    const std::size_t dim = uint_dist(mt);

    auto g1 = random_gspline({0.0, 10.0}, dim, basis::BasisLegendre(6));

    auto g2 = g1.linear_scaling_new_execution_time(20.0);
    EXPECT_TRUE(tools::approx_equal(g1.get_coefficients(),
                                    g2.get_coefficients(), 1.0e-4));
    EXPECT_FLOAT_EQ(g2.get_domain_length(), 20.0);

    auto g3 =
        g1.linear_scaling_new_execution_time_max_velocity_max_acceleration(0.4,
                                                                           0.6);

    EXPECT_TRUE(tools::approx_equal(g3.get_coefficients(),
                                    g2.get_coefficients(), 1.0e-4));
    const Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
        100, g3.get_domain().first, g3.get_domain().second);

    auto g3d = g3.derivate();
    auto g3dd = g3d.derivate();
    const Eigen::MatrixXd gspline_diff_1_evaluated = g3d(time_spam);
    const Eigen::MatrixXd gspline_diff_2_evaluated = g3dd(time_spam);

    EXPECT_LE(gspline_diff_1_evaluated.array().abs().maxCoeff(), 0.401);
    EXPECT_LE(gspline_diff_2_evaluated.array().abs().maxCoeff(), 0.601);
  }
}
TEST(LinearScaling, Basis0101) {
  std::random_device rd;
  std::mt19937 mt(rd());
  std::uniform_int_distribution<std::size_t> uint_dist(2, 10);

  for (int i = 0; i < 2; i++) {
    const std::size_t dim = uint_dist(mt);

    auto g1 = random_gspline({0.0, 10.0}, dim, basis::Basis0101(0.5));

    auto g2 = g1.linear_scaling_new_execution_time(20.0);
    EXPECT_TRUE(tools::approx_equal(g1.get_coefficients(),
                                    g2.get_coefficients(), 1.0e-4));
    EXPECT_FLOAT_EQ(g2.get_domain_length(), 20.0);

    auto g3 =
        g1.linear_scaling_new_execution_time_max_velocity_max_acceleration(0.4,
                                                                           0.6);

    EXPECT_TRUE(tools::approx_equal(g3.get_coefficients(),
                                    g2.get_coefficients(), 1.0e-4));
    const Eigen::VectorXd time_spam = Eigen::VectorXd::LinSpaced(
        100, g3.get_domain().first, g3.get_domain().second);

    auto g3d = g3.derivate();
    auto g3dd = g3d.derivate();
    const Eigen::MatrixXd gspline_diff_1_evaluated = g3d(time_spam);
    const Eigen::MatrixXd gspline_diff_2_evaluated = g3dd(time_spam);

    EXPECT_LE(gspline_diff_1_evaluated.array().abs().maxCoeff(), 0.401);
    EXPECT_LE(gspline_diff_2_evaluated.array().abs().maxCoeff(), 0.601);
  }
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
