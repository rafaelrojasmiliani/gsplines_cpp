

#include <Eigen/Core>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <cstddef>
using namespace gsplines;
TEST(Interpolator, Value) {
  for (std::size_t i = 1; i < 3; i++) {
    const basis::BasisLegendre basis(4);
    const std::size_t dim = 2;
    const std::size_t intervals = 4;
    const Eigen::VectorXd tau = Eigen::VectorXd::Random(4).array() + 1.2;
    const Eigen::MatrixXd wp = Eigen::MatrixXd::Random(intervals + 1, dim);
    Interpolator inter(dim, intervals, basis);
    inter.fill_interpolating_vector(wp);
    const GSpline res = inter.interpolate(tau, wp);

    const Eigen::MatrixXd bp = res.get_waypoints();
    EXPECT_TRUE(tools::approx_equal(bp, wp, 1.0e-9))
        << "actual break points \n"
        << bp << "\n ---\n desired break points \n"
        << wp;
    /*
    for (std::size_t m = 0; m < bp.size(); m++) {
      bp(m);
    }*/
  }
  EXPECT_TRUE(true);
}
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
