#include <eigen3/Eigen/Core>
#include <gtest/gtest.h>
#include <kdl/velocityprofile_trap.hpp>
#include <ruckig/ruckig.hpp>
TEST(KDL, ControlTest) {
  std::vector<KDL::VelocityProfile_Trap> m_velocityProfiles;
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
