#include <eigen3/Eigen/Core>
#include <gtest/gtest.h>
#include <ruckig/ruckig.hpp>

TEST(RUCKIG, ControlTest) {
#ifdef HAS_RUCKIG
  ruckig::Ruckig<3> otg(0.01);  // control cycle
  ruckig::InputParameter<3> input;
  ruckig::OutputParameter<3> output;

  // Set input parameters
  input.current_position = {0.0, 0.0, 0.5};
  input.current_velocity = {0.0, -2.2, -0.5};
  input.current_acceleration = {0.0, 2.5, -0.5};

  input.target_position = {5.0, -2.0, -3.5};
  input.target_velocity = {0.0, -0.5, -2.0};
  input.target_acceleration = {0.0, 0.0, 0.5};

  input.max_velocity = {3.0, 1.0, 3.0};
  input.max_acceleration = {3.0, 2.0, 1.0};
  input.max_jerk = {4.0, 3.0, 2.0};

  // Generate the trajectory within the control loop
  std::cout << "t | position" << std::endl;
  while (otg.update(input, output) == ruckig::Result::Working) {
    std::cout << output.time << " | " << ruckig::join(output.new_position)
              << std::endl;

    output.pass_to_input(input);
  }

  std::cout << "Trajectory duration: " << output.trajectory.get_duration()
            << " [s]." << std::endl;
#endif
}

TEST(RUCKIG, TrajectoryTest) {
#ifdef HAS_RUCKIG
  // Create input parameters
  ruckig::InputParameter<3> input;
  input.current_position = {0.0, 0.0, 0.5};
  input.current_velocity = {0.0, -2.2, -0.5};
  input.current_acceleration = {0.0, 2.5, -0.5};

  input.target_position = {5.0, -2.0, -3.5};
  input.target_velocity = {0.0, -0.5, -2.0};
  input.target_acceleration = {0.0, 0.0, 0.5};

  input.max_velocity = {3.0, 1.0, 3.0};
  input.max_acceleration = {3.0, 2.0, 1.0};
  input.max_jerk = {4.0, 3.0, 2.0};

  // Set different constraints for negative direction
  input.min_velocity = {-2.0, -0.5, -3.0};
  input.min_acceleration = {-2.0, -2.0, -2.0};

  // We don't need to pass the control rate (cycle time) when using only offline
  // features
  ruckig::Ruckig<3> otg;
  ruckig::Trajectory<3> trajectory;

  // Calculate the trajectory in an offline manner (outside of the control loop)
  ruckig::Result result = otg.calculate(input, trajectory);
  if (result == ruckig::Result::ErrorInvalidInput) {
    std::cout << "Invalid input!" << std::endl;
    return;
  }

  // Get duration of the trajectory
  std::cout << "Trajectory duration: " << trajectory.get_duration() << " [s]."
            << std::endl;

  double new_time = 1.0;

  // Then, we can calculate the kinematic state at a given time
  std::array<double, 3> new_position{};
  std::array<double, 3> new_velocity{};
  std::array<double, 3> new_acceleration{};
  trajectory.at_time(new_time, new_position, new_velocity, new_acceleration);

  std::cout << "Position at time " << new_time
            << " [s]: " << ruckig::join(new_position) << std::endl;

  // Get some info about the position extrema of the trajectory
  std::array<ruckig::Bound, 3> position_extrema =
      trajectory.get_position_extrema();
  std::cout << "Position extremas for DoF 4 are " << position_extrema[2].min
            << " (min) to " << position_extrema[2].max << " (max)" << std::endl;
#endif
}

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
