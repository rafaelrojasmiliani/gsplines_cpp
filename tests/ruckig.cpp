#include <eigen3/Eigen/Core>
#include <gtest/gtest.h>
#include <ruckig/ruckig.hpp>
#include <optional>

std::optional<ruckig::InputParameter<ruckig::DynamicDOFs>> paramsFrom(
    Eigen::Ref<const Eigen::MatrixXd> _waypoints,
    const std::vector<double>& _max_abs_vel,
    const std::vector<double>& _max_abs_acc,
    const std::vector<double>& _max_abs_jerk) {
  //
  std::size_t dof = static_cast<std::size_t>(_waypoints.cols());
  ruckig::InputParameter<ruckig::DynamicDOFs> result(dof);

  std::size_t number_of_waypoints = static_cast<std::size_t>(_waypoints.rows());

  if (number_of_waypoints < 2) {
    return std::nullopt;
  }
  if (_waypoints.cols() != dof) {
    return std::nullopt;
  }
  if (_max_abs_vel.size() != dof) {
    return std::nullopt;
  }
  if (_max_abs_acc.size() != dof) {
    return std::nullopt;
  }
  if (_max_abs_jerk.size() != dof) {
    return std::nullopt;
  }

  result.max_velocity = _max_abs_vel;
  result.min_acceleration = _max_abs_acc;
  result.max_jerk = _max_abs_jerk;
  result.current_velocity.assign(dof, 0.0);
  result.current_acceleration.assign(dof, 0.0);

  for (int i = 0; i < dof; ++i) {
    result.current_position[i] = _waypoints.topRows(1)(0, i);
    result.target_position[i] = _waypoints.bottomRows(1)(0, i);
  }

  std::vector<double> vec(dof, 0.0);
  result.intermediate_positions.reserve(number_of_waypoints - 2);
  for (int i = 1; i < number_of_waypoints - 1; ++i) {
    for (int j = 0; j < dof; ++j) {
      vec[j] = _waypoints(i, j);
    }
    result.intermediate_positions.emplace_back(vec);
  }
  ruckig::Ruckig<ruckig::DynamicDOFs> otg(dof, 0.001);
  ruckig::Trajectory<ruckig::DynamicDOFs> trajectory(dof);
  otg.calculate(result, trajectory);
  return result;
}

TEST(RUCKIG, ControlTest) {
#ifdef HAS_RUCKIG
  int dof = 5;
  ruckig::Ruckig<ruckig::DynamicDOFs> otg(dof, 0.001);
  ruckig::InputParameter<ruckig::DynamicDOFs> input(dof);
  ruckig::OutputParameter<ruckig::DynamicDOFs> output(dof);

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
