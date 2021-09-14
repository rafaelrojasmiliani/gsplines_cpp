#ifndef IPOPT_SOLVER_H
#define IPOPT_SOLVER_H

#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>
#include <gsplines/Optimization/ipopt_interface.hpp>
namespace gsplines {

namespace optimization {

gsplines::GSpline
optimal_sobolev_norm(const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
                     const gsplines::basis::Basis &_basis,
                     std::vector<std::pair<std::size_t, double>> _weights,
                     double _exec_time);

gsplines::GSpline
broken_lines_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);

gsplines::GSpline
minimum_acceleration_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);

gsplines::GSpline
minimum_jerk_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);

gsplines::GSpline
minimum_snap_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);

gsplines::GSpline
minimum_crackle_path(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);
} // namespace optimization
} // namespace gsplines
#endif /* IPOPT_SOLVER_H */
