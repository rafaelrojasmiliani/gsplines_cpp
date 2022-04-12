#ifndef IPOPT_SOLVER_H
#define IPOPT_SOLVER_H

#include "gsplines/Functions/FunctionBase.hpp"
#include <cstddef>
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

gsplines::GSpline best_approximation(const functions::FunctionBase &_in,
                                     std::size_t _n_glp, std::size_t _n_inter);
} // namespace optimization
} // namespace gsplines
#endif /* IPOPT_SOLVER_H */
