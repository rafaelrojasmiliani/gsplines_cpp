#ifndef IPOPT_SOLVER_H
#define IPOPT_SOLVER_H

#include <gsplines/Interpolator.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/ipopt_interface.hpp>

namespace gsplines_opt {

gsplines::GSpline
optimal_sobolev_norm(const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
                     const gsplines::basis::Basis &_basis,
                     std::vector<std::pair<std::size_t, double>> _weights,
                     double _exec_time);

} // namespace gsplines_opt
#endif /* IPOPT_SOLVER_H */
