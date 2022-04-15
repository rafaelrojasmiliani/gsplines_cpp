
#include <Eigen/Core>
#include <gsplines/GSpline.hpp>
namespace gsplines {
namespace tools {

bool approx_equal(double _lhs, double _rhs, double _tol);

bool approx_equal(const Eigen::MatrixXd &_lhs, const Eigen::MatrixXd &_rhs,
                  double _tol);

bool approx_equal(const GSplineBase &_lhs, const GSplineBase &_rhs,
                  double _tol);

bool approx_zero(double _rhs, double _tol);

bool approx_zero(const Eigen::MatrixXd &_rhs, double _tol);

bool approx_zero(const GSplineBase &_rhs, double _tol);

} // namespace tools

} // namespace gsplines
