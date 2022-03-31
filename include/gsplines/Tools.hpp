
#include <Eigen/Core>
namespace gsplines {
namespace tools {

bool approx_equal(double _lhs, double _rhs, double _tol);

bool approx_equal(const Eigen::MatrixXd &_lhs, const Eigen::MatrixXd &_rhs,
                  double _tol);

} // namespace tools

} // namespace gsplines
