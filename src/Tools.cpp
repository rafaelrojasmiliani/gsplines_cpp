#include <gsplines/Tools.hpp>

namespace gsplines {
namespace tools {

bool approx_equal(double _lhs, double _rhs, double _tol) {

  double err = std::fabs(_lhs - _rhs);

  if (_lhs < _tol or _rhs < _tol) {
    return err < _tol;
  }
  return err / _lhs < _tol and err / _rhs < _tol;
}

bool approx_equal(const Eigen::MatrixXd &_lhs, const Eigen::MatrixXd &_rhs,
                  double _tol) {

  if (_lhs.rows() != _rhs.rows() or _lhs.cols() != _rhs.cols())
    return false;

  double err = (_lhs - _rhs).array().abs().maxCoeff();
  double lhs_max = (_lhs).array().abs().maxCoeff();
  double rhs_max = (_rhs).array().abs().maxCoeff();

  if (lhs_max < _tol or rhs_max < _tol) {
    return err < _tol;
  }
  return err / lhs_max < _tol and err / rhs_max < _tol;
}

} // namespace tools

} // namespace gsplines
