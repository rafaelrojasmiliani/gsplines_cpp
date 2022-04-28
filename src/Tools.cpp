#include <gsplines/Tools.hpp>

namespace gsplines {
namespace tools {

bool approx_equal(double _lhs, double _rhs, double _tol) {

  double err = std::fabs(_lhs - _rhs);

  double lhs_abs = std::fabs(_lhs);
  double rhs_abs = std::fabs(_rhs);

  if (lhs_abs < _tol or rhs_abs < _tol) {
    return err < _tol;
  }
  return err < lhs_abs * _tol and err < rhs_abs * _tol;
}

bool approx_equal(const Eigen::MatrixXd &_lhs, const Eigen::MatrixXd &_rhs,
                  double _tol) {

  if (_lhs.rows() != _rhs.rows() or _lhs.cols() != _rhs.cols() or
      _lhs.size() != _rhs.size())
    return false;
  if (_lhs.size() == 0)
    return true;

  double err = (_lhs - _rhs).array().abs().maxCoeff();
  double lhs_max = (_lhs).array().abs().maxCoeff();
  double rhs_max = (_rhs).array().abs().maxCoeff();

  if (lhs_max < _tol or rhs_max < _tol) {
    return err < _tol;
  }
  return err / lhs_max < _tol and err / rhs_max < _tol;
}

bool approx_equal(const GSplineBase &_lhs, const GSplineBase &_rhs,
                  double _tol) {
  return approx_equal(_lhs.get_coefficients(), _rhs.get_coefficients(),
                      _tol) and
         approx_equal(_lhs.get_interval_lengths(), _rhs.get_interval_lengths(),
                      _tol) and
         _lhs.get_basis() == _rhs.get_basis();
}

bool approx_zero(double _rhs, double _tol) {

  double rhs_abs = std::fabs(_rhs);

  return (rhs_abs < _tol);
}

bool approx_zero(const Eigen::MatrixXd &_rhs, double _tol) {

  if (_rhs.size() == 0)
    return false;

  double rhs_max = (_rhs).array().abs().maxCoeff();

  return rhs_max < _tol;
}

bool approx_zero(const GSplineBase &_rhs, double _tol) {
  return approx_zero(_rhs.get_coefficients(), _tol);
}
} // namespace tools

} // namespace gsplines
