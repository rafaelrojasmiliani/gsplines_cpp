#include <gsplines/Basis.hpp>
#include <iostream>
#include <math.h>
namespace gsplines {

namespace basis {
/** Returns
 * \int_-1^1 D*D^T d t
 * */
void Basis::init_sobolev_ip_matrix(std::vector<double> _weights) {

  sobolev_ip_matrix_ = _weights[0] * Eigen::MatrixXd::Identity(dim_, dim_);
  Eigen::MatrixXd mat(derivative_matrix_);
  for (unsigned int i = 1; i < _weights.size(); i++) {
    sobolev_ip_matrix_ += _weights[i] * mat * mat.transpose();
    mat *= mat;
  }
}
} // namespace basis
} // namespace gsplines
