#include <eigen3/Eigen/SparseLU>
#include <gsplines++/Interpolator.hpp>
#include <iostream>

namespace gsplines {

Interpolator::Interpolator(std::size_t _codom_dim, std::size_t _num_intervals,
                           basis::Basis &_basis)
    : basis_(&_basis), codom_dim_(_codom_dim), num_intervals_(_num_intervals),
      matrix_size_(_basis.get_dim() * _codom_dim * _num_intervals),
      boundary_buffer_tranposed_(_basis.get_dim(), _basis.get_dim() / 2) {

  interpolating_matrix_.resize(matrix_size_, matrix_size_);
  interpolating_vector_.resize(matrix_size_);
  interpolating_vector_.setZero();
  coefficients_vector_.resize(matrix_size_);

  derivative_buffer_tranposed_.resize(basis_->get_dim(), basis_->get_dim() - 2);

  position_buffer_.resize(basis_->get_dim());

  std::size_t nnz_per_col =
      2 * (basis_->get_dim() / 2 - 1) + 2 + 2 * (basis_->get_dim() - 2);
  interpolating_matrix_.reserve(
      Eigen::VectorXi::Constant(matrix_size_, nnz_per_col));
  interpolating_matrix_.makeCompressed();

  for (int i = 0; i < _num_intervals; i++) {
    {
      Eigen::SparseMatrix<double> mat;
      interpolation_matrix_ders_.push_back(mat);
      interpolation_matrix_ders_.back().resize(matrix_size_, matrix_size_);
      interpolation_matrix_ders_.back().reserve(nnz_per_col *
                                                basis_->get_dim() * codom_dim_);
      interpolation_matrix_ders_.back().makeCompressed();
    }
  }
}

Interpolator::~Interpolator() {}
void Interpolator::fill_interpolating_matrix(
    const Eigen::Ref<const Eigen::VectorXd> _interval_lengths) {

  unsigned int interval_coor = 0;
  unsigned int i0 = 0;
  unsigned int j0 = 0;

  if (_interval_lengths.size() != num_intervals_) {
    fprintf(stderr, "Cannot fill matrix. Vector of interval lenghts is not of "
                    "the requred dimention");
    return;
  }
  fill_buffers(-1.0, _interval_lengths(0));

  fill_position_block(i0, j0, interpolating_matrix_);
  i0 += codom_dim_;

  fill_boundary_derivative_block(i0, j0, interpolating_matrix_);
  i0 += codom_dim_ * (basis_->get_dim() / 2 - 1);

  fill_buffers(1.0, _interval_lengths(0));

  fill_position_block(i0, j0, interpolating_matrix_);
  i0 += codom_dim_;

  if (num_intervals_ == 1) {
    fill_boundary_derivative_block(i0, j0, interpolating_matrix_);
    return;
  }

  fill_continuity_derivative_block(i0, j0, false, interpolating_matrix_);
  j0 += codom_dim_ * basis_->get_dim();

  for (interval_coor = 1; interval_coor < num_intervals_ - 1; interval_coor++) {
    fill_buffers(-1.0, _interval_lengths(interval_coor));

    fill_continuity_derivative_block(i0, j0, true, interpolating_matrix_);
    i0 += codom_dim_ * (basis_->get_dim() - 2);

    fill_position_block(i0, j0, interpolating_matrix_);
    i0 += codom_dim_;

    fill_buffers(1.0, _interval_lengths(interval_coor));

    fill_position_block(i0, j0, interpolating_matrix_);
    i0 += codom_dim_;

    fill_continuity_derivative_block(i0, j0, false, interpolating_matrix_);
    j0 += codom_dim_ * basis_->get_dim();
  }
  fill_buffers(-1.0, _interval_lengths(num_intervals_ - 1));

  fill_continuity_derivative_block(i0, j0, true, interpolating_matrix_);
  i0 += codom_dim_ * (basis_->get_dim() - 2);

  fill_position_block(i0, j0, interpolating_matrix_);
  i0 += codom_dim_;

  fill_buffers(1.0, _interval_lengths(num_intervals_ - 1));
  fill_boundary_derivative_block(i0, j0, interpolating_matrix_);
  i0 += codom_dim_ * (basis_->get_dim() / 2 - 1);

  fill_position_block(i0, j0, interpolating_matrix_);
}
void Interpolator::fill_interpolating_vector(
    const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  unsigned int interval_coor = 0;
  unsigned int codom_coor = 0;
  unsigned int vector_coor = 0;
  for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    interpolating_vector_(vector_coor) = _waypoints(0, codom_coor);
    vector_coor++;
  }
  vector_coor += codom_dim_ * (basis_->get_dim() / 2 - 1);
  for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    interpolating_vector_(vector_coor) = _waypoints(1, codom_coor);
    vector_coor++;
  }
  if (num_intervals_ == 1) {
    return;
  }
  vector_coor += codom_dim_ * (basis_->get_dim() - 2);
  for (interval_coor = 1; interval_coor < num_intervals_ - 1; interval_coor++) {
    for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
      interpolating_vector_(vector_coor) =
          _waypoints(interval_coor, codom_coor);
      vector_coor++;
    }
    for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
      interpolating_vector_(vector_coor) =
          _waypoints(interval_coor + 1, codom_coor);
      vector_coor++;
    }
    vector_coor += codom_dim_ * (basis_->get_dim() - 2);
  }
  for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    interpolating_vector_(vector_coor) =
        _waypoints(num_intervals_ - 1, codom_coor);
    vector_coor++;
  }
  vector_coor += codom_dim_ * (basis_->get_dim() / 2 - 1);
  for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    interpolating_vector_(vector_coor) = _waypoints(num_intervals_, codom_coor);
    vector_coor++;
  }
}

PiecewiseFunction Interpolator::interpolate(
    const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
    const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {
  const Eigen::VectorXd &vector_resut =
      solve_interpolation(_interval_lengths, _waypoints);
  return PiecewiseFunction(codom_dim_, num_intervals_, *basis_, vector_resut,
                           _interval_lengths);
}

void Interpolator::print_interpolating_matrix() {
  Eigen::IOFormat CleanFmt(14, 0, ", ", "\n", "[", "]");
  std::cout << Eigen::MatrixXd(interpolating_matrix_).format(CleanFmt) << '\n';
}

void Interpolator::print_interpolating_vector() {
  Eigen::IOFormat CleanFmt(14, 1, ", ", "\n", "[", "]");
  std::cout << Eigen::MatrixXd(interpolating_vector_).format(CleanFmt) << '\n';
}
void Interpolator::fill_position_block(unsigned int i0, unsigned int j0,
                                       Eigen::SparseMatrix<double> &_mat) {
  unsigned int i, j;
  for (unsigned int codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    for (unsigned int basis_coor = 0; basis_coor < basis_->get_dim();
         basis_coor++) {
      i = i0 + codom_coor;
      j = j0 + basis_->get_dim() * codom_coor + basis_coor;
      _mat.coeffRef(i, j) = position_buffer_(basis_coor);
    }
  }
}
void Interpolator::fill_continuity_derivative_block(
    unsigned int _i0, unsigned int _j0, bool _minus,
    Eigen::SparseMatrix<double> &_mat) {
  unsigned int i, j;
  unsigned int codom_coor;
  unsigned int basis_coor;
  unsigned int der_coor;
  double mutiplier;
  if (_minus) {
    mutiplier = -1.0;
  } else {
    mutiplier = 1.0;
  }
  for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    for (basis_coor = 0; basis_coor < basis_->get_dim(); basis_coor++) {
      for (der_coor = 0; der_coor < basis_->get_dim() - 2; der_coor++) {
        i = _i0 + codom_coor * (basis_->get_dim() - 2) + der_coor;
        j = _j0 + basis_->get_dim() * codom_coor + basis_coor;
        _mat.coeffRef(i, j) =
            mutiplier * derivative_buffer_tranposed_(basis_coor, der_coor);
      }
    }
  }
}

void Interpolator::fill_boundary_derivative_block(
    unsigned int _i0, unsigned int _j0, Eigen::SparseMatrix<double> &_mat) {
  unsigned int i, j;
  unsigned int codom_coor;
  unsigned int basis_coor;
  unsigned int der_coor;
  for (codom_coor = 0; codom_coor < codom_dim_; codom_coor++) {
    for (der_coor = 0; der_coor < basis_->get_dim() / 2 - 1; der_coor++) {
      for (basis_coor = 0; basis_coor < basis_->get_dim(); basis_coor++) {
        i = _i0 + codom_coor * (basis_->get_dim() / 2 - 1) + der_coor;
        j = _j0 + basis_->get_dim() * codom_coor + basis_coor;

        _mat.coeffRef(i, j) =
            derivative_buffer_tranposed_(basis_coor, der_coor);
      }
    }
  }
}
void Interpolator::fill_buffers(double s, double tau) {

  basis_->eval_derivative_on_window(s, tau, 0, position_buffer_);

  for (unsigned int der = 0; der < basis_->get_dim() - 2; der++) {
    basis_->eval_derivative_on_window(s, tau, der + 1,
                                      derivative_buffer_tranposed_.col(der));
  }
}

void Interpolator::fill_buffers_deriv_wrt_tau(double s, double tau) {

  basis_->eval_derivative_wrt_tau_on_window(s, tau, 0, position_buffer_);

  for (unsigned int der = 0; der < basis_->get_dim() - 2; der++) {
    basis_->eval_derivative_wrt_tau_on_window(
        s, tau, der + 1, derivative_buffer_tranposed_.col(der));
  }
}

Eigen::Ref<const Eigen::VectorXd> Interpolator::get_coeff_derivative_wrt_tau(
    Eigen::Ref<const Eigen::VectorXd> _coeff,
    Eigen::Ref<const Eigen::VectorXd> _interval_lengths, std::size_t _tau_idx) {

  if (_interval_lengths.size() != num_intervals_ or _tau_idx < 0 or
      _tau_idx >= num_intervals_) {
    fprintf(stderr, "Cannot fill matrix. Vector of interval lenghts is not of "
                    "the requred dimention");
    return coefficients_vector_;
  }

  Eigen::SparseMatrix<double> &mat = interpolation_matrix_ders_[_tau_idx];
  unsigned int i0 = 0;
  unsigned int j0 = 0;
  if (_tau_idx == 0) {

    fill_buffers_deriv_wrt_tau(-1.0, _interval_lengths(0));

    fill_position_block(i0, j0, mat);
    i0 += codom_dim_;

    fill_boundary_derivative_block(i0, j0, mat);
    i0 += codom_dim_ * (basis_->get_dim() / 2 - 1);

    fill_buffers_deriv_wrt_tau(1.0, _interval_lengths(0));

    fill_position_block(i0, j0, mat);
    i0 += codom_dim_;

    if (num_intervals_ == 1) {
      fill_boundary_derivative_block(i0, j0, mat);
    } else {
      fill_continuity_derivative_block(i0, j0, false, mat);
    }
  } else if (0 < _tau_idx and _tau_idx < num_intervals_ - 1) {

    // rows added in the first interval
    i0 = 2 * codom_dim_ + codom_dim_ * (basis_->get_dim() / 2 - 1);

    i0 += (codom_dim_ * (basis_->get_dim() - 2) + 2 * codom_dim_) *
          (_tau_idx - 1);

    j0 = _tau_idx * basis_->get_dim() * codom_dim_;

    fill_buffers_deriv_wrt_tau(-1.0, _interval_lengths(_tau_idx));

    fill_continuity_derivative_block(i0, j0, true, mat);
    i0 += codom_dim_ * (basis_->get_dim() - 2);

    fill_position_block(i0, j0, mat);
    i0 += codom_dim_;

    fill_buffers_deriv_wrt_tau(1.0, _interval_lengths(_tau_idx));

    fill_position_block(i0, j0, mat);
    i0 += codom_dim_;

    fill_continuity_derivative_block(i0, j0, false, mat);

  } else if (_tau_idx == num_intervals_ - 1) {

    i0 = (basis_->get_dim() / 2 + 1) * codom_dim_ +
         basis_->get_dim() * codom_dim_ * (_tau_idx - 1);

    j0 = _tau_idx * basis_->get_dim() * codom_dim_;

    fill_buffers_deriv_wrt_tau(-1.0, _interval_lengths(num_intervals_ - 1));

    fill_continuity_derivative_block(i0, j0, true, mat);
    i0 += codom_dim_ * (basis_->get_dim() - 2);

    fill_position_block(i0, j0, mat);
    i0 += codom_dim_;

    fill_buffers_deriv_wrt_tau(1.0, _interval_lengths(num_intervals_ - 1));
    fill_boundary_derivative_block(i0, j0, mat);
    i0 += codom_dim_ * (basis_->get_dim() / 2 - 1);

    fill_position_block(i0, j0, mat);
  }
  mat.makeCompressed();
  coefficients_vector_ = mat * _coeff;
  fill_interpolating_matrix(_interval_lengths);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(interpolating_matrix_);
  coefficients_vector_ = -solver.solve(coefficients_vector_);
  return coefficients_vector_;
}
const Eigen::Ref<const Eigen::VectorXd> Interpolator::solve_interpolation(
    const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
    const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  // 1. fill the interpolating matrix
  fill_interpolating_matrix(_interval_lengths);
  // 2. fill the interpolating vector
  fill_interpolating_vector(_waypoints);
  // 3. Solve the interpolation problem
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(interpolating_matrix_);
  sol_buffer_ = solver.solve(interpolating_vector_);
  // 4. Return the interpolating function
  return sol_buffer_;
}
} // namespace gsplines
