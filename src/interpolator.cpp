#include <eigen3/Eigen/SparseLU>
#include <gsplines++/interpolator.hpp>
#include <iostream>

namespace gsplines {

Interpolator::Interpolator(std::size_t _codom_dim, std::size_t _num_intervals,
                           basis::Basis &_basis)
    : basis_(&_basis), num_intervals_(_num_intervals),
      matrix_size_(_basis.get_dim() * _codom_dim * _num_intervals),
      boundary_buffer_tranposed_(_basis.get_dim(), _basis.get_dim() / 2),
      codom_dim_(_codom_dim), nnz_vec_(), interpolating_vector_(matrix_size_) {
  interpolating_matrix_.resize(matrix_size_, matrix_size_);
  interpolating_vector_.resize(matrix_size_);

  printf("rows %zu cols %zu\n", interpolating_matrix_.rows(),
         interpolating_matrix_.cols());
  printf("matriz size %zu\n", matrix_size_);
  std::size_t nnz_per_col =
      2 * (basis_->get_dim() / 2 - 1) + 2 + 2 * (basis_->get_dim() - 2);
  interpolating_matrix_.reserve(
      Eigen::VectorXi::Constant(matrix_size_, nnz_per_col));
  interpolating_matrix_.makeCompressed();
  printf("rows %zu cols %zu\n", interpolating_matrix_.rows(),
         interpolating_matrix_.cols());
}

Interpolator::~Interpolator() {}
void Interpolator::fill_interpolating_matrix(
    const Eigen::Ref<const Eigen::VectorXd> _interval_lengths) {

  Eigen::IOFormat CleanFmt(14, 1, ", ", "\n", "[", "]");
  for (unsigned int der = 0; der < basis_->get_dim() / 2; der++)
    basis_->eval_derivative_on_window(-1.0, _interval_lengths(0), der,
                                      boundary_buffer_tranposed_.col(der));
  unsigned int nnz = 0;
  unsigned int nptr = 0;
  unsigned int i0 = 0;
  for (int j = 0; j < basis_->get_dim() * codom_dim_; j++) {
    i0 = (j / basis_->get_dim()) * basis_->get_dim() / 2;
    for (int i = i0; i < i0 + basis_->get_dim() / 2; i++) {
      printf("i %i j %i\n", i, j);
      fflush(stdout);
      interpolating_matrix_.insert(i, j) =
          boundary_buffer_tranposed_(j % basis_->get_dim(), i - i0);
    }
  }
  for (unsigned int der = 0; der < basis_->get_dim() / 2; der++)
    basis_->eval_derivative_on_window(1.0, _interval_lengths(0), der,
                                      boundary_buffer_tranposed_.col(der));
  int i0_prev = codom_dim_ * basis_->get_dim() / 2;
  if (num_intervals_ == 1) {
    for (int j = 0; j < basis_->get_dim() * codom_dim_; j++) {
      i0 = (j / basis_->get_dim()) * basis_->get_dim() / 2 + i0_prev;
      for (int i = i0; i < i0 + basis_->get_dim() / 2; i++) {
        printf("i %i j %i\n", i, j);
        fflush(stdout);
        interpolating_matrix_.insert(i, j) =
            boundary_buffer_tranposed_(j % basis_->get_dim(), i - i0);
      }
    }
    std::cout << Eigen::MatrixXd(interpolating_matrix_).format(CleanFmt)
              << '\n';
    return;
  }
}
void Interpolator::fill_interpolating_vector(
    const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {

  int i = 0;
  for (int d = 0; d < codom_dim_; d++) {
    interpolating_vector_(i) = _waypoints(0, d);
    i++;
    for (int j = 1; j < basis_->get_dim() / 2; j++) {
      printf("i=%i\n", i);
      interpolating_vector_(i) = 0;
      i++;
    }
  }
  if (num_intervals_ == 1) {
    for (int d = 0; d < codom_dim_; d++) {
      interpolating_vector_(i) = _waypoints(1, d);
      i++;
      for (int j = 1; j < basis_->get_dim() / 2; j++) {
        interpolating_vector_(i) = 0;
        i++;
      }
    }
    return;
  }
}
void Interpolator::interpolate(
    const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
    const Eigen::Ref<const Eigen::MatrixXd> _waypoints) {
  Eigen::VectorXd vector_resut;
  fill_interpolating_matrix(_interval_lengths);
  fill_interpolating_vector(_waypoints);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(interpolating_matrix_);
  vector_resut = solver.solve(interpolating_vector_);
}
} // namespace gsplines
