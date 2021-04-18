#include <gsplines++/interpolator.hpp>

namespace gsplines {

Interpolator::Interpolator(std::size_t _codom_dim, std::size_t _num_intervals,
                           basis::Basis &_basis)
    : basis_(&_basis), num_intervals_(_num_intervals),
      matrix_size_(_basis.get_dim() * _codom_dim * _num_intervals),
      boundary_buffer_tranposed_(_basis.get_dim(), _basis.get_dim() / 2),
      codom_dim_(_codom_dim), nnz_vec_() {
  interpolating_matrix_.resize(matrix_size_, matrix_size_);
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
    Eigen::Ref<Eigen::VectorXd> _interval_lengths) {

  for (unsigned int der = 0; der < basis_->get_dim() / 2; der++)
    basis_->eval_derivative_on_window(-1.0, _interval_lengths(0), der,
                                      boundary_buffer_tranposed_.col(der));
  unsigned int nnz = 0;
  unsigned int nptr = 0;
  unsigned int i0 = 0;
  for (int j = 0; j < basis_->get_dim() * codom_dim_; j++) {
    i0 = (j / basis_->get_dim()) * basis_->get_dim() / 2;
    for (int i = i0; i < i0 + basis_->get_dim() / 2; i++) {
      printf("rows %zu cols %zu\n", interpolating_matrix_.rows(),
             interpolating_matrix_.cols());
      printf("i %zu j %zu\n", i, j);
      printf("j % basis_->get_dim() %zu i - i0 %zu\n", j % basis_->get_dim(),
             i - i0);
      fflush(stdout);
      interpolating_matrix_.insert(i, j) =
          boundary_buffer_tranposed_(j % basis_->get_dim(), i - i0);
    }
    i0 += basis_->get_dim() / 2 * codom_dim_;
  }
  for (unsigned int der = 0; der < basis_->get_dim() / 2; der++)
    basis_->eval_derivative_on_window(1.0, _interval_lengths(0), der,
                                      boundary_buffer_tranposed_.col(der));
  if (num_intervals_ == 1) {
    for (int j = 0; j < basis_->get_dim() * codom_dim_; j++) {
      i0 = (j / basis_->get_dim()) * basis_->get_dim() / 2;
      for (int i = i0; i < i0 + basis_->get_dim() / 2; i++) {
        printf("rows %zu cols %zu\n", interpolating_matrix_.rows(),
               interpolating_matrix_.cols());
        printf("i %zu j %zu\n", i, j);
        printf("j % basis_->get_dim() %zu i - i0 %zu\n", j % basis_->get_dim(),
               i - i0);
        fflush(stdout);
        interpolating_matrix_.insert(i, j) =
            boundary_buffer_tranposed_(j % basis_->get_dim(), i - i0);
      }
      i0 += basis_->get_dim() / 2 * codom_dim_;
    }
    return;
  }
  for (int k = 0; k < interpolating_matrix_.outerSize(); ++k) {
    Eigen::SparseMatrix<double>::InnerIterator it(interpolating_matrix_, k);
    for (; it; ++it)
      printf("%zu, %zu\n", it.row(), it.col());
    printf("-\n");
  }
  fflush(stdout);
  printf("outerSize = %zu innerSize = %zu\n", interpolating_matrix_.outerSize(),
         interpolating_matrix_.innerSize());
  fflush(stdout);
}
} // namespace gsplines
