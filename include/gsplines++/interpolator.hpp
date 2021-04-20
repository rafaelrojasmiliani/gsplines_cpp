#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include <eigen3/Eigen/SparseCore>
#include <gsplines++/Basis.hpp>
#include <gsplines++/PiecewiseFunction.hpp>

namespace gsplines {

class Interpolator {
private:
  Interpolator(const Interpolator &that);
  Interpolator &operator=(const Interpolator &);
  basis::Basis *const basis_;
  std::size_t codom_dim_;
  std::size_t num_intervals_;
  Eigen::SparseMatrix<double> interpolating_matrix_;
  Eigen::VectorXd interpolating_vector_;
  std::size_t matrix_size_;
  std::size_t nnz_size_;
  std::size_t k_factor_; // maxisum derivative present on the cost function
  Eigen::MatrixXd
      boundary_buffer_tranposed_; // this is  basis.dim \times basis.dim/2
  Eigen::VectorXi nnz_vec_;

  void fill_matrix();

public:
  Interpolator(std::size_t _codom_dim, std::size_t _num_intervals,
               basis::Basis &_basis);
  virtual ~Interpolator();
  void fill_interpolating_matrix(
      const Eigen::Ref<const Eigen::VectorXd> _interval_lengths);
  void
  fill_interpolating_vector(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);
  PiecewiseFunction
  interpolate(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
              const Eigen::Ref<const Eigen::MatrixXd> _waypoints);
};

} // namespace gsplines
#endif /* INTERPOLATOR_H */
