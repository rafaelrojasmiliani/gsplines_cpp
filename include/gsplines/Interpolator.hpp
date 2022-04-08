#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include <eigen3/Eigen/SparseCore>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/GSpline.hpp>

namespace gsplines {
class Interpolator {
private:
  Interpolator(const Interpolator &that);
  Interpolator &operator=(const Interpolator &);
  std::unique_ptr<basis::Basis> basis_;
  std::size_t codom_dim_;
  std::size_t num_intervals_;
  std::size_t matrix_size_;
  std::size_t nnz_size_;
  std::size_t k_factor_; // maxisum derivative present on the cost function
  Eigen::SparseMatrix<double> interpolating_matrix_;
  Eigen::VectorXd interpolating_vector_;
  Eigen::VectorXd coefficients_vector_;
  Eigen::MatrixXd
      boundary_buffer_tranposed_; // this is  basis.dim \times basis.dim/2
  Eigen::MatrixXd
      continuity_buffer_tranposed_; // this is  basis.dim \times (basis.dim -2)
  Eigen::MatrixXd
      derivative_buffer_tranposed_; // this is (basis.dim - 2) \times basis.dim
  std::vector<Eigen::SparseMatrix<double>> interpolation_matrix_ders_;
  Eigen::VectorXd position_buffer_; // this is basis.dim vector
  Eigen::VectorXi nnz_vec_;
  Eigen::VectorXd sol_buffer_;

public:
  Interpolator(std::size_t _codom_dim, std::size_t _num_intervals,
               const basis::Basis &_basis);
  virtual ~Interpolator();
  void fill_interpolating_matrix(
      const Eigen::Ref<const Eigen::VectorXd> _interval_lengths);
  void
  fill_interpolating_vector(const Eigen::Ref<const Eigen::MatrixXd> _waypoints);
  GSpline interpolate(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
                      const Eigen::Ref<const Eigen::MatrixXd> _waypoints);
  void print_interpolating_matrix();
  void print_interpolating_vector();

  void print_info();

  void fill_position_block(unsigned int i0, unsigned int j0,
                           Eigen::SparseMatrix<double> &_mat);
  void fill_boundary_derivative_block(unsigned int i0, unsigned int j0,
                                      Eigen::SparseMatrix<double> &_mat);
  void fill_continuity_derivative_block(unsigned int i0, unsigned int j0,
                                        bool _minus,
                                        Eigen::SparseMatrix<double> &_mat);
  void fill_buffers(double s, double tau);

  void fill_buffers_deriv_wrt_tau(double s, double tau);

  Eigen::Ref<const Eigen::VectorXd> get_coeff_derivative_wrt_tau(
      Eigen::Ref<const Eigen::VectorXd> _coeff,
      Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
      std::size_t _tau_idx);

  const Eigen::Ref<const Eigen::VectorXd>
  solve_interpolation(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
                      const Eigen::Ref<const Eigen::MatrixXd> _waypoints);
};

GSpline interpolate(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
                    const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
                    const basis::Basis &_basis);
} // namespace gsplines
#endif /* INTERPOLATOR_H */
