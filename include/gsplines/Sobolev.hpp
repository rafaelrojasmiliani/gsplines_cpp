
#ifndef SOBOLEV_H
#define SOBOLEV_H
#include <cstddef>
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>

#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
class SobolevNorm {
  friend class GSpline;

private:
  SobolevNorm &operator=(const SobolevNorm &);
  std::unique_ptr<basis::Basis> basis_;
  std::size_t num_intervals_;
  std::size_t codom_dim_;
  Interpolator interpolator_;
  SobolevNorm(const SobolevNorm &that);
  std::vector<std::pair<std::size_t, double>> weights_;
  Eigen::MatrixXd waypoints_;

  double inner_prod(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
                    const Eigen::Ref<const Eigen::VectorXd> _v1,
                    const Eigen::Ref<const Eigen::VectorXd> _v2);

protected:
  Eigen::MatrixXd matrix_;
  Eigen::MatrixXd matrix_2_;

public:
  SobolevNorm(const Eigen::Ref<const Eigen::MatrixXd> _waypoints,
              const basis::Basis &_basis,
              std::vector<std::pair<std::size_t, double>> _weights);
  virtual ~SobolevNorm() {}

  double operator()(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths);
  void deriv_wrt_interval_len(
      const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
      Eigen::Ref<Eigen::VectorXd> _buff);
};
} // namespace gsplines

#endif /* SOBOLEV_H */
