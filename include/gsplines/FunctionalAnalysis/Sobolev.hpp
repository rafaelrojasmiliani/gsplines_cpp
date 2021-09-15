
#ifndef SOBOLEV_H
#define SOBOLEV_H
#include <cstddef>
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Functions/FunctionBase.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/GSpline.hpp>
#include <gsplines/Interpolator.hpp>

#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functional_analysis {

class SobolevNorm {
  friend class GSpline;

private:
  SobolevNorm &operator=(const SobolevNorm &) = delete;
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
  virtual ~SobolevNorm() = default;

  double operator()(const Eigen::Ref<const Eigen::VectorXd> _interval_lengths);
  void deriv_wrt_interval_len(
      const Eigen::Ref<const Eigen::VectorXd> _interval_lengths,
      Eigen::Ref<Eigen::VectorXd> _buff);
};

double
sobolev_semi_inner_product(const functions::FunctionBase &_lhs,
                           const functions::FunctionBase &_rhs,
                           std::vector<std::pair<std::size_t, double>> _terms,
                           std::size_t _n_glp = 10, std::size_t _n_int = 1);

double sobolev_semi_norm(const functions::FunctionBase &,
                         std::vector<std::pair<std::size_t, double>> _terms,
                         std::size_t _n_glp = 10, std::size_t _n_int = 1);

double sobolev_seminorm(const functions::FunctionBase &_in,
                        const std::vector<std::pair<std::size_t, double>> &,
                        std::size_t _n_glp = 10, std::size_t _n_int = 1);

double sobolev_norm(const functions::FunctionBase &, std::size_t _deriv_deg,
                    std::size_t _n_glp = 10, std::size_t _n_int = 1);

double l2_norm(const functions::FunctionBase &_in, std::size_t _n_glp = 10,
               std::size_t _n_int = 1);

double l2_inner_product(const functions::FunctionBase &_lhs,
                        const functions::FunctionBase &_rhs,
                        std::size_t _n_glp = 10, std::size_t _n_int = 1);

double
sobolev_semi_inner_product(const GSpline &_lhs, const GSpline &_rhs,
                           std::vector<std::pair<std::size_t, double>> _terms,
                           std::size_t _n_glp = 10, std::size_t _n_int = 1);

double sobolev_semi_norm(const GSpline &,
                         std::vector<std::pair<std::size_t, double>> _terms,
                         std::size_t _n_glp = 10, std::size_t _n_int = 1);

double sobolev_seminorm(const GSpline &_in,
                        const std::vector<std::pair<std::size_t, double>> &,
                        std::size_t _n_glp = 10, std::size_t _n_int = 1);

double sobolev_norm(const GSpline &, std::size_t _deriv_deg,
                    std::size_t _n_glp = 10, std::size_t _n_int = 1);

double l2_norm(const GSpline &_in, std::size_t _n_glp = 10,
               std::size_t _n_int = 1);

double l2_inner_product(const GSpline &_lhs, const GSpline &_rhs,
                        std::size_t _n_glp = 10, std::size_t _n_int = 1);
} // namespace functional_analysis
} // namespace gsplines

#endif /* SOBOLEV_H */
