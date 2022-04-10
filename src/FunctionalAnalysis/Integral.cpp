#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <gsplines/FunctionalAnalysis/Integral.hpp>

namespace gsplines {
namespace functional_analysis {

double integral(const functions::FunctionBase &_in, std::size_t _n_glp,
                std::size_t _n_int) {

  assert(_in.get_codom_dim() == 1);

  Eigen::VectorXd glp;
  Eigen::VectorXd glw;
  double left_bound;
  double right_bound;
  std::tie(left_bound, right_bound) = _in.get_domain();
  double local_interval_length = (right_bound - left_bound) / _n_int;
  std::tie(glp, glw) =
      collocation::legendre_gauss_lobatto_points_and_weights(_n_glp);

  double result = 0.0;
  Eigen::MatrixXd local_value(_n_glp, 1);
  for (std::size_t interval = 0; interval < _n_int; interval++) {
    double local_left_bound = left_bound + interval * local_interval_length;
    Eigen::VectorXd local_points =
        (glp.array() + 1.0) * local_interval_length / 2.0 + local_left_bound;
    _in.value(local_points, local_value);
    result += (local_value.array() * glw.array()).col(0).sum() *
              local_interval_length / 2.0;
  }
  return result;
}
} // namespace functional_analysis

} // namespace gsplines
