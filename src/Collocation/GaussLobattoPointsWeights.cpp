

#include <cmath>
#include <gsplines/Collocation/GaussLobattoPointsWeights.hpp>
#include <tuple>
#include <vector>

namespace gsplines {
namespace collocation {
std::tuple<double, double, double> q_and_evaluation(double _t, std::size_t _n) {
  /*
   *  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 24: qAndLEvaluation
   * */
  double l_n_minus_2 = 1.0;
  double l_n_minus_1 = _t;
  double lprime_n_minus_2 = 0.0;
  double lprime_n_minus_1 = 1.0;
  double l_n;
  double lprime_n;
  double lprime_n_plus_1;
  double l_n_plus_1;
  double q;
  double qprime;

  std::size_t uick;

  for (uick = 2; uick < _n + 1; uick++) {

    l_n = (2.0 * (double)uick - 1.0) / ((double)uick) * _t * l_n_minus_1 -
          ((double)uick - 1.0) / ((double)uick) * l_n_minus_2;

    lprime_n = lprime_n_minus_2 + (2.0 * (double)uick - 1.0) * l_n_minus_1;

    l_n_minus_2 = l_n_minus_1;
    l_n_minus_1 = l_n;

    lprime_n_minus_2 = lprime_n_minus_1;
    lprime_n_minus_1 = lprime_n;
  }

  l_n_plus_1 = (2.0 * (double)uick - 1.0) / ((double)uick) * _t * l_n -
               ((double)uick - 1.0) / ((double)uick) * l_n_minus_2;

  lprime_n_plus_1 = lprime_n_minus_2 + (2.0 * (double)uick - 1.0) * l_n;

  q = l_n_plus_1 - l_n_minus_2;

  qprime = lprime_n_plus_1 - lprime_n_minus_2;

  return std::tuple<double, double, double>(q, qprime, l_n);
}

std::pair<Eigen::VectorXd, Eigen::VectorXd>
legendre_gauss_lobatto_points_and_weights(std::size_t _n) {
  /*
   *  David A. Kopriva
   *  Implementing Spectral
   *  Methods for Partial
   *  Differential Equations
   *  Algorithm 25: LegendreGaussLobattoNodesAndWeights
   * */

  Eigen::VectorXd points(_n);
  Eigen::VectorXd weights(_n);
  double q;
  double q_prime;
  double l_n;
  std::size_t book_n = _n - 1;

  if (_n == 2) {
    points(0) = -1.0;
    points(_n - 1) = 1.0;
    weights(0) = 1.0;
    weights(_n - 1) = 1.0;
  } else {
    points(0) = -1.0;
    points(_n - 1) = 1.0;
    weights(0) = 2.0 / ((double)book_n * ((double)book_n + 1.0));
    weights(_n - 1) = weights(0);
    for (std::size_t uicj = 1;
         uicj <= std::trunc(((double)book_n + 1.0) / 2) - 1; uicj++) {

      double arg =
          (((double)uicj + 0.25) * M_PI / (double)book_n) -
          (3.0 / (8.0 * (double)book_n * M_PI)) * (1.0 / ((double)uicj + 0.25));

      points(uicj) = -std::cos(arg);

      for (std::size_t uick = 0; uick < 100; uick++) {
        std::tie(q, q_prime, l_n) = q_and_evaluation(points(uicj), book_n);

        double delta = -q / q_prime;
        points(uicj) += delta;
        if (std::fabs(delta) < 1.0e-9 * std::fabs(points(uicj)))
          break;
      }

      points(book_n - uicj) = -points(uicj);

      std::tie(q, q_prime, l_n) = q_and_evaluation(points[uicj], book_n);

      weights(uicj) =
          2.0 / ((double)book_n * ((double)book_n + 1.0) * std::pow(l_n, 2.0));

      weights(book_n - uicj) = weights(uicj);
    }
    if (book_n % 2u == 0) {
      std::tie(q, q_prime, l_n) = q_and_evaluation(0.0, book_n);
      points[book_n / 2] = 0.0;
      weights[book_n / 2] =
          2.0 / ((double)book_n * ((double)book_n + 1.0) * std::pow(l_n, 2.0));
    }
  }

  return std::pair<Eigen::VectorXd, Eigen::VectorXd>(std::move(points),
                                                     std::move(weights));
}

Eigen::VectorXd legendre_gauss_lobatto_points(std::size_t _n) {
  Eigen::VectorXd result;
  std::tie(result, std::ignore) = legendre_gauss_lobatto_points_and_weights(_n);
  return result;
}

Eigen::VectorXd legendre_gauss_lobatto_weights(std::size_t _n) {

  Eigen::VectorXd result;
  std::tie(std::ignore, result) = legendre_gauss_lobatto_points_and_weights(_n);
  return result;
}

Eigen::VectorXd
legendre_gauss_lobatto_points(std::tuple<double, double> _domain,
                              std::size_t _nglp, std::size_t _n_intervals) {

  Eigen::VectorXd glp = legendre_gauss_lobatto_points(_nglp);

  double left_bound;
  double right_bound;
  std::tie(left_bound, right_bound) = _domain;

  Eigen::VectorXd result(_nglp * _n_intervals);

  double local_interval_length = (right_bound - left_bound) / _n_intervals;

  for (std::size_t interval = 0; interval < _n_intervals; interval++) {

    double local_left_bound = left_bound + interval * local_interval_length;

    result.segment(interval * _nglp, _nglp) =
        (glp.array() + 1.0) * local_interval_length / 2.0 + local_left_bound;
  }
  return result;
}

Eigen::VectorXd
legendre_gauss_lobatto_weights(std::tuple<double, double> _domain,
                               std::size_t _nglp, std::size_t _n_intervals) {

  Eigen::VectorXd glw = legendre_gauss_lobatto_weights(_nglp);

  double left_bound;
  double right_bound;
  std::tie(left_bound, right_bound) = _domain;

  Eigen::VectorXd result(_nglp * _n_intervals);

  double local_interval_length = (right_bound - left_bound) / _n_intervals;

  for (std::size_t interval = 0; interval < _n_intervals; interval++) {

    result.segment(interval * _nglp, _nglp) = glw * local_interval_length / 2.0;
  }
  return result;
}

Eigen::VectorXd
legendre_gauss_lobatto_points(double _left_bound, std::size_t _nglp,
                              const Eigen::VectorXd &_intervals) {

  Eigen::VectorXd glp = legendre_gauss_lobatto_points(_nglp);

  double left_bound = _left_bound;

  Eigen::VectorXd result(_nglp * _intervals.size());

  for (long interval = 0; interval < _intervals.size(); interval++) {

    double local_left_bound = left_bound + interval * _intervals[interval];

    result.segment(interval * _nglp, _nglp) =
        (glp.array() + 1.0) * _intervals[interval] / 2.0 + local_left_bound;
  }
  return result;
}

Eigen::VectorXd
legendre_gauss_lobatto_weights(std::size_t _nglp,
                               const Eigen::VectorXd &_intervals) {

  Eigen::VectorXd glw = legendre_gauss_lobatto_weights(_nglp);

  Eigen::VectorXd result(_nglp * _intervals.size());

  for (long interval = 0; interval < _intervals.size(); interval++) {

    result.segment(interval * _nglp, _nglp) = glw * _intervals[interval] / 2.0;
  }
  return result;
}
} // namespace collocation
} // namespace gsplines
