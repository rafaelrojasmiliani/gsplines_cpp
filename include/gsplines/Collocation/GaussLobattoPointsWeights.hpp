#ifndef GAUSSLOBATTO
#define GAUSSLOBATTO

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <vector>

namespace gsplines {
namespace collocation {

std::pair<Eigen::VectorXd, Eigen::VectorXd>
legendre_gauss_lobatto_points_and_weights(std::size_t _n);

Eigen::VectorXd legendre_gauss_lobatto_points(std::size_t _n);

Eigen::VectorXd legendre_gauss_lobatto_weights(std::size_t _n);

Eigen::VectorXd legendre_gauss_lobatto_points(std::tuple<double, double>,
                                              std::size_t _nglp,
                                              std::size_t _ninter);

Eigen::VectorXd legendre_gauss_lobatto_weights(std::tuple<double, double>,
                                               std::size_t _nglp,
                                               std::size_t _ninter);
std::tuple<double, double, double> q_and_evaluation(double _t, std::size_t _n);

} // namespace collocation
} // namespace gsplines

#endif /* ifndef GAUSSLOBATTO */
