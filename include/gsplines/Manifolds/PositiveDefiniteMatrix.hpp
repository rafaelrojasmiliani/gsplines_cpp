#ifndef POSITIVEDEFINITEMATRIX_H
#define POSITIVEDEFINITEMATRIX_H

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
/**
 * PDM::Eleement element = PDM::random();
 *
 * PDM::Element el_2 = PDM::CanonicParametrization(vector);
 *
 * Eigen::MatrixXd Q = PDM::inclusionOnMatrices(element);
 *
 * Eigen::VectorXd x_element = PDM::CanonicChart(element);
 * */

namespace gsplines {
namespace manifolds {
namespace pdm {

class PDMatrix {
private:
  Eigen::MatrixXd value_;

public:
  PDMatrix(const Eigen::VectorXd &_x);

  const Eigen::MatrixXd &matrix_inclusion() const { return value_; };
  PDMatrix(const PDMatrix &_that);
  PDMatrix(PDMatrix &&_that);
  virtual ~PDMatrix() = default;
};

class CanonicChart {
public:
  PDMatrix operator()(const Eigen::VectorXd &_that);

  Eigen::VectorXd inverse(const PDMatrix &_that);

  const Eigen::MatrixXd &tangent_basis(PDMatrix &_that);
};

} // namespace pdm
} // namespace manifolds
} // namespace gsplines

#endif /* POSITIVEDEFINITEMATRIX_H */
