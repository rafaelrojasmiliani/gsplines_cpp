#ifndef HANKELMATRIX_H
#define HANKELMATRIX_H

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>
namespace gsplines {
namespace manifolds {
namespace hankel {

class HankelMatrix {
private:
public:
  HankelMatrix();
  virtual ~HankelMatrix() = default;
};

class CanonicChart {
public:
  CanonicChart();
  virtual ~CanonicChart() = default;
};

HankelMatrix frobenius_orthogonal_projector(Eigen::MatrixXd &_that);

} // namespace hankel
} // namespace manifolds
} // namespace gsplines

#endif /* HANKELMATRIX_H */
