#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace gsplines {
namespace manifolds {

namespace symmetric_matrices {

class SymmetricMatrix {
private:
  Eigen::MatrixXd value_;

public:
  SymmetricMatrix(const Eigen::VectorXd &_x);
  virtual ~SymmetricMatrix() = default;
};

} // namespace symmetric_matrices
} // namespace manifolds
} // namespace gsplines
#endif /* SYMMETRICMATRIX_H */
