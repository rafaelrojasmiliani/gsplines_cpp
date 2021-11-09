
#include <iostream>
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

int main() {
  const double pi = std::acos(-1.0);

  MatrixXd A(MatrixXd::Random(30, 30));
  std::cout << "The matrix exponential of A is:\n" << A.exp() << "\n\n";
}
