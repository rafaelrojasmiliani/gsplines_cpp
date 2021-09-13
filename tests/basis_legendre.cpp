

#include <eigen3/Eigen/Core>
#include <gsplines/Basis/BasisLegendre.hpp>
int main() {

  gsplines::basis::BasisLegendre basis(6);
  Eigen::VectorXd buff(6);
  basis.eval_on_window(0, 1, buff);
  return 0;
}
