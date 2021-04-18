

#include <gsplines++/BasisLegendre.hpp>
#include <gsplines++/interpolator.hpp>
int main() {

  gsplines::basis::BasisLegendre basis(6);
  Eigen::VectorXd buff(6);
  basis.eval_on_window(0, 1, buff);
  gsplines::Interpolator inter(6, 1, basis);
  Eigen::VectorXd tau(1);
  tau(0) = 1.0;
  inter.fill_interpolating_matrix(tau);
}
