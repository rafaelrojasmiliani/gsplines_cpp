

#include <gsplines++/BasisLegendre.hpp>
int main() {

  gsplines::basis::BasisLegendre basis(6);
  double buff[6];
  basis.eval_on_window(0, 1, buff);
  return 0;
}
