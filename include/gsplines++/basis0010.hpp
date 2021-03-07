#ifndef BASIS0010_H
#define BASIS0010_H
//#include <eigen3/Eigen/Core>
#include <gsplines++/basis.hpp>

namespace gsplines {
namespace basis {
class Basis0010 : public Basis {
private:
  // Eigen::Matrix<double, 6, 6> dmat_;

  Basis0010 &operator=(const Basis0010 &);
  Basis0010(const Basis0010 &that);

public:
  Basis0010();
  ~Basis0010();
  void eval_on_window(double _s, double _tau, double _buff[]);
  void eval_derivative_on_window(double _s, double _tau, unsigned int _deg,
                                 double _buff[]);

  void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                         unsigned int _deg, double _buff[]);
};
} // namespace basis
} // namespace gsplines
#endif /* BASIS0010_H */
