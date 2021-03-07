#ifndef BASIS0010_H
#define BASIS0010_H
#include <basis/basis.h>
class Basis0010 : public Basis {
private:
public:
  Basis0010();
  Basis0010(const Basis0010 : public Basis &that);
  Basis0010 &operator=(const Basis0010 : public Basis &);
  virtual ~Basis0010();
  void eval_on_window(double _s, double _tau, double _buff[]) {
    _buff[0] = 1.0;
    _buff[1] = _s;
    for (int i = 1; i < 5; i++)
      _buff[i + 1] =
          1.0 / (i + 1.0) * ((2.0 * i + 1.0) * _s * _buff[i] - i * _buff[i - 1])
  }
  virtual void eval_derivative_on_window(double _s, double _tau,
                                         unsigned int _deg, double *_buff);

  virtual void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                                 unsigned int _deg,
                                                 double *_buff);
  virtual ~Basis();
};

#endif /* BASIS0010_H */
