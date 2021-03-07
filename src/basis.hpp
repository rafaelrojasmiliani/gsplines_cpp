#ifndef BASIS_H
#define BASIS_H

class Basis {
private:
  unsigned int dim_;

public:
  Basis(unsigned int _dim) : dim_(_dim) {}
  virtual void eval_on_window(double _s, double _tau, double _buff[]) = 0;
  virtual void eval_derivative_on_window(double _s, double _tau,
                                         unsigned int _deg, double _buff[]) = 0;

  virtual void eval_derivative_wrt_tau_on_window(double _s, double _tau,
                                                 unsigned int _deg,
                                                 double _buff[]) = 0;
  virtual ~Basis(){};
  double get_dim() const { return dim_; }
};

#endif /* BASIS_H */
