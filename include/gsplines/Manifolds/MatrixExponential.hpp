#ifndef MATRIXEXPONENTIAL_H
#define MATRIXEXPONENTIAL_H
#include <Eigen/Core>
#include <unsupported/Eigen/MatrixFunctions>

namespace gsplines {
namespace manifolds {

class MatrixExponential {
private:
  mutable Eigen::MatrixXd A;
  mutable Eigen::MatrixXd E;
  mutable Eigen::MatrixXd A2;
  mutable Eigen::MatrixXd M2;
  mutable Eigen::MatrixXd A4;
  mutable Eigen::MatrixXd M4;
  mutable Eigen::MatrixXd A6;
  mutable Eigen::MatrixXd M6;
  mutable Eigen::MatrixXd A8;
  mutable Eigen::MatrixXd M8;
  mutable Eigen::MatrixXd W1;
  mutable Eigen::MatrixXd W2;
  mutable Eigen::MatrixXd Z1;
  mutable Eigen::MatrixXd Z2;
  mutable Eigen::MatrixXd W;
  mutable Eigen::MatrixXd U;
  mutable Eigen::MatrixXd V;
  mutable Eigen::MatrixXd Lw1;
  mutable Eigen::MatrixXd Lw2;
  mutable Eigen::MatrixXd Lz1;
  mutable Eigen::MatrixXd Lz2;
  mutable Eigen::MatrixXd Lw;
  mutable Eigen::MatrixXd Lu;
  mutable Eigen::MatrixXd Lv;
  mutable Eigen::MatrixXd R;
  mutable Eigen::MatrixXd L;
  const Eigen::MatrixXd ident;

  void compute_A_M(const Eigen::MatrixXd &_mat,
                   const Eigen::MatrixXd &_direction);

public:
  MatrixExponential(std::size_t _n);

  void diff_pade3(Eigen::Ref<const Eigen::MatrixXd> _mat,
                  Eigen::Ref<const Eigen::MatrixXd> _direction,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_4);

  void diff_pade5(Eigen::Ref<const Eigen::MatrixXd> _mat,
                  Eigen::Ref<const Eigen::MatrixXd> _direction,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_4);

  void diff_pade7(Eigen::Ref<const Eigen::MatrixXd> _mat,
                  Eigen::Ref<const Eigen::MatrixXd> _direction,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_4);

  void diff_pade9(Eigen::Ref<const Eigen::MatrixXd> _mat,
                  Eigen::Ref<const Eigen::MatrixXd> _direction,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                  Eigen::Ref<Eigen::MatrixXd> _buffer_4);

  Eigen::MatrixXd operator()(Eigen::Ref<const Eigen::MatrixXd> _mat);
  Eigen::MatrixXd
  gataux_derivative(Eigen::Ref<const Eigen::MatrixXd> _mat,
                    Eigen::Ref<const Eigen::MatrixXd> _direction);

  virtual ~MatrixExponential() = default;
};

} // namespace manifolds
} // namespace gsplines

#endif /* MATRIXEXPONENTIAL_H */
