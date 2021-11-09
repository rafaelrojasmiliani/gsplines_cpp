#include <cmath>
#include <gsplines/Manifolds/MatrixExponential.hpp>
#include <iostream>
#include <tuple>
#include <unsupported/Eigen/MatrixFunctions>
#include <utility>
#include <vector>

namespace gsplines {
namespace manifolds {

MatrixExponential::MatrixExponential(std::size_t _n)
    : A(_n, _n), E(_n, _n), A2(_n, _n), M2(_n, _n), A4(_n, _n), M4(_n, _n),
      A6(_n, _n), M6(_n, _n), A8(_n, _n), M8(_n, _n), W1(_n, _n), W2(_n, _n),
      Z1(_n, _n), Z2(_n, _n), W(_n, _n), U(_n, _n), V(_n, _n), Lw1(_n, _n),
      Lw2(_n, _n), Lz1(_n, _n), Lz2(_n, _n), Lw(_n, _n), Lu(_n, _n), Lv(_n, _n),
      ident(Eigen::MatrixXd::Identity(_n, _n)), R(_n, _n), L(_n, _n) {}

void MatrixExponential::diff_pade3(Eigen::Ref<const Eigen::MatrixXd> _mat,
                                   Eigen::Ref<const Eigen::MatrixXd> _direction,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_4) {

  //  return _buffer_1 == U, _buffer_2 == V, _buffer_3 == Lu, _buffer_4 == Lv
  // b = (120., 60., 12., 1.)
  static std::vector<double> b{120., 60., 12., 1.};
  // A2 = A.dot(A)
  A2.noalias() = _mat * _mat;
  // M2 = np.dot(A, E) + np.dot(E, A)
  M2.noalias() = _mat * _direction + _direction * _mat;
  //  U = A.dot(b[3]*A2 + b[1]*ident)
  _buffer_1.noalias() = _mat * (b[3] * A2 + b[1] * ident);
  // V = b[2]*A2 + b[0]*ident
  _buffer_2.noalias() = b[2] * A2 + b[0] * ident;
  // Lu = A.dot(b[3]*M2) + E.dot(b[3]*A2 + b[1]*ident)
  _buffer_3.noalias() =
      _mat * (b[3] * M2) + _direction * (b[3] * A2 + b[1] * ident);
  // Lv = b[2]*M2
  _buffer_4.noalias() = M2 * b[2];
}

void MatrixExponential::diff_pade5(Eigen::Ref<const Eigen::MatrixXd> _mat,
                                   Eigen::Ref<const Eigen::MatrixXd> _direction,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_4) {

  //  return _buffer_1 == U, _buffer_2 == V, _buffer_3 == Lu, _buffer_4 == Lv
  //    b = (30240., 15120., 3360., 420., 30., 1.)
  static std::vector<double> b{30240., 15120., 3360., 420., 30., 1.};

  const Eigen::MatrixXd &A = _mat;
  const Eigen::MatrixXd &E = _direction;
  //    A2 = A.dot(A)
  A2.noalias() = A * A;
  //    M2 = np.dot(A, E) + np.dot(E, A)
  M2.noalias() = A * E + E * A;
  //    A4 = np.dot(A2, A2)
  A4.noalias() = A2 * A2;
  //    M4 = np.dot(A2, M2) + np.dot(M2, A2)
  M4.noalias() = A2 * M2 + M2 * A2;
  //    U = A.dot(           b[5]*A4   + b[3] * A2 + b[1] * ident)
  _buffer_1.noalias() = A * (b[5] * A4 + b[3] * A2 + b[1] * ident);
  //    V =             b[4] * A4 + b[2] * A2 + b[0] * ident
  _buffer_2.noalias() = b[4] * A4 + b[2] * A2 + b[0] * ident;
  //    Lu = ( A.dot(b[5]*M4 + b[3]*M2) +
  //            E.dot(b[5]*A4 + b[3]*A2 + b[1]*ident))
  _buffer_3.noalias() = (A * (b[5] * M4 + b[3] * M2) +
                         E * (b[5] * A4 + b[3] * A2 + b[1] * ident));
  //    Lv = b[4]*M4 + b[2]*M2
  _buffer_4.noalias() = b[4] * M4 + b[2] * M2;
  //    return U, V, Lu, Lv
}
void MatrixExponential::diff_pade7(Eigen::Ref<const Eigen::MatrixXd> _mat,
                                   Eigen::Ref<const Eigen::MatrixXd> _direction,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_4) {

  const Eigen::MatrixXd &A = _mat;
  const Eigen::MatrixXd &E = _direction;
  Eigen::Ref<Eigen::MatrixXd> U = _buffer_1;
  Eigen::Ref<Eigen::MatrixXd> V = _buffer_2;
  Eigen::Ref<Eigen::MatrixXd> Lu = _buffer_3;
  Eigen::Ref<Eigen::MatrixXd> Lv = _buffer_4;
  //    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
  static std::vector<double> b{17297280., 8648640., 1995840., 277200.,
                               25200.,    1512.,    56.,      1.};
  //    A2 = A.dot(A)
  A2 = A * A;
  //    M2 = np.dot(A, E) + np.dot(E, A)
  M2.noalias() = A * E + E * A;
  //    A4 = np.dot(A2, A2)
  A4.noalias() = A2 * A2;
  //    M4 = np.dot(A2, M2) + np.dot(M2, A2)
  M4.noalias() = A2 * M2 + M2 * A2;
  //    M6 = np.dot(A4, M2) + np.dot(M4, A2)
  M6.noalias() = A4 * M2 + M4 * A2;
  //    A6 = np.dot(A2, A4);
  A6.noalias() = A2 * A4;
  //    U = A.dot(   b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * ident)
  U.noalias() = A * (b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * ident);
  //    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
  V.noalias() = b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * ident;
  //    Lu = (A.dot(b[7]*M6 + b[5]*M4 + b[3]*M2) +
  //            E.dot(b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident))
  Lu.noalias() = (A * (b[7] * M6 + b[5] * M4 + b[3] * M2) +
                  E * (b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * ident));
  //    Lv = b[6]*M6 + b[4]*M4 + b[2]*M2
  Lv.noalias() = b[6] * M6 + b[4] * M4 + b[2] * M2;
  //    return U, V, Lu, Lv
}
void MatrixExponential::diff_pade9(Eigen::Ref<const Eigen::MatrixXd> _mat,
                                   Eigen::Ref<const Eigen::MatrixXd> _direction,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                                   Eigen::Ref<Eigen::MatrixXd> _buffer_4) {

  const Eigen::MatrixXd &A = _mat;
  const Eigen::MatrixXd &E = _direction;
  Eigen::Ref<Eigen::MatrixXd> U = _buffer_1;
  Eigen::Ref<Eigen::MatrixXd> V = _buffer_2;
  Eigen::Ref<Eigen::MatrixXd> Lu = _buffer_3;
  Eigen::Ref<Eigen::MatrixXd> Lv = _buffer_4;
  //    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
  //            2162160., 110880., 3960., 90., 1.)
  static std::vector<double> b{
      17643225600., 8821612800., 2075673600., 302702400., 30270240.,
      2162160.,     110880.,     3960.,       90.,        1.};
  //    A2 = A.dot(A)
  A2 = A * A;
  //    M2 = np.dot(A, E) + np.dot(E, A)
  M2.noalias() = A * E + E * A;
  //    A4 = np.dot(A2, A2)
  A4.noalias() = A2 * A2;
  //    M4 = np.dot(A2, M2) + np.dot(M2, A2)
  M4.noalias() = A2 * M2 + M2 * A2;
  //    M6 = np.dot(A4, M2) + np.dot(M4, A2)
  M6.noalias() = A4 * M2 + M4 * A2;
  //    A6 = np.dot(A2, A4);
  A6.noalias() = A2 * A4;
  //    A8 = np.dot(A4, A4)
  A8.noalias() = A4 * A4;
  //    M8 = np.dot(A4, M4) + np.dot(M4, A4)
  M8.noalias() = A4 * M4 + M4 * A4;
  //    U = A.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
  U.noalias() = A * (b[9] * (A4 * A4) + b[7] * (A2 * A4) + b[5] * A4 +
                     b[3] * A2 + b[1] * ident);
  //    V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
  V.noalias() = b[8] * (A4 * A4) + b[6] * (A2 * A4) + b[4] * A4 + b[2] * A2 +
                b[0] * ident;
  //    Lu = (A.dot(b[9]*M8 + b[7]*M6 + b[5]*M4 + b[3]*M2) +
  //            E.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident))
  Lu.noalias() =
      (A * (b[9] * M8 + b[7] * M6 + b[5] * M4 + b[3] * M2) +
       E * (b[9] * A8 + b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * ident));
  //        Lv = b[8] * M8 + b[6] * M6 + b[4] * M4 + b[2] * M2
  Lv.noalias() = b[8] * M8 + b[6] * M6 + b[4] * M4 + b[2] * M2;
}

void MatrixExponential::compute_A_M(const Eigen::MatrixXd &_mat,
                                    const Eigen::MatrixXd &_direction) {

  const Eigen::MatrixXd &A = _mat;
  const Eigen::MatrixXd &E = _direction;
  A2 = A * A;
  //    M2 = np.dot(A, E) + np.dot(E, A)
  M2.noalias() = A * E + E * A;
  //    A4 = np.dot(A2, A2)
  A4.noalias() = A2 * A2;
  //    M4 = np.dot(A2, M2) + np.dot(M2, A2)
  M4.noalias() = A2 * M2 + M2 * A2;
  //    M6 = np.dot(A4, M2) + np.dot(M4, A2)
  M6.noalias() = A4 * M2 + M4 * A2;
  //    A6 = np.dot(A2, A4);
  A6.noalias() = A2 * A4;
  //    A8 = np.dot(A4, A4)
  A8.noalias() = A4 * A4;
  //    M8 = np.dot(A4, M4) + np.dot(M4, A4)
  M8.noalias() = A4 * M4 + M4 * A4;
}

Eigen::MatrixXd
MatrixExponential::operator()(Eigen::Ref<const Eigen::MatrixXd> _mat) {
  return std::move(_mat.exp());
}
Eigen::MatrixXd MatrixExponential::gataux_derivative(
    Eigen::Ref<const Eigen::MatrixXd> _mat,
    Eigen::Ref<const Eigen::MatrixXd> _direction) {

  static std::vector<double> ell_table_61{
      0.0,     2.11e-8, 3.56e-4, 1.08e-2, 6.49e-2, 2.00e-1, 4.37e-1,
      7.83e-1, 1.23e0,  1.78e0,  2.42e0,  3.13e0,  3.90e0,  4.74e0,
      5.63e0,  6.56e0,  7.52e0,  8.53e0,  9.56e0,  1.06e1,  1.17e1};

  std::size_t number_of_final_multiplications = 0;

  bool flag = true;

  // n = A.shape[0]
  std::size_t n = _mat.rows();

  //    A_norm_1 = scipy.linalg.norm(A, 1)
  double A_norm_1 = _mat.cwiseAbs().colwise().sum().maxCoeff();

  auto _diff_pade3 = [this](Eigen::Ref<const Eigen::MatrixXd> _mat,
                            Eigen::Ref<const Eigen::MatrixXd> _direction,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_4) {
    this->diff_pade3(_mat, _direction, _buffer_1, _buffer_2, _buffer_3,
                     _buffer_4);
  };

  auto _diff_pade5 = [this](Eigen::Ref<const Eigen::MatrixXd> _mat,
                            Eigen::Ref<const Eigen::MatrixXd> _direction,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_4) {
    this->diff_pade5(_mat, _direction, _buffer_1, _buffer_2, _buffer_3,
                     _buffer_4);
  };

  auto _diff_pade7 = [this](Eigen::Ref<const Eigen::MatrixXd> _mat,
                            Eigen::Ref<const Eigen::MatrixXd> _direction,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_4) {
    this->diff_pade7(_mat, _direction, _buffer_1, _buffer_2, _buffer_3,
                     _buffer_4);
  };
  auto _diff_pade9 = [this](Eigen::Ref<const Eigen::MatrixXd> _mat,
                            Eigen::Ref<const Eigen::MatrixXd> _direction,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_1,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_2,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_3,
                            Eigen::Ref<Eigen::MatrixXd> _buffer_4) {
    this->diff_pade9(_mat, _direction, _buffer_1, _buffer_2, _buffer_3,
                     _buffer_4);
  };

  typedef std::function<void(
      const Eigen::MatrixXd &_mat, const Eigen::MatrixXd &_direction,
      Eigen::MatrixXd &_buffer_1, Eigen::MatrixXd &_buffer_2,
      Eigen::MatrixXd &_buffer_3, Eigen::MatrixXd &_buffer_4)>
      new_type;

  //    m_pade_pairs = (
  //            (3, _diff_pade3),
  //            (5, _diff_pade5),
  //            (7, _diff_pade7),
  //            (9, _diff_pade9))
  static std::vector<std::pair<std::size_t, new_type>> m_pade_pairs{
      {3, _diff_pade3}, {5, _diff_pade5}, {7, _diff_pade7}, {9, _diff_pade9}};

  // for m, pade in m_pade_pairs:
  //     if A_norm_1 <= ell_table_61[m]:
  //         U, V, Lu, Lv = pade(A, E, ident)
  //         s = 0
  //         break
  for (auto &it : m_pade_pairs) {
    if (A_norm_1 <= ell_table_61[it.first]) {
      printf("pade small deg = %zu\n", it.first);
      it.second(_mat, _direction, U, V, Lu, Lv);
      flag = false;
      break;
    }
  }

  if (flag) {
    printf("pade large\n");
    // s = max(0, int(np.ceil(np.log2(A_norm_1 / ell_table_61[13]))))
    std::size_t s = std::max(
        0, static_cast<int>(std::ceil(std::log2(A_norm_1 / ell_table_61[13]))));
    number_of_final_multiplications = s;
    // s = max(0, int(np.ceil(np.log2(A_norm_1 / ell_table_61[13]))))
    // A = A * 2.0**-s
    A = _mat * std::pow(2.0, -static_cast<double>(s));
    // E = E * 2.0**-s
    E = _direction * std::pow(2.0, -static_cast<double>(s));
    //# pade order 13
    // A2 = np.dot(A, A)
    A2 = A * A;
    // M2 = np.dot(A, E) + np.dot(E, A)
    M2 = A * E + E * A;
    // A4 = np.dot(A2, A2)
    A4 = A2 * A2;
    // M4 = np.dot(A2, M2) + np.dot(M2, A2)
    M4 = A2 * M2 + M2 * A2;
    // A6 = np.dot(A2, A4)
    A6 = A2 * A4;
    // M6 = np.dot(A4, M2) + np.dot(M4, A2)
    M6 = A4 * M2 + M2 * A4;
    // b = (64764752532480000., 32382376266240000., 7771770303897600.,
    // 1187353796428800., 129060195264000., 10559470521600.,
    // 670442572800., 33522128640., 1323241920., 40840800., 960960.,
    // 16380., 182., 1.)
    static std::vector<double> b{64764752532480000.,
                                 32382376266240000.,
                                 7771770303897600.,
                                 1187353796428800.,
                                 129060195264000.,
                                 10559470521600.,
                                 670442572800.,
                                 33522128640.,
                                 1323241920.,
                                 40840800.,
                                 960960.,
                                 16380.,
                                 182.,
                                 1.};
    // W1=b[13]* A6 + b[11] * A4 + b[9] * A2
    W1 = b[13] * A6 + b[11] * A4 + b[9] * A2;
    // W2=b[7]* A6 + b[5] * A4 + b[3] * A2 + b[1] * ident
    W2 = b[7] * A6 + b[5] * A4 + b[3] * A2 + b[1] * ident;
    // Z1=b[12]* A6 + b[10] * A4 + b[8] * A2;
    Z1 = b[12] * A6 + b[10] * A4 + b[8] * A2;
    // Z2=b[6]* A6 + b[4] * A4 + b[2] * A2 + b[0] * ident;
    Z2 = b[6] * A6 + b[4] * A4 + b[2] * A2 + b[0] * ident;
    // W = np.dot(A6, W1) + W2;
    W = A6 * W1 + W2;
    // U = np.dot(A, W);
    U = A * W;
    // V = np.dot(A6, Z1) + Z2;
    V = A6 * Z1 + Z2;
    // Lw1=b[13]* M6 + b[11] * M4 + b[9] * M2;
    Lw1 = b[13] * M6 + b[11] * M4 + b[9] * M2;
    // Lw2=b[7]* M6 + b[5] * M4 + b[3] * M2;
    Lw2 = b[7] * M6 + b[5] * M4 + b[3] * M2;
    // Lz1=b[12]* M6 + b[10] * M4 + b[8] * M2;
    Lz1 = b[12] * M6 + b[10] * M4 + b[8] * M2;
    // Lz2=b[6]* M6 + b[4] * M4 + b[2] * M2;
    Lz2 = b[6] * M6 + b[4] * M4 + b[2] * M2;
    // Lw = np.dot(A6, Lw1) + np.dot(M6, W1) + Lw2;
    Lw = A6 * Lw1 + M6 * W1 + Lw2;
    // Lu = np.dot(A, Lw) + np.dot(E, W);
    Lu = A * Lw + E * W;
    // Lv = np.dot(A6, Lz1) + np.dot(M6, Z1) + Lz2;
    Lv = A6 * Lz1 + M6 * Z1 + Lz2;
    /*
     * */
  }
  // lu_piv = scipy.linalg.lu_factor(-U + V)

  Eigen::PartialPivLU<Eigen::MatrixXd> lu_piv = (-U + V).lu();
  // R = scipy.linalg.lu_solve(lu_piv, U + V)
  R = lu_piv.solve(U + V);
  // L = scipy.linalg.lu_solve(lu_piv, Lu + Lv + np.dot((Lu - Lv), R))
  L = lu_piv.solve(Lu + Lv + (Lu - Lv) * R);

  //  for k in range(s):
  //      L = np.dot(R, L) + np.dot(L, R)
  //      R = np.dot(R, R)
  for (std::size_t k = 0; k < number_of_final_multiplications; k++) {
    L = R * L + L * R;
    R = R * R;
  }
  return L;
}
} // namespace manifolds
} // namespace gsplines
