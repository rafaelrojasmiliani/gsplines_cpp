
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis0101.hpp>

namespace gsplines::basis {

void compute_Q_block(double _tau, double _alpha,
                     Eigen::Ref<Eigen::MatrixXd> _res);

void compute_Qd1_block(double _tau, double _alpha,
                       Eigen::Ref<Eigen::MatrixXd> _res);
void compute_Qd2_block(double _tau, double _alpha,
                       Eigen::Ref<Eigen::MatrixXd> _res);
void compute_Qd3_block(double _tau, double _alpha,
                       Eigen::Ref<Eigen::MatrixXd> _res);

void compute_Qd1_dtau_block(double _tau, double _alpha,
                            Eigen::Ref<Eigen::MatrixXd> _res);
void compute_Qd3_dtau_block(double _tau, double _alpha,
                            Eigen::Ref<Eigen::MatrixXd> _res);

Basis0101::Basis0101(double _alpha)
    : Basis(6, "basis0101", [_alpha]() {
        Eigen::VectorXd res(1);
        res(0) = _alpha;
        return res;
      }()) {}

std::shared_ptr<Basis0101> Basis0101::get(double k) {
  return std::shared_ptr<Basis0101>(new Basis0101(k));
}

double Basis0101::get_alpha() const { return this->get_parameters()(0); }

void Basis0101::eval_on_window(
    double _s, double _tau,
    Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const {
  double alpha = this->get_parameters()(0);
  double k = std::sqrt(2) / 4.0 * std::pow(alpha, 0.25) /
             std::pow((1.0 - alpha), 0.25);
  double p = _tau * k * _s;
  double expp = std::exp(p);
  double cosp = std::cos(p);
  double sinp = std::sin(p);
  _buff[0] = expp * cosp;
  _buff[1] = expp * sinp;
  _buff[2] = cosp / expp;
  _buff[3] = sinp / expp;
  _buff[4] = p;
  _buff[5] = 1.0;
}
void Basis0101::eval_derivative_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const {
  double alpha = this->get_parameters()(0);
  double k = std::sqrt(2) / 4.0 * std::pow(alpha, 0.25) /
             std::pow((1.0 - alpha), 0.25);
  double p = _tau * k * _s;
  double expp = std::exp(p);
  double cosp = std::cos(p);
  double sinp = std::sin(p);
  _buff[0] = expp * cosp;
  _buff[1] = expp * sinp;
  _buff[2] = cosp / expp;
  _buff[3] = sinp / expp;
  _buff[4] = p;
  _buff[5] = 1.0;

  for (unsigned int i = 1; i <= _deg; i++) {
    double v0 = _buff[0];
    double v1 = _buff[1];
    _buff[0] = v0 - v1;
    _buff[1] = v0 + v1;
    v0 = _buff[2];
    v1 = _buff[3];
    _buff[2] = -v0 - v1;
    _buff[3] = v0 - v1;
    _buff[4] = _buff[5];
    _buff[5] = 0;
    _buff *= k * 2;
  }
}

void Basis0101::eval_derivative_wrt_tau_on_window(
    double _s, double _tau, unsigned int _deg,
    Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> _buff) const {
  double alpha = this->get_parameters()(0);
  double k = std::sqrt(2) / 4.0 * std::pow(alpha, 0.25) /
             std::pow((1.0 - alpha), 0.25);

  this->eval_derivative_on_window(_s, _tau, _deg, _buff);
  double v0 = _buff[0];
  double v1 = _buff[1];
  _buff[0] = v0 - v1;
  _buff[1] = v0 + v1;
  v0 = _buff[2];
  v1 = _buff[3];
  _buff[2] = -v0 - v1;
  _buff[3] = v0 - v1;
  _buff[4] = _buff[5];
  _buff[5] = 0;
  _buff *= _s * k;
}

void Basis0101::add_derivative_matrix_deriv_wrt_tau(
    double tau, std::size_t _deg, Eigen::Ref<Eigen::MatrixXd> _mat) {
  qBuffer_.setZero();
  double alpha = this->get_parameters()(0);
  switch (_deg) {
    case 1:
      compute_Qd1_dtau_block(tau, alpha, qBuffer_);
      break;
    case 3:
      compute_Qd3_dtau_block(tau, alpha, qBuffer_);
      break;
    default:
      throw std::invalid_argument(
          "For basis 1010 this derivative wrt tau is not implemented");
  }
  _mat.noalias() += qBuffer_;
}

void Basis0101::add_derivative_matrix(double tau, std::size_t _deg,
                                      Eigen::Ref<Eigen::MatrixXd> _mat) {
  qBuffer_.setZero();
  double alpha = this->get_parameters()(0);
  switch (_deg) {
    case 0:
      compute_Q_block(tau, alpha, qBuffer_);
      break;
    case 1:
      compute_Qd1_block(tau, alpha, qBuffer_);
      break;
    case 2:
      compute_Qd2_block(tau, alpha, qBuffer_);
      break;
    case 3:
      compute_Qd3_block(tau, alpha, qBuffer_);
      break;
    default:
      throw std::invalid_argument(
          "This derivative matrix has not been implemented");
  }
  _mat.noalias() += qBuffer_;
}

std::unique_ptr<Basis> Basis0101::clone() const {
  return std::unique_ptr<Basis>(new Basis0101(*this));
}
std::unique_ptr<Basis> Basis0101::move_clone() {
  return std::unique_ptr<Basis>(new Basis0101(std::move(*this)));
}

Eigen::MatrixXd Basis0101::derivative_matrix_impl(std::size_t _deg) const {
  static constexpr std::array<double, 36> d1matrix = {
      1, -1, 0, 0,  0, 0, 1, 1, 0, 0, 0, 0, 0, 0, -1, -1, 0, 0,
      0, 0,  1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,  0,  0, 0};

  double alpha = this->get_parameters()(0);
  if (_deg == 0) {
    return Eigen::MatrixXd::Identity(6, 6);
  }

  double k = std::sqrt(2.0) / 4.0 * std::pow(alpha, 0.25) /
             std::pow((1.0 - alpha), 0.25);
  double tauk = k * 2.0;

  Eigen::MatrixXd Q(Eigen::Map<const Eigen::MatrixXd>(d1matrix.data(), 6, 6));

  Q *= tauk;

  for (std::size_t i = 2; i <= _deg; i++) {
    Q *= Q;
  }

  return Q;
}

using namespace std;

void compute_Qd3_block(double _tau, double _alpha,
                       Eigen::Ref<Eigen::MatrixXd> _res) {
  double k =
      0.35355339059327379 * pow(_alpha, 0.25) * pow(1.0 / (1.0 - _alpha), 0.25);
  _res(0, 0) =
      96.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (96.0 * pow(k, 5) * pow(sin(_tau * k), 2) -
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       32.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 1) =
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) -
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      (-32.0 * pow(k, 5) * pow(sin(_tau * k), 2) -
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       32.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 2) = -256.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k);
  _res(0, 3) = -256.0 * _tau * pow(k, 6) * pow(sin(_tau * k), 2) -
               256.0 * _tau * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(0, 4) = 0;
  _res(0, 5) = 0;
  _res(1, 0) =
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) -
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      (-32.0 * pow(k, 5) * pow(sin(_tau * k), 2) -
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       32.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(1, 1) =
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      96.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (32.0 * pow(k, 5) * pow(sin(_tau * k), 2) +
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       96.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(1, 2) = 256.0 * _tau * pow(k, 6) * pow(sin(_tau * k), 2) +
               256.0 * _tau * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(1, 3) = 256.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k);
  _res(1, 4) = 0;
  _res(1, 5) = 0;
  _res(2, 0) = -256.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k);
  _res(2, 1) = 256.0 * _tau * pow(k, 6) * pow(sin(_tau * k), 2) +
               256.0 * _tau * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(2, 2) =
      96.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (96.0 * pow(k, 5) * pow(sin(_tau * k), 2) -
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       32.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 3) =
      -32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (-32.0 * pow(k, 5) * pow(sin(_tau * k), 2) -
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       32.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 4) = 0;
  _res(2, 5) = 0;
  _res(3, 0) = -256.0 * _tau * pow(k, 6) * pow(sin(_tau * k), 2) -
               256.0 * _tau * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(3, 1) = 256.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k);
  _res(3, 2) =
      -32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (-32.0 * pow(k, 5) * pow(sin(_tau * k), 2) -
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       32.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(3, 3) =
      32.0 * pow(k, 5) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      64.0 * pow(k, 5) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      96.0 * pow(k, 5) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (32.0 * pow(k, 5) * pow(sin(_tau * k), 2) +
       64.0 * pow(k, 5) * sin(_tau * k) * cos(_tau * k) +
       96.0 * pow(k, 5) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(3, 4) = 0;
  _res(3, 5) = 0;
  _res(4, 0) = 0;
  _res(4, 1) = 0;
  _res(4, 2) = 0;
  _res(4, 3) = 0;
  _res(4, 4) = 0;
  _res(4, 5) = 0;
  _res(5, 0) = 0;
  _res(5, 1) = 0;
  _res(5, 2) = 0;
  _res(5, 3) = 0;
  _res(5, 4) = 0;
  _res(5, 5) = 0;
}

void compute_Qd3_dtau_block(double _tau, double _alpha,
                            Eigen::Ref<Eigen::MatrixXd> _res) {
  double k =
      0.35355339059327379 * pow(_alpha, 0.25) * pow(1.0 / (1.0 - _alpha), 0.25);
  _res(0, 0) =
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      256.0 * pow(k, 6) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (128.0 * pow(k, 6) * pow(sin(_tau * k), 2) -
           256.0 * pow(k, 6) * sin(_tau * k) * cos(_tau * k) +
           128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 1) = 128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
               128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               1.0 *
                   (-128.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
                    128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(0, 2) = 256.0 * pow(k, 6) * pow(sin(_tau * k), 2) -
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(0, 3) = -256.0 * pow(k, 6) * pow(sin(_tau * k), 2) -
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(0, 4) = 0;
  _res(0, 5) = 0;
  _res(1, 0) = 128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
               128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               1.0 *
                   (-128.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
                    128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(1, 1) =
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      256.0 * pow(k, 6) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (128.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
           256.0 * pow(k, 6) * sin(_tau * k) * cos(_tau * k) +
           128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(1, 2) = 256.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(1, 3) = -256.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(1, 4) = 0;
  _res(1, 5) = 0;
  _res(2, 0) = 256.0 * pow(k, 6) * pow(sin(_tau * k), 2) -
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(2, 1) = 256.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(2, 2) =
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      256.0 * pow(k, 6) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (128.0 * pow(k, 6) * pow(sin(_tau * k), 2) -
           256.0 * pow(k, 6) * sin(_tau * k) * cos(_tau * k) +
           128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 3) = -128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
               1.0 *
                   (-128.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
                    128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(2, 4) = 0;
  _res(2, 5) = 0;
  _res(3, 0) = -256.0 * pow(k, 6) * pow(sin(_tau * k), 2) -
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(3, 1) = -256.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
               256.0 * pow(k, 6) * pow(cos(_tau * k), 2);
  _res(3, 2) = -128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
               1.0 *
                   (-128.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
                    128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(3, 3) =
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      256.0 * pow(k, 6) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      128.0 * pow(k, 6) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (128.0 * pow(k, 6) * pow(sin(_tau * k), 2) +
           256.0 * pow(k, 6) * sin(_tau * k) * cos(_tau * k) +
           128.0 * pow(k, 6) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(3, 4) = 0;
  _res(3, 5) = 0;
  _res(4, 0) = 0;
  _res(4, 1) = 0;
  _res(4, 2) = 0;
  _res(4, 3) = 0;
  _res(4, 4) = 0;
  _res(4, 5) = 0;
  _res(5, 0) = 0;
  _res(5, 1) = 0;
  _res(5, 2) = 0;
  _res(5, 3) = 0;
  _res(5, 4) = 0;
  _res(5, 5) = 0;
}

void compute_Qd1_block(double _tau, double _alpha,
                       Eigen::Ref<Eigen::MatrixXd> _res) {
  double k =
      0.35355339059327379 * pow(_alpha, 0.25) * pow(1.0 / (1.0 - _alpha), 0.25);
  _res(0, 0) =
      0.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      1.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      0.5 *
          (k * pow(sin(_tau * k), 2) + 2.0 * k * sin(_tau * k) * cos(_tau * k) +
           3.0 * k * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 1) = -0.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
               0.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               0.5 *
                   (-k * pow(sin(_tau * k), 2) -
                    2.0 * k * sin(_tau * k) * cos(_tau * k) +
                    k * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(0, 2) = -4.0 * k * sin(_tau * k) * cos(_tau * k);
  _res(0, 3) = 4.0 * _tau * pow(k, 2) * pow(sin(_tau * k), 2) +
               4.0 * _tau * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(0, 4) = 2.0 * k * exp(_tau * k) * cos(_tau * k) -
               2.0 * k * exp(-_tau * k) * cos(_tau * k);
  _res(0, 5) = 0;
  _res(1, 0) = -0.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
               0.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               0.5 *
                   (-k * pow(sin(_tau * k), 2) -
                    2.0 * k * sin(_tau * k) * cos(_tau * k) +
                    k * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(1, 1) = 1.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
               0.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               0.5 *
                   (3.0 * k * pow(sin(_tau * k), 2) -
                    2.0 * k * sin(_tau * k) * cos(_tau * k) +
                    k * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(1, 2) = -4.0 * _tau * pow(k, 2) * pow(sin(_tau * k), 2) -
               4.0 * _tau * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(1, 3) = 4.0 * k * sin(_tau * k) * cos(_tau * k);
  _res(1, 4) = 2.0 * k * exp(_tau * k) * sin(_tau * k) +
               2.0 * k * exp(-_tau * k) * sin(_tau * k);
  _res(1, 5) = 0;
  _res(2, 0) = -4.0 * k * sin(_tau * k) * cos(_tau * k);
  _res(2, 1) = -4.0 * _tau * pow(k, 2) * pow(sin(_tau * k), 2) -
               4.0 * _tau * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(2, 2) =
      0.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      1.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      0.5 *
          (k * pow(sin(_tau * k), 2) + 2.0 * k * sin(_tau * k) * cos(_tau * k) +
           3.0 * k * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 3) = 0.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
               k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) -
               0.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
               0.5 *
                   (-k * pow(sin(_tau * k), 2) -
                    2.0 * k * sin(_tau * k) * cos(_tau * k) +
                    k * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(2, 4) = -2.0 * k * exp(_tau * k) * cos(_tau * k) +
               2.0 * k * exp(-_tau * k) * cos(_tau * k);
  _res(2, 5) = 0;
  _res(3, 0) = 4.0 * _tau * pow(k, 2) * pow(sin(_tau * k), 2) +
               4.0 * _tau * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(3, 1) = 4.0 * k * sin(_tau * k) * cos(_tau * k);
  _res(3, 2) = 0.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
               k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) -
               0.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
               0.5 *
                   (-k * pow(sin(_tau * k), 2) -
                    2.0 * k * sin(_tau * k) * cos(_tau * k) +
                    k * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(3, 3) = 1.5 * k * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               k * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
               0.5 * k * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               0.5 *
                   (3.0 * k * pow(sin(_tau * k), 2) -
                    2.0 * k * sin(_tau * k) * cos(_tau * k) +
                    k * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(3, 4) = 2.0 * k * exp(_tau * k) * sin(_tau * k) +
               2.0 * k * exp(-_tau * k) * sin(_tau * k);
  _res(3, 5) = 0;
  _res(4, 0) = 2.0 * k * exp(_tau * k) * cos(_tau * k) -
               2.0 * k * exp(-_tau * k) * cos(_tau * k);
  _res(4, 1) = 2.0 * k * exp(_tau * k) * sin(_tau * k) +
               2.0 * k * exp(-_tau * k) * sin(_tau * k);
  _res(4, 2) = -2.0 * k * exp(_tau * k) * cos(_tau * k) +
               2.0 * k * exp(-_tau * k) * cos(_tau * k);
  _res(4, 3) = 2.0 * k * exp(_tau * k) * sin(_tau * k) +
               2.0 * k * exp(-_tau * k) * sin(_tau * k);
  _res(4, 4) = 4.0 * _tau * pow(k, 2);
  _res(4, 5) = 0;
  _res(5, 0) = 0;
  _res(5, 1) = 0;
  _res(5, 2) = 0;
  _res(5, 3) = 0;
  _res(5, 4) = 0;
  _res(5, 5) = 0;
}

void compute_Qd1_dtau_block(double _tau, double _alpha,
                            Eigen::Ref<Eigen::MatrixXd> _res) {
  double k =
      0.35355339059327379 * pow(_alpha, 0.25) * pow(1.0 / (1.0 - _alpha), 0.25);
  _res(0, 0) =
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      4.0 * pow(k, 2) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (2.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
           4.0 * pow(k, 2) * sin(_tau * k) * cos(_tau * k) +
           2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 1) = -2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
               1.0 *
                   (-2.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
                    2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(0, 2) = 4.0 * pow(k, 2) * pow(sin(_tau * k), 2) -
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(0, 3) = 4.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(0, 4) =
      -2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) +
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) +
      1.0 *
          (2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(0, 5) = 0;
  _res(1, 0) = -2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
               2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
               1.0 *
                   (-2.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
                    2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(1, 1) =
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      4.0 * pow(k, 2) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (2.0 * pow(k, 2) * pow(sin(_tau * k), 2) -
           4.0 * pow(k, 2) * sin(_tau * k) * cos(_tau * k) +
           2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(1, 2) = -4.0 * pow(k, 2) * pow(sin(_tau * k), 2) -
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(1, 3) = -4.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(1, 4) =
      2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) +
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) +
      1.0 *
          (-2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(1, 5) = 0;
  _res(2, 0) = 4.0 * pow(k, 2) * pow(sin(_tau * k), 2) -
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(2, 1) = -4.0 * pow(k, 2) * pow(sin(_tau * k), 2) -
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(2, 2) =
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      4.0 * pow(k, 2) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (2.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
           4.0 * pow(k, 2) * sin(_tau * k) * cos(_tau * k) +
           2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 3) = 2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
               2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               1.0 *
                   (-2.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
                    2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(2, 4) =
      2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) -
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) -
      1.0 *
          (2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(2, 5) = 0;
  _res(3, 0) = 4.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(3, 1) = -4.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
               4.0 * pow(k, 2) * pow(cos(_tau * k), 2);
  _res(3, 2) = 2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
               2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
               1.0 *
                   (-2.0 * pow(k, 2) * pow(sin(_tau * k), 2) +
                    2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
                   exp(-2 * _tau * k);
  _res(3, 3) =
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      4.0 * pow(k, 2) * exp(2 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
      2.0 * pow(k, 2) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      1.0 *
          (2.0 * pow(k, 2) * pow(sin(_tau * k), 2) -
           4.0 * pow(k, 2) * sin(_tau * k) * cos(_tau * k) +
           2.0 * pow(k, 2) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(3, 4) =
      2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) +
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) +
      1.0 *
          (-2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(3, 5) = 0;
  _res(4, 0) =
      -2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) +
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) +
      1.0 *
          (2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(4, 1) =
      2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) +
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) +
      1.0 *
          (-2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(4, 2) =
      2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) -
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) -
      1.0 *
          (2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(4, 3) =
      2.0 * pow(k, 2) * exp(_tau * k) * sin(_tau * k) +
      2.0 * pow(k, 2) * exp(_tau * k) * cos(_tau * k) +
      1.0 *
          (-2.0 * pow(k, 2) * sin(_tau * k) + 2.0 * pow(k, 2) * cos(_tau * k)) *
          exp(-_tau * k);
  _res(4, 4) = 4.0 * pow(k, 2);
  _res(4, 5) = 0;
  _res(5, 0) = 0;
  _res(5, 1) = 0;
  _res(5, 2) = 0;
  _res(5, 3) = 0;
  _res(5, 4) = 0;
  _res(5, 5) = 0;
}

void compute_Qd2_block(double _tau, double _alpha,
                       Eigen::Ref<Eigen::MatrixXd> _res) {
  double k =
      0.35355339059327379 * pow(_alpha, 0.25) * pow(1.0 / (1.0 - _alpha), 0.25);
  _res(0, 0) =
      3.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) *
          pow(sin(_tau * k), 2) -
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) +
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (3.0 * pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) +
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 1) =
      -pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) +
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (-pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) +
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(0, 2) = -8.0 * pow(_tau, 3) * pow(k, 4) +
               8.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k);
  _res(0, 3) = 0;
  _res(0, 4) = 0;
  _res(0, 5) = 0;
  _res(1, 0) =
      -pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) -
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) +
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (-pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) +
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(1, 1) =
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) +
      3.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) *
          pow(cos(_tau * k), 2) -
      (pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) -
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       3.0 * pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(1, 2) = 0;
  _res(1, 3) = -8.0 * pow(_tau, 3) * pow(k, 4) -
               8.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k);
  _res(1, 4) = 0;
  _res(1, 5) = 0;
  _res(2, 0) = -8.0 * pow(_tau, 3) * pow(k, 4) +
               8.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k);
  _res(2, 1) = 0;
  _res(2, 2) =
      3.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) *
          pow(sin(_tau * k), 2) -
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) +
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) -
      (3.0 * pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) +
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 3) =
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) -
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      (-pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) +
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(2, 4) = 0;
  _res(2, 5) = 0;
  _res(3, 0) = 0;
  _res(3, 1) = -8.0 * pow(_tau, 3) * pow(k, 4) -
               8.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k);
  _res(3, 2) =
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) -
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(cos(_tau * k), 2) +
      (-pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) +
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(3, 3) =
      pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * pow(sin(_tau * k), 2) +
      2.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) * sin(_tau * k) *
          cos(_tau * k) +
      3.0 * pow(_tau, 2) * pow(k, 3) * exp(2 * _tau * k) *
          pow(cos(_tau * k), 2) -
      (pow(_tau, 2) * pow(k, 3) * pow(sin(_tau * k), 2) -
       2.0 * pow(_tau, 2) * pow(k, 3) * sin(_tau * k) * cos(_tau * k) +
       3.0 * pow(_tau, 2) * pow(k, 3) * pow(cos(_tau * k), 2)) *
          exp(-2 * _tau * k);
  _res(3, 4) = 0;
  _res(3, 5) = 0;
  _res(4, 0) = 0;
  _res(4, 1) = 0;
  _res(4, 2) = 0;
  _res(4, 3) = 0;
  _res(4, 4) = 0;
  _res(4, 5) = 0;
  _res(5, 0) = 0;
  _res(5, 1) = 0;
  _res(5, 2) = 0;
  _res(5, 3) = 0;
  _res(5, 4) = 0;
  _res(5, 5) = 0;
}

void compute_Q_block(double _tau, double _alpha,
                     Eigen::Ref<Eigen::MatrixXd> _res) {
  double k =
      0.35355339059327379 * pow(_alpha, 0.25) * pow(1.0 / (1.0 - _alpha), 0.25);
  _res(0, 0) =
      0.0625 *
      (exp(4 * _tau * k) * pow(sin(_tau * k), 2) +
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       3.0 * exp(4 * _tau * k) * pow(cos(_tau * k), 2) - pow(sin(_tau * k), 2) +
       2.0 * sin(_tau * k) * cos(_tau * k) - 3.0 * pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(0, 1) =
      -0.0625 *
      (-exp(4 * _tau * k) * pow(sin(_tau * k), 2) -
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       exp(4 * _tau * k) * pow(cos(_tau * k), 2) + pow(sin(_tau * k), 2) -
       2.0 * sin(_tau * k) * cos(_tau * k) - pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(0, 2) = 0.5 * _tau + 0.5 * sin(_tau * k) * cos(_tau * k) / k;
  _res(0, 3) = 0;
  _res(0, 4) = 0.25 * _tau * exp(_tau * k) * sin(_tau * k) +
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) +
               0.25 *
                   (-_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * sin(_tau * k) - sin(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(0, 5) =
      0.25 *
      (exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) +
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(1, 0) =
      -0.0625 *
      (-exp(4 * _tau * k) * pow(sin(_tau * k), 2) -
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       exp(4 * _tau * k) * pow(cos(_tau * k), 2) + pow(sin(_tau * k), 2) -
       2.0 * sin(_tau * k) * cos(_tau * k) - pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(1, 1) =
      0.0625 *
      (3.0 * exp(4 * _tau * k) * pow(sin(_tau * k), 2) -
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       exp(4 * _tau * k) * pow(cos(_tau * k), 2) - 3.0 * pow(sin(_tau * k), 2) -
       2.0 * sin(_tau * k) * cos(_tau * k) - pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(1, 2) = 0;
  _res(1, 3) = 0.5 * _tau - 0.5 * sin(_tau * k) * cos(_tau * k) / k;
  _res(1, 4) = 0.25 * _tau * exp(_tau * k) * sin(_tau * k) -
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) -
               0.25 *
                   (_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * cos(_tau * k) + cos(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(1, 5) =
      -0.25 *
      (-exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) -
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(2, 0) = 0.5 * _tau + 0.5 * sin(_tau * k) * cos(_tau * k) / k;
  _res(2, 1) = 0;
  _res(2, 2) =
      0.0625 *
      (exp(4 * _tau * k) * pow(sin(_tau * k), 2) +
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       3.0 * exp(4 * _tau * k) * pow(cos(_tau * k), 2) - pow(sin(_tau * k), 2) +
       2.0 * sin(_tau * k) * cos(_tau * k) - 3.0 * pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(2, 3) =
      0.0625 *
      (-exp(4 * _tau * k) * pow(sin(_tau * k), 2) -
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       exp(4 * _tau * k) * pow(cos(_tau * k), 2) + pow(sin(_tau * k), 2) -
       2.0 * sin(_tau * k) * cos(_tau * k) - pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(2, 4) = -0.25 * _tau * exp(_tau * k) * sin(_tau * k) -
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) -
               0.25 *
                   (-_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * sin(_tau * k) - sin(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(2, 5) =
      0.25 *
      (exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) +
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(3, 0) = 0;
  _res(3, 1) = 0.5 * _tau - 0.5 * sin(_tau * k) * cos(_tau * k) / k;
  _res(3, 2) =
      0.0625 *
      (-exp(4 * _tau * k) * pow(sin(_tau * k), 2) -
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       exp(4 * _tau * k) * pow(cos(_tau * k), 2) + pow(sin(_tau * k), 2) -
       2.0 * sin(_tau * k) * cos(_tau * k) - pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(3, 3) =
      0.0625 *
      (3.0 * exp(4 * _tau * k) * pow(sin(_tau * k), 2) -
       2.0 * exp(4 * _tau * k) * sin(_tau * k) * cos(_tau * k) +
       exp(4 * _tau * k) * pow(cos(_tau * k), 2) - 3.0 * pow(sin(_tau * k), 2) -
       2.0 * sin(_tau * k) * cos(_tau * k) - pow(cos(_tau * k), 2)) *
      exp(-2 * _tau * k) / k;
  _res(3, 4) = 0.25 * _tau * exp(_tau * k) * sin(_tau * k) -
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) -
               0.25 *
                   (_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * cos(_tau * k) + cos(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(3, 5) =
      0.25 *
      (-exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) -
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(4, 0) = 0.25 * _tau * exp(_tau * k) * sin(_tau * k) +
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) +
               0.25 *
                   (-_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * sin(_tau * k) - sin(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(4, 1) = 0.25 * _tau * exp(_tau * k) * sin(_tau * k) -
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) -
               0.25 *
                   (_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * cos(_tau * k) + cos(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(4, 2) = -0.25 * _tau * exp(_tau * k) * sin(_tau * k) -
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) -
               0.25 *
                   (-_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * sin(_tau * k) - sin(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(4, 3) = 0.25 * _tau * exp(_tau * k) * sin(_tau * k) -
               0.25 * _tau * exp(_tau * k) * cos(_tau * k) -
               0.25 *
                   (_tau * k * sin(_tau * k) + _tau * k * cos(_tau * k) -
                    exp(2 * _tau * k) * cos(_tau * k) + cos(_tau * k)) *
                   exp(-_tau * k) / k;
  _res(4, 4) = 0.33333333333333331 * pow(_tau, 3) * pow(k, 2);
  _res(4, 5) = 0;
  _res(5, 0) =
      0.25 *
      (exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) +
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(5, 1) =
      -0.25 *
      (-exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) -
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(5, 2) =
      0.25 *
      (exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) +
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(5, 3) =
      0.25 *
      (-exp(2 * _tau * k) * sin(_tau * k) + exp(2 * _tau * k) * cos(_tau * k) -
       sin(_tau * k) - cos(_tau * k)) *
      exp(-_tau * k) / k;
  _res(5, 4) = 0;
  _res(5, 5) = _tau;
}
}  // namespace gsplines::basis
