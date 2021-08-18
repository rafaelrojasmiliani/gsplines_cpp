#include <gsplines/Functions/ElementalFunctions.hpp>
#include <iostream>
namespace gsplines {
namespace functions {

ConstFunction::ConstFunction(std::pair<double, double> _domain,
                             Eigen::Ref<Eigen::VectorXd> _values)
    : Function(_domain, _values.size(), "ConstFunction"), values_(_values) {}

ConstFunction::ConstFunction(std::pair<double, double> _domain,
                             std::size_t _codom_dim, double _value)
    : Function(_domain, _codom_dim, "ConstFunction"),
      values_(Eigen::VectorXd::Ones(_codom_dim) * _value) {}

ConstFunction::ConstFunction(const ConstFunction &_that)
    : Function(_that), values_(_that.values_) {}

void Exponential::value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                        Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = Eigen::exp(_domain_points.array()).matrix();
};

void Cos::value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = Eigen::cos(_domain_points.array()).matrix();
};

std::unique_ptr<FunctionExpression> Cos::deriv(int _deg) const {
  return std::make_unique<FunctionExpression>(
      ConstFunction(get_domain(), 1, -1.0) * Sin(get_domain()));
}

void Sin::value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = Eigen::sin(_domain_points.array()).matrix();
};

std::unique_ptr<FunctionExpression> Sin::deriv(int _deg) const {
  return std::make_unique<Cos>(get_domain());
}

void CanonicPolynomial::value(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {

  _result.setConstant(coefficients_[coefficients_.size() - 1]);
  for (int uici = coefficients_.size() - 2; uici >= 0; uici--) {
    _result.array() *= _domain_points.array();
    _result.array() += coefficients_[uici];
  }
}

std::unique_ptr<FunctionExpression> CanonicPolynomial::deriv(int _deg) const {

  if (_deg == 0)
    return std::make_unique<CanonicPolynomial>(*this);
  Eigen::VectorXd result(coefficients_);

  std::size_t vsize = result.size();

  if ((int)vsize - _deg <= 0) {
    return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(), 0.0);
  }
  for (std::size_t k = 1; k <= _deg; k++) {
    for (std::size_t uici = 0; uici < vsize - 1; uici++) {
      result(uici) = ((double)uici + 1.0) * result(uici + 1);
    }
    vsize--;
    if (vsize == 0) {
      return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(),
                                             result(0));
    }
  }
  return std::make_unique<CanonicPolynomial>(get_domain(),
                                             std::move(result.head(vsize)));
} // namespace functions

} // namespace functions
} // namespace gsplines
