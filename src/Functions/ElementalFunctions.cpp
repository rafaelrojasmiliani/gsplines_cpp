#include <gsplines/Functions/ElementalFunctions.hpp>
#include <iostream>
namespace gsplines {
namespace functions {

// ConstFunction
ConstFunction::ConstFunction(std::pair<double, double> _domain,
                             Eigen::Ref<Eigen::VectorXd> _values)
    : FunctionInheritanceHelper(_domain, _values.size(), "ConstFunction"),
      values_(_values) {}

ConstFunction::ConstFunction(std::pair<double, double> _domain,
                             std::size_t _codom_dim, double _value)
    : FunctionInheritanceHelper(_domain, _codom_dim, "ConstFunction"),
      values_(Eigen::VectorXd::Ones(_codom_dim) * _value) {}

ConstFunction::ConstFunction(const ConstFunction &_that)
    : FunctionInheritanceHelper(_that), values_(_that.values_) {}

void ConstFunction::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {

  _result = Eigen::MatrixXd::Ones(_domain_points.size(), get_codom_dim())
                .array()
                .rowwise() *
            values_.transpose().array();
}

ConstFunction *ConstFunction::deriv_impl(std::size_t _deg) const {
  if (_deg == 0) {
    return new ConstFunction(*this);
  }
  return new ConstFunction(get_domain(), get_codom_dim(), 0.0);
}

// DomainLinearDilation
DomainLinearDilation::DomainLinearDilation(std::pair<double, double> _domain,
                                           double _dilation_factor,
                                           const std::string &_name)
    : FunctionInheritanceHelper(_domain, 1, _name),
      dilation_factor_(_dilation_factor) {}
DomainLinearDilation::DomainLinearDilation(const DomainLinearDilation &that)
    : FunctionInheritanceHelper(that), dilation_factor_(that.dilation_factor_) {
}
void DomainLinearDilation::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result.noalias() = dilation_factor_ * _domain_points;
};

ConstFunction *DomainLinearDilation::deriv_impl(std::size_t _deg) const {
  return new ConstFunction(get_domain(), get_codom_dim(), dilation_factor_);
}
// Identity
Identity::Identity(std::pair<double, double> _domain)
    : DomainLinearDilation(_domain, 1.0, "Identity") {}

Identity::Identity(const Identity &that) : DomainLinearDilation(that) {}

void Identity::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = _domain_points;
};

// Exponential
Exponential::Exponential(std::pair<double, double> _domain)
    : FunctionInheritanceHelper(_domain, 1, "Exponential") {}

Exponential::Exponential(const Exponential &that)
    : FunctionInheritanceHelper(that) {}

void Exponential::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = Eigen::exp(_domain_points.array()).matrix();
};
Exponential *Exponential::deriv_impl(std::size_t _deg) const {
  return new Exponential(*this);
}

// Cos
Cos::Cos(std::pair<double, double> _domain)
    : FunctionInheritanceHelper(_domain, 1, "Cos") {}

Cos::Cos(const Cos &that) : FunctionInheritanceHelper(that) {}

void Cos::value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                     Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = Eigen::cos(_domain_points.array()).matrix();
};

FunctionExpression *Cos::deriv_impl(std::size_t _deg) const {
  switch (_deg % 4) {
  case 0:
    return new FunctionExpression(Cos(*this));
  case 3:
    return new FunctionExpression(-Sin(get_domain()));
  case 2:
    return new FunctionExpression(-Cos(get_domain()));
  case 1:
    return new FunctionExpression(Sin(get_domain()));
  }
  throw std::invalid_argument("This case shoul not happen: Cos:deriv_impl");
}

// Sin
Sin::Sin(std::pair<double, double> _domain)
    : FunctionInheritanceHelper(_domain, 1, "Sin") {}

Sin::Sin(const Sin &that) : FunctionInheritanceHelper(that) {}

void Sin::value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                     Eigen::Ref<Eigen::MatrixXd> _result) const {
  _result = Eigen::sin(_domain_points.array()).matrix();
};

FunctionExpression *Sin::deriv_impl(std::size_t _deg) const {
  switch (_deg % 4) {
  case 0:
    return new FunctionExpression(Sin(*this));
  case 3:
    return new FunctionExpression(Cos(get_domain()));
  case 2:
    return new FunctionExpression(-Sin(get_domain()));
  case 1:
    return new FunctionExpression(-Cos(get_domain()));
  }
  throw std::invalid_argument("This case shoul not happen: Sin:deriv_impl");
}

// CanonicPolynomial
CanonicPolynomial::CanonicPolynomial(std::pair<double, double> _domain,
                                     const Eigen::VectorXd &_coefficients)
    : FunctionInheritanceHelper(_domain, 1, "CanonicPolynomial"),
      coefficients_(_coefficients) {}

CanonicPolynomial::CanonicPolynomial(std::pair<double, double> _domain,
                                     Eigen::VectorXd &&_coefficients)
    : FunctionInheritanceHelper(_domain, 1, "CanonicPolynomial"),
      coefficients_(std::move(_coefficients)) {}

CanonicPolynomial::CanonicPolynomial(const CanonicPolynomial &that)
    : FunctionInheritanceHelper(that), coefficients_(that.coefficients_) {}

CanonicPolynomial::CanonicPolynomial(CanonicPolynomial &&that)
    : FunctionInheritanceHelper(std::move(that)),
      coefficients_(std::move(that.coefficients_)) {}

void CanonicPolynomial::value_impl(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points,
    Eigen::Ref<Eigen::MatrixXd> _result) const {

  _result.setConstant(coefficients_[coefficients_.size() - 1]);
  for (int uici = coefficients_.size() - 2; uici >= 0; uici--) {
    _result.array() *= _domain_points.array();
    _result.array() += coefficients_[uici];
  }
}

CanonicPolynomial *CanonicPolynomial::deriv_impl(std::size_t _deg) const {

  if (_deg == 0)
    return new CanonicPolynomial(*this);
  Eigen::VectorXd result(coefficients_);

  std::size_t vsize = result.size();

  /*
  if ((int)vsize - _deg <= 0) {
    return new ConstFunction(get_domain(), get_codom_dim(), 0.0);
  }*/
  for (std::size_t k = 1; k <= _deg; k++) {
    for (std::size_t uici = 0; uici < vsize - 1; uici++) {
      result(uici) = ((double)uici + 1.0) * result(uici + 1);
    }
    vsize--;
    /*
    if (vsize == 0) {
      return new ConstFunction(get_domain(), get_codom_dim(), result(0));
    }
    */
  }
  return new CanonicPolynomial(get_domain(), std::move(result.head(vsize)));
} // namespace functions

} // namespace functions
} // namespace gsplines
