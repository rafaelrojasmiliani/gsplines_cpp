#ifndef ELEMENTAL_FUNCTIONS
#define ELEMENTAL_FUNCTIONS

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <functional>
#include <gsplines/Functions/Function.hpp>
#include <gsplines/Functions/FunctionExpression.hpp>
#include <gsplines/Functions/FunctionInheritanceHelper.hpp>
#include <iostream>
#include <list>
#include <memory>
#include <utility>
namespace gsplines {
namespace functions {

class DomainLinearDilation;
// ----------------------------
// Const Function
// ---------------------------

class ConstFunction
    : public FunctionInheritanceHelper<ConstFunction, Function, ConstFunction> {

private:
  Eigen::VectorXd values_;

public:
  ConstFunction(std::pair<double, double> _domain,
                Eigen::Ref<Eigen::VectorXd> _values);

  ConstFunction(std::pair<double, double> _domain, std::size_t _codom_dim,
                double _value);

  ConstFunction(const ConstFunction &_that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  ConstFunction *deriv_impl(std::size_t _deg = 1) const override;
};
// ----------------------------
// Domain Linear Dilation
// ---------------------------
class DomainLinearDilation
    : public FunctionInheritanceHelper<DomainLinearDilation, Function,
                                       FunctionExpression> {
private:
  double dilation_factor_;

public:
  DomainLinearDilation(std::pair<double, double> _domain,
                       double _dilation_factor);

  DomainLinearDilation(const DomainLinearDilation &that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  virtual FunctionExpression *deriv_impl(std::size_t _deg) const override;
};

// ----------------------------
// Identity
// ---------------------------
class Identity : public DomainLinearDilation {

public:
  Identity(std::pair<double, double> _domain);

  Identity(const Identity &that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;
};

// ----------------------------
// Exponential
// ---------------------------
class Exponential
    : public FunctionInheritanceHelper<Exponential, Function, Exponential> {
public:
  Exponential(std::pair<double, double> _domain);

  Exponential(const Exponential &that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  virtual Exponential *deriv_impl(std::size_t _deg) const override;
};

// ----------------------------
// Cos
// ---------------------------
class Cos
    : public FunctionInheritanceHelper<Cos, Function, FunctionExpression> {
public:
  Cos(std::pair<double, double> _domain);

  Cos(const Cos &that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  FunctionExpression *deriv_impl(std::size_t _deg) const override;
};

class Sin
    : public FunctionInheritanceHelper<Sin, Function, FunctionExpression> {
public:
  Sin(std::pair<double, double> _domain);

  Sin(const Sin &that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  FunctionExpression *deriv_impl(std::size_t _deg) const override;
};

class CanonicPolynomial
    : public FunctionInheritanceHelper<CanonicPolynomial, Function,
                                       CanonicPolynomial> {
protected:
  Eigen::VectorXd coefficients_;

public:
  CanonicPolynomial(std::pair<double, double> _domain,
                    const Eigen::VectorXd &_coefficients);

  CanonicPolynomial(std::pair<double, double> _domain,
                    Eigen::VectorXd &&_coefficients);

  CanonicPolynomial(const CanonicPolynomial &that);

  CanonicPolynomial(CanonicPolynomial &&that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  CanonicPolynomial *deriv_impl(std::size_t _deg) const override;
};

class DotProduct : public FunctionInheritanceHelper<DotProduct, Function,
                                                    FunctionExpression> {
private:
  FunctionExpression f1_;
  FunctionExpression f2_;

public:
  DotProduct(const FunctionBase &_f1, const FunctionBase &_f2);
  DotProduct(FunctionBase &&_f1, FunctionBase &&_f2);
  DotProduct(const FunctionBase &_f1, FunctionBase &&_f2);

  DotProduct(const DotProduct &that);

  DotProduct(DotProduct &&that);

  void value_impl(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                  Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  FunctionExpression *deriv_impl(std::size_t _deg) const override;
  FunctionExpression *first_deriv_impl(std::size_t _deg) const;
};
} // namespace functions
} // namespace gsplines
#endif
