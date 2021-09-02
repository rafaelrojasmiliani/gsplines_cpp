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

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  ConstFunction *deriv_impl(std::size_t _deg = 1) const override;
};

class DomainLinearDilation
    : public FunctionInheritanceHelper<DomainLinearDilation, Function,
                                       ConstFunction> {
private:
  double dilation_factor_;

public:
  DomainLinearDilation(std::pair<double, double> _domain,
                       double _dilation_factor,
                       const std::string &_name = "DomainLinearDilation");

  DomainLinearDilation(const DomainLinearDilation &that);

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  virtual ConstFunction *deriv_impl(std::size_t _deg) const override;
};

class Identity : public DomainLinearDilation {

public:
  Identity(std::pair<double, double> _domain);

  Identity(const Identity &that);

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;
};

class Exponential
    : public FunctionInheritanceHelper<Exponential, Function, Exponential> {
public:
  Exponential(std::pair<double, double> _domain);

  Exponential(const Exponential &that);

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  virtual Exponential *deriv_impl(std::size_t _deg) const override;
};

class Cos
    : public FunctionInheritanceHelper<Cos, Function, FunctionExpression> {
public:
  Cos(std::pair<double, double> _domain);

  Cos(const Cos &that);

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  FunctionExpression *deriv_impl(std::size_t _deg) const override;
};

class Sin
    : public FunctionInheritanceHelper<Sin, Function, FunctionExpression> {
public:
  Sin(std::pair<double, double> _domain);

  Sin(const Sin &that);

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
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

  void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
             Eigen::Ref<Eigen::MatrixXd> _result) const override;

protected:
  CanonicPolynomial *deriv_impl(std::size_t _deg) const override;
};
} // namespace functions
} // namespace gsplines
#endif
