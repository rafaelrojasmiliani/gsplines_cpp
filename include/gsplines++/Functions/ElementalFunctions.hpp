#ifndef ELEMENTAL_FUNCTIONS
#define ELEMENTAL_FUNCTIONS

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <functional>
#include <gsplines++/Functions/Function.hpp>
#include <iostream>
#include <list>
#include <memory>
#include <utility>
namespace gsplines {
namespace functions {

class DomainLinearDilation;

class ConstFunction : public Function {

private:
  Eigen::VectorXd values_;

public:
  ConstFunction(std::pair<double, double> _domain,
                Eigen::Ref<Eigen::VectorXd> _values);

  ConstFunction(std::pair<double, double> _domain, std::size_t _codom_dim,
                double _value);

  ConstFunction(const ConstFunction &_that);

  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const override {
    Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());

    result = Eigen::MatrixXd::Ones(_domain_points.size(), get_codom_dim())
                 .array()
                 .rowwise() *
             values_.transpose().array();

    return result;
  };

  std::unique_ptr<FunctionExpression> clone() const override {
    return std::make_unique<ConstFunction>(*this);
  }
  std::unique_ptr<FunctionExpression> deriv(int _deg) const override {
    return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(), 0.0);
  }
};

class DomainLinearDilation : public Function {
private:
  double dilation_factor_;

public:
  DomainLinearDilation(std::pair<double, double> _domain,
                       double _dilation_factor,
                       const std::string &_name = "DomainLinearDilation")
      : Function(_domain, 1, _name), dilation_factor_(_dilation_factor) {}

  DomainLinearDilation(const DomainLinearDilation &that)
      : Function(that), dilation_factor_(that.dilation_factor_) {}

  virtual Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const override {
    return dilation_factor_ * _domain_points;
  };

  virtual std::unique_ptr<FunctionExpression> clone() const override {
    return std::make_unique<DomainLinearDilation>(*this);
  }
  virtual std::unique_ptr<FunctionExpression> deriv(int _deg) const override {
    return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(),
                                           dilation_factor_);
  }
};

class Identity : public DomainLinearDilation {

public:
  Identity(std::pair<double, double> _domain)
      : DomainLinearDilation(_domain, 1.0, "Identity") {}

  Identity(const Identity &that) : DomainLinearDilation(that) {}

  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const override {
    return _domain_points;
  };

  std::unique_ptr<FunctionExpression> clone() const override {
    return std::make_unique<Identity>(*this);
  }
  std::unique_ptr<FunctionExpression> deriv(int _deg) const override {
    return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(), 1.0);
  }
};

class Exponential : public Function {
public:
  Exponential(std::pair<double, double> _domain)
      : Function(_domain, 1, "Exponential") {}

  Exponential(const Exponential &that) : Function(that) {}

  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const override;

  std::unique_ptr<FunctionExpression> clone() const override {
    return std::make_unique<Exponential>(*this);
  }
  std::unique_ptr<FunctionExpression> deriv(int _deg) const override {
    return std::make_unique<Exponential>(*this);
  }
};

class Sin;

class Cos : public Function {
public:
  Cos(std::pair<double, double> _domain) : Function(_domain, 1, "Cos") {}

  Cos(const Exponential &that) : Function(that) {}

  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const override;

  std::unique_ptr<FunctionExpression> clone() const override {
    return std::make_unique<Cos>(*this);
  }
  std::unique_ptr<FunctionExpression> deriv(int _deg) const override;
};

class Sin : public Function {
public:
  Sin(std::pair<double, double> _domain) : Function(_domain, 1, "Sin") {}

  Sin(const Exponential &that) : Function(that) {}

  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) const override;

  std::unique_ptr<FunctionExpression> clone() const override {
    return std::make_unique<Sin>(*this);
  }

  std::unique_ptr<FunctionExpression> deriv(int _deg) const override;
};
} // namespace functions
} // namespace gsplines
#endif
