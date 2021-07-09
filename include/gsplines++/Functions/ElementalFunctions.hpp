#ifndef ELEMENTAL_FUNCTIONS
#define ELEMENTAL_FUNCTIONS

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <functional>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <iostream>
#include <list>
#include <memory>
#include <utility>
namespace gsplines {
namespace functions {

class ConstFunction : public Function {

private:
  Eigen::VectorXd values_;

public:
  ConstFunction(std::pair<double, double> _domain,
                Eigen::Ref<Eigen::VectorXd> _values);

  ConstFunction(std::pair<double, double> _domain, std::size_t _codom_dim,
                double _value);

  ConstFunction(const ConstFunction &_that);

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override {
    Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());

    result = Eigen::MatrixXd::Ones(_domain_points.size(), get_codom_dim())
                 .array()
                 .rowwise() *
             values_.transpose().array();

    return result;
  };

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<ConstFunction>(*this);
  }
  std::unique_ptr<Function> deriv(int _deg) override {
    return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(), 0.0);
  }
};

class Identity : public Function {

public:
  Identity(std::pair<double, double> _domain)
      : Function(_domain, 1, "Identity") {}

  Identity(const Identity &that) : Function(that) {}

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override {
    return _domain_points;
  };

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<Identity>(*this);
  }
  std::unique_ptr<Function> deriv(int _deg) override {
    return std::make_unique<ConstFunction>(get_domain(), get_codom_dim(), 1.0);
  }
};

class Exponential : public Function {
public:
  Exponential(std::pair<double, double> _domain)
      : Function(_domain, 1, "Exponential") {}

  Exponential(const Exponential &that) : Function(that) {}

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<Exponential>(*this);
  }
  std::unique_ptr<Function> deriv(int _deg) override {
    return std::make_unique<Exponential>(*this);
  }
};

class Sin;

class Cos : public Function {
public:
  Cos(std::pair<double, double> _domain) : Function(_domain, 1, "Cos") {}

  Cos(const Exponential &that) : Function(that) {}

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<Cos>(*this);
  }
  std::unique_ptr<Function> deriv(int _deg) override;
};

class Sin : public Function {
public:
  Sin(std::pair<double, double> _domain) : Function(_domain, 1, "Sin") {}

  Sin(const Exponential &that) : Function(that) {}

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<Sin>(*this);
  }

  std::unique_ptr<Function> deriv(int _deg) override;
};
} // namespace functions
} // namespace gsplines
#endif
