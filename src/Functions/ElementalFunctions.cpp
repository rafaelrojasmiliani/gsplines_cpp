#include <gsplines++/Functions/ElementalFunctions.hpp>
namespace gsplines {
namespace functions {

ConstFunction::ConstFunction(std::pair<double, double> _domain,
                             Eigen::Ref<Eigen::VectorXd> _values)
    : Function(_domain, _values.size()), values_(_values) {}

ConstFunction::ConstFunction(std::pair<double, double> _domain,
                             std::size_t _codom_dim, double _value)
    : Function(_domain, _codom_dim),
      values_(Eigen::VectorXd::Ones(_codom_dim) * _value) {}

ConstFunction::ConstFunction(const ConstFunction &_that)
    : Function(_that), values_(_that.values_) {}

Eigen::MatrixXd Exponential::operator()(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
  return Eigen::exp(_domain_points.array()).matrix();
};

Eigen::MatrixXd
Cos::operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
  return Eigen::cos(_domain_points.array()).matrix();
};

std::unique_ptr<Function> Cos::deriv(int _deg) {
  return std::make_unique<FunctionExpression>(
      ConstFunction(get_domain(), 1, -1.0) * Sin(get_domain()));
}

Eigen::MatrixXd
Sin::operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
  return Eigen::cos(_domain_points.array()).matrix();
};

std::unique_ptr<Function> Sin::deriv(int _deg) {
  return std::make_unique<Cos>(get_domain());
}
} // namespace functions
} // namespace gsplines
