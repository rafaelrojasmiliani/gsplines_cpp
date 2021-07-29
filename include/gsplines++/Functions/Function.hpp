#ifndef FUNCTION
#define FUNCTION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class FunctionExpression;

class Function : public FunctionExpression {

public:
  Function(std::pair<double, double> _domain, std::size_t _codom_dim,
           const std::string &_name);
  Function(const Function &that);

  virtual std::unique_ptr<FunctionExpression>
  deriv(int _deg = 1) const override = 0;

  virtual void value(const Eigen::Ref<const Eigen::VectorXd> _domain_points,
                     Eigen::Ref<Eigen::MatrixXd> _result) const override = 0;

  virtual std::unique_ptr<FunctionExpression> clone() const override = 0;

  virtual std::unique_ptr<FunctionExpression> move_clone() override = 0;

  virtual ~Function() {}

  FunctionExpression derivate(int _deg = 1) const override;
};

} // namespace functions
} // namespace gsplines
#endif
