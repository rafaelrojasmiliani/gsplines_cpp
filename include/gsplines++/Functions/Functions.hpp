

#ifndef FUNCTION
#define FUNCTION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class FunctionExpression;

class FunctionBase {
public:
  std::size_t codom_dim_;
  std::pair<double, double> window_;
  std::pair<double, double> domain_;
  static const double dom_tollerance_;

public:
  FunctionBase(std::pair<double, double> _domain, std::size_t _codom_dim);
  FunctionBase(const FunctionBase &that);

  std::pair<double, double> get_domain() const { return domain_; };
  std::size_t get_codom_dim() const { return codom_dim_; };

  static bool same_domain(const FunctionBase &_f1, const FunctionBase &_f2);
  static bool same_codomain(const FunctionBase &_f1, const FunctionBase &_f2);

  virtual ~FunctionBase() {}

  bool is_point_in_domain(double _domain_point);
};

class Function : public FunctionBase {
private:
public:
  Function(std::pair<double, double> _domain, std::size_t _codom_dim);
  Function(const Function &that);

  virtual Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) = 0;

  virtual Eigen::MatrixXd
  value(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
    return this->operator()(_domain_points);
  }

  virtual std::unique_ptr<Function> deriv(int _deg = 1) = 0;

  static bool same_set(Function &_f1, Function &_f2);
  static bool same_domain(Function &_f1, Function &_f2);

  virtual std::unique_ptr<Function> clone() const = 0;
  virtual ~Function() {}
  FunctionExpression operator-() const;
};

/*
std::unique_ptr<Function> operator-(Function &, Function &);

class CompositionOfFunctions : public Function {
private:
  std::unique_ptr<Function> f1_;
  std::unique_ptr<Function> f2_;

public:
  CompositionOfFunctions(Function &_f1, Function &_f2);
  CompositionOfFunctions(const CompositionOfFunctions &that);
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg) override;
};
*/
} // namespace functions
} // namespace gsplines
#endif
