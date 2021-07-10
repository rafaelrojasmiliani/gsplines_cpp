

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

  const std::pair<double, double> &get_domain() const { return domain_; };
  std::size_t get_codom_dim() const { return codom_dim_; };

  static bool same_domain(const FunctionBase &_f1, const FunctionBase &_f2);
  static bool same_codomain(const FunctionBase &_f1, const FunctionBase &_f2);

  virtual ~FunctionBase() {}

  bool is_point_in_domain(double _domain_point);
};

class Function : public FunctionBase {
private:
  const std::string name_;

public:
  Function(std::pair<double, double> _domain, std::size_t _codom_dim);
  Function(std::pair<double, double> _domain, std::size_t _codom_dim,
           const std::string &_name);
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
  virtual FunctionExpression operator-() const;

  const std::string &get_name() const { return name_; }

  virtual void print(std::size_t _indent = 0) {
    printf("%*s- %s\n", 4 * (int)_indent, "", get_name().c_str());
  }

  virtual FunctionExpression compose(const Function &that) const;
  virtual FunctionExpression compose(const FunctionExpression &that) const;
  virtual FunctionExpression compose(FunctionExpression &&that) const;
};

} // namespace functions
} // namespace gsplines
#endif
