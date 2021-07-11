#ifndef FUNCTION_BASE
#define FUNCTION_BASE

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class FunctionBase {
private:
  std::size_t codom_dim_;
  std::pair<double, double> window_;
  std::pair<double, double> domain_;
  const std::string name_;

public:
  static const double dom_tollerance_;
  FunctionBase(std::pair<double, double> _domain, std::size_t _codom_dim,
               const std::string &_name);
  FunctionBase(const FunctionBase &that);

  const std::pair<double, double> &get_domain() const { return domain_; };

  std::size_t get_codom_dim() const { return codom_dim_; };

  static bool same_domain(const FunctionBase &_f1, const FunctionBase &_f2);

  static bool same_codomain(const FunctionBase &_f1, const FunctionBase &_f2);

  virtual ~FunctionBase() {}

  bool is_point_in_domain(double _domain_point) const;

  virtual void print(std::size_t _indent = 0) const;

  const std::string &get_name() const { return name_; }
};

} // namespace functions
} // namespace gsplines
#endif
