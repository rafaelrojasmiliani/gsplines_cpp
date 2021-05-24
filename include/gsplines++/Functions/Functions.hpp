

#ifndef FUNCTION
#define FUNCTION

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <utility>
#include <vector>

namespace gsplines {
namespace functions {

class Function {
private:
  const std::size_t codom_dim_;
  const std::pair<double, double> window_;
  const std::pair<double, double> domain_;
  static const double dom_tollerance_;

public:
  Function(std::pair<double, double> _domain, std::size_t _codom_dim);
  Function(const Function &that);

  virtual Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) = 0;

  virtual Eigen::MatrixXd
  value(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {
    return this->operator()(_domain_points);
  }

  virtual std::unique_ptr<Function> deriv(int _deg) = 0;
  std::pair<double, double> get_domain() { return domain_; };
  std::size_t get_codom_dim() { return codom_dim_; };

  static bool same_set(Function &_f1, Function &_f2);
  static bool same_domain(Function &_f1, Function &_f2);

  virtual std::unique_ptr<Function> clone() const = 0;
  virtual ~Function() {}
};

class SumOfFunctions: public Function {
  public:
    SumOfFunctions(std::vector<std::unique_ptr<Function>> &_function_array);
    SumOfFunctions(const SumOfFunctions &that);
    SumOfFunctions(SumOfFunctions &&that);
  std::unique_ptr<Function> clone() const override {return std::make_unique<SumOfFunctions>(*this);}
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg) override;
  private:
    std::vector<std::unique_ptr<Function>> function_array_;


};
SumOfFunctions operator+(const Function &_f1,
                                    const Function &_f2);
/*
std::unique_ptr<Function> operator-(Function &, Function &);
std::unique_ptr<Function> operator*(Function &, Function &);

class MulScalarFunction : public Function {
private:
  std::unique_ptr<Function> vf_;
  std::unique_ptr<Function> sf_;

public:
  MulScalarFunction(Function &_sf, Function &_vf);
  MulScalarFunction(const MulScalarFunction &that);
  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) override final;
  std::unique_ptr<Function> deriv(int _deg) override final;
  ~MulScalarFunction() {}
};


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
