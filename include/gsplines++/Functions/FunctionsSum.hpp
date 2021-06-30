#ifndef FUNCTIONSSUM_H
#define FUNCTIONSSUM_H

#include <gsplines++/Functions/Functions.hpp>

namespace gsplines {
namespace functions {
class SumOfFunctions : public Function {
public:
  SumOfFunctions(std::vector<std::unique_ptr<Function>> &_function_array);
  SumOfFunctions(const SumOfFunctions &that);
  SumOfFunctions(SumOfFunctions &&that);
  std::unique_ptr<Function> clone() const override {
    return std::make_unique<SumOfFunctions>(*this);
  }
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg = 1) override;

private:
  std::vector<std::unique_ptr<Function>> function_array_;
};

SumOfFunctions operator+(const Function &_f1, const Function &_f2);

Eigen::MatrixXd
sum_functions(std::vector<std::unique_ptr<Function>> &_function_array);
} // namespace functions
} // namespace gsplines

#endif /* FUNCTIONSSUM_H */
