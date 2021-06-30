#ifndef SUMOFFUNCTIONS_H
#define SUMOFFUNCTIONS_H
#include <gsplines++/Functions/Functions.hpp>
namespace gsplines {
namespace functions {
class SumOfFunctions : public Function {
private:
  std::vector<std::unique_ptr<Function>> array_;

public:
  SumOfFunctions(std::vector<std::unique_ptr<Function>> _array);
  SumOfFunctions(const SumOfFunctions &that);
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg) override;
  std::unique_ptr<Function> clone() const override {

    return std::make_unique<SumOfFunctions>(*this);
  }
};

std::unique_ptr<Function> operator+(const std::unique_ptr<Function> &_f1,
                                    const std::unique_ptr<Function> &_f2);
} // namespace functions
} // namespace gsplines
#endif /* SUMOFFUNCTIONS_H */
