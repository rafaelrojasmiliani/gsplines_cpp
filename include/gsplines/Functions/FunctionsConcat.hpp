#ifndef FUNCTIONSCONCAT_H
#define FUNCTIONSCONCAT_H

#include <gsplines/Functions/Functions.hpp>
namespace gsplines {
namespace functions {
class FunctionsConcat : public Function {

public:
  FunctionsConcat(std::vector<std::unique_ptr<Function>> &_function_array);
  FunctionsConcat(const FunctionsConcat &that);
  FunctionsConcat(FunctionsConcat &&that);
  std::unique_ptr<Function> clone() const override {
    return std::make_unique<FunctionsConcat>(*this);
  }
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg = 1) override;

private:
  std::vector<std::unique_ptr<Function>> function_array_;
  std::size_t number_of_intervals_;
  std::vector<double> domain_break_points_;
  std::size_t get_interval(double _domain_point) const;
};

FunctionsConcat concatenate(const Function &_f1, const Function &_f2);

} // namespace functions
} // namespace gsplines
#endif /* FUNCTIONSCONCAT_H */
