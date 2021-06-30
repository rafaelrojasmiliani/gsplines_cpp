
#include <gsplines++/Functions/FunctionsSum.hpp>
namespace gsplines {
namespace functions {
SumOfFunctions::SumOfFunctions(std::vector<std::unique_ptr<Function>> &_array)
    : Function(_array.at(0)->get_domain(), _array.at(0)->get_codom_dim()) {
  for (std::unique_ptr<Function> &f : _array) {
    if (not Function::same_set(*(_array.at(0)), *f))
      throw 1;
    function_array_.push_back(f->clone());
  }
}

SumOfFunctions::SumOfFunctions(const SumOfFunctions &that)
    : Function(that.function_array_.at(0)->get_domain(),
               that.function_array_.at(0)->get_codom_dim()) {

  for (const std::unique_ptr<Function> &f : that.function_array_) {
    function_array_.push_back(f->clone());
  }
}

SumOfFunctions::SumOfFunctions(SumOfFunctions &&that)
    : Function(that.function_array_.at(0)->get_domain(),
               that.function_array_.at(0)->get_codom_dim()) {

  function_array_ = std::move(that.function_array_);
}

Eigen::MatrixXd SumOfFunctions::operator()(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());
  result.setZero();
  for (std::unique_ptr<Function> &f : function_array_) {
    result += f->value(_domain_points);
  }
  return result;
}

std::unique_ptr<Function> SumOfFunctions::deriv(int _deg) {
  std::vector<std::unique_ptr<Function>> result_array;
  for (std::unique_ptr<Function> &f : function_array_) {
    std::unique_ptr<Function> dfdt = f->deriv(_deg);
    result_array.push_back(std::move(dfdt));
  }
  return std::make_unique<SumOfFunctions>(result_array);
}

SumOfFunctions operator+(const Function &_f1, const Function &_f2) {
  std::vector<std::unique_ptr<Function>> array;
  array.push_back(_f1.clone());
  array.push_back(_f2.clone());
  return SumOfFunctions(array);
}
Eigen::MatrixXd
sum_functions(std::vector<std::unique_ptr<Function>> &_function_array,
              const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(), get_codom_dim());
  result.setZero();
  for (std::unique_ptr<Function> &f : function_array_) {
    result += f->value(_domain_points);
  }
  return result;
}
} // namespace functions
} // namespace gsplines
