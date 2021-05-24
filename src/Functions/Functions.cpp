
#include <cmath>
#include <gsplines++/Functions/Functions.hpp>
#include <stdexcept>

namespace gsplines {
namespace functions {
const double Function::dom_tollerance_ = 1.0e-5;
Function::Function(std::pair<double, double> _domain, std::size_t _codom_dim)
    : codom_dim_(_codom_dim), window_(_domain), domain_(_domain) {}

Function::Function(const Function &that)
    : codom_dim_(that.codom_dim_), window_(that.window_),
      domain_(that.domain_) {}

bool Function::same_set(Function &_f1, Function &_f2) {
  bool err1 =
      abs(_f1.domain_.first - _f2.domain_.first) < Function::dom_tollerance_;
  bool err2 =
      abs(_f1.domain_.second - _f2.domain_.second) < Function::dom_tollerance_;

  bool err3 = _f1.codom_dim_ == _f2.codom_dim_;

  return err1 and err2 and err3;
}

bool Function::same_domain(Function &_f1, Function &_f2) {
  bool err1 =
      abs(_f1.domain_.first - _f2.domain_.first) < Function::dom_tollerance_;
  bool err2 =
      abs(_f1.domain_.second - _f2.domain_.second) < Function::dom_tollerance_;

  return err1 and err2;
}
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

Eigen::MatrixXd SumOfFunctions::
operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

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

SumOfFunctions operator+(const Function &_f1,
                                    const Function &_f2) {
  std::vector<std::unique_ptr<Function>> array;
  array.push_back(_f1.clone());
  array.push_back(_f2.clone());
  return SumOfFunctions(array);
}
/*
MulScalarFunction::MulScalarFunction(Function &_sf, Function &_vf)
    : Function(_vf.get_domain(), _vf.get_codom_dim()), sf_(_sf.clone()),
      vf_(_vf.clone()) {}

MulScalarFunction::MulScalarFunction(const MulScalarFunction &that)
    : Function(that), sf_(that.sf_->clone()), vf_(that.vf_->clone()) {}

Eigen::MatrixXd MulScalarFunction::
operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd scalar_res = sf_->operator()(_domain_points);
  Eigen::MatrixXd vector_res = vf_->operator()(_domain_points);

  return vector_res.array().colwise() * scalar_res.array();
}

std::unique_ptr<Function> MulScalarFunction::deriv(int _deg) {
  return sf_->deriv() * vf_ + sf_ * vf_->deriv();
};


*/
} // namespace functions
} // namespace gsplines
