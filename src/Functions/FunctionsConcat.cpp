
#include <gsplines++/Functions/FunctionsConcat.hpp>
namespace gsplines {
namespace functions {

FunctionsConcat::FunctionsConcat(
    std::vector<std::unique_ptr<Function>> &_function_array)
    : Function({0, 1}, 8) {}

Eigen::MatrixXd FunctionsConcat::operator()(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(), codom_dim_);
  std::size_t result_cols(_domain_points.size());
  std::size_t current_interval = 0;
  double s, tau;

  int i, j;

  for (i = 0; i < result_cols; i++) {
    current_interval = get_interval(_domain_points(i));
    result.row(i) = function_array_[current_interval]->operator()(
        _domain_points.segment(i, 1));
  }
  return result;
}
std::unique_ptr<Function> FunctionsConcat::deriv(int _deg) {
  std::vector<std::unique_ptr<Function>> result_array;
  for (const std::unique_ptr<Function> &f : function_array_) {
    result_array.push_back(f->deriv(_deg));
  }
  return std::make_unique<FunctionsConcat>(result_array);
}

std::size_t FunctionsConcat::get_interval(double _domain_point) const {
  if (_domain_point <= domain_break_points_[0])
    return 0;
  for (int i = 0; i < number_of_intervals_; i++) {
    if (domain_break_points_[i] < _domain_point and
        _domain_point <= domain_break_points_[i + 1]) {
      return i;
    }
  }
  return number_of_intervals_ - 1;
}
} // namespace functions
} // namespace gsplines
