
#include <gsplines++/Functions/FunctionExpression.hpp>
namespace gsplines {
namespace functions {

/* -----
 *  Function Evaluation
 * -----*/

std::unique_ptr<Function>
deriv_concat_functions(std::list<std::unique_ptr<Function>> &_function_array,
                       std::size_t _deg) {

  std::list<std::unique_ptr<Function>> result_array;
  for (std::unique_ptr<Function> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }
  std::size_t codom_dim = _function_array.front()->get_codom_dim();
  std::pair<double, double> domain = _function_array.front()->get_domain();
  return std::make_unique<FunctionExpression>(
      domain, codom_dim, FunctionExpression::Type::CONCATENATION, result_array);
}
Eigen::MatrixXd
eval_concat_functions(std::list<std::unique_ptr<Function>> &_function_array,
                      const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  return Eigen::MatrixXd();
}
} // namespace functions
} // namespace gsplines
