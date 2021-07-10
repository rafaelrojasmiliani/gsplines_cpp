
#include <algorithm>
#include <gsplines++/Functions/FunctionExpression.hpp>
namespace gsplines {
namespace functions {

/* -----
 *  Function Evaluation
 * -----*/

std::list<std::unique_ptr<Function>>::const_iterator get_interval_function(
    double _domain_point,
    const std::list<std::unique_ptr<Function>> &_function_array) {

  std::list<std::unique_ptr<Function>>::const_iterator result = std::find_if(
      _function_array.begin(), _function_array.end(),
      [_domain_point](const std::unique_ptr<Function> &element) {
        return (element->get_domain().first - 1.0e-6 <= _domain_point) and
               (_domain_point <= element->get_domain().second + 1.0e-6);
      });

  return result;
}

std::unique_ptr<Function>
deriv_concat_functions(std::list<std::unique_ptr<Function>> &_function_array,
                       std::size_t _deg) {

  std::list<std::unique_ptr<Function>> result_array;

  for (std::unique_ptr<Function> &f : _function_array) {
    result_array.push_back(f->deriv(_deg));
  }

  return std::make_unique<FunctionExpression>(
      _function_array.front()->get_domain(),
      _function_array.front()->get_codom_dim(),
      FunctionExpression::Type::CONCATENATION, result_array);
}

Eigen::MatrixXd
eval_concat_functions(std::list<std::unique_ptr<Function>> &_function_array,
                      const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd result(_domain_points.size(),
                         _function_array.front()->get_codom_dim());

  Eigen::VectorXd vector_one(1);

  for (std::size_t i = 0; i < _domain_points.size(); i++) {

    std::list<std::unique_ptr<Function>>::const_iterator f =
        get_interval_function(_domain_points[i], _function_array);

    if (f != _function_array.end()) {
      result.row(i) = (*f)->value(_domain_points.segment(i, 1));
    } else if (_domain_points[i] <
               _function_array.front()->get_domain().first) {
      _function_array.front()->value(Eigen::VectorXd::Constant(
          1, _function_array.front()->get_domain().first));
    } else {

      _function_array.front()->value(Eigen::VectorXd::Constant(
          1, _function_array.back()->get_domain().second));
    }
  }

  return result;
}
} // namespace functions
} // namespace gsplines
