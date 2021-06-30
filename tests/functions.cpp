#include <cmath>
#include <gsplines++/Functions/Functions.hpp>
#include<iostream>
using namespace gsplines::functions;
class Exponential : public gsplines::functions::Function {
public:
  Exponential() : Function({-1, 1}, 1) {}

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override {
    Eigen::MatrixXd result(_domain_points.size(), 1);
    for (int i = 0; i < _domain_points.size(); i++)
      result(i) = sin(_domain_points(i));
    return result;
  };

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<Exponential>(*this);
  }
  std::unique_ptr<Function> deriv(int _deg) override {
    return std::make_unique<Exponential>(*this);
  }
};

int main() {

  Exponential s, g;
  Eigen::VectorXd time_span = Eigen::VectorXd::Random(4);
  Eigen::MatrixXd res = s(time_span);
  
 SumOfFunctions f = s + g + s +g;

  Eigen::MatrixXd res_2 = f(time_span);

  std::cout << (4*res - res_2) << std::endl;


  
  return 0;
}
