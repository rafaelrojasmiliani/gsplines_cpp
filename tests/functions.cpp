#include <cmath>
#include <gsplines++/Functions/FunctionExpression.hpp>
#include <iostream>
using namespace gsplines::functions;
class Exponential : public gsplines::functions::Function {
public:
    static std::size_t num_call_clone_;
  Exponential() : Function({-1, 1}, 1){}

  Exponential(const Exponential& that) : Function(that) {}

  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override {
    Eigen::MatrixXd result(_domain_points.size(), 1);
    for (int i = 0; i < _domain_points.size(); i++)
      result(i) = exp(_domain_points(i));
    return result;
  };

  std::unique_ptr<Function> clone() const override {
    return std::make_unique<Exponential>(*this);
    num_call_clone_++;
  }
  std::unique_ptr<Function> deriv(int _deg) override {
    return std::make_unique<Exponential>(*this);
  }
};
std::size_t Exponential::num_call_clone_ = 0;

int main() {
    Exponential s, g;
    Eigen::VectorXd time_span = Eigen::VectorXd::Random(4);
    Eigen::MatrixXd exp_value = s(time_span);

    gsplines::functions::FunctionExpression f = s + g + s + g;

    printf(".................f = s + g + s +g\n");
    f.print_performace();
    printf(".................\n");
    gsplines::functions::FunctionExpression z = f;
    printf(".................z=f\n");
    f.print_performace();
    printf(".................\n");
    gsplines::functions::FunctionExpression m = f + g + s;
    printf(".................m = f + g + s;\n");
    f.print_performace();
    printf(".................\n");

    printf("|||||||||||||||||||||||||");
    assert((4*exp_value - f(time_span)).norm()<1.0e-9);
    printf("|||||||||||||||||||||||||");
    assert((m(time_span) - 6*exp_value).norm()<1.0e-9);
    printf("|||||||||||||||||||||||||");

    f+=g;
    assert((f(time_span) - 5*exp_value).norm()<1.0e-9);
    printf("|||||||||||||||||||||||||");

    f+=m;
    assert((f(time_span) - 11*exp_value).norm()<1.0e-9);
    printf("|||||||||||||||||||||||||");

    f+=(m +s);
    f.print_performace();

  return 0;
}
