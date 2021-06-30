
#ifndef FUNCTIONSMUL_H
#define FUNCTIONSMUL_H

#include<gsplines++/Functions/Functions.hpp>

namespace gsplines {
namespace functions {
class MulScalarFunction : public Function {
private:
  std::unique_ptr<Function> vf_;
  std::unique_ptr<Function> sf_;

public:
  MulScalarFunction(const Function &_sf, const Function &_vf);
  MulScalarFunction(const MulScalarFunction &that);
  Eigen::MatrixXd operator()(
      const Eigen::Ref<const Eigen::VectorXd> _domain_points) override final;
  std::unique_ptr<Function> deriv(int _deg=1) override final;
  std::unique_ptr<Function> clone() const override {return std::make_unique<MulScalarFunction>(*this);}
  ~MulScalarFunction() {}
};

std::unique_ptr<Function> operator*(Function &, Function &);
}
}
#endif /* FUNCTIONSMUL_H */
