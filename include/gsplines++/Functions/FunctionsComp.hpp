#ifndef FUNCTIONCOMP_H
#define FUNCTIONCOMP_H

#include<gsplines++/Functions/Functions.hpp>

namespace gsplines {
namespace functions {
class FunctionsComp: public Function {
  public:
    FunctionsComp(const Function &_sf, const Function &_vf);
    FunctionsComp(const FunctionsComp &that);
    FunctionsComp(FunctionsComp &&that);
  std::unique_ptr<Function> clone() const override {return std::make_unique<FunctionsComp>(*this);}
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg=1) override;
  private:
  std::unique_ptr<Function> vf_;
  std::unique_ptr<Function> sf_;


};

FunctionsComp compose(const Function &_f1,
                                    const Function &_f2);
}
}


#endif /* FUNCTIONCOMP_H */
