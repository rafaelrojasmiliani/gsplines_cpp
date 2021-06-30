

#include <gsplines++/Functions/FunctionsComp.hpp>
namespace gsplines {
namespace functions {

FunctionsComp::FunctionsComp(const Function &_sf, const Function &_vf)
    : Function(_vf.get_domain(), _vf.get_codom_dim()), sf_(_sf.clone()),
      vf_(_vf.clone()) {}

FunctionsComp::FunctionsComp(const FunctionsComp &that)
    : Function(that), sf_(that.sf_->clone()), vf_(that.vf_->clone()) {}

Eigen::MatrixXd FunctionsComp::
operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd scalar_res = sf_->operator()(_domain_points);

  return  vf_->operator()(scalar_res);
}

FunctionsComp compose(const Function &_f1,
                                    const Function &_f2) {
 assert(_f2.get_codom_dim()==1);
  return FunctionsComp(_f2, _f1);
}

std::unique_ptr<Function> FunctionsComp::deriv(int _deg) {
  return nullptr;
};
}
}
