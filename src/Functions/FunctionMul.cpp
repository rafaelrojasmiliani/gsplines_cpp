
#include <gsplines++/Functions/FunctionsMul.hpp>
namespace gsplines {
namespace functions {

MulScalarFunction::MulScalarFunction(const Function &_sf, const Function &_vf)
    : Function(_vf.get_domain(), _vf.get_codom_dim()), sf_(_sf.clone()),
      vf_(_vf.clone()) {}

MulScalarFunction::MulScalarFunction(const MulScalarFunction &that)
    : Function(that), sf_(that.sf_->clone()), vf_(that.vf_->clone()) {}

Eigen::MatrixXd MulScalarFunction::operator()(
    const Eigen::Ref<const Eigen::VectorXd> _domain_points) {

  Eigen::MatrixXd scalar_res = sf_->operator()(_domain_points);
  Eigen::MatrixXd vector_res = vf_->operator()(_domain_points);

  return vector_res.array().colwise() * scalar_res.col(0).array();
}

MulScalarFunction operator*(const Function &_f1, const Function &_f2) {
  return MulScalarFunction(_f1, _f2);
}

std::unique_ptr<Function> MulScalarFunction::deriv(int _deg) {
  return nullptr;
};
} // namespace functions
} // namespace gsplines
