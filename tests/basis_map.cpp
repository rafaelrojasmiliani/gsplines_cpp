
#include <eigen3/Eigen/Core>
#include <gsplines/Basis/Basis.hpp>
#include <gsplines/Basis/BasisLegendre.hpp>
#include <gsplines/Tools.hpp>
#include <gtest/gtest.h>
#include <unordered_map>
using namespace gsplines::basis;

template <class T> inline void hash_combine(std::size_t &seed, const T &v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct BasisIndexer {
  std::string name_;
  std::vector<int32_t> p1_;
  std::vector<double> p2_;
  BasisIndexer(const std::string &_name, const std::vector<int32_t> &_p1,
               const std::vector<double> &_p2)
      : name_(_name), p1_(_p1), p2_(_p2) {}
  bool operator==(const BasisIndexer & /*_in*/) const { return true; }
};

namespace std {
template <> struct hash<BasisIndexer> {
  size_t operator()(const BasisIndexer &_in) const {

    std::size_t result = 0;

    hash_combine(result, _in.name_);

    for (int int_par : _in.p1_) {
      hash_combine(result, int_par);
    }

    for (double double_par : _in.p2_) {
      // take just 6 decimals
      int par = double_par * 10e6;
      hash_combine(result, par);
    }

    return result;
  }
};
} // namespace std

class BasisTest : public Basis {
public:
  std::vector<double> double_pars_;
  std::vector<int> int_pars_;
  void eval_on_window(double /*_s*/, double /*_tau*/,
                      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>>
                      /*_buff*/) const override {}
  void eval_derivative_on_window(
      double /*_s*/, double /*_tau*/, unsigned int /*_deg*/,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> /*_buff*/)
      const override{};
  void eval_derivative_wrt_tau_on_window(
      double /*_s*/, double /*_tau*/, unsigned int /*_deg*/,
      Eigen::Ref<Eigen::VectorXd, 0, Eigen::InnerStride<>> /*_buff*/)
      const override{};
  void add_derivative_matrix(double /*tau*/, std::size_t /*_deg*/,
                             Eigen::Ref<Eigen::MatrixXd> /*_mat*/) override {}
  void add_derivative_matrix_deriv_wrt_tau(
      double /*tau*/, std::size_t /*_deg*/,
      Eigen::Ref<Eigen::MatrixXd> /*_mat*/) override {}

  std::unique_ptr<Basis> clone() const override { return nullptr; }
  std::unique_ptr<Basis> move_clone() override { return nullptr; }

  Eigen::MatrixXd derivative_matrix_impl(std::size_t /*_deg*/) const override {
    return Eigen::MatrixXd::Zero(1, 1);
  };
};

TEST(BasisMap, Value) {
  //
  std::unordered_map<BasisIndexer, std::shared_ptr<Basis>> map;

  BasisIndexer index("Legendre", {5}, {});
  map[index] = std::make_shared<gsplines::basis::BasisLegendre>(5);
}

int main(int argc, char **argv) {

  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
