
#ifndef PIECEWISEFUNCTION_H
#define PIECEWISEFUNCTION_H

#include <cstddef>
#include <eigen3/Eigen/Core>
#include <memory>
#include <utility>
#include <vector>
#include <gsplines++/Functions/Functions.hpp>

namespace gsplines {
namespace functions {
class PieceWiseFunction: public Function {
  public:
    PieceWiseFunction(std::vector<std::unique_ptr<Function>> &_function_array);
    PieceWiseFunction(const PieceWiseFunction &that);
    PieceWiseFunction(PieceWiseFunction &&that);
  std::unique_ptr<Function> clone() const override {return std::make_unique<PieceWiseFunction>(*this);}
  Eigen::MatrixXd
  operator()(const Eigen::Ref<const Eigen::VectorXd> _domain_points) override;
  std::unique_ptr<Function> deriv(int _deg) override;
  private:
    std::vector<std::unique_ptr<Function>> function_array_;


};


PieceWiseFunction concat(const Function &_f1,
                                    const Function &_f2);

}
}
#endif /* PIECEWISEFUNCTION_H */
