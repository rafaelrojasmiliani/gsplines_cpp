#ifndef FUNCTIONINHERITANCEHELPER_H
#define FUNCTIONINHERITANCEHELPER_H
#include <cstddef>
#include <memory>

namespace gsplines {
namespace functions {

template <typename Current, typename Base, typename DerivativeClass>
class FunctionInheritanceHelper : public Base {
 public:
  using Base::Base;
  FunctionInheritanceHelper(const Base& _that) : Base(_that){};

  FunctionInheritanceHelper(Base&& _that) : Base(std::move(_that)){};

  std::unique_ptr<Current> clone() const& {
    return std::unique_ptr<Current>(static_cast<Current*>(this->clone_impl()));
  }

  std::unique_ptr<Current> clone() && {
    return std::unique_ptr<Current>(
        static_cast<Current*>(this->move_clone_impl()));
  }

  std::unique_ptr<Current> move_clone() {
    return std::unique_ptr<Current>(
        static_cast<Current*>(this->move_clone_impl()));
  }

  std::unique_ptr<DerivativeClass> deriv(std::size_t _deg) const& {
    return std::unique_ptr<DerivativeClass>(
        static_cast<DerivativeClass*>(this->deriv_impl(_deg)));
  }

  DerivativeClass derivate(std::size_t _deg = 1) const {
    return DerivativeClass(*deriv(_deg));
  }

  ~FunctionInheritanceHelper() = default;

 protected:
  FunctionInheritanceHelper* clone_impl() const override {
    return new Current(*static_cast<const Current*>(this));  // NOLINT
  }

  FunctionInheritanceHelper* move_clone_impl() override {
    return new Current(std::move(*static_cast<Current*>(this)));  // NOLINT
  }
};
}  // namespace functions
}  // namespace gsplines
#endif /* FUNCTIONINHERITANCEHELPER_H */
