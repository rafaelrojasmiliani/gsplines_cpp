#ifndef DOMAININHERITANCEHELPER_H
#define DOMAININHERITANCEHELPER_H
#include <memory>

namespace gsplines {
namespace functions {

template <typename Current, typename Base>
class DomainInheritanceHelper : public Base {
public:
  using Base::Base;
  DomainInheritanceHelper(const DomainInheritanceHelper &_that)
      : Base(_that){};

  DomainInheritanceHelper(DomainInheritanceHelper &&_that)
      : Base(std::move(_that)){};

  std::unique_ptr<Current> clone() const & {
    return std::unique_ptr<Current>(static_cast<Current *>(this->clone_impl()));
  }

  std::unique_ptr<Current> clone() && {
    return std::unique_ptr<Current>(
        static_cast<Current *>(this->move_clone_impl()));
  }

  std::unique_ptr<Current> move_clone() {
    return std::unique_ptr<Current>(
        static_cast<Current *>(this->move_clone_impl()));
  }


  virtual ~DomainInheritanceHelper() = default;

protected:
  virtual DomainInheritanceHelper *clone_impl() const override {
    return new Current(*static_cast<const Current *>(this));
  }

  virtual DomainInheritanceHelper *move_clone_impl() override {
    return new Current(std::move(*static_cast<Current *>(this)));
  }
};
} // namespace functions
} // namespace gsplines
#endif /* DOMAININHERITANCEHELPER_H */

