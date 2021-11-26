
#ifndef INTERVAL_H
#define INTERVAL_H

#include <gsplines/Functions/Domain.hpp>
#include <gsplines/Functions/DomainInheritanceHelper.hpp>
#include <vector>
namespace gsplines {

namespace functions {
template <std::size_t DIM>
class Interval : public DomainInheritanceHelper<Interval<DIM>, Domain<DIM>> {
private:
  std::array<std::pair<double, double>, DIM> interval_;

public:
  typedef typename Domain<DIM>::Element Element;
  Interval(std::vector<std::pair<double, double>> &_interval);
  virtual ~Interval() = default;
  virtual bool contains(Element &_domain) {}
};
} // namespace functions
} // namespace gsplines
#endif /* INTERVAL_H */
