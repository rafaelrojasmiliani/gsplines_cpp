#ifndef DOMAINBASE_H
#define DOMAINBASE_H

#include <Eigen/Core>
namespace gsplines {
namespace functions {

template <std::size_t DIM> class Domain {
private:
  const std::size_t dim_ = DIM;

public:
  bool contains(Domain &_domain) = 0;

  virtual ~Domain() = default;

  std::size_t dim() { return dim_; }

  typedef Eigen::Matrix<double, 1, DIM> Element;

  typedef Eigen::Matrix<double, Eigen::Dynamic, DIM> ElementArray;

protected:
  virtual Domain *clone_impl() const = 0;

  virtual Domain *move_clone_impl() = 0;
};

} // namespace functions
} // namespace gsplines
#endif /* DOMAINBASE_H */
