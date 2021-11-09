#ifndef MANIFOLD_H
#define MANIFOLD_H
#include <Eigen/Core>
#include <map>
#include <utility>
#include <vector>

class ManifoldElement {
public:
  ManifoldElement();
  virtual ~ManifoldElement() = default;
  const Eigen::VectorXd &charts();
  const Eigen::VectorXd &default_chart();
};

template <typename X, typename Y> class Map;
template <typename X> class Chart;
class CoordinateTransformation;


template <typename Element> class Manifold {
private:
  std::vector<Chart<Manifold>> atlas_;
  std::map<Chart<Manifold>, std::map<Chart<Manifold>, CoordinateTransformation>>
      smooth_structure_;

public:
  Manifold();
};

class Rn;

template <typename X, typename Y> class Map {
public:
  Map(Map<X, Rn> _domain_chart, Map<Rn, Y> _codomain_chart);
};

#endif /* MANIFOLD_H */
