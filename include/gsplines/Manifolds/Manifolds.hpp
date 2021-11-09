#ifndef MANIFOLDS_H
#define MANIFOLDS_H
#include <eigen3/Eigen/Core>

class Manifold {
private:
public:
  Manifold();
  std::size_t dim() const;
  virtual ~Manifold();
};

class Topology {
public:
  Topology(arguments);
  virtual ~Topology();

private:
  /* data */
};

class Atlas {
public:
  Atlas();
  virtual ~Atlas();

private:
  /* data */
};
class Chart {
public:
  Chart(arguments);
  virtual ~Chart();

private:
  /* data */
};
#endif /* MANIFOLDS_H */
