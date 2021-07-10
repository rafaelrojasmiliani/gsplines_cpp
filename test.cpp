#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <vector>

#include <memory>

struct A {

  virtual int function() { return 0; }
};

struct B : public A {

  virtual int function() override = 0;
};

struct C : public B {

  std::vector<std::unique_ptr<A>> array_;

  C() { array_.push_back(std::unique_ptr<A>(this)); }

  int function() override { return 1; };
};
int main() {

  int a = 1;

  std::unique_ptr<int> m;

  return 0;
}
