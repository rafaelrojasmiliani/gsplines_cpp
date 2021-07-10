#include <stdio.h>
#include <stdlib.h>
#include <time.h>

struct A {

  virtual int function() { return 0; }
};

struct B : public A {

  virtual int function() override = 0;
};

struct C : public B {

  int function() override { return 1; };
};
int main() {

  C b;
  return 0;
}
