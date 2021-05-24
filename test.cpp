
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

class SumOfFunctions{
  public:
    SumOfFunctions(std::vector<std::unique_ptr<int>> array){
        for(const auto &i:array){
            function_array_.push_back(std::unique_ptr<int>(new int(*i)));
        }
    }
    SumOfFunctions(const SumOfFunctions &that);
  std::unique_ptr<int> clone() const {return nullptr;}
  std::unique_ptr<int> deriv(int _deg);
  private:
    std::vector<std::unique_ptr<int>> function_array_;


};

int main() {
std::vector<std::unique_ptr<int>> array;
int a,b,c;
array.push_back(std::unique_ptr<int>(new int(a)));
array.push_back(std::unique_ptr<int>(new int(b)));
SumOfFunctions k(array);
  return 0;
}
