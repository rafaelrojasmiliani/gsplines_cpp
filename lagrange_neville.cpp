
#include <vector>
class LagrangeBase{
struct PointClass {
  double x, y;
};

std::size_t dim_;
std::vector<double> points_;
PointClass point_list_;
LagrangeBase(std::vector<double> &_points):dim_(_points.size()), points_(_points){

}
double neville_interpolation(std::vector<PointClass> pt, const double t) {
  const int n = pt.size();
  double total;
  for (int j = 1; j < n; j++) {
    for (int i = n - 1; i >= j; i--) {
      pt[i].y = ((t - pt[i - j].x) * pt[i].y - (t - pt[i].x) * pt[i - 1].y) /
                (pt[i].x - pt[i - j].x);
    }
  }
  total = pt[n - 1].y;
  return (total);
}
};
