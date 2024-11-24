#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

double findS(double x1, double y1, double r1, double x2, double y2, double r2, double x3, double y3, double r3) {
  double leftBorder = std::max(x1 - r1, std::max(x2 - r2, x3 - r3));
  double rightBorder = std::min(x1 + r1, std::min(x2 + r2, x3 + r3));
  double upBorder = std::min(y1 + r1, std::min(y2 + r2, y3 + r3));
  double downBorder = std::max(y1 - r1, std::max(y2 - r2, y3 - r3));
  double srec = (upBorder - downBorder) * (rightBorder - leftBorder);
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<> distrx(leftBorder, rightBorder);
  std::uniform_real_distribution<> distry(downBorder, upBorder);
  int m = 0;
  long long n = 15000000;
  for (int i = 0; i < n; ++i) {
    double x = distrx(generator);
    double y = distry(generator);
    if ((x - x1) * (x - x1) + (y - y1) * (y - y1) <= r1 * r1 && (x - x2) * (x - x2) + (y - y2) * (y - y2) <= r2 * r2 &&
        (x - x3) * (x - x3) + (y - y3) * (y - y3) <= r3 * r3) {
      m++;
    }
  }
  return srec * m / n;
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);
  double x1 = 0;
  double y1 = 0;
  double r1 = 0;
  double x2 = 0;
  double y2 = 0;
  double r2 = 0;
  double x3 = 0;
  double y3 = 0;
  double r3 = 0;
  std::cin >> x1 >> y1 >> r1 >> x2 >> y2 >> r2 >> x3 >> y3 >> r3;
  double s = findS(x1, y1, r1, x2, y2, r2, x3, y3, r3);
  std::cout << s << std::endl;
  return 0;
}