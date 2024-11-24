#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <iomanip>

double findS(double x1, double y1, double r1, double x2, double y2, double r2, double x3, double y3, double r3,
             long long n, double k) {
  double leftBorderVnesh = std::min(x1 - r1, std::min(x2 - r2, x3 - r3));
  double rightBorderVnesh = std::max(x1 + r1, std::max(x2 + r2, x3 + r3));
  double upBorderVnesh = std::max(y1 + r1, std::max(y2 + r2, y3 + r3));
  double downBorderVnesh = std::min(y1 - r1, std::min(y2 - r2, y3 - r3));
  double leftBorderVnutr = std::max(x1 - r1, std::max(x2 - r2, x3 - r3));
  double rightBorderVnutr = std::min(x1 + r1, std::min(x2 + r2, x3 + r3));
  double upBorderVnutr = std::min(y1 + r1, std::min(y2 + r2, y3 + r3));
  double downBorderVnutr = std::max(y1 - r1, std::max(y2 - r2, y3 - r3));
  double leftBorder = leftBorderVnutr - (leftBorderVnutr - leftBorderVnesh) * k;
  double rightBorder = rightBorderVnutr + (rightBorderVnesh - rightBorderVnutr) * k;
  double upBorder = upBorderVnutr + (upBorderVnesh - upBorderVnutr) * k;
  double downBorder = downBorderVnutr - (downBorderVnutr - downBorderVnesh) * k;
  double srec = (upBorder - downBorder) * (rightBorder - leftBorder);

  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<> distrx(leftBorder, rightBorder);
  std::uniform_real_distribution<> distry(downBorder, upBorder);

  int m = 0;
  for (long long i = 0; i < n; ++i) {
    double x = distrx(generator);
    double y = distry(generator);
    if ((x - x1) * (x - x1) + (y - y1) * (y - y1) <= r1 * r1 && (x - x2) * (x - x2) + (y - y2) * (y - y2) <= r2 * r2 &&
        (x - x3) * (x - x3) + (y - y3) * (y - y3) <= r3 * r3) {
      m++;
    }
  }
  return srec * m / n;
}

template <typename T>
void printVector(const std::vector<T>& vec) {
  for (int i = 0; i < vec.size(); ++i) {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl<<std::endl;
}

int main() {
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  double x1 = 1;
  double y1 = 1;
  double r1 = 1;
  double x2 = 1.5;
  double y2 = 2;
  double r2 = std::sqrt(5) / 2;
  double x3 = 2;
  double y3 = 1.5;
  double r3 = std::sqrt(5) / 2;

  double exact = 0.25 * M_PI + 1.25 * std::asin(0.8) - 1;

  std::vector<long long> ns;
  std::vector<double> results1;
  std::vector<double> deviations1;
  std::vector<double> ks;
  std::vector<double> results2;
  std::vector<double> deviations2;

  for (long long n = 100; n <= 100000; n += 500) {
    ns.push_back(n);
    double estimatedArea = findS(x1, y1, r1, x2, y2, r2, x3, y3, r3, n, 0);
    results1.push_back(estimatedArea);
    double deviation = std::abs((estimatedArea - exact) / exact) * 100.0;  // В процентах
    deviations1.push_back(deviation);
   }
  for (double k = 0; k <= 1; k += 0.01) {
    ks.push_back(k);
    double estimatedArea = findS(x1, y1, r1, x2, y2, r2, x3, y3, r3, 100000, k);
    results2.push_back(estimatedArea);
    double deviation = std::abs((estimatedArea - exact) / exact) * 100.0;
    deviations2.push_back(deviation);
   }
  printVector(ns);
  printVector(results1);
  printVector(deviations1);
  printVector(ks);
  printVector(results2);
  printVector(deviations2);
  return 0;
}
