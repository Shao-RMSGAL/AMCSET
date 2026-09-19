// Quaternion implementations moved out of amcset_common.cpp
#include "amcset_common.h"

#include <cmath>

namespace amcset {
namespace common {

Quaternion &Quaternion::operator+=(const Quaternion &b) {
  r_ += b.r_;
  x_ += b.x_;
  y_ += b.y_;
  z_ += b.z_;
  return *this;
}

Quaternion Quaternion::operator+(const Quaternion &b) const {
  return Quaternion(*this) += b;
}

Quaternion operator+(const Quaternion &a, const Quaternion &b) {
  return Quaternion(a) += b;
}

Quaternion &Quaternion::operator*=(const Quaternion &b) {
  auto r_1 = r_;
  auto r_2 = b.r_;
  auto x_1 = x_;
  auto y_1 = y_;
  auto z_1 = z_;
  auto x_2 = b.x_;
  auto y_2 = b.y_;
  auto z_2 = b.z_;

  r_ = r_1 * r_2 - x_1 * x_2 - y_1 * y_2 - z_1 * z_2;
  x_ = r_1 * x_2 + r_2 * x_1 + y_1 * z_2 - z_1 * y_2;
  y_ = r_1 * y_2 + r_2 * y_1 + z_1 * x_2 - x_1 * z_2;
  z_ = r_1 * z_2 + r_2 * z_1 + x_1 * y_2 - y_1 * x_2;
  return *this;
}

Quaternion Quaternion::operator*(const Quaternion &b) const {
  return Quaternion(*this) *= b;
}

Quaternion &Quaternion::operator*=(double a) {
  r_ *= a;
  x_ *= a;
  y_ *= a;
  z_ *= a;
  return *this;
}

Quaternion Quaternion::operator*(const double b) const {
  return Quaternion(*this) *= b;
}

Quaternion operator*(double a, const Quaternion &b) {
  return Quaternion(b) *= a;
}

Quaternion &Quaternion::operator/=(double a) {
  r_ /= a;
  x_ /= a;
  y_ /= a;
  z_ /= a;
  return *this;
}

Quaternion Quaternion::operator/(const double b) const {
  return Quaternion(*this) /= b;
}

Quaternion operator/(const Quaternion &a, double b) {
  return Quaternion(a) /= b;
}

Quaternion &Quaternion::normalize_inplace() {
  *this /= (this->magnitude());
  return *this;
}

Quaternion &Quaternion::conjugate_inplace() {
  x_ = -x_;
  y_ = -y_;
  z_ = -z_;
  return *this;
}

Quaternion Quaternion::conjugate_copy() const {
  return Quaternion(r_, -x_, -y_, -z_);
}

double Quaternion::magnitude() const {
  return std::sqrt(r_ * r_ + x_ * x_ + y_ * y_ + z_ * z_);
}

} // namespace common
} // namespace amcset
