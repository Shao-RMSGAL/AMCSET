// Copyright 2025, Texas A&M University
//
// This file is part of AMCSET.
//
// AMCSET is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

/*! \file quaternion.cpp
 *  \brief Quaternion arithmetic used for particle direction rotation.
 */

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

Quaternion &Quaternion::normalize_inplace() {
  *this /= (this->magnitude());
  return *this;
}

Quaternion Quaternion::normalize_copy() const {
  return Quaternion(*this).normalize_inplace();
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
