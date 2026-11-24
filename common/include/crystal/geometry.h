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

/*! \file geometry.h
 *  \brief Small 3-vector and 3x3-matrix helpers for crystal geometry.
 *
 * The crystal library works in Angstrom with plain doubles: lattice vectors,
 * Cartesian positions and distances are in Angstrom (members and functions
 * say so with an `_A` suffix or in their comment), fractional coordinates
 * and directions are dimensionless. The public setters of the scenario take
 * Boost.Units quantities and convert at that boundary (docs/CRYSTAL_DESIGN.md).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace amcset {
namespace crystal {

//! A 3-vector (a direction, fractional coordinates, or a position in Angstrom).
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  constexpr double operator[](int i) const noexcept { return i == 0 ? x : i == 1 ? y : z; }
  constexpr double &operator[](int i) noexcept { return i == 0 ? x : i == 1 ? y : z; }
  constexpr Vec3 operator+(const Vec3 &o) const noexcept { return {x + o.x, y + o.y, z + o.z}; }
  constexpr Vec3 operator-(const Vec3 &o) const noexcept { return {x - o.x, y - o.y, z - o.z}; }
  constexpr Vec3 operator-() const noexcept { return {-x, -y, -z}; }
  constexpr Vec3 operator*(double s) const noexcept { return {x * s, y * s, z * s}; }
  constexpr Vec3 operator/(double s) const noexcept { return {x / s, y / s, z / s}; }
  constexpr Vec3 &operator+=(const Vec3 &o) noexcept {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
  constexpr bool operator==(const Vec3 &o) const noexcept = default;
};

constexpr Vec3 operator*(double s, const Vec3 &v) noexcept { return v * s; }
constexpr double dot(const Vec3 &a, const Vec3 &b) noexcept { return a.x * b.x + a.y * b.y + a.z * b.z; }
constexpr Vec3 cross(const Vec3 &a, const Vec3 &b) noexcept {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3 &v) noexcept { return std::sqrt(dot(v, v)); }
//! \p v / |v|; the zero vector stays zero.
inline Vec3 normalized(const Vec3 &v) noexcept {
  const double n = norm(v);
  return n > 0.0 ? v / n : v;
}

/*!
 * \brief A 3x3 matrix, stored by rows.
 *
 * For a lattice the rows are the lattice vectors a, b, c (Angstrom); for a
 * rotation crystal -> laboratory, v_lab = R * v_crystal.
 */
struct Mat3 {
  std::array<Vec3, 3> row{Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};

  static constexpr Mat3 identity() noexcept { return Mat3{}; }
  static constexpr Mat3 from_rows(const Vec3 &a, const Vec3 &b, const Vec3 &c) noexcept {
    Mat3 m;
    m.row = {a, b, c};
    return m;
  }
  constexpr const Vec3 &operator[](int i) const noexcept { return row[static_cast<std::size_t>(i)]; }
  constexpr Vec3 &operator[](int i) noexcept { return row[static_cast<std::size_t>(i)]; }
  constexpr Vec3 column(int j) const noexcept { return {row[0][j], row[1][j], row[2][j]}; }
  //! M * v.
  constexpr Vec3 operator*(const Vec3 &v) const noexcept {
    return {dot(row[0], v), dot(row[1], v), dot(row[2], v)};
  }
  //! M * N.
  constexpr Mat3 operator*(const Mat3 &n) const noexcept {
    Mat3 m;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) m[i][j] = dot(row[static_cast<std::size_t>(i)], n.column(j));
    }
    return m;
  }
  constexpr Mat3 transposed() const noexcept { return from_rows(column(0), column(1), column(2)); }
  constexpr double determinant() const noexcept { return dot(row[0], cross(row[1], row[2])); }
  //! The inverse (the caller checks that the determinant is not zero).
  Mat3 inverse() const noexcept {
    const double d = determinant();
    // Columns of the inverse are the cross products of the rows / det.
    const Vec3 c0 = cross(row[1], row[2]) / d;
    const Vec3 c1 = cross(row[2], row[0]) / d;
    const Vec3 c2 = cross(row[0], row[1]) / d;
    return from_rows(c0, c1, c2).transposed();
  }
};

/*!
 * \brief Rotation by \p angle_rad about the unit vector \p axis
 * (right-handed: counterclockwise when the axis points at the viewer).
 */
inline Mat3 rotation_about(const Vec3 &axis, double angle_rad) noexcept {
  const Vec3 u = normalized(axis);
  const double c = std::cos(angle_rad);
  const double s = std::sin(angle_rad);
  const double t = 1.0 - c;
  return Mat3::from_rows({c + u.x * u.x * t, u.x * u.y * t - u.z * s, u.x * u.z * t + u.y * s},
                         {u.y * u.x * t + u.z * s, c + u.y * u.y * t, u.y * u.z * t - u.x * s},
                         {u.z * u.x * t - u.y * s, u.z * u.y * t + u.x * s, c + u.z * u.z * t});
}

//! Largest |M M^T - I| element and |det M - 1|: how far \p m is from a rotation.
inline double rotation_error(const Mat3 &m) noexcept {
  const Mat3 p = m * m.transposed();
  double e = std::abs(m.determinant() - 1.0);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) e = std::max(e, std::abs(p[i][j] - (i == j ? 1.0 : 0.0)));
  }
  return e;
}

} // namespace crystal
} // namespace amcset
