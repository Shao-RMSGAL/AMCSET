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

/*! \file fast_math.h
 *  \brief Building blocks of SRIM mode's fast approximations
 *  (TrimConfig::fast_math): inline sine and cosine, and tables of smooth
 *  functions indexed by the bits of their argument.
 *
 * std::sin and std::cos accept any argument and spend much of their time on
 * the general range reduction. The transport only needs angles in a quarter
 * turn (theta/2) or a fraction of a full turn (the azimuth), where a short
 * polynomial is accurate to 1e-11 (fast_math_test.cpp), far below the
 * tabulated scattering angle's 4e-5, and faster. OctaveTable replaces
 * pow, cbrt, log and exp of smooth one-dimensional functions (the stopping
 * power, the Kinchin-Pease partition, the flight path's energy dependence)
 * by a lookup without any logarithm. None of this is bit-identical to the
 * C runtime, so the exact (default) transport does not use it.
 * (Polynomial exp and log were tried as well and are slower than the C
 * runtime's table-driven ones, so there are none here.)
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace amcset {
namespace physics {
namespace fast {

/*!
 * \brief A smooth function of x > 0 tabulated on 2^sub_bits bins per factor
 * of two (octave) and interpolated linearly in x within a bin. The bin is
 * read from the bits of x (its exponent and top mantissa bits), so a lookup
 * needs no logarithm: a shift, a subtraction, a conversion and one
 * multiply-add. For f ~ x^b the relative error is at most about
 * |b (b - 1)| / 8 * 2^(-2 sub_bits) = 7.6e-6 |b (b - 1)|: below 2e-6 for
 * 0 <= b <= 1, 1.1e-5 for b = -0.8 (fast_math_test.cpp).
 */
class OctaveTable {
public:
  static constexpr int sub_bits = 7;

  OctaveTable() = default;
  //! Tabulate \p f on [2^lo_exponent, 2^hi_exponent) (|exponents| < 1022).
  template <class F>
  OctaveTable(int lo_exponent, int hi_exponent, F f) {
    key0_ = key(std::ldexp(1.0, lo_exponent));
    n_ = key(std::ldexp(1.0, hi_exponent)) - key0_;
    nodes_.resize(static_cast<std::size_t>(n_));
    double prev = f(edge(key0_));
    for (std::uint64_t k = 0; k < n_; ++k) {
      const double next = f(edge(key0_ + k + 1));
      nodes_[static_cast<std::size_t>(k)] = {prev, next - prev};
      prev = next;
    }
  }
  //! f(x) in \p y if x lies in the table (false for x outside it, x <= 0,
  //! infinities and NaN; y is then untouched).
  bool lookup(double x, double &y) const noexcept {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof bits);
    const std::uint64_t k = (bits >> shift) - key0_; // wraps below the table
    if (k >= n_) return false;
    const double t = static_cast<double>(bits & low_mask) * low_scale;
    const Node &c = nodes_[static_cast<std::size_t>(k)];
    y = c.v + t * c.dv;
    return true;
  }
  bool empty() const noexcept { return nodes_.empty(); }
  std::size_t size() const noexcept { return nodes_.size(); }

private:
  static constexpr int shift = 52 - sub_bits;
  static constexpr std::uint64_t low_mask = (std::uint64_t{1} << shift) - 1;
  static constexpr double low_scale = 1.0 / static_cast<double>(std::uint64_t{1} << shift);
  struct Node {
    double v = 0.0;  //!< f at the bin's lower edge
    double dv = 0.0; //!< f at the upper edge minus v
  };
  static std::uint64_t key(double x) noexcept {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof bits);
    return bits >> shift;
  }
  static double edge(std::uint64_t k) noexcept {
    const std::uint64_t bits = k << shift;
    double x = 0.0;
    std::memcpy(&x, &bits, sizeof x);
    return x;
  }
  std::uint64_t key0_ = 0;
  std::uint64_t n_ = 0;
  std::vector<Node> nodes_;
};

//! TRIM's Kinchin-Pease (Robinson) partition g(eps) = eps + 0.4 eps^(3/4) +
//! 3.4 eps^(1/6) from a table (relative error ~1e-6; the formula outside
//! 2^-40 <= eps < 2^40).
inline double kp_partition(double eps) {
  const auto g = [](double e) { return e + 0.4 * std::pow(e, 0.75) + 3.4 * std::pow(e, 1.0 / 6.0); };
  static const OctaveTable table(-40, 40, g);
  double y = 0.0;
  return table.lookup(eps, y) ? y : g(eps);
}

//! The energy dependence of TRIM's free flight path (0x44ab60): with
//! u = eps epsdg and e = sqrt(u), e + sqrt(e) + 0.125 e^0.1, from a table
//! (relative error ~1e-6; the formula outside 2^-60 <= u < 2^60).
inline double flight_path_denominator(double u) {
  const auto h = [](double v) {
    const double e = std::sqrt(v);
    return e + std::sqrt(e) + 0.125 * std::pow(e, 0.1);
  };
  static const OctaveTable table(-60, 60, h);
  double y = 0.0;
  return table.lookup(u, y) ? y : h(u);
}

//! log10(x) for x >= 1 from a table (absolute error ~3e-6; exact at powers
//! of two; std::log10 outside 1 <= x < 2^64).
inline double log10_at_least_one(double x) {
  static const OctaveTable table(0, 64, [](double v) { return std::log10(v); });
  double y = 0.0;
  return table.lookup(x, y) ? y : std::log10(x);
}

namespace detail {
//! sin y and cos y for |y| <= pi/4 (Taylor to y^11 / y^12: relative
//! error < 1e-11 for sin, < 1e-12 for cos; the scattering table that
//! feeds them is accurate to 4e-5).
inline void sincos_kernel(double y, double &s, double &c) noexcept {
  const double y2 = y * y;
  double ps = -1.0 / 39916800.0;
  ps = ps * y2 + 1.0 / 362880.0;
  ps = ps * y2 - 1.0 / 5040.0;
  ps = ps * y2 + 1.0 / 120.0;
  ps = ps * y2 - 1.0 / 6.0;
  s = y + y * y2 * ps;
  double pc = 1.0 / 479001600.0;
  pc = pc * y2 - 1.0 / 3628800.0;
  pc = pc * y2 + 1.0 / 40320.0;
  pc = pc * y2 - 1.0 / 720.0;
  pc = pc * y2 + 1.0 / 24.0;
  pc = pc * y2 - 0.5;
  c = 1.0 + y2 * pc;
}
} // namespace detail

/*!
 * \brief sin x and cos x for 0 <= x <= pi/2 (x is clamped to that range).
 * Both keep their relative precision near the ends (sin near 0, cos near
 * pi/2).
 */
inline void sincos_quarter(double x, double &s, double &c) noexcept {
  constexpr double pi_4 = 0.78539816339744830962;
  constexpr double pi_2_hi = 1.57079632679489655800; // pi/2 = hi + lo
  constexpr double pi_2_lo = 6.12323399573676603587e-17;
  if (!(x > 0.0)) x = 0.0;
  if (x > pi_2_hi) x = pi_2_hi;
  if (x <= pi_4) {
    detail::sincos_kernel(x, s, c);
  } else {
    // sin x = cos(pi/2 - x), cos x = sin(pi/2 - x).
    detail::sincos_kernel((pi_2_hi - x) + pi_2_lo, c, s);
  }
}

/*!
 * \brief cos(2 pi u) and sin(2 pi u) for a fraction of a turn 0 <= u < 1
 * (the azimuth of a collision drawn from a uniform deviate). The quadrant
 * is applied without branches (a random quadrant defeats branch
 * prediction).
 */
inline void sincos_turn(double u, double &cos_phi, double &sin_phi) noexcept {
  const double q4 = 4.0 * u;
  int q = static_cast<int>(q4);
  q = q < 0 ? 0 : (q > 3 ? 3 : q);
  const double r = q4 - static_cast<double>(q); // exact
  double s = 0.0;
  double c = 1.0;
  // r in [0, 1): r pi/2 in [0, pi/2), one polynomial branch at pi/4 (taken
  // half the time; cheaper than a second kernel evaluation).
  sincos_quarter(r * 1.57079632679489661923, s, c);
  // cos(q pi/2 + a) = c, -s, -c, s and sin(q pi/2 + a) = s, c, -s, -c.
  const double v[4] = {c, -s, -c, s};
  cos_phi = v[q];
  sin_phi = v[(q + 3) & 3];
}

} // namespace fast
} // namespace physics
} // namespace amcset
