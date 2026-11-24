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

/*! \file random.h
 *  \brief Random-number generation for SRIM mode.
 *
 * TRIM uses P. L'Ecuyer's combined multiplicative congruential generator
 * ("Efficient and Portable Combined Random Number Generators", Comm. ACM 31,
 * 742 (1988)). AMCSET uses the same generator, but gives every incident ion
 * two independent streams (one for the ion, one for its recoils) derived
 * from the user seed and the ion number. Results are therefore reproducible
 * for a given seed and independent of the number of threads, and the ion
 * trajectories do not depend on how the recoils are treated (quick vs. full
 * cascade), exactly as in TRIM.
 */

#pragma once

#include <cmath>
#include <cstdint>

namespace amcset {
namespace trim {

/*!
 * \brief L'Ecuyer (1988) combined MLCG with period ~2.3e18.
 *
 * Returns uniform deviates in the open interval (0, 1).
 */
class LecuyerRandom {
public:
  static constexpr std::int64_t m1 = 2147483563;
  static constexpr std::int64_t m2 = 2147483399;
  static constexpr std::int64_t a1 = 40014;
  static constexpr std::int64_t a2 = 40692;

  //! Seed both components (values are reduced into their valid ranges).
  explicit LecuyerRandom(std::uint64_t s1 = 12345, std::uint64_t s2 = 67890) {
    seed(s1, s2);
  }

  void seed(std::uint64_t s1, std::uint64_t s2) {
    s1_ = static_cast<std::int64_t>(s1 % static_cast<std::uint64_t>(m1 - 1)) + 1;
    s2_ = static_cast<std::int64_t>(s2 % static_cast<std::uint64_t>(m2 - 1)) + 1;
  }

  /*!
   * \brief (a s) mod m for m = 2^31 - c without a division: with
   * a s = h 2^31 + l and 2^31 = c (mod m), (a s) mod m = (l + c h) mod m,
   * and l + c h < 2 m (a < 2^16, s < 2^31), so one conditional subtraction
   * remains. The same value as (a * s) % m for 0 <= s < 2^31
   * (TrimRandom.ModularReductionMatchesTheRemainder).
   */
  static constexpr std::int64_t mul_mod(std::int64_t a, std::int64_t s, std::int64_t m) noexcept {
    const std::int64_t x = a * s;
    const std::int64_t c = (std::int64_t{1} << 31) - m;
    std::int64_t r = (x & 0x7fffffff) + c * (x >> 31);
    if (r >= m) r -= m;
    return r;
  }

  //! Next uniform deviate in (0, 1).
  double operator()() noexcept {
    s1_ = mul_mod(a1, s1_, m1);
    s2_ = mul_mod(a2, s2_, m2);
    std::int64_t z = s1_ - s2_;
    if (z < 1) z += m1 - 1;
    return static_cast<double>(z) * (1.0 / static_cast<double>(m1));
  }

  //! Standard normal deviate (Box-Muller; consumes two uniforms).
  double gaussian() noexcept {
    const double u1 = (*this)();
    const double u2 = (*this)();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
  }

  std::int64_t state1() const noexcept { return s1_; }
  std::int64_t state2() const noexcept { return s2_; }

private:
  std::int64_t s1_ = 1;
  std::int64_t s2_ = 1;
};

//! SplitMix64 mixing function (Steele, Lea & Flood 2014).
constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}

/*!
 * \brief The random stream \p stream (0 = ion, 1 = recoils) of ion number
 * \p ion for the user seed \p seed.
 */
inline LecuyerRandom make_stream(std::uint32_t seed, std::uint64_t ion,
                                 unsigned stream) {
  const std::uint64_t base =
      splitmix64(static_cast<std::uint64_t>(seed) * 0x100000001B3ULL ^
                 splitmix64(ion * 4 + stream));
  return LecuyerRandom(splitmix64(base), splitmix64(base ^ 0xD1B54A32D192ED03ULL));
}

} // namespace trim
} // namespace amcset
