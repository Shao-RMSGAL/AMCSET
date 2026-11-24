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

#include "trim/single_precision.h"

#include <algorithm>
#include <cmath>

namespace amcset {
namespace trim {

namespace {

// One single-precision operation: the exact sum of two floats fits a double
// here (their exponents differ by far less than 29), so this rounds once.
float add(float s, float d) { return static_cast<float>(static_cast<double>(s) + static_cast<double>(d)); }

} // namespace

double single_precision_sum(std::uint64_t increments, std::uint64_t decrements, float step) {
  const std::uint64_t total = increments + decrements;
  if (total == 0 || !(step > 0.0f)) return 0.0;
  const double ratio = static_cast<double>(decrements) / static_cast<double>(total);
  // Decrements among the first k operations (Bresenham-like spreading).
  auto before = [&](std::uint64_t k) -> std::uint64_t {
    if (k >= total) return decrements;
    return std::min<std::uint64_t>(decrements, static_cast<std::uint64_t>(std::floor(static_cast<double>(k) * ratio)));
  };
  float s = 0.0f;
  std::uint64_t i = 0;
  auto single = [&](std::uint64_t k) {
    const bool down = before(k + 1) > before(k);
    s = add(s, down ? -step : step);
  };
  while (i < total) {
    std::uint64_t chunk = 0;
    if (s > 64.0f * step) {
      int e = 0;
      std::frexp(s, &e); // s in [2^(e-1), 2^e)
      const double lo = std::ldexp(1.0, e - 1);
      const double hi = std::ldexp(1.0, e);
      const double ulp = std::ldexp(1.0, e - 24);
      const double up = static_cast<double>(add(s, step)) - s;
      const double dn = s - static_cast<double>(add(s, -step));
      // A step of exactly half an ulp beyond a multiple rounds to even and
      // alternates: no fast path then.
      const double rem = std::fmod(static_cast<double>(step), ulp);
      const bool tie = rem == 0.5 * ulp;
      if (!tie && up > 0.0 && dn > 0.0) {
        const double room_up = hi - static_cast<double>(step) - 2.0 * ulp - s;
        const double room_dn = s - lo - static_cast<double>(step) - 2.0 * ulp;
        if (room_up > 0.0 && room_dn > 0.0) {
          const double c = std::min(std::floor(room_up / up), std::floor(room_dn / dn));
          if (c >= 16.0) chunk = std::min<std::uint64_t>(total - i, static_cast<std::uint64_t>(c));
        }
      } else if (!tie && up == 0.0 && dn == 0.0) {
        chunk = total - i; // the step is lost against the sum: it stays put
      }
      if (chunk > 0) {
        const std::uint64_t d = before(i + chunk) - before(i);
        const std::uint64_t u = chunk - d;
        s = static_cast<float>(s + static_cast<double>(u) * up - static_cast<double>(d) * dn);
        i += chunk;
        continue;
      }
    }
    single(i);
    ++i;
  }
  return s;
}

double trim_table_count(double count, double removed, double bin_width_A) {
  if (!(bin_width_A > 0.0)) return count;
  const double r = std::max(0.0, std::round(removed));
  const double n = std::round(count + r);
  if (!(n >= 0.0) || n > 9.0e15) return count;
  // TRIM: bin width = (window end - start) / 100 and weight = 1 / width, as
  // Singles.
  const float width = static_cast<float>(bin_width_A);
  const float step = static_cast<float>(1.0 / static_cast<double>(width));
  const double sum = single_precision_sum(static_cast<std::uint64_t>(n), static_cast<std::uint64_t>(r), step);
  return sum * static_cast<double>(width);
}

} // namespace trim
} // namespace amcset
