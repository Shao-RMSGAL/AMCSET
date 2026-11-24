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

/*! \file fast_math_test.cpp
 *  \brief The exact and approximate fast paths of the scattering code:
 *  shared exponentials (bit-identical), the kinematics from the half angle,
 *  MagicTable against the magic formula, the fast sine/cosine, and the
 *  bit-indexed tables (OctaveTable: stopping power, Kinchin-Pease partition,
 *  flight path, log10).
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstring>
#include <random>
#include <string>

#include "amcset_units.h"
#include "physics/fast_math.h"
#include "physics/scattering.h"
#include "physics/stopping.h"

using namespace amcset;
using namespace amcset::common;
using namespace amcset::physics;

namespace {

bool same_bits(double a, double b) { return std::memcmp(&a, &b, sizeof a) == 0; }

} // namespace

// The exact transport's speed-up rests on this: phi and dphi from one set of
// exponentials must be the very same numbers as from phi() and dphi().
TEST(Scattering, PhiAndDerivativeAreBitIdenticalToSeparateCalls) {
  const ScreeningFunction f = ScreeningFunction::universal();
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int i = 0; i < 20000; ++i) {
    const double x = 200.0 * u(rng) * u(rng);
    double p = 0.0;
    double dp = 0.0;
    f.phi_and_derivative(x, p, dp);
    ASSERT_TRUE(same_bits(p, f.phi(x))) << "x = " << x;
    ASSERT_TRUE(same_bits(dp, f.dphi(x))) << "x = " << x;
  }
}

TEST(Scattering, KinematicsFromTheHalfAngleMatchSin2Form) {
  const CollisionPair pair(26, 55.935 * atomic_mass_unit, 26, 55.847 * atomic_mass_unit);
  const CollisionPair light(1, 1.008 * atomic_mass_unit, 74, 183.84 * atomic_mass_unit);
  const auto e = 1000.0 * electron_volt;
  for (const auto *p : {&pair, &light}) {
    for (double s2 : {1e-9, 1e-4, 0.03, 0.4, 0.7, 0.999, 1.0}) {
      const auto a = collision_kinematics(*p, e, s2);
      const auto b = collision_kinematics_half_angle(*p, e, std::sqrt(s2), std::sqrt(1.0 - s2));
      EXPECT_NEAR(in_units(b.transferred, electron_volt), in_units(a.transferred, electron_volt),
                  1e-12 * in_units(e, electron_volt));
      EXPECT_NEAR(b.cos_theta_cm, a.cos_theta_cm, 1e-12);
      EXPECT_NEAR(b.cos_projectile, a.cos_projectile, 1e-12);
      EXPECT_NEAR(b.sin_projectile, a.sin_projectile, 1e-12);
      EXPECT_NEAR(b.cos_recoil, a.cos_recoil, 1e-12);
      EXPECT_NEAR(b.sin_recoil, a.sin_recoil, 1e-12);
    }
  }
}

// MagicTable stands in for the magic formula with TrimConfig::fast_math. The
// documented accuracy (relative error of sin^2 and cos^2 of theta/2 below
// 4e-5 wherever the value is above 1e-8; measured 3.5e-5 and 1e-6, recorded
// as test properties) over the energies and impact parameters of a
// transport (log-uniform eps, area-uniform b).
TEST(Scattering, MagicTableReproducesTheMagicFormula) {
  const ScreeningFunction f = ScreeningFunction::universal();
  const MagicTable &table = MagicTable::universal();
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  double worst_s = 0.0;
  double worst_c = 0.0;
  std::string where_s;
  std::string where_c;
  int inside = 0;
  for (int i = 0; i < 20000; ++i) {
    const double eps = std::exp(std::log(1e-8) + u(rng) * (std::log(10.0) - std::log(1e-8)));
    // Solids (b up to ~30) and, for every fourth sample, gases (up to 250).
    const double b = (i % 4 == 0 ? 250.0 : 30.0) * std::sqrt(u(rng));
    double sh = 0.0;
    double ch = 1.0;
    if (!table.half_angle(eps, b, sh, ch)) continue;
    ++inside;
    const double s2 = magic_sin2_half_theta(f, eps, b);
    const double c2 = 1.0 - s2;
    const auto at = [&] { return "eps " + std::to_string(eps) + " b " + std::to_string(b); };
    if (s2 > 1e-8 && std::abs(sh * sh - s2) / s2 > worst_s) {
      worst_s = std::abs(sh * sh - s2) / s2;
      where_s = at();
    }
    if (c2 > 1e-6 && std::abs(ch * ch - c2) / c2 > worst_c) {
      worst_c = std::abs(ch * ch - c2) / c2;
      where_c = at() + " c2 " + std::to_string(c2) + " table " + std::to_string(ch * ch);
    }
  }
  EXPECT_EQ(inside, 20000);
  RecordProperty("worst_sin2", std::to_string(worst_s));
  RecordProperty("worst_cos2", std::to_string(worst_c));
  EXPECT_LT(worst_s, 5e-5) << where_s;
  EXPECT_LT(worst_c, 5e-5) << where_c;
}

TEST(Scattering, MagicTableDeclinesOutsideItsRange) {
  const MagicTable &table = MagicTable::universal();
  double sh = -1.0;
  double ch = -1.0;
  EXPECT_FALSE(table.half_angle(1e-12, 1.0, sh, ch));
  EXPECT_FALSE(table.half_angle(1.0, 1000.0, sh, ch));
  EXPECT_FALSE(table.half_angle(std::nan(""), 1.0, sh, ch));
  EXPECT_EQ(sh, -1.0); // outputs untouched
  // fast_half_angle then falls back to the formula itself.
  const CollisionPair pair(26, 55.935 * atomic_mass_unit, 26, 55.847 * atomic_mass_unit);
  fast_half_angle(pair, 1.0, 1000.0, ScatteringMethod::magic, sh, ch);
  const double s2 = trim_sin2_half_theta(pair.screening, 1.0, 1000.0);
  EXPECT_NEAR(sh * sh, s2, 1e-15);
  // Screened Rutherford above eps = 10, with both halves precise.
  fast_half_angle(pair, 50.0, 0.3, ScatteringMethod::magic, sh, ch);
  EXPECT_NEAR(sh * sh, rutherford_sin2_half_theta(50.0, 0.3), 1e-14);
  EXPECT_NEAR(sh * sh + ch * ch, 1.0, 1e-14);
}

TEST(FastMath, SinCosOfAQuarterAndOfAFullTurn) {
  double worst = 0.0;
  for (int i = 0; i <= 100000; ++i) {
    const double x = 1.5707963267948966 * i / 100000.0;
    double s = 0.0;
    double c = 0.0;
    fast::sincos_quarter(x, s, c);
    worst = std::max(worst, std::abs(s - std::sin(x)));
    worst = std::max(worst, std::abs(c - std::cos(x)));
  }
  EXPECT_LT(worst, 1e-11);
  worst = 0.0;
  for (int i = 0; i < 100000; ++i) {
    const double u = i / 100000.0;
    double c = 0.0;
    double s = 0.0;
    fast::sincos_turn(u, c, s);
    worst = std::max(worst, std::abs(c - std::cos(6.283185307179586 * u)));
    worst = std::max(worst, std::abs(s - std::sin(6.283185307179586 * u)));
  }
  EXPECT_LT(worst, 2e-11);
  // Relative precision near the ends of the quarter turn.
  double s = 0.0;
  double c = 0.0;
  fast::sincos_quarter(1e-9, s, c);
  EXPECT_NEAR(s / 1e-9, 1.0, 1e-15);
  fast::sincos_quarter(1.5707963267948966 - 1e-9, s, c);
  EXPECT_NEAR(c / std::cos(1.5707963267948966 - 1e-9), 1.0, 1e-7);
}

// OctaveTable: the documented error 7.6e-6 |b (b - 1)| for powers x^b, and
// no answer outside the table, for x <= 0, infinities and NaN.
TEST(FastMath, OctaveTableOfPowerLaws) {
  std::mt19937_64 rng(5);
  std::uniform_real_distribution<double> u(-30.0, 30.0);
  for (const double b : {0.5, -0.8, 1.0 / 6.0, 0.75, 2.0}) {
    const fast::OctaveTable t(-40, 40, [b](double x) { return std::pow(x, b); });
    EXPECT_EQ(t.size(), 80u * 128u);
    double worst = 0.0;
    for (int i = 0; i < 200000; ++i) {
      const double x = std::exp2(u(rng));
      double y = 0.0;
      ASSERT_TRUE(t.lookup(x, y));
      worst = std::max(worst, std::abs(y / std::pow(x, b) - 1.0));
    }
    EXPECT_LT(worst, 7.7e-6 * std::abs(b * (b - 1.0)) + 1e-15) << "b " << b;
    // At a bin edge the table is exact.
    double y = 0.0;
    ASSERT_TRUE(t.lookup(1.0, y));
    EXPECT_EQ(y, 1.0);
    double untouched = -7.0;
    for (const double x : {0.0, -1.0, std::ldexp(1.0, -41), std::ldexp(1.0, 40), HUGE_VAL, std::nan("")}) {
      EXPECT_FALSE(t.lookup(x, untouched)) << x;
    }
    EXPECT_EQ(untouched, -7.0);
  }
}

TEST(FastMath, TabulatedTrimFunctionsMatchTheirFormulas) {
  std::mt19937_64 rng(9);
  std::uniform_real_distribution<double> u(-12.0, 12.0);
  double worst_g = 0.0;
  double worst_h = 0.0;
  double worst_l = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const double x = std::pow(10.0, u(rng));
    const double g = x + 0.4 * std::pow(x, 0.75) + 3.4 * std::pow(x, 1.0 / 6.0);
    worst_g = std::max(worst_g, std::abs(fast::kp_partition(x) / g - 1.0));
    const double e = std::sqrt(x);
    const double h = e + std::sqrt(e) + 0.125 * std::pow(e, 0.1);
    worst_h = std::max(worst_h, std::abs(fast::flight_path_denominator(x) / h - 1.0));
    const double big = std::pow(10.0, std::abs(u(rng)) * 1.25); // 1 to 1e15
    worst_l = std::max(worst_l, std::abs(fast::log10_at_least_one(big) - std::log10(big)));
  }
  EXPECT_LT(worst_g, 2e-6);
  EXPECT_LT(worst_h, 3e-6);
  EXPECT_LT(worst_l, 5e-6);
  EXPECT_EQ(fast::log10_at_least_one(1.0), 0.0); // the spectra's first bin starts there
  // Outside the tables: the formulas themselves.
  EXPECT_EQ(fast::kp_partition(1e-20), 1e-20 + 0.4 * std::pow(1e-20, 0.75) + 3.4 * std::pow(1e-20, 1.0 / 6.0));
  EXPECT_EQ(fast::log10_at_least_one(1e30), 30.0);
}

TEST(FastMath, TabulatedStoppingPowerMatchesTheCurve) {
  using namespace amcset::common;
  const std::vector<StoppingComponent> fe = {{26, 1.0}};
  const number_density_quantity n = 0.0849 * per_cubic_angstrom;
  StoppingCurve curve(ElectronicStoppingModel::builtin, nullptr, 26, 55.935 * atomic_mass_unit, fe, false, 1.0, n);
  const auto e = 123456.0 * electron_volt;
  // Before enable_fast_lookup() the fast lookup is the curve itself.
  EXPECT_EQ(in_units(curve.stopping_power_fast(e), electron_volt_per_angstrom),
            in_units(curve.stopping_power(e), electron_volt_per_angstrom));
  curve.enable_fast_lookup();
  std::mt19937_64 rng(3);
  std::uniform_real_distribution<double> u(-3.0, 11.0);
  double worst = 0.0;
  for (int i = 0; i < 200000; ++i) {
    const auto en = std::pow(10.0, u(rng)) * electron_volt;
    const double exact = in_units(curve.stopping_power(en), electron_volt_per_angstrom);
    const double fast_value = in_units(curve.stopping_power_fast(en), electron_volt_per_angstrom);
    worst = std::max(worst, std::abs(fast_value / exact - 1.0));
  }
  // (Bins that straddle a node of the curve, where its log-log slope
  // changes, are the least accurate.)
  EXPECT_LT(worst, 4e-5);
  EXPECT_EQ(in_units(curve.stopping_power_fast(0.0 * electron_volt), electron_volt_per_angstrom), 0.0);
}
