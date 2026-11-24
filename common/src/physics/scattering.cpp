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

/*! \file scattering.cpp
 *  \brief Interatomic potentials, magic formula, scattering integral and
 *  collision kinematics.
 */

#include "physics/scattering.h"

#include "physics/fast_math.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <vector>

namespace amcset {
namespace physics {

using namespace common;

// ---------------------------------------------------------------------------
// Screening functions
// ---------------------------------------------------------------------------

ScreeningFunction ScreeningFunction::universal() noexcept {
  // Ziegler, Biersack & Littmark (1985); the coefficients TRIM-85 uses.
  ScreeningFunction f;
  f.c = {0.18175, 0.50986, 0.28022, 0.028171};
  f.d = {3.1998, 0.94229, 0.4029, 0.20162};
  return f;
}

length_quantity universal_screening_length(unsigned z1, unsigned z2) {
  // 0.8854 is the Thomas-Fermi constant (9 pi^2 / 128)^(1/3).
  return 0.8854 * bohr_radius /
         (std::pow(static_cast<double>(z1), 0.23) +
          std::pow(static_cast<double>(z2), 0.23));
}

// ---------------------------------------------------------------------------
// Pair potentials (SRIM-2013 data files)
// ---------------------------------------------------------------------------

namespace {
std::pair<unsigned, unsigned> pair_key(unsigned z1, unsigned z2) {
  return {std::min(z1, z2), std::max(z1, z2)};
}
} // namespace

void PairPotentialTable::add_potential(unsigned z1, unsigned z2,
                                       const ScreeningFunction &f) {
  potentials_[pair_key(z1, z2)] = f;
}

void PairPotentialTable::add_nuclear_fit(unsigned z1, unsigned z2,
                                         const std::array<double, 4> &f) {
  nuclear_[pair_key(z1, z2)] = f;
}

std::optional<ScreeningFunction> PairPotentialTable::screening(unsigned z1,
                                                               unsigned z2) const {
  const auto it = potentials_.find(pair_key(z1, z2));
  if (it == potentials_.end()) return std::nullopt;
  return it->second;
}

std::optional<std::array<double, 4>>
PairPotentialTable::nuclear_stopping_fit(unsigned z1, unsigned z2) const {
  const auto it = nuclear_.find(pair_key(z1, z2));
  if (it == nuclear_.end()) return std::nullopt;
  return it->second;
}

std::shared_ptr<const PairPotentialTable>
PairPotentialTable::load_srim(const std::filesystem::path &srim_data_dir) {
  std::filesystem::path dir = srim_data_dir;
  if (!std::filesystem::exists(dir / "POTCOEF03.DAT") &&
      std::filesystem::exists(dir / "Data" / "POTCOEF03.DAT")) {
    dir /= "Data";
  }
  const auto potfile = dir / "POTCOEF03.DAT";
  std::ifstream pin(potfile);
  if (!pin) {
    throw std::runtime_error("Cannot open " + potfile.string() +
                             " (expected the Data directory of a SRIM-2013 installation)");
  }
  auto table = std::make_shared<PairPotentialTable>();
  std::string line;
  while (std::getline(pin, line)) {
    std::istringstream ss(line);
    unsigned z1 = 0;
    unsigned z2 = 0;
    ScreeningFunction f;
    if (!(ss >> z1 >> z2)) continue;
    bool ok = true;
    for (int i = 0; i < 4 && ok; ++i) {
      ok = static_cast<bool>(ss >> f.c[static_cast<std::size_t>(i)] >>
                             f.d[static_cast<std::size_t>(i)]);
    }
    if (!ok) continue;
    const double sum = f.c[0] + f.c[1] + f.c[2] + f.c[3];
    if (std::abs(sum - 1.0) > 1e-3) {
      throw std::runtime_error("POTCOEF03.DAT: amplitudes for Z1=" + std::to_string(z1) +
                               " Z2=" + std::to_string(z2) + " do not sum to 1");
    }
    table->add_potential(z1, z2, f);
  }
  const auto snucfile = dir / "SNUC03.dat";
  std::ifstream nin(snucfile);
  if (nin) {
    while (std::getline(nin, line)) {
      std::istringstream ss(line);
      unsigned z1 = 0;
      unsigned z2 = 0;
      std::array<double, 4> f{};
      if (ss >> z1 >> z2 >> f[0] >> f[1] >> f[2] >> f[3]) {
        table->add_nuclear_fit(z1, z2, f);
      }
    }
  }
  if (table->size() == 0) {
    throw std::runtime_error("No potentials found in " + potfile.string());
  }
  return table;
}

double reduced_nuclear_stopping(double epsilon, const std::array<double, 4> *fit) {
  if (epsilon <= 0.0) return 0.0;
  if (epsilon > 30.0) return std::log(epsilon) / (2.0 * epsilon);
  if (fit != nullptr) {
    const auto &f = *fit;
    return std::log(1.0 + f[0] * epsilon) /
           (2.0 * (epsilon + f[1] * std::pow(epsilon, f[2]) + f[3] * std::sqrt(epsilon)));
  }
  return std::log(1.0 + 1.1383 * epsilon) /
         (2.0 * (epsilon + 0.01321 * std::pow(epsilon, 0.21226) +
                 0.19593 * std::sqrt(epsilon)));
}

// ---------------------------------------------------------------------------
// CollisionPair
// ---------------------------------------------------------------------------

CollisionPair::CollisionPair(unsigned z1_, mass_quantity m1_, unsigned z2_,
                             mass_quantity m2_,
                             const PairPotentialTable *potentials)
    : z1(z1_), z2(z2_), m1(m1_), m2(m2_),
      screening_length(universal_screening_length(z1_, z2_)) {
  if (z1 == 0 || z2 == 0) throw std::invalid_argument("CollisionPair: Z must be positive");
  const double Z1 = static_cast<double>(z1);
  const double Z2 = static_cast<double>(z2);
  // eps = a E_cm / (Z1 Z2 e^2) with E_cm = E M2 / (M1 + M2).
  reduced_energy_per_energy = screening_length * (m2 / (m1 + m2)).value() /
                              (Z1 * Z2 * coulomb_energy_length);
  mass_ratio = (m1 / m2).value();
  gamma = (4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2))).value();
  if (potentials != nullptr) {
    if (auto f = potentials->screening(z1, z2)) {
      screening = *f;
      universal_screening = false;
    }
    nuclear_fit = potentials->nuclear_stopping_fit(z1, z2);
  }
}

stopping_cross_section_quantity
CollisionPair::nuclear_stopping(energy_quantity e) const {
  const double eps = reduced_energy(e);
  const double sn =
      reduced_nuclear_stopping(eps, nuclear_fit ? &*nuclear_fit : nullptr);
  // S_n = 4 pi a Z1 Z2 e^2 M1/(M1+M2) s_n(eps)
  return 4.0 * pi_v * screening_length * static_cast<double>(z1) *
         static_cast<double>(z2) * coulomb_energy_length *
         (m1 / (m1 + m2)).value() * sn;
}

// ---------------------------------------------------------------------------
// Scattering angle
// ---------------------------------------------------------------------------

double closest_approach(const ScreeningFunction &f, double eps, double b) {
  // Root of g(R) = 1 - phi(R)/(eps R) - b^2/R^2, which increases
  // monotonically in R for a repulsive potential.
  // TRIM-85 starting value.
  double r = b;
  if (eps * b > 0.0) {
    double rr = -2.7 * std::log(eps * b);
    if (rr >= b) {
      rr = -2.7 * std::log(eps * rr);
      if (rr >= b) r = rr;
    }
  }
  if (!(r > 0.0)) r = std::max(b, 1e-6);
  // Bracket for the bisection fallback.
  double lo = 0.0;
  // (Nested std::max: the initializer-list overload compiles to an
  // out-of-line call with MSVC, measurable in this hot loop.)
  double hi = std::max(std::max(2.0 / eps, 2.0 * b), 1e-6);
  for (int iter = 0; iter < 100; ++iter) {
    // Newton on F(R) = R g(R) = R - phi(R)/eps - b^2/R. (This loop is the
    // hottest code of SRIM mode: phi and phi' share their exponentials.)
    double ph = 0.0;
    double dph = 0.0;
    f.phi_and_derivative(r, ph, dph);
    const double F = r - ph / eps - b * b / r;
    const double dF = 1.0 - dph / eps + b * b / (r * r);
    if (F < 0.0) lo = std::max(lo, r); else hi = std::min(hi, r);
    double next = r - F / dF;
    if (!(next > lo) || !(next < hi) || !std::isfinite(next)) {
      next = 0.5 * (lo + hi); // bisection step
    }
    if (std::abs(next - r) <= 1e-9 * next) return next;
    r = next;
  }
  return r;
}

double magic_sin2_half_theta(const ScreeningFunction &f, double eps, double b) {
  constexpr double C1 = 0.99229;
  constexpr double C2 = 0.011615;
  constexpr double C3 = 0.007122;
  constexpr double C4 = 9.3066;
  constexpr double C5 = 14.813;
  if (b <= 0.0) return 1.0;
  const double r0 = closest_approach(f, eps, b);
  double ph0 = 0.0;
  double dph0 = 0.0;
  f.phi_and_derivative(r0, ph0, dph0);
  const double v = ph0 / r0;                        // reduced potential V(R0)
  const double v1 = dph0 / r0 - v / r0;             // dV/dR at R0
  const double rho = -2.0 * (eps - v) / v1;         // radius of curvature
  // Far outside the screening radius (R0 above ~3500, impact parameters of
  // hundreds of Angstrom: only in targets of gas density marked as solids)
  // the potential is denormal or zero, rho overflows and cos(theta/2) below
  // would be inf/inf = NaN. The limit is no deflection. (Never reached for
  // the impact parameters of condensed targets: their results are unchanged.)
  if (!(std::abs(rho) <= std::numeric_limits<double>::max())) return 0.0;
  const double sqe = std::sqrt(eps);
  const double alpha = 1.0 + C1 / sqe;
  const double beta = (C2 + sqe) / (C3 + sqe);
  const double gam = (C4 + eps) / (C5 + eps);
  const double A = 2.0 * alpha * eps * std::pow(b, beta);
  // G = 1 / (gamma ((1 + A^2)^(1/2) - A)), as in the TRIM-85 code. (The SRIM
  // book prints Eq. 7-11 as G = gamma / ((1 + A^2)^(1/2) - A), which only
  // agrees for gamma -> 1 and underestimates sin^2(theta/2) by up to ~40 %
  // at eps ~ 0.1-1; the code form reproduces the scattering integral.)
  const double G = 1.0 / (gam * (std::sqrt(1.0 + A * A) - A));
  const double delta = A * (r0 - b) / (1.0 + G);
  const double cos_half = (b + rho + delta) / (r0 + rho);
  const double s2 = 1.0 - cos_half * cos_half;
  return std::clamp(s2, 0.0, 1.0);
}

double rutherford_sin2_half_theta(double eps, double b) noexcept {
  const double x = 2.0 * eps * b;
  return 1.0 / (1.0 + (1.0 + b * (1.0 + b)) * x * x);
}

double integral_sin2_half_theta(const ScreeningFunction &f, double eps,
                                double b, int nodes) {
  if (b <= 0.0) return 1.0;
  const double r0 = closest_approach(f, eps, b);
  const int n = std::max(4, nodes + (nodes % 2));
  double sum = 0.0;
  for (int j = 1; j <= n / 2; ++j) {
    const double u = std::cos((2.0 * j - 1.0) * pi_v / (2.0 * n));
    const double r = r0 / u;
    const double g = 1.0 - f.phi(r) / (eps * r) - b * b / (r * r);
    sum += std::sqrt((1.0 - u * u) / std::max(g, 1e-300));
  }
  const double theta = pi_v * (1.0 - 2.0 * b / (n * r0) * sum);
  const double s = std::sin(0.5 * theta);
  return std::clamp(s * s, 0.0, 1.0);
}

double trim_sin2_half_theta(const ScreeningFunction &f, double eps, double b,
                            ScatteringMethod method) {
  if (eps > 10.0) return rutherford_sin2_half_theta(eps, b);
  if (method == ScatteringMethod::integral) return integral_sin2_half_theta(f, eps, b);
  return magic_sin2_half_theta(f, eps, b);
}

CollisionKinematics collision_kinematics(const CollisionPair &pair,
                                         energy_quantity energy,
                                         double sin2_half_theta) {
  CollisionKinematics k;
  const double s2 = std::clamp(sin2_half_theta, 0.0, 1.0);
  const double c2 = 1.0 - s2;
  k.sin2_half_theta = s2;
  k.cos_theta_cm = 1.0 - 2.0 * s2;
  const double sin_theta = 2.0 * std::sqrt(s2 * c2);
  k.transferred = pair.gamma * s2 * energy;
  const double mu = pair.mass_ratio;
  const double den = std::sqrt(1.0 + 2.0 * mu * k.cos_theta_cm + mu * mu);
  if (den < 1e-12) {
    k.cos_projectile = 0.0;
    k.sin_projectile = 1.0;
  } else {
    k.cos_projectile = (k.cos_theta_cm + mu) / den;
    k.sin_projectile = sin_theta / den;
  }
  k.cos_recoil = std::sqrt(s2);
  k.sin_recoil = std::sqrt(c2);
  return k;
}

double trim_max_reduced_impact_parameter(double eps, double eps_min,
                                         double gamma) noexcept {
  const double xi = std::sqrt(eps * eps_min / gamma);
  return 1.0 / (xi + std::sqrt(xi) + 0.125 * std::pow(xi, 0.1));
}

CollisionKinematics collision_kinematics_half_angle(const CollisionPair &pair,
                                                    energy_quantity energy,
                                                    double sin_half, double cos_half) {
  CollisionKinematics k;
  const double s2 = sin_half * sin_half;
  const double c2 = cos_half * cos_half;
  k.sin2_half_theta = s2;
  k.cos_theta_cm = c2 - s2;
  const double sin_theta = 2.0 * sin_half * cos_half;
  k.transferred = pair.gamma * s2 * energy;
  const double mu = pair.mass_ratio;
  const double den = std::sqrt(1.0 + 2.0 * mu * k.cos_theta_cm + mu * mu);
  if (den < 1e-12) {
    k.cos_projectile = 0.0;
    k.sin_projectile = 1.0;
  } else {
    const double inv = 1.0 / den;
    k.cos_projectile = (k.cos_theta_cm + mu) * inv;
    k.sin_projectile = sin_theta * inv;
  }
  k.cos_recoil = sin_half;
  k.sin_recoil = cos_half;
  return k;
}

// ---------------------------------------------------------------------------
// Tabulated magic formula (fast approximations)
// ---------------------------------------------------------------------------

namespace {

/*!
 * The scattering angle of the magic formula, computed for the table: the
 * same formula as magic_sin2_half_theta, but 1 - cos(theta/2) is formed
 * without cancellation (small angles keep their relative precision).
 */
double magic_theta(const ScreeningFunction &f, double eps, double b) {
  constexpr double C1 = 0.99229;
  constexpr double C2 = 0.011615;
  constexpr double C3 = 0.007122;
  constexpr double C4 = 9.3066;
  constexpr double C5 = 14.813;
  if (b <= 0.0) return pi_v;
  const double r0 = closest_approach(f, eps, b);
  double ph0 = 0.0;
  double dph0 = 0.0;
  f.phi_and_derivative(r0, ph0, dph0);
  const double v = ph0 / r0;
  const double v1 = dph0 / r0 - v / r0;
  const double rho = -2.0 * (eps - v) / v1;
  if (!(std::abs(rho) <= std::numeric_limits<double>::max())) return 0.0; // as magic_sin2_half_theta
  const double sqe = std::sqrt(eps);
  const double alpha = 1.0 + C1 / sqe;
  const double beta = (C2 + sqe) / (C3 + sqe);
  const double gam = (C4 + eps) / (C5 + eps);
  const double A = 2.0 * alpha * eps * std::pow(b, beta);
  const double G = 1.0 / (gam * (std::sqrt(1.0 + A * A) - A));
  const double delta = A * (r0 - b) / (1.0 + G);
  const double cos_half = (b + rho + delta) / (r0 + rho);
  const double one_minus_cos_half = (r0 - b - delta) / (r0 + rho);
  const double s2 = std::clamp(one_minus_cos_half * (1.0 + cos_half), 0.0, 1.0);
  return 2.0 * std::asin(std::sqrt(s2));
}

//! Catmull-Rom weights at fraction t of the interval between nodes 1 and 2.
inline void catmull_rom(double t, double w[4]) noexcept {
  const double t2 = t * t;
  const double t3 = t2 * t;
  w[0] = -0.5 * t3 + t2 - 0.5 * t;
  w[1] = 1.5 * t3 - 2.5 * t2 + 1.0;
  w[2] = -1.5 * t3 + 2.0 * t2 + 0.5 * t;
  w[3] = 0.5 * t3 - 0.5 * t2;
}

//! Cubic Lagrange weights for nodes 0..3 at position t (0 <= t < 1): the
//! one-sided stencil of the first interval.
inline void lagrange_first(double t, double w[4]) noexcept {
  w[0] = (t - 1.0) * (t - 2.0) * (t - 3.0) / -6.0;
  w[1] = t * (t - 2.0) * (t - 3.0) / 2.0;
  w[2] = t * (t - 1.0) * (t - 3.0) / -2.0;
  w[3] = t * (t - 1.0) * (t - 2.0) / 6.0;
}

} // namespace

MagicTable::MagicTable() {
  constexpr int nu = energy_points;
  constexpr int nt = impact_points;
  u0_ = std::log(min_eps);
  const double du = (std::log(max_eps) - u0_) / (nu - 1);
  const double dt = std::sqrt(max_b) / (nt - 1);
  inv_du_ = 1.0 / du;
  inv_dt_ = 1.0 / dt;
  ln_theta_.resize(static_cast<std::size_t>(nu) * nt);
  const ScreeningFunction f = ScreeningFunction::universal();
  // 640 000 evaluations of the formula: rows in parallel.
  const auto fill = [&](int first, int step) {
    for (int i = first; i < nu; i += step) {
      const double eps = std::exp(u0_ + du * i);
      for (int j = 0; j < nt; ++j) {
        const double t = dt * j;
        const double theta = magic_theta(f, eps, t * t);
        // Far out the formula's angle underflows to 0 in double precision.
        // A floor of 1e-22 rad (sin^2 ~ 1e-44: no energy is transferred)
        // keeps ln(theta) smooth; a floor near the smallest double would
        // make the interpolation ring between -690 and the last real value.
        ln_theta_[static_cast<std::size_t>(i) * nt + j] = std::log(std::max(theta, 1e-22));
      }
    }
  };
  const unsigned hw = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
  std::vector<std::thread> workers;
  try {
    for (unsigned k = 1; k < hw; ++k) workers.emplace_back(fill, static_cast<int>(k), static_cast<int>(hw));
  } catch (const std::system_error &) {
    // No more threads: the rows of the missing workers are filled below.
  }
  const int step = static_cast<int>(workers.size()) + 1;
  fill(0, step);
  for (auto &w : workers) w.join();
}

const MagicTable &MagicTable::universal() {
  static const MagicTable table;
  return table;
}

bool MagicTable::half_angle(double eps, double b, double &sin_half, double &cos_half) const noexcept {
  if (!(eps >= min_eps && eps <= max_eps) || !(b >= 0.0 && b < max_b)) return false;
  // Node i - 1 .. i + 2 in energy must exist: 1 <= x < nu - 2.
  const double x = (std::log(eps) - u0_) * inv_du_;
  const double y = std::sqrt(b) * inv_dt_;
  if (!(x >= 1.0 && x < energy_points - 2) || !(y < impact_points - 2)) return false;
  const int i = static_cast<int>(x);
  const int j = static_cast<int>(y);
  double wx[4];
  double wy[4];
  catmull_rom(x - i, wx);
  const double *base = nullptr;
  if (j == 0) {
    lagrange_first(y, wy);
    base = &ln_theta_[static_cast<std::size_t>(i - 1) * impact_points];
  } else {
    catmull_rom(y - j, wy);
    base = &ln_theta_[static_cast<std::size_t>(i - 1) * impact_points + (j - 1)];
  }
  double lt = 0.0;
  for (int a = 0; a < 4; ++a) {
    const double *row = base + static_cast<std::size_t>(a) * impact_points;
    lt += wx[a] * (wy[0] * row[0] + wy[1] * row[1] + wy[2] * row[2] + wy[3] * row[3]);
  }
  // theta <= pi (interpolation may overshoot slightly near head-on
  // collisions; sincos_quarter also clamps theta/2 to [0, pi/2]).
  constexpr double ln_pi = 1.1447298858494002;
  if (lt > ln_pi) lt = ln_pi;
  fast::sincos_quarter(0.5 * std::exp(lt), sin_half, cos_half);
  return true;
}

void fast_half_angle(const CollisionPair &pair, double eps, double b, ScatteringMethod method,
                     double &sin_half, double &cos_half) {
  if (eps > 10.0) {
    // Screened Rutherford: s2 = 1/(1 + q), 1 - s2 = q/(1 + q).
    const double x = 2.0 * eps * b;
    const double q = (1.0 + b * (1.0 + b)) * x * x;
    sin_half = std::sqrt(1.0 / (1.0 + q));
    cos_half = std::sqrt(q / (1.0 + q));
    return;
  }
  if (method == ScatteringMethod::magic && pair.universal_screening &&
      MagicTable::universal().half_angle(eps, b, sin_half, cos_half)) {
    return;
  }
  const double s2 = trim_sin2_half_theta(pair.screening, eps, b, method);
  sin_half = std::sqrt(s2);
  cos_half = std::sqrt(1.0 - s2);
}

} // namespace physics
} // namespace amcset
