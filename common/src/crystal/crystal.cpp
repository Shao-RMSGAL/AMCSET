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

#include "crystal/crystal.h"

#include <cmath>
#include <stdexcept>

namespace amcset {
namespace crystal {

using namespace common;

std::shared_ptr<const GrainMap> make_grain_map(const Crystal &crystal, double depth_A) {
  if (!crystal.polycrystal) throw std::invalid_argument("crystal: not a polycrystal");
  const auto &p = *crystal.polycrystal;
  if (p.source == Polycrystal::Source::voronoi) {
    return std::make_shared<const GrainMap>(generate_voronoi(p, crystal.rotation(), depth_A));
  }
  if (!p.map) {
    throw std::invalid_argument("polycrystal: the grains of '" + p.file + "' have not been read");
  }
  return p.map;
}

double debye_phi(double x) {
  if (!(x >= 0.0)) throw std::invalid_argument("debye_phi: x must be non-negative");
  if (x < 1e-6) return 1.0 - 0.25 * x; // t/(e^t - 1) = 1 - t/2 + ...
  if (x > 60.0) return pi_v * pi_v / (6.0 * x); // the rest of the integral is below 1e-24
  // Composite Simpson's rule; the integrand t/(e^t - 1) is smooth (1 at t = 0).
  const int n = 2000;
  const double h = x / n;
  const auto f = [](double t) { return t > 0.0 ? t / std::expm1(t) : 1.0; };
  double s = f(0.0) + f(x);
  for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0 : 2.0) * f(i * h);
  return s * h / 3.0 / x;
}

length_quantity debye_rms_displacement(mass_quantity mass, temperature_quantity debye,
                                       temperature_quantity temperature) {
  if (!(mass > 0.0 * kilogram)) throw std::invalid_argument("Debye amplitude: the mass must be positive");
  if (!(debye > 0.0 * kelvin)) throw std::invalid_argument("Debye amplitude: the Debye temperature must be positive");
  if (!(temperature >= 0.0 * kelvin)) throw std::invalid_argument("Debye amplitude: the temperature must not be negative");
  // <u_x^2> = 3 hbar^2 / (M k_B Theta) * (Phi(x)/x + 1/4), x = Theta / T.
  double bracket = 0.25;
  if (temperature > 0.0 * kelvin) {
    const double x = (debye / temperature).value();
    bracket += debye_phi(x) / x;
  }
  const auto u2 = 3.0 * reduced_planck_constant * reduced_planck_constant / (mass * boltzmann_constant * debye) * bracket;
  return boost::units::sqrt(u2);
}

std::optional<length_quantity> Crystal::rms_displacement(unsigned z, mass_quantity mass) const {
  for (const auto &v : vibrations) {
    if (v.z != z) continue;
    if (v.rms) return *v.rms;
    if (v.debye_temperature) return debye_rms_displacement(mass, *v.debye_temperature, temperature);
  }
  return std::nullopt;
}

} // namespace crystal
} // namespace amcset
