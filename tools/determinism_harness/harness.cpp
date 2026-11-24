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

/*! \file harness.cpp
 *  \brief Determinism harness for the AMCSET simulation core.
 *
 *  Fires particles with fixed RNG seeds (bypassing EntropySeeder) and prints
 *  every recorded CollisionPoint with 17 significant digits. Run it on two
 *  builds and compare the outputs byte-for-byte (cmp/diff) to prove that a
 *  refactor did not change numerical results. See docs/TESTING.md.
 *
 *  Run from the repository root (or any directory containing common/assets)
 *  so the electron parameter tables resolve.
 */

#include "amcset_common.h"

#include <absl/log/initialize.h>

#include <cstdio>
#include <memory>

using namespace amcset::common;

namespace {

Layer water(length_quantity depth) {
  return Layer({{2.0, Particle::Properties(1, 1)},
                {1.0, Particle::Properties(8, 16)}},
               depth, 1000.0 * kg_per_cubic_meter);
}

Layer iron(length_quantity depth) {
  Layer::material_vector v{{0.05845, Particle::Properties(26, 54)},
                           {0.91754, Particle::Properties(26, 56)},
                           {0.02119, Particle::Properties(26, 57)},
                           {0.00282, Particle::Properties(26, 58)}};
  return Layer(std::move(v), depth, 7874.0 * kg_per_cubic_meter);
}

void run(const char *name, size_t z, size_t a, bool electron, double keV,
         Simulation::Settings::ScreeningModel screening,
         Simulation::Settings::ElectronicStoppingModel stopping, int seeds,
         double angle_deg) {
  const auto range = 20000.0 * angstrom;
  // Settings overload (3); custom-model parameters are deliberately not the
  // defaults so that the Custom models are exercised.
  Simulation::Settings settings(
      1.0 * kilo_electron_volt, z, a, true, 0.04 * kilo_electron_volt,
      0.04 * kilo_electron_volt, 64, 64, range, 1, electron,
      keV * kilo_electron_volt, 1, angle_deg * amcset::common::degree, 50, 50,
      50, 50, 50, "harness_output/d.csv", "harness_output/e.csv",
      "harness_output/n.csv", "harness_output/t.csv", false,
      "harness_output/b/", "common/assets/electron_screening_parameters.csv",
      "common/assets/mott_scattering_parameters.csv", screening,
      "harness_output/r.csv", stopping, 2.0, 120.0, 1.5, 110.0, 1.3);
  Simulation simulation(
      settings,
      Volume({water(10.0 * angstrom), iron(5000.0 * angstrom), water(range)}));
  for (int seed = 1; seed <= seeds; ++seed) {
    boost::random::uniform_01<double> uniform;
    boost::random::mt19937 rng(static_cast<unsigned>(seed));
    ParticleDataLeaf leaf;
    Bombardment bombardment(simulation, 0);
    std::unique_ptr<Particle> particle;
    if (electron) {
      particle = std::make_unique<Electron>(z, a, simulation, uniform, rng,
                                            leaf, bombardment);
    } else {
      particle = std::make_unique<Ion>(z, a, simulation, uniform, rng, leaf,
                                       bombardment);
    }
    particle->fire();
    const auto points = bombardment.take_collision_points();
    std::printf("%s seed=%d npts=%zu\n", name, seed, points.size());
    for (const auto &c : points) {
      std::printf("%.17g %.17g %.17g %.17g %.17g %.17g %.17g %zu %zu %.17g "
                  "%.17g\n",
                  c.coordinate.x_.value(), c.coordinate.y_.value(),
                  c.coordinate.z_.value(), c.vx, c.vy, c.vz, c.energy_eV,
                  c.generation, c.z_number, c.atomic_mass_number, c.time);
    }
  }
}

} // namespace

int main() {
  absl::InitializeLog();
  using SM = Simulation::Settings::ScreeningModel;
  using EM = Simulation::Settings::ElectronicStoppingModel;
  run("Si50_TF_ZBL", 14, 28, false, 50, SM::ThomasFermi, EM::ZBL, 15, 0);
  run("Fe200_ZBL_LSS", 26, 56, false, 200, SM::ZBL, EM::LSS, 8, 7);
  run("H10_Mol_Custom", 1, 1, false, 10, SM::Moliere, EM::CustomPowerLaw, 8,
      30);
  run("He5_Firsov_ZBL", 2, 4, false, 5, SM::QuickBasic, EM::ZBL, 8, 0);
  run("Cu30_Custom_ZBL", 29, 63, false, 30, SM::Custom, EM::ZBL, 8, 45);
  run("e20", 0, 0, true, 20, SM::ThomasFermi, EM::ZBL, 5, 0);
  return 0;
}
