// Simple integration test for tallies
#include <gtest/gtest.h>

#include <numeric>

#include "amcset_common.h"

using namespace amcset::common;

TEST(TalliesIntegration, MergeAfterBombardments) {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 26;
  size_t mass_number = 56;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  size_t divisor_angle_number = 100;
  size_t flying_distance_number = 100;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 2;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 1;

  auto volume =
      Volume({
          Layer({{1.0, Particle::Properties(26,56)}}, length_quantity(500.0*angstrom), 7874.0*kg_per_cubic_meter),
      });

  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
    ion_stopping_energy, ion_displacement_energy,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count);

  Simulation sim(settings, std::move(volume));
  ASSERT_EQ(sim.get_settings(&Simulation::Settings::bombardment_count_), bombardment_count);
  sim.run_simulation();

  const auto &t = sim.get_tallies();
  EXPECT_EQ(t.incidents, bombardment_count);
  EXPECT_GT(t.range_stats.n, 0u);
  EXPECT_EQ(t.displacement_histogram.bins.size(), settings.displacement_bin_count_);
  EXPECT_GT(t.total_displacements, 0u);
  EXPECT_GT(t.displacement_histogram.bins.size(), 0u);
}

TEST(TalliesIntegration, EnergyLossHistogramsRecordNonZeroData) {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 26;
  size_t mass_number = 56;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  size_t divisor_angle_number = 100;
  size_t flying_distance_number = 100;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 1;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 1;

  auto volume = Volume({
      Layer({{1.0, Particle::Properties(26, 56)}}, length_quantity(500.0 * angstrom),
            7874.0 * kg_per_cubic_meter),
  });

  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
    ion_stopping_energy, ion_displacement_energy,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count);

  Simulation sim(settings, std::move(volume));
  sim.run_simulation();

  const auto &t = sim.get_tallies();
  EXPECT_GT(t.energy_loss_stats.n, 0u);
  EXPECT_GT(std::accumulate(t.electronic_energy_loss_histogram.bins.begin(),
                           t.electronic_energy_loss_histogram.bins.end(), 0.0),
            0.0);
  EXPECT_GT(std::accumulate(t.nuclear_energy_loss_histogram.bins.begin(),
                           t.nuclear_energy_loss_histogram.bins.end(), 0.0),
            0.0);
  EXPECT_GT(std::accumulate(t.total_energy_loss_histogram.bins.begin(),
                           t.total_energy_loss_histogram.bins.end(), 0.0),
            0.0);
}
