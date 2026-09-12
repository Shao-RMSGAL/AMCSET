// Additional tallies tests
#include <gtest/gtest.h>
#include "amcset_common.h"

using namespace amcset::common;

TEST(TalliesUnit, MergeBehavior) {
  Tallies a(10, 100.0);
  Tallies b(10, 100.0);

  a.incidents = 2;
  a.total_displacements = 5;
  a.sputtered = 1;
  a.range_stats.add(10.0L);
  a.range_stats.add(20.0L);
  a.energy_loss_stats.add(1.0L);
  a.depth_histogram.add(10.0);

  b.incidents = 3;
  b.total_displacements = 2;
  b.sputtered = 4;
  b.range_stats.add(30.0L);
  b.energy_loss_stats.add(2.0L);
  b.depth_histogram.add(10.0);

  a.merge(b);

  EXPECT_EQ(a.incidents, 5u);
  EXPECT_EQ(a.total_displacements, 7u);
  EXPECT_EQ(a.sputtered, 5u);
  ASSERT_EQ(a.depth_histogram.bins.size(), 10u);

  uint64_t sum = 0;
  for (auto v : a.depth_histogram.bins) sum += v;
  EXPECT_EQ(sum, 2u);

  EXPECT_EQ(a.range_stats.n, 3u);
  EXPECT_GT((double)a.range_stats.mean, 0.0);
  EXPECT_EQ(a.energy_loss_stats.n, 2u);
}

TEST(TalliesUnit, EnergyLossHistogramStoresEnergyPerDepthBin) {
  Histogram energy_loss(10, 0.0, 10.0);
  energy_loss.add_value(1.5, 3.0);
  energy_loss.add_value(1.5, 2.5);
  energy_loss.add_value(9.0, 4.0);

  EXPECT_DOUBLE_EQ(energy_loss.bins[1], 5.5);
  EXPECT_DOUBLE_EQ(energy_loss.bins[9], 4.0);
}

TEST(TalliesIntegration, SimulationAggregatesTallies) {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 26;
  size_t mass_number = 56;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 100;
  size_t flying_distance_number = 100;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 3;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 1;

  auto volume =
      Volume({Layer({{1.0, Particle::Properties(26,56)}}, length_quantity(500.0*angstrom), 7874.0*kg_per_cubic_meter),});

  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count);

  Simulation sim(settings, std::move(volume));
  ASSERT_EQ(sim.get_settings(&Simulation::Settings::bombardment_count_), bombardment_count);
  sim.run_simulation();

  const auto &t = sim.get_tallies();
  EXPECT_EQ(t.incidents, bombardment_count);
  EXPECT_GT(t.range_stats.n, 0u);
  EXPECT_GT(t.energy_loss_stats.n, 0u);
}
