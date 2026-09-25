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

/*!
 * \file common_test.cpp
 *
 * \brief The testing code for the amcset_common library.
 *
 * This testing suite verifies the functionality and correctness of the AMCSET
 * library. It is designed for close to 100% code coverage.
 */

#include <any>
#include <sstream>

#define private public
#include "amcset_common.h"
#undef private

#include <gtest/gtest.h>

#include <boost/random/mersenne_twister.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <glog/logging.h>

#include "amcset_utilities.h"

using namespace amcset::common;

/*!
 * \brief The testing fixture class for the amcset_common tests.
 *
 * This test fixture allows tests which utilize this fixture to run the same
 * setup code (contained within the member functions defined in CommonTest).
 */
class CommonTest : public testing::Test {
protected:
  //! The constructor for CommonTest
  CommonTest() {}
  //! The override for the SetUp() function
  void SetUp() override {
    static bool logging_initialized = false;
    if (!logging_initialized) {
      google::InitGoogleLogging("common-test");
      FLAGS_logtostderr = 1;
      ASSERT_NE(std::freopen("log.txt", "w", stderr), nullptr);
      logging_initialized = true;
    }
  }
  //! The override for the TearDown() function
  void TearDown() override {}

  //! Convenience function for creating a water layer for testing
  Layer create_water_layer(length_quantity depth) const {
    constexpr Particle::Properties hydrogen(1, 1);
    constexpr Particle::Properties oxygen(8, 16);
    constexpr auto mass_density = 1000.0 * kg_per_cubic_meter;
    return Layer({{2.0, hydrogen}, {1.0, oxygen}}, depth, mass_density);
  }

  //! Convenience function for creating an iron layer for testing
  Layer create_iron_layer(length_quantity depth) const try {
    // TODO: Replace with natural-abundance function
    constexpr double natural_abundance_54 = 0.05845;
    constexpr double natural_abundance_56 = 0.91754;
    constexpr double natural_abundance_57 = 0.02119;
    constexpr double natural_abundance_58 = 0.00282;
    Layer::material_vector iron{
        {natural_abundance_54, Particle::Properties(26, 54)},
        {natural_abundance_56, Particle::Properties(26, 56)},
        {natural_abundance_57, Particle::Properties(26, 57)},
        {natural_abundance_58, Particle::Properties(26, 58)},
    };
    constexpr auto mass_density = 7874.0 * kg_per_cubic_meter;
    return Layer(std::move(iron), depth, mass_density);
  } catch (const std::exception &e) {
    print_exception(e);
    throw e;
  }

  /*!
   * \brief A convenience function for creating a simulation
   */
  Simulation create_simulation() const try {
    auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
    size_t z_number = 14; // Silicon
    size_t mass_number = 28;
    bool enable_damage_cascade = true;
    auto ion_stopping_energy =
        energy_quantity(0.04 * kilo_electron_volt); // Stop at 4 eV
    auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
    bool log_single_displacement = true;
    size_t divisor_angle_number = 1000;
    size_t flying_distance_number = 1000;
    auto range = length_quantity(80000000.0 * angstrom);
    size_t bombardment_count = 1;
    bool is_electron = false;
    auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
    size_t thread_count = 8;
    angle_quantity incident_angle = 0;

    auto volume =
        Volume({create_water_layer(length_quantity(2.0 * angstrom)),
                create_iron_layer(length_quantity(500.0 * angstrom)),
                create_water_layer(length_quantity(1000.0 * angstrom))});
    Simulation::Settings settings(
        electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
        ion_stopping_energy, ion_displacement_energy, log_single_displacement,
        divisor_angle_number, flying_distance_number, range, bombardment_count,
        is_electron, incident_energy, thread_count, incident_angle);

    return Simulation(settings, std::move(volume));
  } catch (const std::exception &e) {
    print_exception(e);
    throw e;
  }
};

//! Template test. Does not test the library, just a simple demonstration.
TEST_F(CommonTest, common_assertions) try {
  EXPECT_STRNE("hello", "world");
  EXPECT_EQ(7 * 6, 42);
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

TEST_F(CommonTest, simulation_tracks_completed_bombardments) try {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 14;
  size_t mass_number = 28;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 1000;
  size_t flying_distance_number = 1000;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 3;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 2;
  angle_quantity incident_angle = 0;

  auto volume = Volume({create_water_layer(length_quantity(2.0 * angstrom)),
                        create_iron_layer(length_quantity(500.0 * angstrom)),
                        create_water_layer(length_quantity(1000.0 * angstrom))});
  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle);

  Simulation simulation(settings, std::move(volume));
  EXPECT_EQ(simulation.get_completed_bombardment_count(), size_t{0});
  EXPECT_EQ(simulation.run_simulation(), 0);
  EXPECT_EQ(simulation.get_completed_bombardment_count(), bombardment_count);
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

TEST_F(CommonTest, electron_bombardment_runs) try {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 26;
  size_t mass_number = 56;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 64;
  size_t flying_distance_number = 64;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 1;
  bool is_electron = true;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 1;
  angle_quantity incident_angle = angle_quantity(0.0 * si::radian);

  auto volume = Volume({create_iron_layer(range)});
  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle);

  Simulation simulation(settings, std::move(volume));
  EXPECT_EQ(simulation.run_simulation(), 0);
  EXPECT_EQ(simulation.get_completed_bombardment_count(), bombardment_count);
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Coordinate struct test.
TEST_F(CommonTest, coordinate_constructor_test) try {
  constexpr Coordinate coordinate(length_quantity(-1.0 * angstrom),
                                  length_quantity(0.0 * angstrom),
                                  length_quantity(1.0 * angstrom));

  ASSERT_EQ(to_string(coordinate.x_), "-1e-10 m");
  ASSERT_EQ(to_string(coordinate.y_), "0 m");
  ASSERT_EQ(to_string(coordinate.z_), "1e-10 m");
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Velocity struct test.
TEST_F(CommonTest, velocity_constructor_test) try {
  // constexpr Velocity velocity(angle_quantity(-pi * si::radian),
  //                             angle_quantity(pi * si::radian),
  //                             energy_quantity(10000.0 * kilo_electron_volt));
  Velocity velocity(Quaternion(0.0, 0.0, 0.0, 2.0),
                    energy_quantity(10000.0 * kilo_electron_volt));

  // ASSERT_EQ(to_string(velocity.x_angle_), "-3.14159 rad");
  ASSERT_EQ(velocity.quaternion_.r_, 0.0);
  ASSERT_EQ(velocity.quaternion_.x_, 0.0);
  ASSERT_EQ(velocity.quaternion_.y_, 0.0);
  ASSERT_EQ(velocity.quaternion_.z_, 1.0);
  ASSERT_EQ(to_string(velocity.energy_), "1.60218e-12 m^2 kg s^-2");
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Isotope data test
TEST_F(CommonTest, isotope_data_get_isotope_mass_test) try {
  ASSERT_EQ(IsotopeData::get_isotope_mass(50, 132),
            double(131.9178267) * atomic_mass_unit);
  ASSERT_EQ(IsotopeData::get_isotope_mass(2, 4),
            double(4.00260325413) * atomic_mass_unit);
  ASSERT_EQ(IsotopeData::get_isotope_mass(118, 295),
            double(295.21624) * atomic_mass_unit);

  ASSERT_EQ(IsotopeData::get_isotope_mass(0, 0), double(1.0) * constants::m_e);

  ASSERT_THROW(IsotopeData::get_isotope_mass(0, 1), std::out_of_range);
  ASSERT_THROW(IsotopeData::get_isotope_mass(1, 0), std::out_of_range);
  // Looking up a non-existent mass number for a valid element should throw
  // out_of_range (no matching isotope in the database).
  ASSERT_THROW(IsotopeData::get_isotope_mass(1, 100), std::out_of_range);
} catch (const std::exception &e) {
  print_exception(e);
  throw e;
}

TEST_F(CommonTest, screening_length_quickbasic_equivalence) try {
  // Construct a minimal simulation with QuickBasic screening model
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 14; // Silicon
  size_t mass_number = 28;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 1000;
  size_t flying_distance_number = 1000;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 1;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 1;
  angle_quantity incident_angle = 0;

  auto volume = Volume({create_water_layer(length_quantity(2.0 * angstrom)),
                        create_iron_layer(length_quantity(500.0 * angstrom)),
                        create_water_layer(length_quantity(1000.0 * angstrom))});

  // Use QuickBasic screening model explicitly
  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle,
      100, 100, 100, 100, std::string("output/displacement_histogram.csv"),
      std::string("output/electronic_loss_histogram.csv"),
      std::string("output/nuclear_loss_histogram.csv"),
      std::string("output/total_loss_histogram.csv"), false,
      std::string("output/particle_tree.json"),
      Simulation::Settings::ScreeningModel::QuickBasic);

  Simulation simulation(settings, std::move(volume));

  const std::vector<size_t> tests = {1, 14, 29};
  for (auto z : tests) {
    auto z1 = amcset::common::dimensionless_quantity(static_cast<double>(z));
    auto z2 = amcset::common::dimensionless_quantity(static_cast<double>(z));
    auto got = Ion::screening_length_for_model(z1, z2, true);
    auto qb = 0.8854 * 0.529 * angstrom / pow(pow(z1, 0.5) + pow(z2, 0.5), 2.0 / 3.0);
    const double rel = std::abs(((got - qb) / qb).value());
    EXPECT_LE(rel, 0.01);
  }
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
  print_exception(e);
  FAIL();
}

//! Particle::Properties test
TEST_F(CommonTest, properties_test) try {
  Particle::Properties electron(0, 0);

  ASSERT_EQ(electron.mass_, double(1.0) * constants::m_e);
  ASSERT_EQ(electron.charge_, -1);
  ASSERT_THROW(Particle::Properties(0, 1), std::out_of_range);
  ASSERT_THROW(Particle::Properties(1, 0), std::out_of_range);

  Particle::Properties iron(26, 56);
  ASSERT_EQ(iron.mass_, IsotopeData::get_isotope_mass(26, 56));
  ASSERT_EQ(iron.charge_, 26);

} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Layer class test
TEST_F(CommonTest, layer_test) try {
  const auto water_layer = create_water_layer(length_quantity(1.0 * angstrom));
  const auto depth = length_quantity(1.0 * angstrom);

  // Testing get_relative_compositions()
  const auto composition_view = water_layer.get_relative_compositions();
  std::vector<double> compare_doubles = {2.0 / 3.0, 1.0 / 3.0};

  ASSERT_TRUE(std::equal(composition_view.begin(), composition_view.end(),
                         compare_doubles.begin(), compare_doubles.end()));
  compare_doubles.at(0) = 4.0;
  ASSERT_FALSE(std::equal(composition_view.begin(), composition_view.end(),
                          compare_doubles.begin(), compare_doubles.end()));

  // Testing get_property
  ASSERT_EQ(water_layer.get_property(0).mass_,
            IsotopeData::get_isotope_mass(1, 1));
  ASSERT_EQ(water_layer.get_property(1).mass_,
            IsotopeData::get_isotope_mass(8, 16));

  // Testing get_relative_composition
  ASSERT_EQ(water_layer.get_relative_composition(0), 2.0 / 3.0);
  ASSERT_EQ(water_layer.get_relative_composition(1), 1.0 / 3.0);

  constexpr auto mass_density = 1000.0 * kg_per_cubic_meter;

  // Testing contructor
  Layer::material_vector empty_vec{};
  ASSERT_THROW(Layer(std::move(empty_vec), depth, mass_density),
               std::invalid_argument);
  ASSERT_THROW(Layer(std::move(empty_vec), double(-1.0) * depth, mass_density),
               std::invalid_argument);

  // Testing get_depth
  ASSERT_EQ(water_layer.get_depth(), length_quantity(1.0 * angstrom));

  // Testing get_number_of_components
  ASSERT_EQ(water_layer.get_number_of_components(), 2);

  // Testing get_mass_density
  ASSERT_EQ(water_layer.get_mass_density(), mass_density);

  // Testing get_number_density
  ASSERT_EQ(water_layer.get_number_density(),
            mass_density /
                (((2.0 / 3.0) * IsotopeData::get_isotope_mass(1, 1) +
                  (1.0 / 3.0) * IsotopeData::get_isotope_mass(8, 16)) *
                 constants::N_A));

  // Testing discrete distibution
  auto iron_layer = create_iron_layer(500.0 * angstrom);
  boost::random::mt19937 rng;
  auto distribution = iron_layer.get_discrete_distribution();
  size_t random_number;
  for (int i = 0; i < 1000; i++) {
    random_number = distribution(rng);
    ASSERT_TRUE(random_number < iron_layer.get_number_of_components());
  }

} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Volume class test
TEST_F(CommonTest, volume_test) try {
  const auto water_layer_1 =
      create_water_layer(length_quantity(2.0 * angstrom));
  const auto iron_layer = create_iron_layer(length_quantity(500.0 * angstrom));
  const auto water_layer_2 =
      create_water_layer(length_quantity(1000.0 * angstrom));
  auto volume = Volume({water_layer_1, iron_layer, water_layer_2});

  // Test constructor
  ASSERT_THROW(Volume({water_layer_2, iron_layer, water_layer_1}),
               std::invalid_argument);
  // ASSERT_THROW(Volume(std::vector<Layer>()), std::invalid_argument);

  // Test get_layer
  ASSERT_EQ(
      volume.get_layer(length_quantity(1.0 * angstrom)).first->get_depth(),
      water_layer_1.get_depth());
  ASSERT_EQ(
      volume.get_layer(length_quantity(2.0 * angstrom)).first->get_depth(),
      water_layer_1.get_depth());
  ASSERT_EQ(
      volume.get_layer(length_quantity(250.0 * angstrom)).first->get_depth(),
      iron_layer.get_depth());
  ASSERT_EQ(
      volume.get_layer(length_quantity(501.0 * angstrom)).first->get_depth(),
      water_layer_2.get_depth());
  ASSERT_THROW(volume.get_layer(length_quantity(1001.0 * angstrom)),
               std::out_of_range);

  ASSERT_EQ(volume.get_layer(length_quantity(1.0 * angstrom)).second, 0);
  ASSERT_EQ(volume.get_layer(length_quantity(250.0 * angstrom)).second, 1);
  ASSERT_EQ(volume.get_layer(length_quantity(501.0 * angstrom)).second, 2);

  ASSERT_EQ(
      volume.get_layer(length_quantity(1.0 * angstrom)).first->get_depth(),
      volume.get_layer(0).get_depth());
  ASSERT_EQ(
      volume.get_layer(length_quantity(250.0 * angstrom)).first->get_depth(),
      volume.get_layer(1).get_depth());
  ASSERT_EQ(
      volume.get_layer(length_quantity(501.0 * angstrom)).first->get_depth(),
      volume.get_layer(2).get_depth());
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Simulation and Simulation::Settings class test
TEST_F(CommonTest, simulation_and_settings_test) try {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 14;
  size_t mass_number = 28;
  bool enable_damage_cascade = false;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 1000;
  size_t flying_distance_number = 1000;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 1;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 8;

  auto volume =
      Volume({create_water_layer(length_quantity(2.0 * angstrom)),
              create_iron_layer(range),
              });

  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count,
      0.0 * amcset::common::radian, 64, 128, 256, 512,
      "output/displacement_histogram.csv",
      "output/electronic_loss_histogram.csv",
      "output/nuclear_loss_histogram.csv",
      "output/total_loss_histogram.csv", false,
      "output/particle_tree.json");

  ASSERT_EQ(settings.electron_stopping_energy_, electron_stopping_energy);
  ASSERT_EQ(settings.incident_particle_properties_.mass_,
            IsotopeData::get_isotope_mass(z_number, mass_number));
  ASSERT_EQ(settings.incident_particle_properties_.charge_, z_number);
  ASSERT_EQ(settings.enable_damage_cascade_, enable_damage_cascade);
  ASSERT_EQ(settings.ion_stopping_energy_, ion_stopping_energy);
  ASSERT_EQ(settings.ion_displacement_energy_, ion_displacement_energy);
  ASSERT_EQ(settings.log_single_displacement_, log_single_displacement);
  ASSERT_EQ(settings.divisor_angle_number_, divisor_angle_number);
  ASSERT_EQ(settings.flying_distance_number_, flying_distance_number);
  ASSERT_EQ(settings.range_, range);
  ASSERT_EQ(settings.bombardment_count_, bombardment_count);
  ASSERT_EQ(settings.is_electron_, is_electron);
  ASSERT_EQ(settings.incident_energy_, incident_energy);
  ASSERT_EQ(settings.thread_count_, thread_count);
  ASSERT_EQ(settings.displacement_bin_count_, 64u);
  ASSERT_EQ(settings.electronic_energy_loss_bin_count_, 128u);
  ASSERT_EQ(settings.nuclear_energy_loss_bin_count_, 256u);
  ASSERT_EQ(settings.total_energy_loss_bin_count_, 512u);
  ASSERT_EQ(settings.displacement_histogram_output_file_,
            std::string("output/displacement_histogram.csv"));
  ASSERT_FALSE(settings.log_particle_tree_);
  ASSERT_EQ(settings.particle_tree_output_file_,
            std::string("output/particle_tree.json"));

  Simulation simulation(settings, std::move(volume));

  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::electron_stopping_energy_),
      electron_stopping_energy);
  auto simulation_incident_particle_properties = simulation.get_settings(
      &Simulation::Settings::incident_particle_properties_);
  ASSERT_EQ(simulation_incident_particle_properties.charge_, z_number);
  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::enable_damage_cascade_),
      enable_damage_cascade);
  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::ion_stopping_energy_),
      ion_stopping_energy);
  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::ion_displacement_energy_),
      ion_displacement_energy);
  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::log_single_displacement_),
      log_single_displacement);
  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::divisor_angle_number_),
      divisor_angle_number);
  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::flying_distance_number_),
      flying_distance_number);
  ASSERT_EQ(simulation.get_settings(&Simulation::Settings::range_), range);
  ASSERT_EQ(simulation.get_settings(&Simulation::Settings::bombardment_count_),
            bombardment_count);
  ASSERT_EQ(simulation.get_settings(&Simulation::Settings::is_electron_),
            is_electron);
  ASSERT_EQ(simulation.get_settings(&Simulation::Settings::incident_energy_),
            incident_energy);
  ASSERT_EQ(simulation.get_settings(&Simulation::Settings::thread_count_),
            thread_count);

  EXPECT_NO_THROW(simulation.run_simulation());
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

TEST_F(CommonTest, simulation_writes_all_output_files) try {
  std::filesystem::remove_all("output/test_regression");
  std::filesystem::remove("output/collisions.csv");

  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 14;
  size_t mass_number = 28;
  bool enable_damage_cascade = false;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 1000;
  size_t flying_distance_number = 1000;
  auto range = length_quantity(5000000.0 * angstrom);
  size_t bombardment_count = 1;
  bool is_electron = false;
  auto incident_energy = energy_quantity(10.0 * kilo_electron_volt);
  size_t thread_count = 1;
  angle_quantity incident_angle = 0.0 * amcset::common::radian;

  auto volume = Volume({create_water_layer(length_quantity(2.0 * angstrom)),
                        create_iron_layer(length_quantity(100.0 * angstrom))});

  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle, 64, 64, 64,
      64, "output/test_regression/displacement_histogram.csv",
      "output/test_regression/electronic_loss_histogram.csv",
      "output/test_regression/nuclear_loss_histogram.csv",
      "output/test_regression/total_loss_histogram.csv", false,
      "output/test_regression/particle_tree.json");

  Simulation simulation(settings, std::move(volume));
  EXPECT_NO_THROW(simulation.run_simulation());

  EXPECT_TRUE(std::filesystem::exists("output/test_regression/displacement_histogram.csv"));
  EXPECT_TRUE(std::filesystem::exists("output/test_regression/electronic_loss_histogram.csv"));
  EXPECT_TRUE(std::filesystem::exists("output/test_regression/nuclear_loss_histogram.csv"));
  EXPECT_TRUE(std::filesystem::exists("output/test_regression/total_loss_histogram.csv"));
  EXPECT_FALSE(std::filesystem::exists("output/collisions.csv"));
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Bombardment test
TEST_F(CommonTest, bombardment_test) try {
  // Construct a minimal simulation that uses an iron incident particle
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  size_t z_number = 26; // Iron
  size_t mass_number = 56;
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 1000;
  size_t flying_distance_number = 1000;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 1;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 1;
  angle_quantity incident_angle = 0;

  auto volume = Volume({create_iron_layer(range),create_water_layer(range + length_quantity(1000.0 * angstrom))}); // Single iron layer

  Simulation::Settings settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle);

  Simulation simulation(settings, std::move(volume));

  // Run a single bombardment and assert it does not throw
  Bombardment bombardment(simulation, 0);
  EXPECT_NO_THROW(bombardment.run_bombardment());
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Ion test
TEST_F(CommonTest, ion_test) try {
  Simulation simulation = create_simulation();
  // TODO: Implement Ion test — skip for now until fully specified
  GTEST_SKIP() << "Ion test not implemented";
} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}

//! Electronic stopping energy test
TEST_F(CommonTest, electronic_stopping_energy_test) try {
  auto simulation = create_simulation();
  std::cout << "Simulation created" << std::endl; // TODO: Remove

  // Populate discrete distribution with iron layer composition
  boost::random::discrete_distribution<size_t, double> discrete_distribution(
      simulation.get_volume().get_layer(1).get_relative_compositions().begin(),
      simulation.get_volume().get_layer(1).get_relative_compositions().end());
  boost::random::uniform_01<double> uniform_distribution;
  boost::random::mt19937 random_number_generator;

  const size_t z_number = 14;
  const size_t mass_number = 28;

  ParticleDataLeaf test_data_leaf; // Required for Ion constructor
  Bombardment bombardment(simulation, 0);

  Ion test(z_number, mass_number, simulation, uniform_distribution,
           random_number_generator, test_data_leaf, bombardment);

  // auto silicon_28_mass = IsotopeData::get_isotope_mass(14, 28);

  // auto charge = 14;
  // auto mass = silicon_28_mass;
  // auto energy = energy_quantity(50.0 * kilo_electron_volt);
  auto layer = simulation.get_volume().get_layer(1);
  std::cout << "In-test substrate charge value: "
            << layer.get_property(0).charge_ << std::endl;
  // auto final_energy = test.electronic_stopping_energy(
  //     charge, mass, energy, layer.get_number_density(),
  //     layer.get_property(0));
  // std::cout << std::setprecision(17) << "Energy: " << final_energy <<
  // std::endl;

  // ASSERT_EQ(7.9853010195075595e-15 * si::joule, final_energy);

} catch (const std::exception &e) {
  print_exception(e);
  FAIL();
}
