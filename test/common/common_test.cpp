// Copyright 2024, Texas A&M University
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

/*!
 * \file common_test.cpp
 *
 * \brief The testing code for the amcset_common library.
 *
 * This testing suite verifies the functionality and correctness of the AMCSET
 * library. It is designed for close to 100% code coverage.
 */

#include <gtest/gtest.h>

#include "amcset_common.h"
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
  void SetUp() override {}
  //! The override for the TearDown() function
  void TearDown() override {}

  //! Convenience function for creating a water layer for testing
  Volume::Layer create_water_layer(length_quantity depth) const {
    constexpr Particle::Properties hydrogen(1, 1);
    constexpr Particle::Properties oxygen(8, 16);
    return Volume::Layer({{2.0, hydrogen}, {1.0, oxygen}}, depth);
  }

  //! Convenience function for creating an iron layer for testing
  Volume::Layer create_iron_layer(length_quantity depth) const try {
    // TODO: Replace with natural-abundance function
    constexpr double natural_abundance_54 = 0.05845;
    constexpr double natural_abundance_56 = 0.91754;
    constexpr double natural_abundance_57 = 0.02119;
    constexpr double natural_abundance_58 = 0.00282;
    Volume::Layer::material_vector iron{
        {natural_abundance_54, Particle::Properties(26, 54)},
        {natural_abundance_56, Particle::Properties(26, 56)},
        {natural_abundance_57, Particle::Properties(26, 57)},
        {natural_abundance_58, Particle::Properties(26, 58)},
    };
    return Volume::Layer(std::move(iron), depth);
  } catch (const std::exception& e) {
    print_exception(e);
    throw e;
  }
};

//! Template test. Does not test the library, just a simple demonstration.
TEST_F(CommonTest, common_assertions) try {
  EXPECT_STRNE("hello", "world");
  EXPECT_EQ(7 * 6, 42);
  std::cout << "Hi!" << std::endl;
} catch (const std::exception& e) {
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
} catch (const std::exception& e) {
  print_exception(e);
  FAIL();
}

//! Velocity struct test.
TEST_F(CommonTest, velocity_constructor_test) try {
  constexpr Velocity velocity(angle_quantity(-pi * si::radian),
                              angle_quantity(pi * si::radian),
                              energy_quantity(10000.0 * kilo_electron_volt));

  ASSERT_EQ(to_string(velocity.x_angle_), "-3.14159 rad");
  ASSERT_EQ(to_string(velocity.z_angle_), "3.14159 rad");
  ASSERT_EQ(to_string(velocity.energy_), "1.60218e-12 m^2 kg s^-2");
} catch (const std::exception& e) {
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
  ASSERT_THROW(IsotopeData::get_isotope_mass(1, 100), std::invalid_argument);
} catch (const std::exception& e) {
  print_exception(e);
  FAIL();
}

//! Particle::Properties test
TEST_F(CommonTest, properties_test) try {
  Particle::Properties electron(0, 0);

  ASSERT_EQ(electron.mass_, double(1.0) * constants::m_e);
  ASSERT_EQ(electron.charge_, double(-1.0) * elementary_charge);
  ASSERT_THROW(Particle::Properties(0, 1), std::out_of_range);
  ASSERT_THROW(Particle::Properties(1, 0), std::out_of_range);

  Particle::Properties iron(26, 56);
  ASSERT_EQ(iron.mass_, IsotopeData::get_isotope_mass(26, 56));
  ASSERT_EQ(iron.charge_, double(26) * elementary_charge);

} catch (const std::exception& e) {
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

  // Testing contructor
  Volume::Layer::material_vector empty_vec{};
  ASSERT_THROW(Volume::Layer(std::move(empty_vec), depth),
               std::invalid_argument);
  ASSERT_THROW(Volume::Layer(std::move(empty_vec), double(-1.0) * depth),
               std::invalid_argument);

  // Testing get_depth
  ASSERT_EQ(water_layer.get_depth(), length_quantity(1.0 * angstrom));
} catch (const std::exception& e) {
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
  ASSERT_THROW(Volume(std::vector<Volume::Layer>()), std::invalid_argument);

  // Test get_layer
  ASSERT_EQ(volume.get_layer(length_quantity(1.0 * angstrom)).get_depth(),
            water_layer_1.get_depth());
  ASSERT_EQ(volume.get_layer(length_quantity(2.0 * angstrom)).get_depth(),
            water_layer_1.get_depth());
  ASSERT_EQ(volume.get_layer(length_quantity(250.0 * angstrom)).get_depth(),
            iron_layer.get_depth());
  ASSERT_EQ(volume.get_layer(length_quantity(501.0 * angstrom)).get_depth(),
            water_layer_2.get_depth());
  ASSERT_THROW(volume.get_layer(length_quantity(1001.0 * angstrom)),
               std::out_of_range);
} catch (const std::exception& e) {
  print_exception(e);
  FAIL();
}

//! Simulation class test
TEST_F(CommonTest, simulation_test) try {
  auto electron_stopping_energy = energy_quantity(1.0 * kilo_electron_volt);
  auto incident_particle_properties = Particle::Properties(14, 28);
  bool enable_damage_cascade = true;
  auto ion_stopping_energy = energy_quantity(0.04 * kilo_electron_volt);
  auto ion_displacement_energy = energy_quantity(0.04 * kilo_electron_volt);
  bool log_single_displacement = true;
  size_t divisor_angle_number = 1000;
  size_t flying_distance_number = 1000;
  auto range = length_quantity(80000000.0 * angstrom);
  size_t bombardment_count = 1000;
  bool is_electron = false;
  auto incident_energy = energy_quantity(50.0 * kilo_electron_volt);
  size_t thread_count = 8;

  auto volume =
      Volume({create_water_layer(length_quantity(2.0 * angstrom)),
              create_iron_layer(length_quantity(500.0 * angstrom)),
              create_water_layer(length_quantity(1000.0 * angstrom))});

  Simulation::Settings settings(
      electron_stopping_energy, incident_particle_properties,
      enable_damage_cascade, ion_stopping_energy, ion_displacement_energy,
      log_single_displacement, divisor_angle_number, flying_distance_number,
      range, bombardment_count, is_electron, incident_energy, thread_count);

  ASSERT_EQ(settings.electron_stopping_energy_, electron_stopping_energy);
  ASSERT_EQ(settings.incident_particle_properties_.mass_,
            incident_particle_properties.mass_);
  ASSERT_EQ(settings.incident_particle_properties_.charge_,
            incident_particle_properties.charge_);
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

  Simulation simulation(settings, std::move(volume));

  ASSERT_EQ(
      simulation.get_settings(&Simulation::Settings::electron_stopping_energy_),
      electron_stopping_energy);
  auto simulation_incident_particle_properties = simulation.get_settings(
      &Simulation::Settings::incident_particle_properties_);
  ASSERT_EQ(simulation_incident_particle_properties.charge_,
            incident_particle_properties.charge_);
  ASSERT_EQ(simulation_incident_particle_properties.charge_,
            incident_particle_properties.charge_);
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
} catch (const std::exception& e) {
  print_exception(e);
  FAIL();
}
