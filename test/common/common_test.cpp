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

/bin/bash: line 1: The: command not found
 * \file common_test.cpp
 *
 * \brief The testing code for the amcset_common library.
 *
 * This testing suite verifies the functionality and correctness of the AMCSET
 * library. It is designed for close to 100% code coverage.
 */

#include <gtest/gtest.h>

#include "amcset_common.h"

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
};

//! Template test. Does not test the library, just a simple demonstration.
TEST_F(CommonTest, common_assertions) {
  EXPECT_STRNE("hello", "world");
  EXPECT_EQ(7 * 6, 42);
  std::cout << "Hi!" << std::endl;
}

//! Coordinate struct test.
TEST_F(CommonTest, coordinate_constructor_test) {
  constexpr Coordinate coordinate(-1.0, 0, 1.0);

  ASSERT_EQ(to_string(coordinate.x_), "-1e-10 m");
  ASSERT_EQ(to_string(coordinate.y_), "0 m");
  ASSERT_EQ(to_string(coordinate.z_), "1e-10 m");
}

//! Velocity struct test.
TEST_F(CommonTest, velocity_constructor_test) {
  constexpr Velocity velocity(-pi, pi, 10000);

  ASSERT_EQ(to_string(velocity.x_angle_), "-3.14159 rad");
  ASSERT_EQ(to_string(velocity.z_angle_), "3.14159 rad");
  ASSERT_EQ(to_string(velocity.energy_), "1.60218e-12 m^2 kg s^-2");
}

//! TODO: Simulation class test
TEST_F(CommonTest, simulation_test) {
  // Simulation simulation({
  //     .bombardment_count = 100,
  //     .incident_ion = 5.0,
  //     .substrate = "Silicon"
  //     });

  // ASSERT_EQ(simulation.getSettings(&Simulation::Settings::bombardment_count),
  // 100);
}

//! Isotope data test
TEST_F(CommonTest, isotope_data_get_isotope_mass_test) {
  ASSERT_EQ(IsotopeData::getIsotopeMass(50, 132),
            double(131.9178267) * atomic_mass_unit);
  ASSERT_EQ(IsotopeData::getIsotopeMass(2, 4),
            double(4.00260325413) * atomic_mass_unit);
  ASSERT_EQ(IsotopeData::getIsotopeMass(118, 295),
            double(295.21624) * atomic_mass_unit);
  ASSERT_THROW(IsotopeData::getIsotopeMass(0, 1), std::out_of_range);
  ASSERT_THROW(IsotopeData::getIsotopeMass(1, 0), std::out_of_range);
  ASSERT_THROW(IsotopeData::getIsotopeMass(1, 100), std::invalid_argument);
}

//! Layer class test
TEST_F(CommonTest, layer_constructor_test) {
  constexpr Particle::Properties hydrogen(1, 1);
  constexpr Particle::Properties oxygen(8, 16);
  const Volume::Layer::material_vector water{{2.0, hydrogen}, {1.0, oxygen}};
  constexpr auto depth = length_quantity(1.0 * angstrom);
  const Volume::Layer water_layer(std::move(water), depth);

  ASSERT_EQ(water_layer.get_relative_compositions(),
            std::vector<double>({2.0 / 3.0, 1.0 / 3.0}));
  ASSERT_EQ(water_layer.get_property(0).mass_, hydrogen.mass_);
  ASSERT_EQ(water_layer.get_relative_composition(0), 2.0 / 3.0);
  ASSERT_EQ(water_layer.get_property(1).mass_, oxygen.mass_);
  ASSERT_EQ(water_layer.get_relative_composition(1), 1.0 / 3.0);

  const Volume::Layer::material_vector empty_vec{};
  ASSERT_THROW(Volume::Layer(std::move(empty_vec), depth),
               std::invalid_argument);
  ASSERT_THROW(Volume::Layer(std::move(empty_vec), double(-1.0) * depth),
               std::invalid_argument);
}

