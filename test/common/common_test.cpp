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

#include <gtest/gtest.h>

#include "amcset_common.h"

using namespace amcset::common;

class CommonTest : public testing::Test {
 protected:
  CommonTest() {}

  void SetUp() override {}

  void TearDown() override {}
};

// Template test
TEST_F(CommonTest, common_assertions) {
  EXPECT_STRNE("hello", "world");
  EXPECT_EQ(7 * 6, 42);
  std::cout << "Hi!" << std::endl;
}

// Coordinate struct test
TEST_F(CommonTest, coordinate_constructor_test) {
  Coordinate coordinate(-1.0, 0, 1.0);

  ASSERT_EQ(to_string(coordinate.x_), "-1e-10 m");
  ASSERT_EQ(to_string(coordinate.y_), "0 m");
  ASSERT_EQ(to_string(coordinate.z_), "1e-10 m");
}

// Velocity struct test
TEST_F(CommonTest, velocity_constructor_test) {
  Velocity velocity(-pi, pi, 10000);

  ASSERT_EQ(to_string(velocity.x_angle_), "-3.14159 rad");
  ASSERT_EQ(to_string(velocity.z_angle_), "3.14159 rad");
  ASSERT_EQ(to_string(velocity.energy_), "1.60218e-12 m^2 kg s^-2");
}

// TODO: Simulation class test
TEST_F(CommonTest, simulation_test) {
  // Simulation simulation({
  //     .bombardment_count = 100,
  //     .incident_ion = 5.0,
  //     .substrate = "Silicon"
  //     });

  // ASSERT_EQ(simulation.getSettings(&Simulation::Settings::bombardment_count),
  // 100);
}

// Isotope data test
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
