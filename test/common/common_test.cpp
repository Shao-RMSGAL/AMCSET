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

TEST_F(CommonTest, common_assertions) {
  EXPECT_STRNE("hello", "world");
  EXPECT_EQ(7 * 6, 42);
  std::cout << "Hi!" << std::endl;
}

TEST_F(CommonTest, coordinate_test) {
  Coordinate coordinate(-1.0, 0, 1.0);

  ASSERT_EQ(to_string(coordinate.x_), "-1e-10 m");
  ASSERT_EQ(to_string(coordinate.y_), "0 m");
  ASSERT_EQ(to_string(coordinate.z_), "1e-10 m");
}

TEST_F(CommonTest, velocity_test) {
  Velocity velocity(-pi, pi, 10000);

  ASSERT_EQ(to_string(velocity.x_angle_), "-3.14159 rad");
  ASSERT_EQ(to_string(velocity.z_angle_), "3.14159 rad");
  ASSERT_EQ(to_string(velocity.energy_), "1.60218e-12 m^2 kg s^-2");
}
