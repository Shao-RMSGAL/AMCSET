// Copyright 2024, Texas A&M University
//
// This file is part of AMCSET.
//
// AMCSET is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with AMCSET.
// If not, see <https://www.gnu.org/licenses/>.

#include <gtest/gtest.h>

class CommonTest : public testing::Test {
    protected: 

    CommonTest() {

    }

    void SetUp() override {
    
    }

    void TearDown() override {

    }

};

TEST_F(CommonTest, common_assertions) {
    EXPECT_STRNE("hello", "world");
    EXPECT_EQ( 7 * 6, 42);
    std::cout << "Hi!" << std::endl;
}
