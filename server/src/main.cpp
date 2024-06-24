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

#include "amcset_server.h"
#include "amcset_common.h"

#include <iostream>
#include <boost/math/constants/constants.hpp>

int main(const int argc, const char* argv[]) {
    std::cout << amcset::server::server_greeting() << std::endl;
    std::cout << amcset::common::common_greeting() << std::endl;
    std::cout << "Value of pi: " << boost::math::constants::pi<double>() 
              << std::endl;
    return 0;
}