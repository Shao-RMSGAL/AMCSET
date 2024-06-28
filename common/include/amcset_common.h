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


#if defined(_WIN32)
#  if defined(EXPORTING_AMCSET)
#    define DECLSPEC __declspec(dllexport)
#  else
#    define DECLSPEC __declspec(dllimport)
#  endif
#else // non windows
#  define DECLSPEC
#endif

#pragma once

#include <string>
#include <boost/random/uniform_01.hpp>
#include <boost/random/niederreiter_base2.hpp> 

namespace amcset {
    namespace common {
        std::string common_greeting();

        struct Coordinate {
            double x;
            double y;
            double z;

            Coordinate(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
        };

        struct Velocity {
            double z_angle;
            double x_angle;
            double energy;
        };        

    }
}
