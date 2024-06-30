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
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

/*! \file amcset_common.h
 *  \brief The main include file for AMCSET
 *
 *  This file contains all of the major structs and classes in AMCSET.
 *  This includes objects for simulation, as well as objects used for
 *  interfacing with the front-end gui.
 *
 */

#include <boost/units/base_units/angle/radian.hpp>
#if defined(_WIN32)
#if defined(EXPORTING_AMCSET)
#define DECLSPEC __declspec(dllexport)
#else
#define DECLSPEC __declspec(dllimport)
#endif
#else  // non windows
#define DECLSPEC
#endif

#pragma once

// Boost random numbers
#include <boost/random/niederreiter_base2.hpp>
#include <boost/random/uniform_01.hpp>

// Boost units
#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/units/io.hpp>
#include <boost/units/physical_dimensions/energy.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/electric_potential.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/plane_angle.hpp>

// Boost math
#include <boost/math/constants/constants.hpp>

/*!
 * \brief The namespace used for all library code related to AMCSET.
 */
namespace amcset {
/*!
 * \brief The namespace specific to AMCSET common.
 *
 * \detail
 * Only code within the AMCSET common library is included in this namespace.
 * For example, all code within amcset_common.h is contained within this
 * namespace.
 */
namespace common {

using namespace boost::units;
namespace si = boost::units::si;
namespace metric = boost::units::metric;
namespace constants = boost::units::si::constants::codata;

using namespace boost::math::double_constants;

using length_quantity = quantity<si::length>;      // Meter
using energy_quantity = quantity<si::energy>;      // Joule
using angle_quantity = quantity<si::plane_angle>;  // Radian

constexpr auto kilo_electron_volt = double(1000) * constants::e * si::volt;
constexpr auto angstrom = metric::angstrom_base_unit::unit_type();
constexpr auto radian = angle::radian_base_unit::unit_type();

std::string common_greeting();  // TODO: Remove this

//! A struct to store 3D coordinate data
struct Coordinate {
  length_quantity x_;  //* X position in angstroms
  length_quantity y_;  //* Y position in angstroms
  length_quantity z_;  //* Z position in angstroms

  /*!
   * \brief Construct a Coordinate struct with \c double parameters.
   *
   * Pass three parameters to describe a 3D coordinate in terms of
   * angstroms.
   *
   * \param x The value to set to x_
   * \param y The value to set to y_
   * \param z The value to set to z_
   */
  constexpr Coordinate(double x, double y, double z)
      : x_(x * angstrom), y_(y * angstrom), z_(z * angstrom) {}
};

//! \brief A struct for storing velocity information.
struct Velocity {
  angle_quantity x_angle_;  //* Angle relative to the x-axis in radians
  angle_quantity z_angle_;  //* Angle relative to the z-axis in  radians
  energy_quantity energy_;  //* Energy of the particle in keV

  /*!
   * \brief Construct a Velocity struct with \c double parameters
   *
   * Pass three parameters, two angle parameters (in radians) and one
   * energy (in kiloelectron volts) parameter, to describe velocity.
   *
   * \param x_angle The value to set to x_angle
   * \param z_angle The value to set to z_angle
   * \param energy The value to set to energy_
   */
  constexpr Velocity(double x_angle, double z_angle, double energy)
      : x_angle_(x_angle * radian),
        z_angle_(z_angle * radian),
        energy_(energy * kilo_electron_volt) {};
};

class Simulation;
class Particle;
class Ion;
class Electron;

}  // namespace common
}  // namespace amcset
