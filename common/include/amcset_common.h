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
#include <boost/units/systems/si/mass.hpp>
#include <string_view>
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

#include "amcset_utilities.h"
#include "isotope_data.h"

/*!
 * \brief The namespace used for all library code related to AMCSET.
 */
namespace amcset {
/*!
 * \brief The namespace specific to AMCSET common.
 *
 * Only code within the AMCSET common library is included in this namespace.
 * For example, all code within amcset_common.h is contained within this
 * namespace.
 */
namespace common {

//! A struct to store 3D coordinate data
struct Coordinate {
  length_quantity x_;  //!< X position in angstroms
  length_quantity y_;  //!< Y position in angstroms
  length_quantity z_;  //!< Z position in angstroms

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

//! A struct for storing velocity information.
struct Velocity {
  angle_quantity x_angle_;  //!< Angle relative to the x-axis in radians
  angle_quantity z_angle_;  //!< Angle relative to the z-axis in  radians
  energy_quantity energy_;  //!< Energy of the particle in keV

  /*!
   * \brief Construct a Velocity struct with \c double parameters
   *
   * Pass three parameters, two angle parameters (in radians) and one
   * energy (in kiloelectron volts) parameter, to describe velocity.
   *
   * \param x_angle The value to set to x_angle_
   * \param z_angle The value to set to z_angle_
   * \param energy The value to set to energy_
   */
  constexpr Velocity(double x_angle, double z_angle, double energy)
      : x_angle_(x_angle * radian),
        z_angle_(z_angle * radian),
        energy_(energy * kilo_electron_volt) {};
};

/*!
 * \brief A class for representing simulation particles.
 *
 * A Particle will store its velocity, current coordinate, and a vector with
 * all the coordinates that it has occupied throughout a simulation. In
 * addition, it also contains a Properties struct that stores information
 * about the properties of the particle, including the charge and the mass of
 * the particle.
 */
class Particle {
 public:
  /*!
   * \struct Properties
   * \brief A struct for storing property information associated with a
   * particle.
   *
   * A partcle has two pieces of information of interest, a charge (Z number)
   * and a mass. Both of these pieces of information are stored as
   */
  struct Properties {
    const charge_quantity
        charge_;  //!< The charge of the particle (Elementary charge)
    const mass_quantity mass_;  //!< The exact mass of the particle (amu)

    /*!
     * \brief Constructs a Properties struct using a z number and mass number.
     *
     * This constructor uses the z number and mass number to determine the
     * charge number and exact isotopic mass to the  charge_ and mass_ fields in
     * the Properties struct.
     *
     * \param z_number The atomic (Z) number of the particle.
     * \param mass_number The mass number of the particle
     */
    constexpr Properties(unsigned int z_number, unsigned int mass_number)
        : charge_(double(z_number) * elementary_charge),
          mass_(IsotopeData::getIsotopeMass(z_number, mass_number)) {};
  };

 private:
  // const struct Properties properties;
  // const struct Coordinate coordinate;
  // const struct Velocity veloicty;
};

/*!
 * \brief A class for initating, running, and terminating a multi-bombardment
 * simulation.
 *
 * This simulation class mangages particle simulations. It accepts a number of
 * settings for a simulation, and then performs those simulations, before
 * reporting the results to the front end. For example, the front end may
 * request a 1,000 electron bombardment in Iron, which will cause the back-end
 * to construct a new Simulation with the requested settings (1,000 electrons
 * into Iron) and run it, before returning the results to the front end.
 */
class Simulation final {
 public:
  /*!
   * \brief A struct for storing settings information for instantiation of a
   * Simulation.
   *
   * The Settings struct is used to determine a variety of properties of a
   * simulation, including settings such as the incident particle properties,
   * the incident energy, the substrate composition, cascade information, and
   * many other settings.
   */
  struct Settings {};

 private:
  const Settings settings;

 public:
  /*!
   * \brief The constructor for the Simulation class. Accepts a Settings struct
   * to describe simulation settings.
   *
   * This constructor is used to construct a Simulation object. It accepts a
   * Settings struct to configure the simulation.
   *
   * \param s This is the struct that describes the configuration of the
   * Simulation.
   */
  constexpr Simulation(Settings s) : settings(std::move(s)) {};

  /*!
   * \brief A function to retrieve a \c const reference to the
   * settings member variable.
   */
  template <typename T>
  constexpr const T& getSettings(T Settings::* member) const noexcept {
    return settings.*member;
  };
};

}  // namespace common
}  // namespace amcset
