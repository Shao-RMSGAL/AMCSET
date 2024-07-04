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

#include <boost/units/base_units/angle/radian.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <string_view>
#include <vector>

// Boost random numbers
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
  const length_quantity x_;  //!< X position in angstroms
  const length_quantity y_;  //!< Y position in angstroms
  const length_quantity z_;  //!< Z position in angstroms

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
  const angle_quantity x_angle_;  //!< Angle relative to the x-axis in radians
  const angle_quantity z_angle_;  //!< Angle relative to the z-axis in  radians
  const energy_quantity energy_;  //!< Energy of the particle in keV

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
    constexpr Properties(size_t z_number, size_t mass_number)
        : charge_(double(z_number) * elementary_charge),
          mass_(IsotopeData::getIsotopeMass(z_number, mass_number)) {};
  };

 private:
  // const struct Properties properties;
  // const struct Coordinate coordinate;
  // const struct Velocity veloicty;
};

/*!
 * \brief A class for representing the simulation environment.
 *
 * This class allows for one or more layers of materials to be
 * simulated. A material is composed of reltiave fractions of one or more
 * isotopes of any element.
 */
class Volume final {
 public:
  /*!
   * \brief A class for representing a single layer of material in the Volume
   * class.
   *
   * A Layer can consist of a single isotope, or many isotopes present at a
   * relative fraction. The layer is asssumed to be homogenous, so the
   * probability of interaction with a particular isotope solely depends on the
   * relative fraction.
   */
  class Layer final {
   public:
    using material_vector =
        std::vector<std::pair<double,
                              Particle::Properties>>;  //!< Type for storing
                                                       //!< material information

    /*!
     * \brief Constructor for a Layer that accepts a material_vector
     *
     * Accepts a vector and efficiently moves it into the material member.
     *
     * \param material the material_vector containing the Layer information of
     * interest.
     *
     */
    Layer(const material_vector&& material, length_quantity depth);

    /*!
     * \brief Return relative compositions of isotopes in the layer
     *
     * This function returns a vector which is a copy of the first elements in
     * each std::pair which makes up the elements in material_.
     *
     * \return A vector of doubles representing the relative compositions of the
     * layer
     *u
     * TODO: Reimplement to return references to the data inside of material_
     * instead of copying it. This can be done with custom iterators, but
     * requires working with raw pointers, as well as other unpleasant low level
     * C++ stuff. However, it would be more efficient.
     */
    const std::vector<double> get_relative_compositions() const;

    /*!
     *  \brief Returns a property at the corresponding index
     *
     *  \param index The index of the property to return
     */
    const Particle::Properties& get_property(size_t index) const {
      return material_.at(index).second;
    };

    /*!
     * \brief Returns a relative composition at the corresponding index
     *
     * \param index The index of the property to return
     */
    double get_relative_composition(size_t index) const {
      return material_.at(index).first;
    }

   private:
    material_vector material_;
    const length_quantity depth_;
  };

 private:
  const std::vector<Layer> layers;
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
  struct Settings {
    const energy_quantity
        electron_stopping_energy;  //!< Threshold electron stopping energy
                                   /*!< The energy at which precision of
                                    * electron stopping calculations should be increased
                                    */
    const struct Particle::Properties
        incident_particle_properties;  //!< Properties of the incident particle
    const bool
        enable_damage_cascade;  //!< Control whether damange cascade is enabled
    const energy_quantity ion_stopping_energy;  //!< Energy at which ions stop
    const energy_quantity ion_displacement_energy;  //!< Energy to displace ions
    const bool log_single_displacement;  //!< Track single displacements
    /*!< Determine whether to log when a substrate atom is displaced only once.
     */
    const size_t divisor_angle_number;  //!< Angle segment count
    /*!< used for cross-section calculations, higher value means finer
     * granularity for cross-section integration, but also more computation.
     */
    const size_t
        flying_distance_number;  //!< Number of flying distances per group
    /*!< Used for electron bombardment simulation optimization. Allows for the
     * use of cross-section calculations on many electron interactions at
     * similar energies.
     */
    const length_quantity range;     //!< Maximum range of particles
    const size_t bombardment_count;  //!< Number of bombardments to simulate
  };

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
   * \brief A function to access settings.
   *
   * This templated function can be used to access member variables of the
   * Settings settings member variable.
   *
   * \param member The settings field to be accessed.
   * \return A reference to the requested setting.
   *
   * \code{.cpp}
   * Simulation simulation;
   * Simulation::Settings settings;
   *
   * // Configure settings
   *
   * auto electron_stopping_energy =
   * simulation.getSettings(&Simulation::Settings::electron_stopping_energy);
   * \endcode
   */
  template <typename T>
  constexpr const T& getSettings(T Settings::* member) const noexcept {
    return settings.*member;
  };
};

}  // namespace common
}  // namespace amcset
