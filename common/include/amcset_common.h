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
 */

#include <boost/units/static_rational.hpp>
#include <functional>
#include <memory>
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

// Standard library
#include <ranges>

// Boost random numbers
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

// Boost asio
#include <boost/asio.hpp>

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

  Coordinate() = delete;

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
  constexpr Coordinate(length_quantity x, length_quantity y,
                       length_quantity z) noexcept
      : x_(x), y_(y), z_(z) {}
};

//! A struct for storing velocity information.
struct Velocity {
  angle_quantity x_angle_;  //!< Angle relative to the x-axis in radians
  angle_quantity z_angle_;  //!< Angle relative to the z-axis in  radians
  energy_quantity energy_;  //!< Energy of the particle in keV

  Velocity() = delete;

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
  constexpr Velocity(angle_quantity x_angle, angle_quantity z_angle,
                     energy_quantity energy) noexcept
      : x_angle_(x_angle), z_angle_(z_angle), energy_(energy) {};
};

class Simulation;
class Volume;
class Layer;

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

    Properties() = delete;

    /*!
     * \brief Constructs a Properties struct using a z number and mass number.
     *
     * This constructor uses the z number and mass number to determine the
     * charge number and exact isotopic mass to the  charge_ and mass_ fields in
     * the Properties struct. This can also represent an electron if the
     * z_number and mass_number as specified as 0.
     *
     * \param z_number The atomic (Z) number of the particle. If z_number is 0,
     * an electron is specified.
     * \param mass_number The mass number of the particle. If mass_number is 0,
     * an electron is specified.
     */
    constexpr Properties(size_t z_number, size_t mass_number) try
        : charge_(z_number > 0 ? double(z_number) * elementary_charge
                               : double(-1.0) * elementary_charge),
          mass_(IsotopeData::get_isotope_mass(z_number, mass_number)) {
    } catch (...) {
      rethrow(EXCEPTION_MESSAGE(""));
    };
  };

  /*!
   * \brief Deleted default constructor for Particle.
   */
  Particle() = delete;

  /*!
   * \brief construct a Particle by passing a reference to a Simulation object.
   *
   * All the information needed to create the particle is stored in the
   * Simulation object, hence no other parameters are needed.
   */
  Particle(const Simulation& simulation,
           boost::random::discrete_distribution<size_t, double> discrete_distribution,
           const boost::random::uniform_01<double>& uniform_distribution,
           boost::random::mt19937& random_number_generator);

  /*!
   * \brief A pure virtual firing function for Particle.
   *
   * This function prevents an instance of Particle from being directly
   * instantiated. Either an Electron or Ion must be instantiated for a Particle
   * type.
   */
  virtual void fire() = 0;

  /*!
   * \brief Default destructor for Particle
   */
  virtual ~Particle() = default;

 protected:
  std::vector<std::unique_ptr<Particle>>
      particles_;  //!< List of cascade particles
  std::vector<Coordinate>
      coordinates_;  //!< List of coordinates of the current particle
  const Simulation& simulation_;  //!< Reference to the simulation object
  Velocity velocity_;             //!< Current velocity of the particle
  Coordinate coordinate_;         //!< Current coordinate of the particle
  boost::random::discrete_distribution<size_t, double>
      discrete_distribution_;  //!< Discrete distribution for Volume particle
                               //!< interaction calculations.
  const boost::random::uniform_01<double>&
      uniform_distribution_;  //!< Reference to uniform distribution between 0
                              //!< and 1 for various uses.
  boost::random::mt19937&
      random_number_generator_;  //!< Reference to random number generator
  std::reference_wrapper<const Layer>
      current_layer_;  //!< Current layer of the particle
};

/*!
 *  \brief A class to represent electrons.
 *
 *  This class contains functions specific to electron bombardment simulation.
 */
class Electron : public Particle {
 public:
  /*!
   * \brief Construct an electron.
   *
   * Accepts a Simulation reference and calls the Particle default constructor.
   * Immediately fires Electron once created.
   *
   * \param simulation The reference to the Simulation object
   */
  Electron(const Simulation& simulation,
           boost::random::discrete_distribution<size_t, double> discrete_distribution,
           const boost::random::uniform_01<double>& uniform_distribution,
           boost::random::mt19937& random_number_generator)
      : Particle(simulation, discrete_distribution, uniform_distribution,
                 random_number_generator) {};

  /*!
   * \brief Fire an electron
   *
   * Fires an electron using electron-specific energy loss functions.
   */
  void fire() override;
};

class Volume;

/*!
 * \brief A class to represent ions.
 *
 * This class contains functions specific to proton bombardment simulation.
 */
class Ion : public Particle {
 public:
  /*!
   * \brief Construct an electron.
   *
   * Accepts a Simulation reference and calls the Particle default constructor.
   * Immediately fires Ion once created.
   *
   * \param simulation The reference to the Simulation object
   */
     Ion(const Simulation& simulation,
      boost::random::discrete_distribution<size_t, double> discrete_distribution,
      const boost::random::uniform_01<double>& uniform_distribution,
      boost::random::mt19937& random_number_generator)
      : Particle(simulation, discrete_distribution, uniform_distribution,
                 random_number_generator) {};

  /*!
   * \brief Fire an ion
   *
   * Fires an ion using ion-specific energy loss functions
   */
  void fire() override;

  /*!
   * \brief Calculate the electronic stopping energy of ions in matter.
   *
   * Using the charge, mass, depth, and substrate volume, this function
   * calculates the electronic stopping energy.
   *
   * \param charge The charge of the particle
   * \param mass The mass of the particle
   * \param depth The depth of the particle
   * \param volume The volume in which the particle is traveling.
   * \return The energy lost from electronic ion interaction with matter.
   */
  energy_quantity electronic_stopping_energy(charge_quantity charge,
                                             mass_quantity mass,
                                             energy_quantity energy,
                                             const Layer& layer) const;

  /*!
   * \brief Returns the interactomic spacing, depending on the current layer.
   */
  length_quantity interatomic_spacing() const;
};

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

  Layer() = delete;

  /*!
   * \brief Constructor for a Layer that accepts a material_vector
   *
   * Accepts a vector and efficiently moves it into the material member.
   * Maximum depth of the layer is also specified.
   *
   * \param material the material_vector containing the Layer information
   * of interest.
   *
   * \param depth The maximum depth the layer extends into the material.
   * This is not the thickness of the material, but rather the maximum
   * depth that the layer extends from the bombardment surface of the
   * volume.
   *
   */
  Layer(material_vector&& material, length_quantity depth,
        mass_density_quantity density);

  /*!
   * \brief Return relative compositions of isotopes in the layer
   *
   * This function returns a transform which is view of the first elements in
   * each std::pair which makes up the elements in material_.
   *
   * \return A transform of the doubles representing the relative compositions
   * of the layer
   */
  auto get_relative_compositions() const
      -> decltype(std::declval<const material_vector&>() |
                  std::views::transform(
                      &std::pair<double, Particle::Properties>::first));

  /*!
   *  \brief Returns a property at the corresponding index
   *
   *  \param index The index of the property to return
   *  \return A reference to the Properties struct at the specified index
   */
  const Particle::Properties& get_property(size_t index) const try {
    return material_.at(index).second;
  } catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
  };

  /*!
   * \brief Returns a relative composition at the corresponding index.
   *
   * \param index The index of the property to return.
   *
   * \return The relative composition of the component at the index.
   */
  double get_relative_composition(size_t index) const try {
    return material_.at(index).first;
  } catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
  }

  /*!
   * \brief Get the maximum depth of the layer.
   *
   * \return The depth of the layer.
   */
  length_quantity get_depth() const noexcept { return depth_; };

  /*!
   * \brief Get the number of components in a layer
   *
   * \return The number of components in a layer
   */
  size_t get_number_of_components() const noexcept { return material_.size(); };

  /*!
   * \brief Get the mass density of the layer
   *
   * \return Return the mass density of the layer
   */
  mass_density_quantity get_mass_density() const noexcept {
    return mass_density_;
  }

  /*!
   * \brief Get the number density of the layer
   *
   * \return Return the number density of the layer
   */
  number_density_quantity get_number_density() const noexcept {
    return number_density_;
  }

 private:
  material_vector material_;
  const length_quantity depth_;
  const mass_density_quantity mass_density_;
  const number_density_quantity number_density_;
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
  Volume() = delete;

  /*!
   * \brief Constructor for a Volume.
   *
   * Builds a Volume from a std::vector<Layer>. The constructor ensures that the
   * layers are provided in order. i.e. the depth of each subsequent layer is
   * strictly greater than the previous one.
   *
   * \param layers An rvalue reference to a vector of layers. Each layer must
   * have a  greater depth value than the previous one.
   */
  Volume(std::vector<Layer>&& layers);

  /*!
   * \brief A function to retrieve the layer corresponding to the depth
   * provided.
   *
   * Provides a reference to the layer corresponding to the provided depth.
   *
   * \param depth The depth at which to return the corresponding layer.
   *
   * \return A reference to the Layer at which depth falls into.
   */
  std::pair<const Layer&, size_t> get_layer(length_quantity depth) const;

  /*!
   * \brief Get a layer using an index
   *
   * \param index The index of the layer to return
   *
   * \return A constant reference to the Layer of interest
   */
  const Layer& get_layer(size_t index) const;

 private:
  const std::vector<Layer> layers_;
};

/*!
 * \brief Representation of a single particle bombardment.
 *
 * This class represents a bombardment of a single particle incident on a
 * volume. It initiates a particle bobardment and manages IO of the bombardment
 * information, including what information is sent to the front end, and what is
 * saved to a file on disk.
 */
class Bombardment {
 public:
  /*!
   * \brief Construct a Bombardment.
   *
   * The Simulation reference contains all the needed information for a
   * bombardment run, so it is the only required parameter.
   *
   * \param simulation The simulation which the Bombardment belongs to.
   * Simulation parameters are extracted from this.
   */
  Bombardment(const Simulation& simulation);

  /*!
   * \brief Run a bombardment;
   *
   * Depending on the simulation parameters, runs a ion or electron bombardment
   * in the Volume contained in the Simulation object.
   */
  void run_bombardment();

 private:
  const Simulation& simulation_;
  std::unique_ptr<Particle> incident_particle_;
  boost::random::mt19937 random_number_generator_;
  boost::random::uniform_01<double> uniform_distribution_;
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
        electron_stopping_energy_;  //!< Threshold electron stopping energy
                                    /*!< The energy at which precision of
                                     * electron stopping calculations should be increased
                                     */
    const struct Particle::Properties
        incident_particle_properties_;  //!< Properties of the incident
                                        //!< particle
    const bool enable_damage_cascade_;  //!< Control whether damange cascade is
                                        //!< enabled
    const energy_quantity ion_stopping_energy_;  //!< Energy at which ions stop
    const energy_quantity
        ion_displacement_energy_;         //!< Energy to displace ions
    const bool log_single_displacement_;  //!< Track single displacements
    /*!< Determine whether to log when a substrate atom is displaced only
     * once.
     */
    const size_t divisor_angle_number_;  //!< Angle segment count
    /*!< used for cross-section calculations, higher value means finer
     * granularity for cross-section integration, but also more computation.
     */
    const size_t
        flying_distance_number_;  //!< Number of flying distances per group
    /*!< Used for electron bombardment simulation optimization. Allows for the
     * use of cross-section calculations on many electron interactions at
     * similar energies.
     */
    const length_quantity range_;     //!< Maximum range of particles
    const size_t bombardment_count_;  //!< Number of bombardments to simulate

    const bool is_electron_;  //!< Determines if bombardment is for an electron

    const energy_quantity incident_energy_;  //!< Incident particle energy

    const size_t
        thread_count_;  //!< The number of threads to use in the simulation

    Settings() = delete;  //!< Default constructor ensures intialization

    /*!
     * \brief Constructor for a Settings struct.
     *
     * By deleting the default constructor and allowing construction using this
     * struct, initialization of all fields is guarenteed. All members of
     * Settings are represented as parameters in this constructor.
     */
    Settings(energy_quantity electron_stopping_energy,
             Particle::Properties incident_particle_properties,
             bool enable_damage_cascade, energy_quantity ion_stopping_energy,
             energy_quantity ion_displacement_energy,
             bool log_single_displacement, size_t divisor_angle_number,
             size_t flying_distance_number, length_quantity range,
             size_t bombardment_count, bool is_electron,
             energy_quantity incident_energy, size_t thread_count);
  };

  Simulation() = delete;

  /*!
   * \brief The constructor for the Simulation class. Accepts a Settings
   * struct to describe simulation settings.
   *
   * This constructor is used to construct a Simulation object. It accepts a
   * Settings struct to configure the simulation.
   *
   * \param settings This is the struct that describes the configuration of the
   * Simulation.
   * \param volume The volume that the Simulation will use to simulate
   * collisions
   */
  Simulation(Settings settings, Volume&& volume);

  /*!
   * \brief A function to access settings.
   *
   * This templated function can be used to access member variables of the
   * Settings settings member variable.
   *
   * \param member The settings field to be accessed.
   * \return A reference to the requested setting.
   * \code{.cpp}
   * Simulation simulation;
   * Simulation::Settings settings;
   *
   * // Configure settings
   *
   * auto electron_stopping_energy =
   * simulation.getSettings(&Simulation::Settings::electron_stopping_energy_);
   * \endcode
   */
  template <typename T>
  constexpr const T& get_settings(T Settings::* member) const noexcept {
    return settings_.*member;
  };

  /*!
   * \brief Run the simulation.
   *
   * This function creates a thread pool with Settings::thread_count_ workers
   * before posting Settings::bombardment_count_ Bombardment run_bombardment()
   * commands.
   */
  void run_simulation();

  /*!
   * \brief Return a string of settings.
   *
   * This function provides the Simulation settings in a textual format.
   *
   * \return A string of settings for the simulation
   */
  std::string print_settings() const;

  /*!
   * \brief Return a constant reference to the simulation volume.
   *
   * \return The volume of the simulation environment.
   */
  const Volume& get_volume() const { return volume_; };

 private:
  const Settings settings_;
  Volume volume_;
  std::vector<Bombardment> bombardments_;
  boost::asio::thread_pool thread_pool_;
};

}  // namespace common
}  // namespace amcset
