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

#pragma once

// Includes{{{

#include <boost/units/systems/si/length.hpp>
#include <memory>

#if defined(_WIN32)
#if defined(EXPORTING_AMCSET)
#define DECLSPEC __declspec(dllexport)
#else
#define DECLSPEC __declspec(dllimport)
#endif
#else // non windows
#define DECLSPEC
#endif

// Standard library
#include <mutex>
#include <ranges>
#include <shared_mutex>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

// Boost random numbers
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

// Boost units
#include <boost/random/discrete_distribution.hpp>
#include <boost/units/static_rational.hpp>

// Boost asio
#include <boost/asio.hpp>

// Boost asio
#include <boost/lockfree/spsc_queue.hpp>

// Local includes
#include "amcset_utilities.h"
#include "isotope_data.h" // }}}
// Tallies for statistics
#include "tallies.h"

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

using namespace amcset::common;
using namespace amcset::common;

// Struct declarations and definitions{{{

//! A struct to store 3D coordinate data
struct Coordinate {
  length_quantity x_; //!< X position in angstroms
  length_quantity y_; //!< Y position in angstroms
  length_quantity z_; //!< Z position in angstroms

  Coordinate() = delete;

  // Assignment copy operator
  // Coordinate operator=(const Coordinate &a) { return Coordinate(a); };

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

struct Quaternion {
  double r_; //!< The scalar portion of a quaterion.
  double x_; //!< The x portion of a quaterion.
  double y_; //!< The y portion of a quaterion.
  double z_; //!< The z portion of a quaterion.

  /*!
   * \brief Uninitialized constructor. Use with caution.
   */
  Quaternion() {};

  /*!
   * \brief Constructor for a quaternion using direct values
   */
  Quaternion(double r, double x, double y, double z) {
    double sum = std::sqrt(r * r + x * x + y * y + z * z);
    r_ = r / sum;
    x_ = x / sum;
    y_ = y / sum;
    z_ = z / sum;
  }

  // Compiler-generated copy and move constructors is sufficient.

  /*!
   * \brief += operator.
   */
  Quaternion &operator+=(const Quaternion &b);

  /*!
   * \brief + operator.
   */
  Quaternion operator+(const Quaternion &b) const;

  /*!
   * \brief *= operator.
   */
  Quaternion &operator*=(const Quaternion &b);

  /*!
   * \brief * operator.
   */
  Quaternion operator*(const Quaternion &b) const;

  /*!
   * \brief *= operator for scalar multiplication.
   */
  Quaternion &operator*=(double b);

  /*!
   * \brief * operator for scalar multiplication.
   */
  Quaternion operator*(double b) const;

  /*!
   * \brief /= operator for scalar division.
   */
  Quaternion &operator/=(double b);

  /*!
   * \brief / operator for scalar division.
   */
  Quaternion operator/(double b) const;

  /*!
   * \brief Noramlize the Quaternion.
   */
  Quaternion &normalize_inplace();

  /*!
   * \brief Create and return a normalized Quaternion.
   */
  Quaternion normalize_copy() const;

  /*!
   * \brief Conjugate the Quaternion.
   */
  Quaternion &conjugate_inplace();

  /*!
   * \brief Create and return a conjugated Quaternion.
   */
  Quaternion conjugate_copy() const;

  /*!
   * \brief Retrieve the magnitude of the quaternion
   */
  double magnitude() const;
};

Quaternion operator+(const Quaternion &a, const Quaternion &b);
// Quaternion operator*(const Quaternion &a, const Quaternion &b);
Quaternion operator*(double a, const Quaternion &b);
Quaternion operator/(const Quaternion &a, double b);

/*! \brief A struct for storing velocity information.
 *
 * A Quaternion is used to represent the orientation information
 * of the velocity:
 *
 * (scalar, x, y, z)
 *
 * where the scalar is related to the rotation of the vector around
 * the x,y,z vector. See (https://en.wikipedia.org/wiki/Quaternion) for more
 * information.
 *
 * The quaternion should be normalized, or close to normalized, while the energy
 * is used to represent the speed of the particle.
 */
struct Velocity {
  Quaternion quaternion_;
  energy_quantity energy_; //!< Energy of the particle in keV

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
  Velocity(Quaternion quaternion, energy_quantity energy) noexcept
      : quaternion_(quaternion), energy_(energy) {};
};

/*!
 * \brief A class for storing data generated during a history. This includes
 * coordinate data. TODO: Add other types of raw data to track.
 *
 * This class holds coordinate data, as well as pointers to ParticleDataLeafs
 * of generated knock-on particles.
 */
template <typename... Types> class DataFrame {
public:
  //!< For batch read operations
  class ReadLock {
  private:
    std::shared_lock<std::shared_mutex> lock;

  public:
    /**
     * @brief Construct a read guard for a DataFrame.
     *
     * Locks the underlying mutex for shared (read) access for the lifetime
     * of this object.
     * @param df The DataFrame being read.
     */
    explicit ReadLock(const DataFrame &df) : lock(df.frame_mutex) {}
  };

  //!< For batch write operations
  class WriteLock {
  private:
    std::unique_lock<std::shared_mutex> lock;

  public:
    /**
     * @brief Construct a write guard for a DataFrame.
     *
     * Locks the underlying mutex for exclusive (write) access for the
     * lifetime of this object.
     * @param df The DataFrame being written to.
     */
    explicit WriteLock(const DataFrame &df) : lock(df.frame_mutex) {}
  };

  /**
   * @brief Read-only access to a named column.
   *
   * The caller must hold a `ReadLock` while accessing the returned
   * reference.
   * @tparam T Column element type.
   * @param column_name Name of the column to access.
   * @param _lock A `ReadLock` instance that keeps the frame locked.
   * @return Const reference to the column vector.
   */
  template <typename T>
  const std::vector<T> &get_column(const std::string &column_name,
                                   const ReadLock &) const {
    return std::get<std::vector<T>>(columns.at(column_name));
  }

  /**
   * @brief Read-write access to a named column.
   *
   * The caller must hold a `WriteLock` while accessing the returned
   * reference.
   * @tparam T Column element type.
   * @param column_name Name of the column to access.
   * @param _lock A `WriteLock` instance that keeps the frame locked.
   * @return Reference to the column vector.
   */
  template <typename T>
  std::vector<T> &get_column(const std::string &column_name,
                             const WriteLock &) {
    return std::get<std::vector<T>>(columns.at(column_name));
  }

  /**
   * @brief Add a new column to the DataFrame.
   *
   * The new column will be default-initialized for existing rows.
   * @tparam T Column element type.
   * @param name Column name.
   */
  template <typename T> void add_column(const std::string &name) {
    std::unique_lock lock(frame_mutex);
    columns[name] = std::vector<T>(num_rows);
  }

  /**
   * @brief Read a single cell value.
   *
   * Acquires a shared lock for the duration of the operation.
   * @tparam T Cell type.
   * @param column_name Column name.
   * @param row Row index.
   * @return Copy of the requested cell value.
   */
  template <typename T>
  T read_cell(const std::string &column_name, size_t row) {
    std::shared_lock lock(frame_mutex);
    auto &column = std::get<std::vector<T>>(columns[column_name]);
    return column[row];
  }

  /**
   * @brief Write a value to a single cell.
   *
   * Acquires an exclusive lock for the duration of the operation.
   * @tparam T Cell type.
   * @param column_name Column name.
   * @param row Row index.
   * @param value Value to store.
   */
  template <typename T>
  void write_cell(const std::string &column_name, size_t row, T value) {
    std::unique_lock lock(frame_mutex);
    auto &column = std::get<std::vector<T>>(columns[column_name]);
    column[row] = value;
  }

  /**
   * @brief Append a new (default-initialized) row to the DataFrame.
   */
  void add_row() {
    std::unique_lock lock(frame_mutex);
    for (auto &[_, column] : columns) {
      std::visit([](auto &vec) { vec.push_back({}); }, column);
    }
    num_rows++;
  }

  /**
   * @brief Return the number of rows in the DataFrame.
   * @return Row count.
   */
  size_t size() const {
    std::shared_lock lock(frame_mutex);
    return num_rows;
  }

  /**
   * @brief Increment the per-frame sputtered atom counter.
   * @param count Number of sputtered atoms to add.
   */
  void add_sputtered(size_t count) {
    std::unique_lock lock(frame_mutex);
    sputtered_count += count;
  }

private:
  using ColumnData = std::variant<std::vector<Types>...>;
  std::unordered_map<std::string, ColumnData> columns;
  mutable std::shared_mutex frame_mutex;
  size_t num_rows = 0;
  size_t sputtered_count = 0;
};

class ParticleDataLeaf {
public:
  //!< Construct a ParticleDataLeaf using a move constructor
  ParticleDataLeaf() = default;

  /**
   * @brief Append a coordinate to this leaf.
   * @param coordinate Coordinate to append (typically a position visited by a
   * particle).
   */
  void add_coordinate(const Coordinate &coordinate) {
    data.push_back(coordinate);
  };

  /**
   * @brief Add a child data leaf (used for storing cascade knock-on data).
   * @param child A child ParticleDataLeaf to be moved into this node.
   */
  void add_child(ParticleDataLeaf &&child) {
    children.emplace_back(
        std::move(child)); // Efficient construction of ParticleDataLeaf right
                           // in the vector.
  };

  /**
   * @brief Return stored coordinates for this leaf.
   * @return Vector of coordinates recorded for the particle.
   */
  const std::vector<Coordinate> &get_data() const { return data; }

  /**
   * @brief Return the child leaves of this node (mutable).
   * @return Vector of child ParticleDataLeaf instances.
   */
  std::vector<ParticleDataLeaf> &get_children() { return children; }

  /**
   * @brief Mark this leaf as representing a sputtered particle.
   */
  void sputter() { sputtered = true; }

  /**
   * @brief Mark this leaf as representing a particle that escaped the
   * volume.
   */
  void escape() { escaped = true; }

private:
  std::vector<ParticleDataLeaf> children;
  std::vector<Coordinate> data;
  bool sputtered = false;
  bool escaped = false;
};

/*!
 * \brief A class for storing various tallies that are collected across
 * simulations.
 *
 * A Simulation object will have several
 */

// }}}

class Simulation;
class Volume;
class Layer;

// Particle class declarations{{{

class ParticleDataLeaf;
class Bombardment;

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
   * A particle has two pieces of information of interest, a charge (Z number)
   * and a mass. Both of these pieces of information are stored as
   */
  struct Properties { // TODO: Store mass number as well.
    const dimensionless_quantity
        charge_; //!< The atomic number of the particle (dimensionless)
    const mass_quantity mass_; //!< The exact mass of the particle (amu)
    const dimensionless_quantity
        mass_number_; //!< The mass number of the particle (dimensionless)

    Properties() = delete;

    /*!
     * \brief Constructs a Properties struct using a z number and mass
     * number.
     *
     * This constructor uses the z number and mass number to determine the
     * charge number and exact isotopic mass to the  charge_ and mass_
     * fields in the Properties struct. This can also represent an electron
     * if the z_number and mass_number as specified as 0.
     *
     * \param z_number The atomic (Z) number of the particle. If z_number is
     * 0, an electron is specified. \param mass_number The mass number of
     * the particle. If mass_number is 0, an electron is specified.
     */
    constexpr Properties(size_t z_number, size_t mass_number) try
        : charge_(z_number > 0
                      ? dimensionless_quantity(static_cast<double>(z_number))
                      : dimensionless_quantity(-1)),
          mass_(IsotopeData::get_isotope_mass(z_number, mass_number)),
          mass_number_(
              dimensionless_quantity(static_cast<double>(mass_number))) {
    } catch (...) {
      rethrow(EXCEPTION_MESSAGE(""));
    };
  };

  /*!
   * \brief Deleted default constructor for Particle.
   */
  Particle() = delete;

  Particle(const Particle &) = default;
  Particle(Particle &&) = delete;
  Particle &operator=(const Particle &) = delete;
  Particle &operator=(Particle &&) = delete;
  /*!
   * \brief construct a Particle by passing a reference to a Simulation
   * object.
   *
   * All the information needed to create the particle is stored in the
   * Simulation object, hence no other parameters are needed.
   */
  Particle(size_t z_number, size_t mass_number, const Simulation &simulation,
           boost::random::uniform_01<double> &uniform_distribution,
           boost::random::mt19937 &random_number_generator,
           ParticleDataLeaf &data_leaf, Bombardment &bomardment);

  /*!
   * \brief A pure virtual firing function for Particle.
   *
   * This function prevents an instance of Particle from being directly
   * instantiated. Either an Electron or Ion must be instantiated for a
   * Particle type.
   */
  virtual void fire() = 0;

  /*!
   * \brief Returns a const reference to the particle coordinate data
   */
  virtual const std::vector<Coordinate> &get_coordinates() const {
    return coordinates_;
  };

  /*!
   * \brief Returns a reference to the ParticleDataLeaf.
   *
   * The ParticleDataLeaf contains data that is generated during a particle
   * history.
   */
  virtual ParticleDataLeaf &get_data_leaf() { return data_leaf_; };

  /*!
   * \brief Default destructor for Particle
   */
  virtual ~Particle() = default;

protected:
  /*!
   * \brief Calculate velocity from energy.
   */
  /**
   * @brief Convert kinetic energy to a velocity quantity for a particle of
   * the given mass.
   * @param energy Particle kinetic energy.
   * @param mass Particle mass.
   * @return Velocity quantity corresponding to the supplied energy.
   */
  virtual velocity_quantity speed_from_energy(energy_quantity energy,
                                              mass_quantity mass) const;
  /*!
   * \brief Calculate reduced mass.
   */
  /**
   * @brief Compute the reduced mass of two masses.
   * @param mass_1 First mass.
   * @param mass_2 Second mass.
   * @return Reduced mass = m1*m2/(m1+m2).
   */
  virtual mass_quantity reduced_mass(mass_quantity mass_1,
                                     mass_quantity mass_2) const;
  /*!
   * \brief Calculate kinetic energy for center of mass coordinates.
   */
  /**
   * @brief Compute center-of-mass kinetic energy from reduced mass and
   * velocity.
   * @param reduced_mass Reduced mass of the interacting pair.
   * @param velocity Relative velocity.
   * @return Center-of-mass kinetic energy.
   */
  virtual energy_quantity cm_energy(mass_quantity reduced_mass,
                                    velocity_quantity velocity) const;

  /*!
   * \brief Rotate a particle.
   *
   * Rotates a particle around an axis perpendicular to its velocity
   * vector to apply a scattering angle, as well as a randomized azimuthal
   * rotation between 0 and 2pi.
   */
  /**
   * @brief Apply a scattering rotation to the particle in-place.
   * @param scattering_angle The scattering angle to apply.
   * @return Reference to this Particle after rotation.
   */
  virtual Particle &rotate_inplace(angle_quantity scattering_angle);

  /*!
   * \brief Set the velocity of a particle
   */
  virtual void set_velocity(Velocity velocity) { velocity_ = velocity; };

  /*!
   * \brief Set the coordinate of a particle
   */
  virtual void set_coordinate(Coordinate coordinate) {
    coordinate_ = coordinate;
  };

  Properties properties_; //!< Properties of the particle
  // std::vector<std::unique_ptr<Particle>>
  //     particles_; //!< List of cascade particles // TODO: Remove permanently
  std::vector<Coordinate>
      coordinates_; //!< List of coordinates of the current particle
  const Simulation &simulation_; //!< Reference to the simulation object
  Velocity velocity_;            //!< Current velocity of the particle
  Coordinate coordinate_;        //!< Current coordinate of the particle
  boost::random::uniform_01<double>
      &uniform_distribution_; //!< Reference to uniform distribution between 0
                              //!< and 1 for various uses.
  boost::random::mt19937
      &random_number_generator_; //!< Reference to random number generator
  ParticleDataLeaf &data_leaf_;  //!< Reference to a ParticleDataLeaf
  Bombardment &bombardment_;     //!< Particle bombardment_

}; // }}}

// Electron class declarations{{{

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
   * Accepts a Simulation reference and calls the Particle default
   * constructor. Immediately fires Electron once created.
   *
   * \param simulation The reference to the Simulation object
   */
  Electron(size_t z_number, size_t mass_number, const Simulation &simulation,
           boost::random::uniform_01<double> &uniform_distribution,
           boost::random::mt19937 &random_number_generator,
           ParticleDataLeaf &data_leaf, Bombardment &bombardment)
      : Particle(z_number, mass_number, simulation, uniform_distribution,
                 random_number_generator, data_leaf, bombardment) {};

  /*!
   * \brief Fire an electron
   *
   * Fires an electron using electron-specific energy loss functions.
   */
  void fire() override;
}; // }}}

class Volume;

// Ion class declarations{{{

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
   * Accepts a Simulation reference and calls the Particle default
   * constructor. Immediately fires Ion once created.
   *
   * \param simulation The reference to the Simulation object
   */
  Ion(size_t z_number, size_t mass_number, const Simulation &simulation,
      boost::random::uniform_01<double> &uniform_distribution,
      boost::random::mt19937 &random_number_generator,
      ParticleDataLeaf &data_leaf, Bombardment &bombardment)
      : Particle(z_number, mass_number, simulation, uniform_distribution,
                 random_number_generator, data_leaf, bombardment) {};

  /*!
   * \brief Fire an ion
   *
   * Fires an ion using ion-specific energy loss functions
   */
  void fire() override;

private:
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
  energy_quantity
  electronic_stopping_energy(dimensionless_quantity charge, mass_quantity mass,
                             energy_quantity energy,
                             number_density_quantity number_density,
                             const Properties &properties) const;

  /*!
   * \brief Returns the interatomic spacing, depending on the current layer.
   *
   * \param The number density of the substance in moles of atoms/m^3.
   *
   * \return The average interatomic spacing
   */
  length_quantity
  interatomic_spacing(number_density_quantity number_density) const;

  /*!
   * \brief Screening potential.
   */
  energy_quantity screening_potential(length_quantity radius,
                                      dimensionless_quantity z_1,
                                      dimensionless_quantity z_2) const;

  /*!
   * \brief Derivative of screening potential with respect to distance.
   */
  voltage_quantity
  screening_potential_derivative(length_quantity radius,
                                 dimensionless_quantity z_1,
                                 dimensionless_quantity z_2) const;

  /*!
   * \brief Universal screening function.
   */
  dimensionless_quantity
  screening_function(dimensionless_quantity reduced_radius) const;

  /*!
   * \brief Universal screening length.
   */
  length_quantity screening_length(dimensionless_quantity z_1,
                                   dimensionless_quantity z_2) const;

  /*!
   * \brief Reduced energy.
   */
  dimensionless_quantity reduced_energy(length_quantity screening_length,
                                        energy_quantity cm_energy,
                                        dimensionless_quantity z_1,
                                        dimensionless_quantity z_2) const;

  /*!
   * \brief free_flying_path_length
   */
  length_quantity
  free_flying_path_length(mass_quantity m_1, mass_quantity m_2,
                          dimensionless_quantity reduced_energy,
                          length_quantity screening_length,
                          number_density_quantity number_density) const;
  /*!
   * \brief Impact parameter
   */
  length_quantity impact_parameter(number_density_quantity number_density,
                                   length_quantity length,
                                   dimensionless_quantity reduced_energy) const;

  /*!
   * \brief Collision diameter.
   */
  length_quantity collision_diameter(dimensionless_quantity z_1,
                                     dimensionless_quantity z_2,
                                     mass_quantity m_c,
                                     velocity_quantity v_0) const;

  /*!
   * \brief Calculate the closest approach distance of an ion.
   */
  length_quantity closest_approach(dimensionless_quantity z_1,
                                   dimensionless_quantity z_2,
                                   energy_quantity cm_energy,
                                   length_quantity impact_param) const;

  /*!
   * \brief Radius of curvature
   */
  length_quantity radius_of_curvature(length_quantity radius,
                                      energy_quantity cm_energy,
                                      dimensionless_quantity z_1,
                                      dimensionless_quantity z_2) const;

  /*!
   * \brief Magic formula correction parameter.
   */
  dimensionless_quantity magic_formula_correction_parameter(
      dimensionless_quantity reduced_energy, length_quantity impact_parameter,
      length_quantity screening_length, length_quantity closest_approach) const;

  /*!
   * \brief Magic formula scattering angle.
   */
  angle_quantity
  magic_formula_scattering_angle(length_quantity impact_parameter,
                                 length_quantity closest_approach,
                                 length_quantity radius_of_curvature,
                                 dimensionless_quantity correction_factor,
                                 length_quantity screening_distance) const;

  /*!
   * \brief Nuclear energy loss.
   */
  energy_quantity nuclear_energy_loss(mass_quantity atom_mass,
                                      mass_quantity target_mass,
                                      energy_quantity incident_energy,
                                      angle_quantity scattering_angle) const;

  /*!
   * \brief Laboratory scattering angle
   */
  angle_quantity laboratory_scattering_angle(angle_quantity cm_scattering_angle,
                                             mass_quantity atom_mass,
                                             mass_quantity target_mass) const;

  // /*! TODO: Likely to remove
  //  * \brief Convert relative angle to absolute angle.
  //  */
  // std::array<angle_quantity, 4>
  // relative_to_absolute_angle(angle_quantity initial_altitude,
  //                            angle_quantity initial_azimuth,
  //                            angle_quantity incident_deflection,
  //                            angle_quantity target_deflection) const;
}; // }}}

// Layer class declarations{{{

/*!
 * \brief A class for representing a single layer of material in the Volume
 * class.
 *
 * A Layer can consist of a single isotope, or many isotopes present at a
 * relative fraction. The layer is assumed to be homogenous, so the
 * probability of interaction with a particular isotope solely depends on the
 * relative fraction.
 */
class Layer final {
public:
  using material_vector =
      std::vector<std::pair<double,
                            Particle::Properties>>; //!< Type for storing
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
  Layer(material_vector &&material, length_quantity depth,
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
      -> decltype(std::declval<const material_vector &>() | std::views::keys);

  /*!
   *  \brief Returns a property at the corresponding index
   *
   *  \param index The index of the property to return
   *  \return A reference to the Properties struct at the specified index
   */
  const Particle::Properties &get_property(size_t index) const try {
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

  // TODO: Return by reference
  /*!
   * \brief Return the discrete distribution used for selecting material
   * interaction
   *
   * \return Return the discrete distribution corresponding to the relative
   * composition of atoms in the layer
   */
  boost::random::discrete_distribution<size_t, double> &
  get_discrete_distribution() const noexcept {
    // Cannot return a const reference because calling a distribution is not
    // const. However, this program does not really care what happens internally
    // to discrete_distribution, so this non-const behavior is accepted, and
    // distribution is marked 'mutable'
    return discrete_distribution_;
  }

  // TODO: Allow for the use of any random number generator
  /*!
   * \brief Return a random index using the random number generator reference
   * passed as an argument.
   *
   * \param random_number_generator The random number generator to select the
   * index. Currently only mt19937.
   */
  size_t
  get_random_index(boost::random::mt19937 &random_number_generator) const {
    return discrete_distribution_(random_number_generator);
  }

private:
  material_vector material_;
  const length_quantity depth_;
  const mass_density_quantity mass_density_;
  const number_density_quantity number_density_;
  mutable boost::random::discrete_distribution<size_t, double>
      discrete_distribution_;
}; // }}}

// Volume class declarations{{{

/*!
 * \brief A class for representing the simulation environment.
 *
 * This class allows for one or more layers of materials to be
 * simulated. A material is composed of relatiave fractions of one or more
 * isotopes of any element.
 */
class Volume final {
public:
  Volume() = delete;

  /*!
   * \brief Constructor for a Volume.
   *
   * Builds a Volume from a std::vector<Layer>. The constructor ensures that
   * the layers are provided in order. i.e. the depth of each subsequent layer
   * is strictly greater than the previous one.
   *
   * \param layers An rvalue reference to a vector of layers. Each layer must
   * have a  greater depth value than the previous one.
   */
  Volume(std::vector<Layer> &&layers);

  // TODO: Improve performance by returning a reference/pointer to Layer
  // instead
  /*!
   * \brief A function to retrieve the layer corresponding to the depth
   * provided.
   *
   * Provides a reference to the layer corresponding to the provided depth.
   *
   * \param depth The depth at which to return the corresponding layer.
   *
   * \return A pair with the first item being a reference to the Layer at
   *          which depth falls into, and the second being the index of the
   *          layer in the layers_ vector.
   */
  const std::pair<const Layer *, size_t> get_layer(length_quantity depth) const;

  /*!
   * \brief Get a layer using an index
   *
   * \param index The index of the layer to return
   *
   * \return A constant reference to the Layer of interest
   */
  const Layer &get_layer(size_t index) const;

  /*!
   * \brief Get the total depth of a volume
   */
  length_quantity get_depth() const;

private:
  const std::vector<Layer> layers_;
}; // }}}

// Bombardment class declarations{{{

/*!
 * \brief Representation of a single particle bombardment.
 *
 * This class represents a bombardment of a single particle incident on a
 * volume. It initiates a particle bombardment and manages IO of the
 * bombardment information, including what information is sent to the front
 * end, and what is saved to a file on disk.
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
   * \param id The ID of the bombardment
   */
  Bombardment(Simulation &simulation, size_t id);

  /*!
   * \brief Run a bombardment;
   *
   * Depending on the simulation parameters, runs a ion or electron
   * bombardment in the Volume contained in the Simulation object.
   */
  void run_bombardment();

  // Tally recording helpers (used by Particle implementations)
  void record_incident(double range_meters);
  void record_displacement(double depth_meters);
  void record_sputtered();
  void record_energy_loss(long double energy_joules);

  /*!
   * \brief Return the ID of the bombardment
   *
   * \return The ID of the bombardment
   */
  size_t get_id() const noexcept { return id_; }

  /*!
   * \brief Returns a reference to the ParticleDataLeaf.
   *
   * The ParticleDataLeaf contains data that is generated during a particle
   * history.
   */
  ParticleDataLeaf &get_data_leaf() { return data_leaf_; };

  /*!
   * \brief Returns a reference to the simulation
   */
  Simulation &get_simulation() { return simulation_; };

private:
  Simulation &simulation_;
  amcset::common::Tallies local_tallies_;
  std::unique_ptr<Particle> incident_particle_;
  boost::random::mt19937 random_number_generator_;
  boost::random::uniform_01<double> uniform_distribution_;
  size_t id_;
  ParticleDataLeaf data_leaf_; //!< A structure for holding data generated
                               //!< during a particle history.

}; // }}}

// Simulation class declarations{{{

/*!
 * \brief A class for initiating, running, and terminating a multi-bombardment
 * simulation.
 *
 * This simulation class manages particle simulations. It accepts a number of
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
        electron_stopping_energy_; //!< Threshold electron stopping energy
                                   /*!< The energy at which precision of
                                    * electron stopping calculations should be increased
                                    */
    const size_t z_number_;        //!< The z number of the incident particle
    const size_t mass_number_;     //! The mass number of the incident particle
    const struct Particle::Properties
        incident_particle_properties_; //!< Properties of the incident
                                       //!< particle
    const bool enable_damage_cascade_; //!< Control whether damage cascade
                                       //!< is enabled
    const energy_quantity ion_stopping_energy_; //!< Energy at which ions stop
    const energy_quantity ion_displacement_energy_; //!< Energy to displace ions
    const bool log_single_displacement_; //!< Track single displacements
    /*!< Determine whether to log when a substrate atom is displaced only
     * once.
     */
    const size_t divisor_angle_number_; //!< Angle segment count
    /*!< used for cross-section calculations, higher value means finer
     * granularity for cross-section integration, but also more computation.
     */
    const size_t
        flying_distance_number_; //!< Number of flying distances per group
    /*!< Used for electron bombardment simulation optimization. Allows for
     * the use of cross-section calculations on many electron interactions
     * at similar energies.
     */
    const length_quantity range_;    //!< Maximum range of particles
    const size_t bombardment_count_; //!< Number of bombardments to simulate

    const bool is_electron_; //!< Determines if bombardment is for an electron

    const energy_quantity incident_energy_; //!< Incident particle energy

    const size_t
        thread_count_; //!< The number of threads to use in the simulation

    const angle_quantity
        incident_angle_; //!< The incident anlge of the particle

    Settings() = delete; //!< Default constructor ensures initialization

    /*!
     * \brief Constructor for a Settings struct.
     *
     * By deleting the default constructor and allowing construction using
     * this struct, initialization of all fields is guaranteed. All members
     * of Settings are represented as parameters in this constructor.
     *
     * TODO: Modify to accept a map instead
     */
    Settings(energy_quantity electron_stopping_energy, size_t z_number,
             size_t mass_number, bool enable_damage_cascade,
             energy_quantity ion_stopping_energy,
             energy_quantity ion_displacement_energy,
             bool log_single_displacement, size_t divisor_angle_number,
             size_t flying_distance_number, length_quantity range,
             size_t bombardment_count, bool is_electron,
             energy_quantity incident_energy, size_t thread_count,
             angle_quantity incident_angle = 0.0 * radian);
  };

  Simulation() = delete;

  using DrawFunction = std::function<void(const std::vector<Coordinate> &,
                                          size_t bombardment_id)>;

  /*!
   * \brief The constructor for the Simulation class. Accepts a Settings
   * struct to describe simulation settings.
   *
   * This constructor is used to construct a Simulation object. It accepts a
   * Settings struct to configure the simulation.
   *
   * \param settings This is the struct that describes the configuration of
   * the Simulation. \param volume The volume that the Simulation will use to
   * simulate collisions
   */
  Simulation(
      Settings settings, Volume &&volume,
      DrawFunction draw_function =
          [](const std::vector<Coordinate> & /*coords*/,
             size_t /*bombardment_id*/) {});

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
  constexpr const T &get_settings(T Settings::*member) const noexcept {
    return settings_.*member;
  };

  /*!
   * \brief Run the simulation.
   *
   * This function creates a thread pool with Settings::thread_count_ workers
   * before posting Settings::bombardment_count_ Bombardment run_bombardment()
   * commands.
   */
  int run_simulation();

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
  const Volume &get_volume() const { return volume_; };

  /*!
   * \brief Return a mutable reference to the Simulation dataframe.
   *
   * This allows bombardment threads to access the simulation tallies and
   * modify them. It can store either integer types or data types.
   */
  DataFrame<double, int> &get_dataframe() { return dataframe_; };

  using DataVariant =
      std::variant<std::vector<double>>; //!< A type for storing data that
                                         //!< will be communicated to the
                                         //!< frontend asynchronously
  /*!
   * \brief Return a mutable reference to the Simulation dataframe.
   *
   * This allows bombardment threads to access the simulation tallies and
   * modify them. It can store either integer types or data types.
   */
  boost::lockfree::spsc_queue<DataVariant *> &get_spsc_queue() {
    return queue_;
  };

  /*!
   * \brief Return the function used for point drawing
   */
  DrawFunction &get_draw_function() { return draw_function_; };

  // Tallies merge API
  void merge_tallies(const Tallies &tallies);
  const Tallies &get_tallies() const noexcept;

private:
  const Settings settings_;
  Volume volume_;
  Tallies tallies_;
  mutable std::mutex tallies_mutex_;
  std::vector<Bombardment> bombardments_;
  boost::asio::thread_pool thread_pool_;
  DataFrame<double, int> dataframe_; //!< Stores simulation tallies contributed
                                     //!< by multiple bombardments.
  boost::lockfree::spsc_queue<DataVariant *>
      queue_;                  //!< A queue for sending data to the frontend
  DrawFunction draw_function_; //!< Function used for drawing to the canvas
}; // }}}

/*!
 * \brief A class that provides a function to generate a random seed.
 *
 * The generate() function in this class attempts to use a hardware random,
 * number generator to create a seed. If this is not possible, then the code
 * falls back on the system time and thread information instead.
 */
class EntropySeeder {
public:
  /**
   * @brief Generate a reasonably unpredictable seed suitable for
   *        initializing RNGs.
   *
   * The implementation prefers hardware entropy sources when available and
   * falls back to a time/thread-based seed when necessary.
   * @return A 32-bit seed value.
   */
  static boost::mt19937::result_type generate();
};

} // namespace common
} // namespace amcset
