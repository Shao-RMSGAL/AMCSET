// Copyright 2025, Texas A&M University
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

/*! \file amcset_common.h
 *  \brief The main include file for the AMCSET common (simulation) library.
 *
 *  This header declares every public type of the simulation core:
 *
 *  - Plain data: Coordinate, CollisionPoint, Quaternion, Velocity.
 *  - Per-history output: ParticleDataLeaf (a tree of particle tracks).
 *  - Transport: Particle (abstract), Electron, Ion.
 *  - Geometry / materials: Layer, Volume.
 *  - Orchestration: Bombardment (one incident particle history) and
 *    Simulation (many bombardments run on a thread pool, plus tallies and
 *    output files).
 *  - EntropySeeder (RNG seeding helper).
 *
 *  Units: all physical quantities are Boost.Units SI quantities (see
 *  amcset_utilities.h). Note that "number density" is stored in mol/m^3 and
 *  is multiplied by Avogadro's number wherever an atom density is needed.
 */

#pragma once

// Includes{{{

// Standard library
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <tuple>
#include <memory>
#include <mutex>
#include <optional>
#include <ranges>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Boost random numbers
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>

// Boost units
#include <boost/units/static_rational.hpp>
#include <boost/units/systems/si/length.hpp>

// Boost asio (thread pool used by Simulation)
#include <boost/asio/post.hpp>
#include <boost/asio/thread_pool.hpp>

// Local includes
#include "amcset_utilities.h"
#include "isotope_data.h"
#include "physics/stopping.h"
#include "simulation_engine.h"
#include "tallies.h"

// }}}

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

// Forward declarations
class Simulation;
class Bombardment;

// Struct declarations and definitions{{{

/*!
 * \brief A 3D position.
 *
 * The coordinate system has its origin at the bombarded surface; +z points
 * into the target (depth), x and y are lateral. Values are SI lengths
 * (meters); divide by \c angstrom to obtain angstroms.
 */
struct Coordinate {
  length_quantity x_; //!< Lateral x position
  length_quantity y_; //!< Lateral y position
  length_quantity z_; //!< Depth below the surface

  Coordinate() = delete;

  /*!
   * \brief Construct a Coordinate from three lengths.
   *
   * \param x The value to set to x_
   * \param y The value to set to y_
   * \param z The value to set to z_
   */
  constexpr Coordinate(length_quantity x, length_quantity y,
                       length_quantity z) noexcept
      : x_(x), y_(y), z_(z) {}
};

/*!
 * \brief One recorded collision event of one particle in a bombardment.
 *
 * Collision points are accumulated by the owning Bombardment, optionally
 * written to a per-bombardment CSV, and pushed to the Simulation's result
 * queue for the GUI to draw.
 */
struct CollisionPoint {
  Coordinate coordinate;    //!< Position of the collision
  double vx = 0.0;          //!< Velocity x-component (m/s, after the collision's energy loss)
  double vy = 0.0;          //!< Velocity y-component (m/s)
  double vz = 0.0;          //!< Velocity z-component (m/s)
  size_t generation = 0;    //!< Generation number: 0 = incident, 1 = primary knock-on, ...
  size_t z_number = 0;      //!< Atomic number (Z) of the particle that generated this point
  double atomic_mass_number = 0.0; //!< Mass number of the particle (dimensionless)
  /*!
   * \brief Time since the incident particle entered the target (s). SRIM
   * mode (trim2013 transport) fills it from the flight lengths and speeds
   * (an estimate: TRIM has no time); otherwise 0.
   */
  double time = 0.0;
  double energy_eV = 0.0;   //!< Particle energy after the collision (eV)
  /*!
   * \brief Path the point belongs to, within one incident particle's
   * trajectory: 0 for the incident particle, one number per moving recoil.
   * Consecutive points of the same track are joined by straight segments.
   */
  std::uint32_t track = 0;

  // --- SRIM mode, trim2013 transport (0 / defaults elsewhere) ---------------
  //! What the point is.
  enum class Kind : std::uint8_t {
    collision = 0, //!< End of a flight (a collision there, unless the next point says otherwise)
    birth = 1,     //!< First point of a track: the ion's entry, a recoil's site
    stop = 2,      //!< The particle came to rest here (last point of its track)
    exit = 3       //!< The particle left the target here (last point of its track)
  };
  Kind kind = Kind::collision;
  static constexpr std::uint32_t no_parent = 0xffffffffu; //!< parent_track/parent_point of a root
  std::uint64_t ion = 0;          //!< 1-based number of the incident particle (0: not set)
  std::uint32_t point_index = 0;  //!< Index of this point in its ion's trajectory
  //! Track and point (index) of the collision that set this track's atom in
  //! motion (no_parent for the incident particle or a source recoil).
  std::uint32_t parent_track = no_parent;
  std::uint32_t parent_point = no_parent;
  //! At a collision: energy just before it, energy given to the struck
  //! atom, and electronic loss over the flight that led to it (eV).
  double energy_before_eV = 0.0;
  double energy_transferred_eV = 0.0;
  double electronic_loss_eV = 0.0;
};

/*!
 * \brief A minimal quaternion type used to store and rotate direction vectors.
 *
 * Directions are stored as pure quaternions (r_ == 0, (x_, y_, z_) a unit
 * vector). Rotations are applied as q * v * q^-1 in Particle::rotate_inplace.
 *
 * \warning The four-argument constructor *normalizes* its input, so it
 * cannot be used to build a non-unit quaternion, and passing all zeros yields
 * NaNs. conjugate_copy() goes through that constructor and therefore also
 * returns a normalized quaternion.
 */
struct Quaternion {
  double r_; //!< The scalar portion of a quaternion.
  double x_; //!< The x portion of a quaternion.
  double y_; //!< The y portion of a quaternion.
  double z_; //!< The z portion of a quaternion.

  /*!
   * \brief Uninitialized constructor. Members are left indeterminate; assign
   * before reading.
   */
  Quaternion() {};

  /*!
   * \brief Construct a quaternion from components and normalize it to unit
   * magnitude.
   */
  Quaternion(double r, double x, double y, double z) {
    double sum = std::sqrt(r * r + x * x + y * y + z * z);
    r_ = r / sum;
    x_ = x / sum;
    y_ = y / sum;
    z_ = z / sum;
  }

  // Compiler-generated copy and move constructors are sufficient.

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
   * \brief Normalize the Quaternion in place.
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

//! Scalar-on-the-left multiplication (the member overload handles q * s).
Quaternion operator*(double a, const Quaternion &b);

/*! \brief Direction of travel plus kinetic energy of a particle.
 *
 * The direction is stored as a pure unit quaternion (0, x, y, z); the
 * kinetic energy determines the speed (see Particle::speed_from_energy).
 */
struct Velocity {
  Quaternion quaternion_;  //!< Direction of travel as (0, x, y, z), unit length
  energy_quantity energy_; //!< Kinetic energy (SI joules; divide by electron_volt for eV)

  Velocity() = delete;

  /*!
   * \brief Construct a Velocity from a direction and a kinetic energy.
   *
   * \param quaternion Direction of travel (should be a unit pure quaternion).
   * \param energy Kinetic energy of the particle.
   */
  Velocity(Quaternion quaternion, energy_quantity energy) noexcept
      : quaternion_(quaternion), energy_(energy) {};
};

/*!
 * \brief Tree of per-particle track data for one bombardment.
 *
 * Each node stores the coordinates visited by one particle and (optionally)
 * labelled state samples; knock-on particles are stored as children. The
 * tree can be serialized with to_json().
 *
 * \note Not thread-safe; each Bombardment owns its own tree and fills it
 * from a single worker thread.
 */
class ParticleDataLeaf {
public:
  //! A labelled snapshot of position, direction and energy.
  struct Sample {
    double x_m = 0.0;
    double y_m = 0.0;
    double z_m = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double vz = 0.0;
    double energy_joules = 0.0;
    std::string label;
  };

  ParticleDataLeaf() = default;

  //! Set identifying metadata (currently unused by the transport code).
  void set_identity(std::string id, std::string parent_id,
                    std::string particle_type, size_t z_number,
                    size_t mass_number) {
    id_ = std::move(id);
    parent_id_ = std::move(parent_id);
    particle_type_ = std::move(particle_type);
    z_number_ = z_number;
    mass_number_ = mass_number;
  }

  const std::string &get_id() const noexcept { return id_; }
  const std::string &get_parent_id() const noexcept { return parent_id_; }
  const std::string &get_particle_type() const noexcept {
    return particle_type_;
  }

  //! Append a visited coordinate to this particle's track.
  void add_coordinate(const Coordinate &coordinate) { data.push_back(coordinate); }

  //! Append a knock-on particle's subtree.
  void add_child(ParticleDataLeaf &&child) { children.emplace_back(std::move(child)); }

  const std::vector<Coordinate> &get_data() const { return data; }
  std::vector<ParticleDataLeaf> &get_children() { return children; }

  //! Record a labelled sample; the stored velocity is the unit direction.
  void add_sample(const Coordinate &coordinate, const Velocity &velocity,
                  std::string label = "state") {
    Sample sample;
    sample.x_m = coordinate.x_.value();
    sample.y_m = coordinate.y_.value();
    sample.z_m = coordinate.z_.value();
    sample.energy_joules = velocity.energy_.value();

    const auto &q = velocity.quaternion_;
    const auto norm = std::sqrt(q.x_ * q.x_ + q.y_ * q.y_ + q.z_ * q.z_);
    if (norm > 0.0) {
      sample.vx = q.x_ / norm;
      sample.vy = q.y_ / norm;
      sample.vz = q.z_ / norm;
    }
    sample.label = std::move(label);
    samples.push_back(sample);
  }

  void sputter() { sputtered = true; } //!< Mark: particle left through the surface
  void escape() { escaped = true; }    //!< Mark: particle passed the deepest layer

  //! Serialize this node and all children to (pretty-printed) JSON.
  std::string to_json(size_t indent = 0) const {
    std::ostringstream ss;
    const auto pad = std::string(indent, ' ');
    ss << pad << "{\n";
    ss << pad << "  \"id\": \"" << id_ << "\",\n";
    ss << pad << "  \"parent_id\": \"" << parent_id_ << "\",\n";
    ss << pad << "  \"type\": \"" << particle_type_ << "\",\n";
    ss << pad << "  \"z\": " << z_number_ << ",\n";
    ss << pad << "  \"mass_number\": " << mass_number_ << ",\n";
    ss << pad << "  \"sputtered\": " << (sputtered ? "true" : "false")
       << ",\n";
    ss << pad << "  \"escaped\": " << (escaped ? "true" : "false")
       << ",\n";
    ss << pad << "  \"samples\": [";
    for (size_t i = 0; i < samples.size(); ++i) {
      if (i != 0) ss << ", ";
      ss << "{\"label\":\"" << samples[i].label << "\",\"x_m\":"
         << samples[i].x_m << ",\"y_m\":" << samples[i].y_m
         << ",\"z_m\":" << samples[i].z_m << ",\"vx\":"
         << samples[i].vx << ",\"vy\":" << samples[i].vy
         << ",\"vz\":" << samples[i].vz << ",\"energy_joules\":"
         << samples[i].energy_joules << "}";
    }
    ss << "],\n";
    ss << pad << "  \"coordinates\": [";
    for (size_t i = 0; i < data.size(); ++i) {
      if (i != 0) ss << ", ";
      ss << "{\"x_m\":" << data[i].x_.value() << ",\"y_m\":"
         << data[i].y_.value() << ",\"z_m\":" << data[i].z_.value()
         << "}";
    }
    ss << "],\n";
    ss << pad << "  \"children\": [";
    for (size_t i = 0; i < children.size(); ++i) {
      if (i != 0) ss << ", ";
      ss << children[i].to_json(indent + 2);
    }
    ss << "]\n";
    ss << pad << "}";
    return ss.str();
  }

private:
  std::string id_;
  std::string parent_id_;
  std::string particle_type_ = "incident";
  size_t z_number_ = 0;
  size_t mass_number_ = 0;
  std::vector<ParticleDataLeaf> children;
  std::vector<Coordinate> data;
  std::vector<Sample> samples;
  bool sputtered = false;
  bool escaped = false;
};

// }}}

// Particle class declarations{{{

/*!
 * \brief Abstract base class for a transported particle.
 *
 * A Particle stores its current position and direction/energy, the list of
 * coordinates it has visited, and its Properties (charge, mass, mass
 * number). It borrows (by reference) the RNG, the owning Bombardment and a
 * ParticleDataLeaf from the caller; those must outlive the Particle.
 *
 * Particles start at the origin travelling at Settings::incident_angle_
 * from the surface normal (in the x-z plane) with Settings::incident_energy_.
 * Knock-on ions override position/velocity after construction.
 */
class Particle {
public:
  /*!
   * \struct Properties
   * \brief Charge, exact isotopic mass and mass number of a particle species.
   *
   * Used both for the projectile and for each isotope in a Layer.
   */
  struct Properties {
    const dimensionless_quantity
        charge_; //!< The atomic number (Z) of the particle; -1 for an electron
    const mass_quantity mass_; //!< The exact isotopic mass (SI kg; NIST table)
    const dimensionless_quantity
        mass_number_; //!< The mass number of the particle (dimensionless)

    Properties() = delete;

    /*!
     * \brief Constructs a Properties struct using a z number and mass
     * number.
     *
     * The exact mass is looked up in IsotopeData. Passing (0, 0) creates an
     * electron (charge -1, electron rest mass).
     *
     * \param z_number The atomic (Z) number of the particle (0 for electron).
     * \param mass_number The mass number of the particle (0 for electron).
     * \throws std::out_of_range if the isotope is not in the database.
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
   * \brief Construct a Particle.
   *
   * \param z_number Atomic number of the particle (0 for an electron).
   * \param mass_number Mass number of the particle (0 for an electron).
   * \param simulation Source of settings, volume and model functions.
   * \param uniform_distribution U(0,1) distribution shared within a bombardment.
   * \param random_number_generator RNG shared within a bombardment.
   * \param data_leaf Track-data node this particle writes into.
   * \param bombardment Owning bombardment (receives tallies and collision points).
   * \param generation 0 for the incident particle, n+1 for knock-ons of generation n.
   * \param birth_time_seconds Simulation time at which this particle was created.
   */
  Particle(size_t z_number, size_t mass_number, const Simulation &simulation,
           boost::random::uniform_01<double> &uniform_distribution,
           boost::random::mt19937 &random_number_generator,
           ParticleDataLeaf &data_leaf, Bombardment &bombardment,
           size_t generation = 0, double birth_time_seconds = 0.0);

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

  /*!
   * \brief Get the coordinate of a particle
   */
  virtual Coordinate get_coordinate() const { return coordinate_; };
  
  /*!
   * \brief Get the velocity of a particle
   */
  virtual Velocity get_velocity() const { return velocity_; };
  
  /**
   * @brief Convert kinetic energy to a (non-relativistic) speed for a
   * particle of the given mass: v = sqrt(2E/m).
   * @param energy Particle kinetic energy.
   * @param mass Particle mass.
   * @return Velocity quantity corresponding to the supplied energy.
   */
  virtual velocity_quantity speed_from_energy(energy_quantity energy,
                                              mass_quantity mass) const;

  // Kinematics helpers (public so that the physics can be unit tested).

  /**
   * @brief Compute the reduced mass of two masses.
   * @param mass_1 First mass.
   * @param mass_2 Second mass.
   * @return Reduced mass = m1*m2/(m1+m2).
   */
  virtual mass_quantity reduced_mass(mass_quantity mass_1,
                                     mass_quantity mass_2) const;
  /**
   * @brief Compute center-of-mass kinetic energy from reduced mass and
   * velocity.
   * @param reduced_mass Reduced mass of the interacting pair.
   * @param velocity Relative velocity.
   * @return Center-of-mass kinetic energy.
   */
  virtual energy_quantity cm_energy(mass_quantity reduced_mass,
                                    velocity_quantity velocity) const;

  /**
   * @brief Apply a scattering rotation to the particle's direction in place.
   *
   * Deflects the direction by \p scattering_angle about an axis
   * perpendicular to it, then rotates about the original direction by a
   * uniformly random azimuth in (-pi, pi]. Consumes one random number.
   * (Before v0.1.0 the polar axis was not perpendicular to the direction;
   * physics fix P1.)
   *
   * @param scattering_angle The (laboratory-frame) scattering angle to apply.
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

protected:
  Properties properties_; //!< Properties of the particle
  std::vector<Coordinate>
      coordinates_; //!< List of coordinates of the current particle
  const Simulation &simulation_; //!< Reference to the simulation object
  Velocity velocity_;            //!< Current velocity of the particle
  Coordinate coordinate_;        //!< Current coordinate of the particle
  size_t generation_ = 0;        //!< Generation (0 incident, >0 knock-ons)
  double birth_time_seconds_ = 0.0;   //!< Simulation time this particle was created
  double elapsed_time_seconds_ = 0.0; //!< Time accumulated along this particle's own path
  std::uint32_t track_id_ = 0;        //!< Path number of this particle within its history
  boost::random::uniform_01<double>
      &uniform_distribution_; //!< Reference to uniform distribution between 0
                              //!< and 1 for various uses.
  boost::random::mt19937
      &random_number_generator_; //!< Reference to random number generator
  ParticleDataLeaf &data_leaf_;  //!< Reference to a ParticleDataLeaf
  Bombardment &bombardment_;     //!< Owning bombardment

  /*!
   * \brief Advance this particle's clock by \p delta_time_seconds, stamp
   * \p point with the resulting time and hand it to the Bombardment.
   */
  void record_collision_point(CollisionPoint point, double delta_time_seconds);
  //! Absolute simulation time of this particle (birth + elapsed).
  double current_time_seconds() const noexcept {
    return birth_time_seconds_ + elapsed_time_seconds_;
  }

}; // }}}

// Electron class declarations{{{

/*!
 *  \brief An electron transported with a Mott elastic-scattering /
 *  continuous-slowing-down model.
 *
 *  Elastic scattering uses tabulated Mott cross-section parameters and
 *  electron screening parameters (CSV files, see
 *  Settings::mott_scattering_parameters_path_ and
 *  Settings::electron_screening_parameters_path_). Inelastic losses
 *  (ionization, bremsstrahlung) are applied continuously along each step.
 *  Electrons do not generate knock-on cascades.
 *
 *  The parameter tables are process-global and loaded once (std::call_once)
 *  on the first Electron::fire(); later simulations with different paths
 *  reuse the first tables.
 */
class Electron : public Particle {
public:
  /*!
   * \brief Construct an electron. See Particle::Particle for parameters.
   * The particle is not fired until fire() is called.
   */
  Electron(size_t z_number, size_t mass_number, const Simulation &simulation,
           boost::random::uniform_01<double> &uniform_distribution,
           boost::random::mt19937 &random_number_generator,
           ParticleDataLeaf &data_leaf, Bombardment &bombardment,
           size_t generation = 0, double birth_time_seconds = 0.0)
      : Particle(z_number, mass_number, simulation, uniform_distribution,
                 random_number_generator, data_leaf, bombardment, generation,
                 birth_time_seconds) {};

  /*!
   * \brief Transport the electron until its energy drops to
   * Settings::electron_stopping_energy_ or it leaves the target.
   *
   * \throws std::runtime_error if the parameter tables cannot be opened.
   */
  void fire() override;

  // Physics helpers (public so that the physics can be unit tested).

  //! Load the Mott/screening tables once per process (relative paths are
  //! resolved against the current working directory).
  static void initialize_tables(const std::string &screening_path,
                               const std::string &mott_path);
  //! Angle of the i-th node of the logarithmic angular grid in (0, pi].
  static double find_angle(size_t i, size_t divisor_count);

  //! Mott differential elastic cross section dsigma/dOmega (per steradian).
  area_quantity mott_differential_cross_section(double theta_rad,
                                                double substrate_charge,
                                                energy_quantity energy) const;
  //! Total elastic cross section integrated over the angular grid.
  area_quantity mott_total_cross_section(energy_quantity energy,
                                         double substrate_charge,
                                         size_t divisor_count) const;
  //! Sample a polar scattering angle (radians) from the Mott distribution.
  double sample_scattering_angle(area_quantity total_cross_section,
                                 double substrate_charge, energy_quantity energy,
                                 size_t divisor_count) const;
  //! Relativistic recoil energy given to a target nucleus of mass \p target_mass.
  energy_quantity elastic_energy_loss(double scattering_angle_rad,
                                      energy_quantity energy,
                                      mass_quantity target_mass) const;
  //! Collision (ionization) stopping power, relativistic Bethe form with the
  //! Joy-Luo low-energy correction of the mean ionization energy.
  stopping_power_quantity ionization_stopping_power(
      energy_quantity energy, double substrate_charge,
      number_density_quantity atomic_density) const;
  //! Radiative (bremsstrahlung) stopping power, high-energy Bethe-Heitler form.
  stopping_power_quantity bremsstrahlung_stopping_power(
      energy_quantity energy, double substrate_charge,
      number_density_quantity atomic_density) const;
  //! Flight length: exponential with mean 1/(N sigma) plus one interatomic spacing.
  length_quantity generate_flying_distance(area_quantity total_cross_section,
                                           number_density_quantity atomic_density) const;
}; // }}}

// Ion class declarations{{{

/*!
 * \brief An ion transported with a binary-collision (BCA) model.
 *
 * Each step: continuous electronic loss over one interatomic spacing, a
 * binary nuclear collision with a randomly chosen target isotope (scattering
 * angle from the ZBL "magic formula" with the selected screening model),
 * optional recursive knock-on (recoil) transport, then a free flight.
 *
 * Generation-0 ions use Settings::custom_electronic_coefficient_/exponent_;
 * recoils (generation >= 1) use the custom_recoil_* values (only relevant for
 * the CustomPowerLaw electronic stopping model).
 */
class Ion : public Particle {
public:
  /*!
   * \brief Construct an ion. See Particle::Particle for parameters.
   * The particle is not fired until fire() is called.
   */
  Ion(size_t z_number, size_t mass_number, const Simulation &simulation,
      boost::random::uniform_01<double> &uniform_distribution,
      boost::random::mt19937 &random_number_generator,
      ParticleDataLeaf &data_leaf, Bombardment &bombardment,
      size_t generation = 0, double birth_time_seconds = 0.0);

  /*!
   * \brief Transport the ion (and, if enabled, its recoil cascade) until it
   * stops (Settings::ion_stopping_energy_), is sputtered (z < 0) or escapes
   * below the deepest layer.
   */
  void fire() override;

  /*!
   * \brief Screening length for testing: Firsov/QuickBasic form when
   * \p use_quickbasic is true, otherwise the "Thomas-Fermi" form used by
   * Simulation::thomas_fermi_screening_length.
   */
  static length_quantity screening_length_for_model(dimensionless_quantity z_1,
                                                   dimensionless_quantity z_2,
                                                   bool use_quickbasic);

private:
  const double custom_electronic_coefficient_; //!< C for the custom power law (primary or recoil)
  const double custom_electronic_exponent_;    //!< n for the custom power law (primary or recoil)

  /*!
   * \brief Electronic energy loss over one interatomic spacing, using the
   * Simulation's selected electronic stopping model.
   *
   * \param charge Projectile atomic number.
   * \param mass Projectile mass.
   * \param energy Projectile kinetic energy.
   * \param number_density Layer number density (mol/m^3).
   * \param properties Target atom properties.
   * \return The energy lost to electrons over the step.
   */
  energy_quantity
  electronic_stopping_energy(dimensionless_quantity charge, mass_quantity mass,
                             energy_quantity energy,
                             number_density_quantity number_density,
                             const Properties &properties) const;

  /*!
   * \brief Returns the interatomic spacing, depending on the current layer.
   *
   * \param number_density Layer number density in mol/m^3 (converted to
   * atoms/m^3 internally).
   *
   * \return The mean interatomic spacing N^(-1/3).
   */
  length_quantity
  interatomic_spacing(number_density_quantity number_density) const;

  /*!
   * \brief Screened Coulomb potential V(r) = Z1 Z2 e^2/(4 pi eps0 r) * phi(r/a)
   * using the Simulation's selected screening function and length.
   */
  energy_quantity screening_potential(length_quantity radius,
                                      dimensionless_quantity z_1,
                                      dimensionless_quantity z_2) const;

  /*!
   * \brief dV/dr of the screened Coulomb potential.
   */
  force_quantity
  screening_potential_derivative(length_quantity radius,
                                 dimensionless_quantity z_1,
                                 dimensionless_quantity z_2) const;

  /*!
   * \brief Screening length for the Simulation's selected screening model.
   */
  length_quantity screening_length(dimensionless_quantity z_1,
                                   dimensionless_quantity z_2) const;

  /*!
   * \brief Reduced (dimensionless) energy eps = a * E_cm / (Z1 Z2 e^2).
   */
  dimensionless_quantity reduced_energy(length_quantity screening_length,
                                        energy_quantity cm_energy,
                                        dimensionless_quantity z_1,
                                        dimensionless_quantity z_2) const;

  /*!
   * \brief Free-flight path length between collisions.
   *
   * Currently always the mean interatomic spacing N^(-1/3); the mass,
   * reduced-energy and screening-length arguments are unused (a
   * high-energy path-length model is not implemented).
   */
  length_quantity
  free_flying_path_length(mass_quantity m_1, mass_quantity m_2,
                          dimensionless_quantity reduced_energy,
                          length_quantity screening_length,
                          number_density_quantity number_density) const;
  /*!
   * \brief Sample an impact parameter (consumes one random number).
   *
   * For reduced energy > 10: b = sqrt(-ln(R) / (pi N L)); otherwise
   * b = sqrt(R / (pi N^(2/3))), with N in atoms/m^3.
   */
  length_quantity impact_parameter(number_density_quantity number_density,
                                   length_quantity length,
                                   dimensionless_quantity reduced_energy) const;

  /*!
   * \brief Classical collision diameter 4 Z1 Z2 e^2 / (m_c v^2).
   * \note Currently unused.
   */
  length_quantity collision_diameter(dimensionless_quantity z_1,
                                     dimensionless_quantity z_2,
                                     mass_quantity m_c,
                                     velocity_quantity v_0) const;

  /*!
   * \brief Distance of closest approach r0, solving
   * 1 - V(r0)/E_cm - (b/r0)^2 = 0 by Newton iteration with a bisection
   * fallback when a sign change is detected.
   */
  length_quantity closest_approach(dimensionless_quantity z_1,
                                   dimensionless_quantity z_2,
                                   energy_quantity cm_energy,
                                   length_quantity impact_param) const;

  /*!
   * \brief Radius of curvature of the trajectory at \p radius (magic formula).
   */
  length_quantity radius_of_curvature(length_quantity radius,
                                      energy_quantity cm_energy,
                                      dimensionless_quantity z_1,
                                      dimensionless_quantity z_2) const;

  /*!
   * \brief Correction term Delta of the ZBL magic formula.
   */
  dimensionless_quantity magic_formula_correction_parameter(
      dimensionless_quantity reduced_energy, length_quantity impact_parameter,
      length_quantity screening_length, length_quantity closest_approach) const;

  /*!
   * \brief ZBL magic-formula angle: returns acos((B + Rc + Delta) / (R0 + Rc)).
   *
   * \note In the ZBL magic formula this quantity is theta_cm / 2 (half the
   * center-of-mass angle). The energy-transfer formula used downstream
   * (Tmax * sin^2(angle)) is consistent with that, but
   * laboratory_scattering_angle() treats the value as the full theta_cm.
   * See docs/KNOWN_ISSUES.md.
   */
  angle_quantity
  magic_formula_scattering_angle(length_quantity impact_parameter,
                                 length_quantity closest_approach,
                                 length_quantity radius_of_curvature,
                                 dimensionless_quantity correction_factor,
                                 length_quantity screening_distance) const;

  /*!
   * \brief Energy transferred to the target atom (binary elastic collision).
   */
  energy_quantity nuclear_energy_loss(mass_quantity atom_mass,
                                      mass_quantity target_mass,
                                      energy_quantity incident_energy,
                                      angle_quantity scattering_angle) const;

public:
  /*!
   * \brief Convert a center-of-mass scattering angle to the laboratory frame.
   */
  angle_quantity laboratory_scattering_angle(angle_quantity cm_scattering_angle,
                                             mass_quantity atom_mass,
                                             mass_quantity target_mass) const;
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
   * \brief Construct a Layer from (relative fraction, isotope) pairs.
   *
   * Fractions are normalized to sum to 1. The number density (mol/m^3) is
   * derived from the mass density and the fraction-weighted isotopic mass.
   *
   * \param material The (fraction, Properties) pairs; moved from.
   * \param depth The depth of the *bottom* of the layer measured from the
   * surface (cumulative), not the layer thickness.
   * \param density Mass density of the layer.
   * \throws std::invalid_argument if \p material is empty, a fraction is
   * negative, the fractions sum to zero, or \p depth is negative.
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

  /*!
   * \brief Return the discrete distribution used for selecting the target
   * isotope of a collision.
   *
   * boost::random::discrete_distribution::operator()(URNG&) is const, so a
   * const reference is sufficient for sampling.
   */
  const boost::random::discrete_distribution<size_t, double> &
  get_discrete_distribution() const noexcept {
    return discrete_distribution_;
  }

  /*!
   * \brief Sample a component index weighted by relative composition.
   *
   * \param random_number_generator The RNG used for sampling.
   *
   * \note Thread-safe: sampling only reads the distribution's alias table
   * (operator() is const), so one Layer may be shared by all worker threads
   * as long as each thread passes its own RNG.
   */
  size_t
  get_random_index(boost::random::mt19937 &random_number_generator) const {
    return discrete_distribution_(random_number_generator);
  }

  /*!
   * \brief Set the displacement energy of each component (same order as the
   * material). A negative value means "use the simulation-wide value".
   * \throws std::invalid_argument if the size does not match.
   */
  void set_displacement_energies(std::vector<energy_quantity> energies);

  //! Displacement energy of component \p index, if one was set.
  std::optional<energy_quantity> get_displacement_energy(size_t index) const {
    if (index < displacement_energies_.size() &&
        displacement_energies_[index] >= 0.0 * electron_volt) {
      return displacement_energies_[index];
    }
    return std::nullopt;
  }

private:
  material_vector material_;
  const length_quantity depth_;
  const mass_density_quantity mass_density_;
  const number_density_quantity number_density_;
  boost::random::discrete_distribution<size_t, double> discrete_distribution_;
  std::vector<energy_quantity> displacement_energies_; //!< Per component; empty = global
}; // }}}

// Volume class declarations{{{

/*!
 * \brief The simulated target: an ordered stack of Layers.
 *
 * Layers are stored top (surface) to bottom; each Layer's depth is the
 * cumulative depth of its lower boundary. A particle deeper than the last
 * layer's depth has escaped the target.
 */
class Volume final {
public:
  Volume() = delete;

  /*!
   * \brief Constructor for a Volume.
   *
   * \param layers Layers ordered from the surface down. Each layer's depth
   * must be greater than or equal to the previous layer's depth.
   * \throws std::invalid_argument if \p layers is empty or out of order.
   */
  Volume(std::vector<Layer> &&layers);

  /*!
   * \brief Find the layer containing \p depth (linear search).
   *
   * A depth exactly on a boundary belongs to the upper layer.
   *
   * \return (pointer to the Layer, index of the layer).
   * \throws std::invalid_argument if depth < 0.
   * \throws std::out_of_range if depth is below the deepest layer.
   */
  const std::pair<const Layer *, size_t> get_layer(length_quantity depth) const;

  /*!
   * \brief Get a layer using an index
   *
   * \param index The index of the layer to return
   *
   * \return A constant reference to the Layer of interest
   * \throws std::out_of_range if \p index >= number of layers.
   */
  const Layer &get_layer(size_t index) const;

  //! Number of layers.
  size_t get_layer_count() const noexcept { return layers_.size(); }

  /*!
   * \brief Get the total depth of a volume
   */
  length_quantity get_depth() const;

private:
  const std::vector<Layer> layers_;
}; // }}}

// Bombardment class declarations{{{

/*!
 * \brief One incident-particle history (and its cascade).
 *
 * A Bombardment owns its own RNG (seeded by EntropySeeder), its own
 * thread-local Tallies and the collision points produced during the
 * history. run_bombardment() fires the incident particle, merges the local
 * tallies into the Simulation, optionally writes a per-bombardment CSV, and
 * pushes the collision points onto the Simulation's result queue.
 *
 * Each Bombardment is run by exactly one worker thread.
 */
class Bombardment {
public:
  /*!
   * \brief Construct a Bombardment.
   *
   * \param simulation The simulation which the Bombardment belongs to.
   * Simulation parameters are extracted from this.
   * \param id The ID of the bombardment
   */
  Bombardment(Simulation &simulation, size_t id);

  /*!
   * \brief Run the bombardment: create and fire an Ion or Electron
   * (Settings::is_electron_), then publish tallies and collision points.
   */
  void run_bombardment();

  //! Which energy-loss histogram a record_energy_loss() call feeds.
  enum class EnergyLossType { electronic, nuclear, total };

  // Tally recording helpers (called by Particle implementations).

  //! Record the final depth of the incident particle (range statistics).
  void record_incident(double range_meters);
  //! Record one displacement (knock-on above threshold) at a depth.
  void record_displacement(double depth_meters);
  //! Record one particle leaving through the surface.
  void record_sputtered();
  /*!
   * \brief Record an energy deposit at a depth.
   *
   * Every call (of any \p type) also feeds Tallies::energy_loss_stats, so
   * that running mean mixes electronic, nuclear and total entries.
   */
  void record_energy_loss(long double energy_joules,
                         double depth_meters = 0.0,
                         EnergyLossType type = EnergyLossType::total);

  /*!
   * \brief Return the ID of the bombardment
   *
   * \return The ID of the bombardment
   */
  size_t get_id() const noexcept { return id_; }

  //! Add a collision point (already stamped with its particle-specific time).
  void add_collision_point(const CollisionPoint &pt) {
    collision_points_.push_back(pt);
  }

  //! A new path number for a particle of this history (0 for the first).
  std::uint32_t next_track_id() noexcept { return next_track_++; }

  //! Move out all collision points recorded so far and clear the buffer.
  std::vector<CollisionPoint> take_collision_points() {
    std::vector<CollisionPoint> tmp = std::move(collision_points_);
    collision_points_.clear();
    return tmp;
  }

  /*!
   * \brief Returns a reference to the ParticleDataLeaf.
   *
   * The ParticleDataLeaf contains data that is generated during a particle
   * history.
   */
  ParticleDataLeaf &get_data_leaf() { return data_leaf_; };
  const ParticleDataLeaf &get_data_leaf() const { return data_leaf_; };

  /*!
   * \brief Returns a reference to the simulation
   */
  Simulation &get_simulation() { return simulation_; };

private:
  Simulation &simulation_;
  std::vector<CollisionPoint> collision_points_; //!< Points recorded during the history
  std::uint32_t next_track_ = 0;                 //!< Path numbers handed out so far
  amcset::common::Tallies local_tallies_;        //!< Merged into the Simulation at the end
  //! The incident particle. Kept alive after the run (together with its
  //! coordinate list), so memory grows with the total number of steps.
  std::unique_ptr<Particle> incident_particle_;
  boost::random::mt19937 random_number_generator_; //!< Seeded by EntropySeeder
  boost::random::uniform_01<double> uniform_distribution_;
  size_t id_;
  ParticleDataLeaf data_leaf_; //!< Track tree of the incident particle and its cascade
}; // }}}

// Simulation class declarations{{{

/*!
 * \brief A multi-bombardment simulation: settings, target, model selection,
 * thread pool, aggregated tallies and output files.
 *
 * Lifecycle:
 *  1. Construct from a Settings and a Volume. All Bombardment objects
 *     (Settings::bombardment_count_) are created up front.
 *  2. run_simulation() posts every bombardment to a Boost.Asio thread pool
 *     with Settings::thread_count_ workers and blocks until all finish (or
 *     until request_cancel() causes not-yet-started bombardments to be
 *     skipped), then writes the histogram CSVs and \c output/tallies.json.
 *  3. While running, another thread (the GUI) may poll
 *     get_completed_bombardment_count(), snapshot_tallies() and
 *     try_pop_bombardment_results().
 *
 * Physics models are chosen once at construction and stored as plain
 * function pointers (see the *_function_ members), so the hot loop does not
 * branch on the model enums.
 */
class Simulation final : public SimulationEngine {
public:
  /*!
   * \brief Immutable configuration of a Simulation.
   *
   * Three constructors are provided for backward compatibility; they differ
   * only in which trailing parameters are accepted:
   *  - (1) ..., histogram paths, [save_per_bombardment_info,
   *    per_bombardment_output_dir], [screening_model, ...] -- every parameter
   *    after \c thread_count has a default.
   *  - (2) ..., histogram paths, screening_model, ... -- per-bombardment
   *    output disabled.
   *  - (3) ..., histogram paths, save_per_bombardment_info,
   *    per_bombardment_output_dir, electron_screening_parameters_path,
   *    mott_scattering_parameters_path, [screening_model, ...] -- the full
   *    form; (1) and (2) delegate to it with the default asset paths
   *    "common/assets/electron_screening_parameters.csv" and
   *    "common/assets/mott_scattering_parameters.csv".
   *
   * All constructors validate their arguments and throw
   * std::invalid_argument (wrapped by rethrow()) on bad input.
   *
   * Relative paths are resolved against the process working directory.
   */
  struct Settings {
    //! Nuclear screening model (function + screening length).
    enum class ScreeningModel {
      ThomasFermi = 0, //!< "Thomas-Fermi" function and length (default)
      QuickBasic = 1,  //!< Thomas-Fermi function with the Firsov length
      ZBL = 2,         //!< ZBL universal function and length
      Moliere = 3,     //!< Moliere function and length
      Custom = 4       //!< Thomas-Fermi function, length scaled by custom_nuclear_screening_scale_
    };
    //! Electronic stopping model.
    enum class ElectronicStoppingModel {
      ZBL = 0,            //!< k_L sqrt(E) (ZBL low-energy form; default)
      LSS = 1,            //!< Lindhard-Scharff velocity-proportional form
      CustomPowerLaw = 2, //!< S_e/N = C * E[keV]^(n/100) * 1e-18 keV cm^2
      BetheBloch = 3,     //!< Relativistic Bethe-Bloch (high-energy particles)
      Builtin = 4,        //!< Lindhard-Scharff / Bethe-Bloch interpolation (all energies)
      SrimTable = 5       //!< SRIM-2013 tables (srim_stopping_table_path_), Bethe-Bloch above
    };

    const energy_quantity electron_stopping_energy_; //!< Electrons stop below this energy
    const size_t z_number_;        //!< Atomic number of the incident particle (0 = electron)
    const size_t mass_number_;     //!< Mass number of the incident particle (0 = electron)
    const struct Particle::Properties incident_particle_properties_; //!< Derived from z/mass number
    const bool enable_damage_cascade_; //!< Follow recoils (ions only)
    const energy_quantity ion_stopping_energy_; //!< Ions stop below this energy
    const energy_quantity ion_displacement_energy_; //!< Recoil-energy threshold to create a knock-on
    const size_t divisor_angle_number_; //!< Number of angular grid nodes for the electron Mott integrals
    const size_t flying_distance_number_; //!< Validated and displayed but currently unused by the physics
    const length_quantity range_;    //!< Upper bound of all depth histograms (also GUI plot extent)
    const size_t bombardment_count_; //!< Number of incident particles to simulate
    const bool is_electron_; //!< true: electron bombardment, false: ion bombardment
    const energy_quantity incident_energy_; //!< Incident kinetic energy
    const size_t thread_count_; //!< Worker threads in the pool
    const angle_quantity incident_angle_; //!< Angle from the surface normal, in the x-z plane
    const size_t displacement_bin_count_;           //!< Bins of the displacement histogram
    const size_t electronic_energy_loss_bin_count_; //!< Bins of the electronic loss histogram
    const size_t nuclear_energy_loss_bin_count_;    //!< Bins of the nuclear loss histogram
    const size_t total_energy_loss_bin_count_;      //!< Bins of the total loss histogram
    const size_t range_bin_count_;                  //!< Bins of the range (and depth) histogram
    const std::string displacement_histogram_output_file_;           //!< CSV path
    const std::string electronic_energy_loss_histogram_output_file_; //!< CSV path
    const std::string nuclear_energy_loss_histogram_output_file_;    //!< CSV path
    const std::string total_energy_loss_histogram_output_file_;      //!< CSV path
    const std::string range_histogram_output_file_; //!< CSV path; empty => "<dir of displacement file>/range_histogram.csv"
    const bool save_per_bombardment_info_;          //!< Write bombardment_<id>.csv per history
    const std::string per_bombardment_output_dir_;  //!< Directory for per-bombardment CSVs
    const std::string electron_screening_parameters_path_; //!< Electron screening table (CSV)
    const std::string mott_scattering_parameters_path_;    //!< Mott parameter table (CSV)
    const ScreeningModel screening_model_;                   //!< Selected nuclear screening model
    const ElectronicStoppingModel electronic_stopping_model_; //!< Selected electronic stopping model
    const double custom_electronic_coefficient_;        //!< Custom power law C (primary ion)
    const double custom_electronic_exponent_;           //!< Custom power law n (primary ion)
    const double custom_recoil_electronic_coefficient_; //!< Custom power law C (recoils)
    const double custom_recoil_electronic_exponent_;    //!< Custom power law n (recoils)
    const double custom_nuclear_screening_scale_;       //!< Screening-length multiplier (Custom model)
    const std::string srim_stopping_table_path_; //!< SRIM stopping tables (ElectronicStoppingModel::SrimTable)

    Settings() = delete; //!< Every field must be provided explicitly

    //! Overload (1); see the struct description.
    Settings(energy_quantity electron_stopping_energy, size_t z_number,
             size_t mass_number, bool enable_damage_cascade,
             energy_quantity ion_stopping_energy,
             energy_quantity ion_displacement_energy,
             size_t divisor_angle_number, size_t flying_distance_number,
             length_quantity range, size_t bombardment_count, bool is_electron,
             energy_quantity incident_energy, size_t thread_count,
             angle_quantity incident_angle = 0.0 * radian,
             size_t displacement_bin_count = 100,
             size_t electronic_energy_loss_bin_count = 100,
             size_t nuclear_energy_loss_bin_count = 100,
             size_t total_energy_loss_bin_count = 100,
             size_t range_bin_count = 100,
             std::string displacement_histogram_output_file =
                 "output/displacement_histogram.csv",
             std::string electronic_energy_loss_histogram_output_file =
                 "output/electronic_loss_histogram.csv",
             std::string nuclear_energy_loss_histogram_output_file =
                 "output/nuclear_loss_histogram.csv",
             std::string total_energy_loss_histogram_output_file =
                 "output/total_loss_histogram.csv",
             bool save_per_bombardment_info = false,
             std::string per_bombardment_output_dir = "output/bombardments/",
             ScreeningModel screening_model = ScreeningModel::ThomasFermi,
             std::string range_histogram_output_file = {},
             ElectronicStoppingModel electronic_stopping_model =
                 ElectronicStoppingModel::ZBL,
             double custom_electronic_coefficient = 1.0,
             double custom_electronic_exponent = 100.0,
             double custom_recoil_electronic_coefficient = 1.0,
             double custom_recoil_electronic_exponent = 100.0,
             double custom_nuclear_screening_scale = 1.0,
             std::string srim_stopping_table_path = {});

    //! Overload (2); see the struct description.
    Settings(energy_quantity electron_stopping_energy, size_t z_number,
             size_t mass_number, bool enable_damage_cascade,
             energy_quantity ion_stopping_energy,
             energy_quantity ion_displacement_energy,
             size_t divisor_angle_number, size_t flying_distance_number,
             length_quantity range, size_t bombardment_count, bool is_electron,
             energy_quantity incident_energy, size_t thread_count,
             angle_quantity incident_angle, size_t displacement_bin_count,
             size_t electronic_energy_loss_bin_count,
             size_t nuclear_energy_loss_bin_count,
             size_t total_energy_loss_bin_count, size_t range_bin_count,
             std::string displacement_histogram_output_file,
             std::string electronic_energy_loss_histogram_output_file,
             std::string nuclear_energy_loss_histogram_output_file,
             std::string total_energy_loss_histogram_output_file,
             ScreeningModel screening_model,
             std::string range_histogram_output_file = {},
             ElectronicStoppingModel electronic_stopping_model =
                 ElectronicStoppingModel::ZBL,
             double custom_electronic_coefficient = 1.0,
             double custom_electronic_exponent = 100.0,
             double custom_recoil_electronic_coefficient = 1.0,
             double custom_recoil_electronic_exponent = 100.0,
             double custom_nuclear_screening_scale = 1.0,
             std::string srim_stopping_table_path = {});

    //! Overload (3), the full form; see the struct description.
    Settings(energy_quantity electron_stopping_energy, size_t z_number,
             size_t mass_number, bool enable_damage_cascade,
             energy_quantity ion_stopping_energy,
             energy_quantity ion_displacement_energy,
             size_t divisor_angle_number, size_t flying_distance_number,
             length_quantity range, size_t bombardment_count, bool is_electron,
             energy_quantity incident_energy, size_t thread_count,
             angle_quantity incident_angle, size_t displacement_bin_count,
             size_t electronic_energy_loss_bin_count,
             size_t nuclear_energy_loss_bin_count,
             size_t total_energy_loss_bin_count, size_t range_bin_count,
             std::string displacement_histogram_output_file,
             std::string electronic_energy_loss_histogram_output_file,
             std::string nuclear_energy_loss_histogram_output_file,
             std::string total_energy_loss_histogram_output_file,
             bool save_per_bombardment_info,
             std::string per_bombardment_output_dir,
             std::string electron_screening_parameters_path,
             std::string mott_scattering_parameters_path,
             ScreeningModel screening_model = ScreeningModel::ThomasFermi,
             std::string range_histogram_output_file = {},
             ElectronicStoppingModel electronic_stopping_model =
                 ElectronicStoppingModel::ZBL,
             double custom_electronic_coefficient = 1.0,
             double custom_electronic_exponent = 100.0,
             double custom_recoil_electronic_coefficient = 1.0,
             double custom_recoil_electronic_exponent = 100.0,
             double custom_nuclear_screening_scale = 1.0,
             std::string srim_stopping_table_path = {});
  };

  // Physics model function signatures {{{

  /*!
   * Electronic energy loss over one step:
   * (Z1, M1, E, N [mol/m^3], step length, target Properties, C, n) -> dE.
   * C and n are only used by the custom power law.
   */
  using ElectronicStoppingFunction = energy_quantity (*)(
      dimensionless_quantity, mass_quantity, energy_quantity,
      number_density_quantity, length_quantity, const Particle::Properties &,
      double, double);
  //! (M1, M2, E, angle) -> energy transferred to the target atom.
  using NuclearStoppingFunction = energy_quantity (*)(
      mass_quantity, mass_quantity, energy_quantity, angle_quantity);
  //! Screening function phi(x) or its derivative dphi/dx, x = r / a.
  using NuclearScreeningFunction = dimensionless_quantity (*)(
      dimensionless_quantity);
  //! (Z1, Z2, scale) -> screening length a. \c scale is only used by the Custom model.
  using NuclearScreeningLengthFunction = length_quantity (*)(
      dimensionless_quantity, dimensionless_quantity, double);

  // Electronic stopping models (implemented in ion.cpp).
  static energy_quantity zbl_electronic_stopping(
      dimensionless_quantity, mass_quantity, energy_quantity,
      number_density_quantity, length_quantity,
      const Particle::Properties &, double, double);
  static energy_quantity lss_electronic_stopping(
      dimensionless_quantity, mass_quantity, energy_quantity,
      number_density_quantity, length_quantity,
      const Particle::Properties &, double, double);
  static energy_quantity custom_power_law_electronic_stopping(
      dimensionless_quantity, mass_quantity, energy_quantity,
      number_density_quantity, length_quantity,
      const Particle::Properties &, double, double);

  // Nuclear energy transfer (implemented in ion.cpp).
  static energy_quantity binary_collision_nuclear_stopping(
      mass_quantity, mass_quantity, energy_quantity, angle_quantity);

  // Screening functions and derivatives (implemented in ion.cpp).
  static dimensionless_quantity thomas_fermi_screening(
      dimensionless_quantity);
  static dimensionless_quantity zbl_screening(dimensionless_quantity);
  static dimensionless_quantity moliere_screening(dimensionless_quantity);
  static dimensionless_quantity thomas_fermi_screening_derivative(
      dimensionless_quantity);
  static dimensionless_quantity zbl_screening_derivative(
      dimensionless_quantity);
  static dimensionless_quantity moliere_screening_derivative(
      dimensionless_quantity);

  // Screening lengths (implemented in ion.cpp).
  static length_quantity thomas_fermi_screening_length(
      dimensionless_quantity, dimensionless_quantity, double);
  static length_quantity custom_screening_length(
      dimensionless_quantity, dimensionless_quantity, double);
  static length_quantity firsov_screening_length(
      dimensionless_quantity, dimensionless_quantity, double);
  static length_quantity zbl_screening_length(
      dimensionless_quantity, dimensionless_quantity, double);
  static length_quantity moliere_screening_length(
      dimensionless_quantity, dimensionless_quantity, double);

  // }}}

  Simulation() = delete;

  /*!
   * \brief Construct a Simulation and pre-allocate all of its Bombardments.
   *
   * \param settings The configuration of the Simulation.
   * \param volume The target the Simulation will use (moved from).
   */
  Simulation(Settings settings, Volume &&volume);

  /*!
   * \brief Read a field of the Settings via a pointer-to-member.
   *
   * \code{.cpp}
   * auto energy = simulation.get_settings(
   *     &Simulation::Settings::electron_stopping_energy_);
   * \endcode
   */
  template <typename T>
  constexpr const T &get_settings(T Settings::*member) const noexcept {
    return settings_.*member;
  };

  /*!
   * \brief Run all bombardments on the thread pool and write output files.
   *
   * Blocks until every posted bombardment has finished. Output-file errors
   * are logged and swallowed (they do not throw).
   *
   * \return Always 0 on normal completion (including cancellation).
   * \warning Call at most once per Simulation: the thread pool is joined at
   * the end of the first call.
   */
  int run_simulation();

  // SimulationEngine interface (thin adapters over the native API).
  int run() override { return run_simulation(); }
  size_t completed_count() const noexcept override {
    return get_completed_bombardment_count();
  }
  size_t total_count() const noexcept override {
    return settings_.bombardment_count_;
  }
  bool try_pop_trajectory(std::vector<CollisionPoint> &out) override {
    return try_pop_bombardment_results(out);
  }
  void clear_trajectories() override { clear_bombardment_results(); }
  length_quantity depth_extent() const override { return settings_.range_; }
  std::string describe() const override { return print_settings(); }

  /*!
   * \brief Electronic stopping cross section of projectile (z1, m1) in
   * element z2 for the curve-based models (BetheBloch, Builtin, SrimTable).
   * \returns nullptr for the analytic models.
   */
  const physics::StoppingCurve *stopping_curve(unsigned z1, mass_quantity m1,
                                               unsigned z2) const;

  //! Number of bombardments that have finished. Thread-safe.
  size_t get_completed_bombardment_count() const noexcept {
    return completed_bombardments_.load(std::memory_order_acquire);
  }

  //! Ask the pool to skip bombardments that have not started yet. Thread-safe.
  void request_cancel() noexcept override {
    cancellation_requested_.store(true, std::memory_order_release);
  }

  //! Whether request_cancel() has been called. Thread-safe.
  bool is_cancelled() const noexcept override {
    return cancellation_requested_.load(std::memory_order_acquire);
  }

  //! Called by the worker after each bombardment. Thread-safe.
  void increment_completed_bombardment_count() noexcept {
    completed_bombardments_.fetch_add(1, std::memory_order_acq_rel);
  }

  //! Human-readable dump of the Settings (used in the log).
  std::string print_settings() const;

  //! The target volume.
  const Volume &get_volume() const { return volume_; };

  // Tallies {{{

  //! Merge one bombardment's tallies into the simulation totals. Thread-safe.
  void merge_tallies(const Tallies &tallies);
  //! Copy of the current totals, taken under the lock. Thread-safe.
  Tallies snapshot_tallies() const override;
  /*!
   * \brief Unsynchronized access to the totals.
   * \warning Only safe when no bombardment is running (e.g. after
   * run_simulation() returns). Use snapshot_tallies() otherwise.
   */
  const Tallies &get_tallies() const noexcept;

  /*!
   * \brief Write a histogram as CSV with columns
   * depth_min_angstrom, depth_max_angstrom, \p value_label.
   *
   * Parent directories are created as needed.
   * \throws std::invalid_argument if \p output_file is empty.
   * \throws std::runtime_error if the file cannot be opened.
   */
  void write_histogram_csv(const std::string &output_file,
                           const Histogram &histogram,
                           const std::string &value_label = "count") const;

  // }}}

  // Result queue for the GUI {{{

  //! Push a completed bombardment's collision points. Thread-safe.
  void push_bombardment_results(std::vector<CollisionPoint> &&vec);

  //! Pop the oldest queued result into \p out; false if empty. Thread-safe.
  bool try_pop_bombardment_results(std::vector<CollisionPoint> &out);

  //! Drop all queued results. Thread-safe.
  void clear_bombardment_results();

  // }}}

private:
  friend class Ion; // Ion calls the selected model function pointers below.

  const Settings settings_;
  const ElectronicStoppingFunction electronic_stopping_function_;
  const NuclearStoppingFunction nuclear_stopping_function_;
  const NuclearScreeningFunction nuclear_screening_function_;
  const NuclearScreeningFunction nuclear_screening_derivative_function_;
  const NuclearScreeningLengthFunction nuclear_screening_length_function_;
  Volume volume_;
  //! Curve-based electronic stopping, keyed by (Z1, M1 in micro-amu, Z2).
  std::map<std::tuple<unsigned, long long, unsigned>, physics::StoppingCurve>
      stopping_curves_;
  Tallies tallies_;                    //!< Aggregated totals (guarded by tallies_mutex_)
  mutable std::mutex tallies_mutex_;
  std::vector<Bombardment> bombardments_; //!< One per incident particle, created up front
  boost::asio::thread_pool thread_pool_;
  std::atomic<size_t> completed_bombardments_{0};
  std::atomic<bool> cancellation_requested_{false};
  //! Completed bombardment collision vectors awaiting the GUI (guarded by
  //! bombardment_queue_mutex_). Unbounded.
  std::deque<std::vector<CollisionPoint>> bombardment_queue_;
  mutable std::mutex bombardment_queue_mutex_;
}; // }}}

/*!
 * \brief Generates RNG seeds.
 *
 * generate() XORs ten values from std::random_device when it reports
 * non-zero entropy; otherwise it falls back to the high-resolution clock
 * XOR a hash of the calling thread id (and logs a warning). Because every
 * Bombardment is seeded in quick succession on one thread, the fallback can
 * produce correlated or identical seeds; see docs/KNOWN_ISSUES.md.
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
