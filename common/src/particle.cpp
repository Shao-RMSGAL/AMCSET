/*! \file particle.cpp
 *  \brief Implements particle motion, rotation, and ion/electron interaction logic.
 *
 *  This module captures the runtime behavior of Particle and its derived types,
 *  including scattering, collision bookkeeping, and ion interaction helpers.
 */

#include "amcset_common.h"

#include <glog/logging.h>

#include <boost/math/constants/constants.hpp> 
#include <boost/random/discrete_distribution.hpp> 
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/systems/si/amount.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/systems/si/energy.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass.hpp>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace common {

Particle::Particle(size_t z_number, size_t mass_number,
                   const Simulation &simulation,
                   boost::random::uniform_01<double> &uniform_distribution,
                   boost::random::mt19937 &random_number_generator,
                   ParticleDataLeaf &data_leaf, Bombardment &bombardment,
                   size_t generation, double birth_time_seconds) try
    : properties_(Particle::Properties(z_number, mass_number)),
      coordinates_(std::vector<Coordinate>()), simulation_(simulation),
      velocity_(Velocity(
          Quaternion(0.0,
                     1.0 * sin(simulation_.get_settings(
                               &Simulation::Settings::incident_angle_)),
                     0.0,
                     1.0 * cos(simulation_.get_settings(
                               &Simulation::Settings::incident_angle_))),
          simulation.get_settings(&Simulation::Settings::incident_energy_))),
      coordinate_(Coordinate(angstrom * 0.0, angstrom * 0.0, angstrom * 0.0)),
      generation_(generation), birth_time_seconds_(birth_time_seconds),
      uniform_distribution_(uniform_distribution),
      random_number_generator_(random_number_generator), data_leaf_(data_leaf),
      bombardment_(bombardment) {
  VLOG(1) << "Creating particle with velocity (" << velocity_.quaternion_.r_
          << ", " << velocity_.quaternion_.z_ << ", "
          << velocity_.quaternion_.y_ << ", " << velocity_.quaternion_.z_
          << ").";
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

void Particle::record_collision_point(CollisionPoint point,
                                      double delta_time_seconds) {
  elapsed_time_seconds_ += delta_time_seconds;
  point.time = current_time_seconds();
  bombardment_.add_collision_point(point);
}

velocity_quantity Particle::speed_from_energy(energy_quantity energy,
                                              mass_quantity mass) const {
  return boost::units::sqrt(2.0 * energy / mass);
}

mass_quantity Particle::reduced_mass(mass_quantity mass_1,
                                     mass_quantity mass_2) const {
  return mass_1 * mass_2 / (mass_1 + mass_2);
}

energy_quantity Particle::cm_energy(mass_quantity reduced_mass,
                                    velocity_quantity velocity) const {
  return reduced_mass * boost::units::pow<2>(velocity) / 2.0;
}

Particle &Particle::rotate_inplace(angle_quantity scattering_angle) {
  VLOG(2) << "Beginning rotation with scattering_angle " << scattering_angle
          << " velocity (" << velocity_.quaternion_.r_ << ", "
          << velocity_.quaternion_.x_ << ", " << velocity_.quaternion_.y_
          << ", " << velocity_.quaternion_.z_ << ")";

  auto x = velocity_.quaternion_.x_;
  auto y = velocity_.quaternion_.y_;
  auto z = velocity_.quaternion_.z_;

  if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
    LOG(ERROR) << "Nans detected in bombardment. Position (" << x << ", " << y
               << ", " << z
               << "). Energy: " << velocity_.energy_ / electron_volt << " eV"
               << ". Bombardment ID: " << bombardment_.get_id();
    x = 0.0;
    y = 0.0;
    z = 1.0;
  }

  auto length = std::sqrt(x * x + y * y + z * z);
  if (length == 0.0) {
    throw std::invalid_argument("Input vector has zero length");
  }

  Quaternion rotation_quaternion;
  if (x == 0 && y == 0) {
    rotation_quaternion = Quaternion(0.0, 0.0, 1.0, 0.0);
  } else if (x == 0 && z == 0) {
    rotation_quaternion = Quaternion(0.0, 1.0, 0.0, 0.0);
  } else if (y == 0.0 && z == 0.0) {
    rotation_quaternion = Quaternion(0.0, 0.0, 0.0, 1.0);
  } else {
    rotation_quaternion = Quaternion(0.0, y, x, 0.0);
  }

  const auto alpha_angle =
      radian * 2.0 * pi * (0.5 - uniform_distribution_(random_number_generator_));

  auto alpha_scalar = cos(alpha_angle / 2.0);
  auto alpha_pos = sin(alpha_angle / 2.0);
  auto quaternion_alpha =
      Quaternion(alpha_scalar, x * alpha_pos, y * alpha_pos, z * alpha_pos);

  auto beta_angle = scattering_angle;
  auto beta_scalar = cos(beta_angle / 2.0);
  auto beta_pos = sin(beta_angle / 2.0);
  auto quaternion_beta = Quaternion(
      beta_scalar, rotation_quaternion.x_ * beta_pos,
      rotation_quaternion.y_ * beta_pos, rotation_quaternion.z_ * beta_pos);

  velocity_.quaternion_ =
      (quaternion_alpha * (quaternion_beta * velocity_.quaternion_) *
       quaternion_beta.conjugate_copy()) *
      quaternion_alpha.conjugate_copy();
  velocity_.quaternion_.r_ = 0.0;
  velocity_.quaternion_.normalize_inplace();

  VLOG(2) << "Post-rotation position-only normalized result quaternion ("
          << velocity_.quaternion_.r_ << ", " << velocity_.quaternion_.x_
          << ", " << velocity_.quaternion_.y_ << ", "
          << velocity_.quaternion_.z_ << ")";

  return *this;
}

} // namespace common
} // namespace amcset

