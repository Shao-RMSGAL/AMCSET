// Copyright 2024, Texas A&M University
//
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

// Includes{{{

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

// #include <array>
#include <cmath>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "amcset_utilities.h" // }}}

namespace amcset {
namespace common {

// Layer function implementations{{{
Layer::Layer(material_vector &&material, length_quantity depth,
             mass_density_quantity mass_density) try
    : material_([&material, depth, mass_density]() {
        VLOG(1) << "Creating Layer material with vector of length "
                << material.size() << ", depth " << depth << ", and density "
                << mass_density;
        auto relative_compositions =
            material | std::views::transform(
                           &std::pair<double, Particle::Properties>::first);
        double sum = std::accumulate(relative_compositions.begin(),
                                     relative_compositions.end(), 0.0);
        VLOG(1) << "Sum of compositions is " << sum;

        for (auto &pair : material) {
          if (pair.first < 0) {
            throw std::invalid_argument(EXCEPTION_MESSAGE(
                "Relative composition: " + std::to_string(pair.first) +
                " is less than zero."));
          }
          VLOG(2) << "Initial composition is " << pair.first;
          pair.first /= sum;
          VLOG(2) << "Final composition is " << pair.first;
        }

        return std::move(material);
      }()),
      depth_(depth), mass_density_(mass_density), number_density_([&]() {
        VLOG(1) << "Calculating mass density using " << mass_density;
        return mass_density /
               (std::accumulate(material_.begin(), material_.end(),
                                mass_quantity(0),
                                [](const mass_quantity &sum, const auto &mat) {
                                  VLOG(2)
                                      << "Number density calculation. Sum of "
                                      << sum << ", fraction " << mat.first
                                      << ", and mass of " << mat.second.mass_;
                                  return sum + mat.first * mat.second.mass_;
                                }) *
                constants::N_A);
      }()) {
  VLOG(1) << "Created Layer with " << material_.size()
          << " components, depth of " << depth << ", density of "
          << mass_density << ", and number density of " << number_density_;

  if (depth < length_quantity(0)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE("Depth: " + to_string(depth) +
                                                  " is less than zero."));
  }

  if (material_.size() == 0) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Cannot pass a vector of size 0 to Layer constructor"));
  }

  VLOG(1) << "Passing relative compositions to discrete distribution.";
  auto composition_view = get_relative_compositions();
  discrete_distribution_ = boost::random::discrete_distribution<size_t, double>(
      composition_view.begin(), composition_view.end());

} catch (...) {
  rethrow();
}

auto Layer::get_relative_compositions() const
    -> decltype(std::declval<const material_vector &>() |
                std::views::keys) try {
  VLOG(1) << "Getting relative compositions for layer with " << material_.size()
          << " components";
  // return material_ |
  //        std::views::transform(&std::pair<double,
  //        Particle::Properties>::first);
  return material_ | std::views::keys;
} catch (...) { //!<
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Volume function implementations{{{
Volume::Volume(std::vector<Layer> &&layers) try : layers_(std::move(layers)) {
  VLOG(1) << "Created volume with " << layers_.size() << " layers.";
  if (layers_.size() == 0) {
    throw std::invalid_argument(
        "Vector cannot be empty to Volume constructor.");
  }

  auto prev_depth = layers_.at(0).get_depth();
  for (const auto &layer : layers_) {
    if (layer.get_depth() < prev_depth) {
      throw std::invalid_argument(EXCEPTION_MESSAGE(
          "Provided vector of layers contains a layer which has a "
          "smaller "
          "depth than the previous layer. Depth: " +
          to_string(layer.get_depth()) +
          ". Previous depth: " + to_string(prev_depth)));
    }
  }
} catch (...) {
  rethrow();
}

const Layer &Volume::get_layer(size_t index) const try {
  if (index > layers_.size()) {
    throw std::out_of_range(EXCEPTION_MESSAGE(
        "Index out of range. Index: " + std::to_string(index) +
        ". Vector size: " + std::to_string(layers_.size())));
  }
  return layers_[index];
} catch (...) {
  rethrow();
}

std::pair<const Layer *, size_t> Volume::get_layer(length_quantity depth) const
    try {
  VLOG(1) << "Getting layer at depth " << depth;
  if (depth < length_quantity(0.0 * angstrom)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Depth: " + to_string(depth) + " cannot be less than zero."));
  }

  if (depth > layers_.back().get_depth()) {
    throw std::out_of_range(EXCEPTION_MESSAGE(
        "Depth: " + to_string(depth) + " is deeper than the deepest layer."));
  }

  // Rely on the ordering of the layers (depths are strictly increasing)
  for (size_t i = 0; i < layers_.size(); i++) {
    VLOG(2) << "Checking layer " << i << " with depth "
            << layers_[i].get_depth();
    if (depth <= layers_[i].get_depth()) {
      VLOG(1) << "Found match at " << i << " with depth "
              << layers_[i].get_depth() << " greater than " << depth;
      return std::make_pair(&layers_[i], i);
    }
  }

  // There is a logic error if this is reached
  throw std::logic_error(
      EXCEPTION_MESSAGE("End of layer vector vector reached. Maximum depth: " +
                        to_string(layers_.back().get_depth()) +
                        ". Provided depth: " + to_string(depth)));
} catch (...) {
  rethrow();
} // }}}

// Simulation function implementations{{{
Simulation::Settings::Settings(
    energy_quantity electron_stopping_energy, size_t z_number,
    size_t mass_number, bool enable_damage_cascade,
    energy_quantity ion_stopping_energy,
    energy_quantity ion_displacement_energy, bool log_single_displacement,
    size_t divisor_angle_number, size_t flying_distance_number,
    length_quantity range, size_t bombardment_count, bool is_electron,
    energy_quantity incident_energy, size_t thread_count) try
    : electron_stopping_energy_(electron_stopping_energy), z_number_(z_number),
      mass_number_(mass_number),
      incident_particle_properties_(
          Particle::Properties(z_number, mass_number)),
      enable_damage_cascade_(enable_damage_cascade),
      ion_stopping_energy_(ion_stopping_energy),
      ion_displacement_energy_(ion_displacement_energy),
      log_single_displacement_(log_single_displacement),
      divisor_angle_number_(divisor_angle_number),
      flying_distance_number_(flying_distance_number), range_(range),
      bombardment_count_(bombardment_count), is_electron_(is_electron),
      incident_energy_(incident_energy), thread_count_(thread_count) {
  if (incident_energy <= energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Incident energy should be greater than 0. Value: " +
                          to_string(incident_energy)));
  }
  if (bombardment_count == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Bombardment count must be greater than 0"));
  }
  if (divisor_angle_number == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Divisor angle number must be greater than 0"));
  }
  if (ion_displacement_energy < energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Ion displacement energy cannot be less than 0. Value: " +
        to_string(ion_displacement_energy)));
  }
  if (range < length_quantity(0.0 * angstrom)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Range cannot be less than 0. Value: " + to_string(range)));
  }
  if (flying_distance_number == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Flying distance number must be greater than 0"));
  }
  if (electron_stopping_energy < energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Electron stopping energy cannot be less than 0. Value: " +
        to_string(electron_stopping_energy)));
  }
  if (thread_count == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Thread count must be greater than 0"));
  }
} catch (...) {
  rethrow();
}

Simulation::Simulation(const Settings settings, Volume &&volume) try
    : settings_(settings), volume_(std::move(volume)),
      thread_pool_(settings.thread_count_) {
  LOG(INFO) << "Creating simulation with settings:\n"
            << Simulation::print_settings(); // TODO: Print volume
                                             // information as well.
  std::vector<Bombardment> bombardments;
  bombardments.reserve(settings.bombardment_count_);

  for (size_t i = 0; i < settings.bombardment_count_; i++) {
    VLOG(1) << "Creating bombardment " << i;
    bombardments.push_back(Bombardment(*this, i));
  }
  bombardments_ = std::move(bombardments);
  LOG(INFO) << "Bombardments Created";

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

int Simulation::run_simulation() try {
  LOG(INFO) << "Starting Simulation...";
  for (auto &bombardment : bombardments_) {
    // TODO: Enable multithreading
    // boost::asio::post(thread_pool_, [&bombardment]() {
    LOG(INFO) << "Starting Bombardment...";
    bombardment.run_bombardment();
    LOG(INFO) << "...Bombardment Complete!";
    // });
  }
  // thread_pool_.join();
  LOG(INFO) << "...Simulation complete";
  return 0;
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

std::string Simulation::print_settings() const try {
  std::ostringstream ss;
  ss << "Electron stopping energy: " << settings_.electron_stopping_energy_
     << std::endl;
  ss << "Incident particle z number: " << settings_.z_number_ << std::endl;
  ss << "Incident particle mass number: " << settings_.mass_number_
     << std::endl;
  ss << "Incident particle properties: " << std::endl;
  ss << "\tMass: " << settings_.incident_particle_properties_.mass_
     << std::endl;
  ss << "\tCharge: " << settings_.incident_particle_properties_.charge_
     << std::endl;
  ss << "Enable damage cascade? "
     << (settings_.enable_damage_cascade_ ? "Yes" : "No") << std::endl;
  ss << "Ion stopping energy: " << settings_.ion_stopping_energy_ << std::endl;
  ss << "Ion displacement energy: " << settings_.ion_displacement_energy_
     << std::endl;
  ss << "Log single displacement? "
     << (settings_.log_single_displacement_ ? "Yes" : "No") << std::endl;
  ss << "Divisor angle number: " << settings_.divisor_angle_number_
     << std::endl;
  ss << "Flying distance number: " << settings_.flying_distance_number_
     << std::endl;
  ss << "Range: " << settings_.range_ << std::endl;
  ss << "Bombardment count: " << settings_.bombardment_count_ << std::endl;
  ss << "Is electron? " << (settings_.is_electron_ ? "Yes" : "No") << std::endl;
  ss << "Incident energy: " << settings_.incident_energy_ << std::endl;
  ss << "Thread count: " << settings_.thread_count_ << std::endl;
  return ss.str();
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Bombardment function implementations{{{
Bombardment::Bombardment(const Simulation &simulation, size_t id) try
    : simulation_(simulation), random_number_generator_(),
      uniform_distribution_(), id_(id) {
  VLOG(1) << "Attempting random number: "
          << uniform_distribution_(random_number_generator_);

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

void Bombardment::run_bombardment() try {
  if (simulation_.get_settings(&Simulation::Settings::is_electron_)) {
    VLOG(1) << "Simulating electron bombardment...";
    incident_particle_ = std::make_unique<Electron>(
        simulation_.get_settings(&Simulation::Settings::z_number_),
        simulation_.get_settings(&Simulation::Settings::mass_number_),
        simulation_, uniform_distribution_, random_number_generator_);
    VLOG(1) << "...Electron bombardment complete.";
  } else {
    VLOG(1) << "Simulating ion bombardment";
    incident_particle_ = std::make_unique<Ion>(
        simulation_.get_settings(&Simulation::Settings::z_number_),
        simulation_.get_settings(&Simulation::Settings::mass_number_),
        simulation_, uniform_distribution_, random_number_generator_);
    VLOG(1) << "...Ion bombardment complete.";
  }

  incident_particle_->fire();

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

Particle::Particle(size_t z_number, size_t mass_number,
                   const Simulation &simulation,
                   boost::random::uniform_01<double> &uniform_distribution,
                   boost::random::mt19937 &random_number_generator) try
    : properties_(Particle::Properties(z_number, mass_number)),
      particles_(std::vector<std::unique_ptr<Particle>>()),
      coordinates_(std::vector<Coordinate>()), simulation_(simulation),
      velocity_(Velocity(
          Quaternion(1.0, 1.0, 1.0, 1.0), // Point straight towards z-axis
          simulation.get_settings(&Simulation::Settings::incident_energy_))),
      coordinate_(Coordinate(0, 0, 0)),
      uniform_distribution_(uniform_distribution),
      random_number_generator_(random_number_generator) {
  VLOG(1) << "Creating particle with velocity (" << velocity_.quaternion_.r_
          << ", " << velocity_.quaternion_.z_ << ", "
          << velocity_.quaternion_.y_ << ", " << velocity_.quaternion_.z_
          << ").";
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
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

// Ensure this does not fail with quaternions aligned with axes, nan components,
// and zero-length vector.
Particle &Particle::rotate_inplace(angle_quantity scattering_angle) {
  VLOG(2) << "Beginning rotation with scattering_angle " << scattering_angle
          << " velocity (" << velocity_.quaternion_.r_ << ", "
          << velocity_.quaternion_.x_ << ", " << velocity_.quaternion_.y_
          << ", " << velocity_.quaternion_.z_ << "). ";

  auto x = velocity_.quaternion_.x_;
  auto y = velocity_.quaternion_.y_;
  auto z = velocity_.quaternion_.z_;

  if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
    throw std::invalid_argument("Input vector contains NaNs");
  }

  auto length = std::sqrt(x * x + y * y + z * z);
  if (length == 0.0) {
    throw std::invalid_argument("Input vector has zero length");
  }

  Quaternion rotation_quaternion;

  // Special cases: TODO: Add [[likely]] and [[unlikely]] optimizations
  if (x == 0 && y == 0) { // z-axis aligned
    rotation_quaternion = Quaternion(0.0, 0.0, 1.0, 0.0);
  } else if (x == 0 && z == 0) { // y-axis aligned
    rotation_quaternion = Quaternion(0.0, 1.0, 0.0, 0.0);
  } else if (y == 0.0 && z == 0.0) { // x-axis aligned
    rotation_quaternion = Quaternion(0.0, 0.0, 0.0, 1.0);
  } else { // Non-axis aligned, manual alignment instead.
    rotation_quaternion = Quaternion(0.0, y, x,
                                     0.0); // Direction is not important,
                                           // rotation has random alpha tilt.
  }

  const auto alpha_angle =
      radian * 2.0 * pi *
      (0.5 - uniform_distribution_(random_number_generator_));

  auto alpha_scalar = cos(alpha_angle / 2.0);
  auto alpha_pos = sin(alpha_angle / 2.0);
  auto quaternion_alpha =
      Quaternion(alpha_scalar, x * alpha_pos, y * alpha_pos, z * alpha_pos);

  // TODO: Try pythagorean identity to make more efficient.
  auto beta_angle = scattering_angle;

  // beta_angle = radian * 0.001; // TODO: Temp, remove

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

dimensionless_quantity
Ion::screening_function(dimensionless_quantity reduced_radius) const {
  const auto result = 0.1818 * exp(-3.2 * reduced_radius) +
                      0.5099 * exp(-0.9423 * reduced_radius) +
                      0.2802 * exp(-0.4028 * reduced_radius) +
                      0.2817 * exp(-0.2016 * reduced_radius);
  VLOG(2) << "Screening function with reduced radius " << reduced_radius
          << " is " << result;
  return result;
}

// TODO: Find a way to cache the results of this function for speedup. Consider
// calculating ahead of time and storing the results in a data structure
// constructed at runtime for each Layer material.
length_quantity Ion::screening_length(dimensionless_quantity z_1,
                                      dimensionless_quantity z_2) const {
  const auto result = 0.25 * cbrt(9.0 * pi * pi / 2.0) * bohr_radius /
                      (pow(z_1, 0.23) + pow(z_2, 0.23));
  VLOG(2) << "Screening length with z_1 " << z_1 << " and z_2 " << z_2 << " is "
          << result;
  return result;
}

energy_quantity Ion::screening_potential(length_quantity radius,
                                         dimensionless_quantity z_1,
                                         dimensionless_quantity z_2) const {
  auto reduced_length = radius / screening_length(z_1, z_2);
  auto result = screening_function(reduced_length) * z_1 * z_2 * e_statcoulomb *
                e_statcoulomb / radius;
  VLOG(2) << "Screening potential with radius " << radius << ", z_1 " << z_1
          << ", and z_2 " << z_2 << " is " << result;
  return result;
}

voltage_quantity
Ion::screening_potential_derivative(length_quantity radius,
                                    dimensionless_quantity z_1,
                                    dimensionless_quantity z_2) const {
  // TODO: Consolidate calls to screening_length to above-mentioned method.
  // Optimize this function as well.
  auto s_length = screening_length(z_1, z_2);
  auto reduced_radius = radius / s_length;
  auto result =
      (0.1818 * exp(-3.2 * reduced_radius) * (-3.2 / s_length) +
       0.5099 * exp(-0.9423 * reduced_radius) * (-0.9423 / s_length) +
       0.2802 * exp(-0.4028 * reduced_radius) * (-0.4028 / s_length) +
       0.2817 * exp(-0.2016 * reduced_radius) * (-0.2016 / s_length)) *
          z_1 * z_2 * e_statcoulomb * e_statcoulomb / radius -
      (0.1818 * exp(-3.2 * reduced_radius) +
       0.5099 * exp(-0.9423 * reduced_radius) +
       0.2802 * exp(-0.4028 * reduced_radius) +
       0.2817 * exp(-0.2016 * reduced_radius)) *
          z_1 * z_2 * e_statcoulomb * e_statcoulomb / (radius * radius);
  VLOG(2) << "Screening potential derivative with radius " << radius << ", z_1 "
          << z_1 << ", and z_2 " << z_2 << " is " << result;
  return result;
}

dimensionless_quantity Ion::reduced_energy(length_quantity screening_length,
                                           energy_quantity cm_energy,
                                           dimensionless_quantity z_1,
                                           dimensionless_quantity z_2) const {
  auto result = screening_length * cm_energy /
                (z_1 * z_2 * e_statcoulomb * e_statcoulomb);
  VLOG(2) << "Reduced energy with screening length " << screening_length
          << ", center of mass energy " << cm_energy << ", z_1 " << z_1
          << ", z_2 " << z_2 << ", result is " << result;
  return result;
}

// TODO: Verify that dividing by Avogadro's number is correct.
length_quantity
Ion::free_flying_path_length(mass_quantity m_1, mass_quantity m_2,
                             dimensionless_quantity reduced_energy,
                             length_quantity screening_length,
                             number_density_quantity number_density) const {
  length_quantity L;
  if (reduced_energy > 100) {
    L = (0.02 * pow(1.0 + (m_1 + m_2) / atomic_mass_unit, 2) * reduced_energy *
             reduced_energy +
         0.1 * pow(reduced_energy, 1.38)) /
        (4.0 * pi * screening_length * screening_length * number_density *
         log(1.0 + reduced_energy)) /
        si::constants::codata::N_A;
    VLOG(3) << "High energy free flying path length with m_1 " << m_1 << ", m_2"
            << m_2 << ", reduced_energy " << reduced_energy
            << ", screening_length " << screening_length
            << ", and number_density " << number_density << " with result "
            << L;
  } else {
    L = pow<static_rational<1, 3>>(1.0 / number_density /
                                   si::constants::codata::N_A);
    VLOG(3) << "Low energy free flying path length with m_1 " << m_1 << ", m_2 "
            << m_2 << ", reduced_energy " << reduced_energy
            << ", screening_length " << screening_length
            << ", and number_density " << number_density << " with result "
            << L;
  }
  return L;
}

length_quantity Ion::collision_diameter(dimensionless_quantity z_1,
                                        dimensionless_quantity z_2,
                                        mass_quantity rd_mass,
                                        velocity_quantity speed) const {
  auto result =
      4 * z_1 * z_2 * e_statcoulomb * e_statcoulomb / (rd_mass * speed * speed);
  VLOG(2) << "Collision diameter given z_1 " << z_1 << ", z_2 " << z_2
          << " reduced mass " << rd_mass << ", and speed " << speed << " is "
          << result;
  return result;
}

// TODO: Fix magnitude error. This code reports closest approach on the
// order of meters. Find out why.
length_quantity Ion::closest_approach(dimensionless_quantity z_1,
                                      dimensionless_quantity z_2,
                                      energy_quantity cm_energy,
                                      length_quantity impact_param) const {
  // Potential function
  auto F = [&](length_quantity closest_approach) {
    return 1 - screening_potential(closest_approach, z_1, z_2) / cm_energy -
           pow<2>(impact_param / closest_approach);
  };
  // Derivative of potential function
  auto dF = [&](length_quantity closest_approach) {
    return -screening_potential_derivative(closest_approach, z_1, z_2) /
               cm_energy +
           2.0 * impact_param * impact_param /
               pow<static_rational<3>>(closest_approach);
  };
  auto old_result = impact_param; // Initial guess, adjust if encountering
                                  // convergence issues
  VLOG(3) << "Closest approach initial guess " << old_result;
  auto old_pot = F(old_result);
  auto result = old_result - old_pot / dF(old_result); // First iteration
  auto pot = F(result);
  // Newtonian iteration
  size_t iterations = 0;
  constexpr size_t max_iterations = 100; // If encountering convergence errors,
                                         // increase. Will decrease performance.
  constexpr double rel_tolerance =
      1e-3; // Threshold for relative tolerance. Increase if encountering
            // convergence issues.
  while (iterations < max_iterations) {
    auto rel_diff = abs(result - old_result) / result;

    if (rel_diff < rel_tolerance) {
      VLOG(2) << "Closest approach found after " << iterations
              << " iterations. Result " << result;
      return result;
    }

    auto d_pot = dF(old_result);
    old_pot = pot;
    pot = F(result);

    // Based on verbose log inspection, this is unlikely to occur.
    if ((old_pot < 0.0 && pot > 0.0) || (old_pot > 0.0 && pot < 0.0)) {
      // Go into binary search (Solves oscillation problem)
      size_t binary_search_iterations = 0;

      while (true) {
        if (iterations > max_iterations) {
          break;
        }
        if (rel_diff < rel_tolerance) {
          VLOG(2) << "Closest approach found after " << iterations << " and "
                  << binary_search_iterations << " binary search iterations."
                  << " iterations. Result " << result;
          rel_diff = abs(result - old_result) / result;
          return result;
        }
        auto mid_result = old_result + (result - old_result) / 2.0;
        auto mid_pot = F(mid_result);

        VLOG(3) << "Binary search old_result " << old_result << " old_pot "
                << old_pot << " mid_result " << mid_result << " mid_pot "
                << mid_pot << " result " << result << " pot " << pot;
        if (copysign(old_pot, mid_pot) == old_pot) {
          // If old_pot and mid_pot have the same sign
          old_result = mid_result; // Search the "right" side
          old_pot = mid_pot;
          VLOG(3) << "Searching pot side";
        } else {
          result = mid_result; // Search the "left" side
          pot = mid_pot;
          VLOG(3) << "Searching old_pot side";
        }
        binary_search_iterations++;
      }
    }

    old_result = result;
    result = old_result - pot / d_pot;
    iterations++;

    VLOG(3) << "Zeroing function value " << pot << ", derivative " << d_pot;
    VLOG(3) << "Newtonian iteration " << iterations << ". Previous result "
            << old_result << ", current result " << result << ". Tolerance "
            << rel_tolerance << ", relative difference " << rel_diff;
  }

  // TODO: Write more graceful code to handle this. And throw an exception to
  // be handled by the GUI if convergence doesn't occur.
  LOG(ERROR) << "Convergence error. Please adjust closest approach parameters. "
                "Previous result "
             << old_result << ", current result " << result;

  return result;
}

length_quantity Ion::radius_of_curvature(length_quantity radius,
                                         energy_quantity cm_energy,
                                         dimensionless_quantity z_1,
                                         dimensionless_quantity z_2) const {
  auto result = 2.0 * abs(cm_energy - screening_potential(radius, z_1, z_2)) /
                (-screening_potential_derivative(radius, z_1, z_2));
  VLOG(2) << "Radius of curvature for radius " << radius << " cm_energy "
          << cm_energy << " z_1 " << z_1 << " z_2 " << z_2 << " has result "
          << result;
  return result;
}

dimensionless_quantity Ion::magic_formula_correction_parameter(
    dimensionless_quantity reduced_energy, length_quantity impact_parameter,
    length_quantity screening_length, length_quantity closest_approach) const {
  auto C_1 = 0.99229;
  auto C_2 = 0.011615;
  auto C_3 = 0.007122;
  auto C_4 = 9.3066;
  auto C_5 = 14.813;

  auto alpha = 1 + C_1 / sqrt(reduced_energy);
  auto beta = (C_2 + sqrt(reduced_energy)) / (C_3 + sqrt(reduced_energy));
  auto gamma = (C_4 + reduced_energy) / (C_5 + reduced_energy);
  auto B = impact_parameter / screening_length;
  auto A = 2 * alpha * reduced_energy * pow(B, beta);
  auto G = gamma / (sqrt(1 + A * A) - A);
  auto R_o = closest_approach / screening_length;
  auto result = A * (R_o - B) / (1 + G);
  VLOG(2) << "Magic formula correction parameter with reduced_energy "
          << reduced_energy << ", impact_parameter " << impact_parameter
          << ", screening_length " << screening_length << " result is "
          << result;
  return result;
}

angle_quantity
Ion::magic_formula_scattering_angle(length_quantity impact_parameter,
                                    length_quantity closest_approach,
                                    length_quantity radius_of_curvature,
                                    dimensionless_quantity correction_factor,
                                    length_quantity screening_distance) const {
  // TODO: Correct radius of curvature naming scheme
  const auto B = impact_parameter / screening_distance;
  const auto R_o = closest_approach / screening_distance;
  const auto R_c = radius_of_curvature / screening_distance;
  const auto result = acos((B + R_c + correction_factor) / (R_o + R_c));
  VLOG(2) << "Magic scattering angle with impact_parameter " << impact_parameter
          << " closest_approach " << closest_approach << " radius_of_curvature "
          << radius_of_curvature << " correction_factor " << correction_factor
          << " screening_distance " << screening_distance << " with result "
          << result;
  return result;
}

energy_quantity
Ion::nuclear_energy_loss(mass_quantity atom_mass, mass_quantity target_mass,
                         energy_quantity incident_energy,
                         angle_quantity scattering_angle) const {
  auto result = 4.0 * atom_mass * target_mass /
                pow<2>(atom_mass * target_mass) * incident_energy *
                pow<2>(sin(scattering_angle));
  VLOG(2) << "Nuclear scattering energy loss atom mass " << atom_mass
          << " target_mass " << target_mass << " incident_energy "
          << incident_energy << " scattering_angle " << scattering_angle
          << " has result " << result;
  return electron_volt;
}

angle_quantity
Ion::laboratory_scattering_angle(angle_quantity cm_scattering_angle,
                                 mass_quantity atom_mass,
                                 mass_quantity target_mass) const {
  auto result = atan(sin(cm_scattering_angle) /
                     (cos(cm_scattering_angle) + (atom_mass / target_mass)));

  VLOG(2) << "Laboratory scattering angle with cm_scattering_angle "
          << cm_scattering_angle << " atom_mass " << atom_mass
          << " target_mass " << target_mass << " with result " << result;
  return result;
}

// TODO: Switched to quaternions, likely to remove.{{{
// std::array<angle_quantity, 4>
// Ion::relative_to_absolute_angle(angle_quantity initial_altitude,
//                                 angle_quantity initial_azimuth,
//                                 angle_quantity incident_deflection,
//                                 angle_quantity target_deflection) const {
//   auto THETA1RELATIVE = incident_deflection;
//   auto THETA2RELATIVE = target_deflection;
//   auto ALPHA1RELATIVE =
//       uniform_distribution_(random_number_generator_) * 2.0 * pi;
//   auto ALPHA2RELATIVE = ALPHA1RELATIVE + pi;

//   auto THETAO = initial_altitude;
//   auto ALPHAO = initial_azimuth;

//   // Angle
//   auto X1 = sin(THETA1RELATIVE) * cos(ALPHA1RELATIVE);
//   auto Y1 = sin(THETA1RELATIVE) * sin(ALPHA1RELATIVE);
//   auto Z1 = cos(THETA1RELATIVE);

//   auto Y0 = Y1 * cos(THETAO) + Z1 * sin(THETAO);
//   auto Z0 = -Y1 * sin(THETAO) + Z1 * cos(THETAO);
//   auto X0 = X1;

//   auto Z = Z0;
//   auto X = X0 * sin(ALPHAO) + Y0 * cos(ALPHAO);
//   auto Y = -X0 * cos(ALPHAO) + Y0 * sin(ALPHAO);

//   angle_quantity THETA1;

//   if (Z > 0.0) {
//     THETA1 = atan(sqrt(X * X + Y * Y) / Z);
//   } else if (Z < 0.0) {
//     THETA1 = radian * pi + atan(sqrt(X * X + Y * Y) / Z);
//   } else {
//     THETA1 = radian * pi / 2.0;
//   }

//   angle_quantity ALPHA1;

//   if (sin(THETA1) != 0.0) {
//     if (X > 0.0) {
//       ALPHA1 = atan(Y / X);
//     } else if (X == 0.0) {
//       ALPHA1 = radian * (pi - (Y > 0.0 ? 1.0 : -1.0) * pi / 2.0);
//     } else {
//       ALPHA1 = radian * pi + atan(Y / X);
//     }
//   } else {
//     ALPHA1 = 0.0 * radian;
//   }

//   // Target
//   auto X3 = sin(THETA2RELATIVE) * cos(ALPHA2RELATIVE);
//   auto Y3 = sin(THETA2RELATIVE) * sin(ALPHA2RELATIVE);
//   auto Z3 = cos(THETA2RELATIVE);

//   auto Y2 = Y3 * cos(THETAO) + Z3 * sin(THETAO);
//   auto Z2 = -Y3 * sin(THETAO) + Z3 * cos(THETAO);
//   auto X2 = X3;

//   auto Z5 = Z2;
//   auto X5 = X2 * sin(ALPHAO) + Y2 * cos(ALPHAO);
//   auto Y5 = -X2 * cos(ALPHAO) + Y2 * sin(ALPHAO);

//   angle_quantity THETA2;
//   if (Z5 < 0.0) {
//     THETA2 = radian * pi + atan(sqrt(X5 * X5 + Y5 * Y5) / Z5);
//   } else if (Z5 > 0.0) {
//     THETA2 = atan(sqrt(X5 * X5 + Y5 * Y5) / Z5);
//   } else {
//     THETA2 = radian * pi / 2.0;
//   }

//   angle_quantity ALPHA2;

//   if (sin(THETA2) != 0.0) {
//     if (X5 < 0.0) {
//       ALPHA2 = radian * pi + atan(Y5 / X5);
//     } else if (X5 > 0.0) {
//       atan(Y5 / X5);
//     } else {
//       ALPHA2 = radian * (pi - 0.5 * ((Y5 > 0.0) ? 1 : -1) * pi);
//     }
//   } else {
//     ALPHA2 = radian * 0.0;
//   }
//   // newVelocity.zAngle = THETA1;
//   // newVelocity.xAngle = ALPHA1;
//   // targetVelocity.zAngle = THETA2;
//   // targetVelocity.xAngle = ALPHA2;
//   return {{THETA1, ALPHA1, THETA2, ALPHA2}};
// } // }}}

// Electron function implementations{{{
void Electron::fire() try {
  // TODO: Implement Electron simulation
  throw std::logic_error("Electron simulation is not implemented");
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Ion function implementations // {{{
void Ion::fire() try {
  std::ostringstream string_stream;
  LOG(INFO) << "Firing ion with initial energy " << velocity_.energy_;
  const auto ion_stopping_energy =
      simulation_.get_settings(&Simulation::Settings::ion_stopping_energy_);
  auto current_layer_pair = simulation_.get_volume().get_layer(coordinate_.z_);
  const Layer &current_layer = *(current_layer_pair.first);

#ifdef NDEBUG /*{{{*/
#else
  // Grad an iterator to a view of the first
  auto comp_view = current_layer.get_relative_compositions().begin();
  for (size_t i = 0; i < current_layer.get_number_of_components(); i++) {
    auto direct_composition = current_layer.get_relative_composition(i);
    auto viewed_composition = *(comp_view++);
    VLOG(2) << "Index: " << i << " has directly measured composition of "
            << direct_composition << " view-based composition of "
            << "Relative (view): " << viewed_composition;
    assert(direct_composition == viewed_composition);
  }
#endif /*}}}*/

  for (int i = 0; i < 100; i++) {
    coordinate_ = Coordinate(angstrom * 0.0, angstrom * 0.0, angstrom * 0.0);
    velocity_.energy_ = 100.0 * kilo_electron_volt;
    velocity_.quaternion_ = Quaternion(0.0, 0.0, 0.0, 1.0);

    auto &energy = velocity_.energy_;
    // auto &coordinate = coordinate_;
    while (velocity_.energy_ > ion_stopping_energy) {
      // TODO:
      const size_t atom_index =
          current_layer.get_discrete_distribution()(random_number_generator_);
      const auto layer_properties = current_layer.get_property(atom_index);

      auto z_1 = properties_.charge_;
      auto z_2 = layer_properties.charge_;
      auto m_1 = properties_.mass_;
      auto m_2 = layer_properties.mass_;
      auto number_density = current_layer.get_number_density();

      // TODO: Implement simulation calculations
      //
      // 1. Subtract electronic stopping energy
      velocity_.energy_ = electronic_stopping_energy(
          z_1, m_1, energy, number_density, layer_properties);

      VLOG(2) << "New energy " << velocity_.energy_;
      auto speed = speed_from_energy(energy, m_1);
      // auto rd_mass = reduced_mass(m_1, m_2);
      auto c_mass_energy = cm_energy(m_1, speed);
      auto screen_length = screening_length(z_1, z_2);
      auto rd_energy = reduced_energy(screen_length, c_mass_energy, z_1, z_2);
      auto path_length = free_flying_path_length(m_1, m_2, rd_energy,
                                                 screen_length, number_density);
      auto i_param = impact_parameter(number_density, path_length, rd_energy);

      // auto coll_diameter = collision_diameter(z_1, z_2, rd_mass, speed);
      auto closest_appr = closest_approach(z_1, z_2, c_mass_energy, i_param);
      auto rad_curv =
          radius_of_curvature(closest_appr, c_mass_energy, z_1, z_2);
      auto magic_correction = magic_formula_correction_parameter(
          rd_energy, i_param, screen_length, closest_appr);
      auto cm_scatter_angle = magic_formula_scattering_angle(
          i_param, closest_appr, rad_curv, magic_correction, screen_length);
      auto nuclear_e_loss =
          nuclear_energy_loss(m_1, m_2, energy, cm_scatter_angle);
      energy -= nuclear_e_loss; // TODO: Check for knock-ons
      auto laboratory_scatter_angle =
          laboratory_scattering_angle(cm_scatter_angle, m_1, m_2);
      rotate_inplace(laboratory_scatter_angle); // Changes velocity in-place
      auto &velocity = velocity_.quaternion_;

      VLOG(2) << "Original coordinate (" << coordinate_.x_ << ", "
              << coordinate_.y_ << ", " << coordinate_.z_ << ")";

      VLOG(2) << "Original velocity (" << velocity.x_ << ", " << velocity.y_
              << ", " << velocity.z_ << ")";
      // TODO: Replace with more elegant code
      auto new_coordinate =
          Coordinate{coordinate_.x_ + path_length * velocity.x_,
                     coordinate_.y_ + path_length * velocity.y_,
                     coordinate_.z_ + path_length * velocity.z_};

      coordinate_.operator=(new_coordinate);
      string_stream << velocity.r_ << "," << coordinate_.x_.value() << ","
                    << coordinate_.y_.value() << "," << coordinate_.z_.value()
                    << "\n";

      VLOG(2) << "Final coordinate moved " << path_length << " away: ("
              << coordinate_.x_ << ", " << coordinate_.y_ << ", "
              << coordinate_.z_ << ")";

      // TODO: Implement derivative of screening potential

      // 2. Calculate recoil energy and velocity
      // 3. Calculate new coordinate (push old coordinate to vector)
      // 4. Determine if sputtering occured
      // 5. Check if damage cascade occurs
      //    a. If so, if energy difference is larger than displacement
      //    energy, create a cascade particle and push it to the particles_
      //    vector
    }
  }

  // TODO: FIX
  std::ofstream file_stream("output.csv");
  file_stream << "scalar,x,y,z\n";
  file_stream << string_stream.str();

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

length_quantity
Ion::impact_parameter(number_density_quantity number_density,
                      length_quantity length,
                      dimensionless_quantity reduced_energy) const {
  length_quantity result;

  auto random_number = uniform_distribution_(random_number_generator_);
  // TODO: Fix this. Fix proper random numbers and calculation. High energy
  // gives ridiculous results. Research e >> 10 on 7-8 in SRIM book.
  if (reduced_energy > 10.0) {
    result = sqrt(-log(random_number) /
                  (pi * number_density * length * si::constants::codata::N_A));
    VLOG(2) << "Impact parameter at high energy for number density "
            << number_density << ", length " << length
            << ", and reduced energy " << reduced_energy << " is " << result;
  } else {
    result = sqrt(random_number /
                  (pi * pow<static_rational<2, 3>>(
                            number_density * si::constants::codata::N_A)));
    VLOG(2) << "Impact parameter at low energy for number density "
            << number_density << ", length " << length
            << ", and reduced energy " << reduced_energy << " is " << result;
  }
  return result;
}

length_quantity
Ion::interatomic_spacing(number_density_quantity number_density) const {
  return root<-3>(number_density * constants::N_A);
};

energy_quantity
Ion::electronic_stopping_energy(dimensionless_quantity charge,
                                mass_quantity mass, energy_quantity energy,
                                number_density_quantity number_density,
                                const Properties &properties) const try {
  VLOG(2) << "Calculating ion electronic stopping energy with charge " << charge
          << ", mass " << mass << ", and energy " << energy;

  const auto c_1 =
      1.212 * root<2>(constants::m_u) * root<2>(electron_volt) / angstrom;
  VLOG(3) << "Value of c_1 " << c_1;
  const auto k_l = c_1 * pow<static_rational<7, 6>>(charge) *
                   properties.charge_ /
                   (pow<static_rational<3, 2>>(
                        pow<static_rational<2, 3>>(charge) +
                        pow<static_rational<2, 3>>(properties.charge_)) *
                    root<2>(mass));
  VLOG(3) << "Value of k_l " << k_l;
  const auto stopping_energy = k_l * root<2>(energy);
  VLOG(3) << "Value of stopping_energy " << stopping_energy;
  const auto c_2 = 1.59 * pow<3>(angstrom);
  VLOG(3) << "Value of c_2 " << c_2;
  const auto travel_length = interatomic_spacing(number_density);
  VLOG(3) << "Value of travel_length " << travel_length;
  const auto energy_loss =
      c_2 * travel_length * number_density * constants::N_A * stopping_energy;
  VLOG(3) << "Value of energy_loss " << energy_loss;
  return energy - energy_loss;
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Quaternion function implementations{{{

Quaternion &Quaternion::operator+=(const Quaternion &b) {
  r_ += b.r_;
  x_ += b.x_;
  y_ += b.y_;
  z_ += b.z_;
  return *this;
}

Quaternion Quaternion::operator+(const Quaternion &b) const {
  return Quaternion(*this) += b;
}

Quaternion operator+(const Quaternion &a, const Quaternion &b) {
  return Quaternion(a) += b;
}

Quaternion &Quaternion::operator*=(const Quaternion &b) {
  auto r_1 = r_;
  auto r_2 = b.r_;
  auto x_1 = x_;
  auto y_1 = y_;
  auto z_1 = z_;
  auto x_2 = b.x_;
  auto y_2 = b.y_;
  auto z_2 = b.z_;

  r_ = r_1 * r_2 - x_1 * x_2 - y_1 * y_2 - z_1 * z_2;
  x_ = r_1 * x_2 + r_2 * x_1 + y_1 * z_2 - z_1 * y_2;
  y_ = r_1 * y_2 + r_2 * y_1 + z_1 * x_2 - x_1 * z_2;
  z_ = r_1 * z_2 + r_2 * z_1 + x_1 * y_2 - y_1 * x_2;
  return *this;
}

Quaternion Quaternion::operator*(const Quaternion &b) const {
  return Quaternion(*this) *= b;
}

// Quaternion operator*(const Quaternion &a, const Quaternion &b) {
//   return Quaternion(a) *= b;
// }

Quaternion &Quaternion::operator*=(double a) {
  r_ *= a;
  x_ *= a;
  y_ *= a;
  z_ *= a;
  return *this;
}

Quaternion Quaternion::operator*(const double b) const {
  return Quaternion(*this) *= b;
}

Quaternion operator*(double a, const Quaternion &b) {
  return Quaternion(b) *= a;
}

Quaternion &Quaternion::operator/=(double a) {
  r_ /= a;
  x_ /= a;
  y_ /= a;
  z_ /= a;
  return *this;
}

Quaternion Quaternion::operator/(const double b) const {
  return Quaternion(*this) /= b;
}

Quaternion operator/(const Quaternion &a, double b) {
  return Quaternion(a) /= b;
}

Quaternion &Quaternion::normalize_inplace() {
  *this /= (this->magnitude());
  return *this;
}

Quaternion &Quaternion::conjugate_inplace() {
  x_ = -x_;
  y_ = -y_;
  z_ = -z_;
  return *this;
}

Quaternion Quaternion::conjugate_copy() const {
  return Quaternion(r_, -x_, -y_, -z_);
}

double Quaternion::magnitude() const {
  return std::sqrt(r_ * r_ + x_ * x_ + y_ * y_ + z_ * z_);
} // }}}

} // namespace common
} // namespace amcset
