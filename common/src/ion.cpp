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

/*! \file ion.cpp
 *  \brief Ion transport (binary-collision approximation) and the physics
 *  model functions selected through Simulation::Settings: screening
 *  functions/lengths, nuclear energy transfer and electronic stopping.
 */

#include "amcset_common.h"

#include <absl/log/log.h>

#include <boost/units/cmath.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>

#include <algorithm>
#include <cmath>

namespace amcset {
namespace common {

Ion::Ion(size_t z_number, size_t mass_number, const Simulation &simulation,
         boost::random::uniform_01<double> &uniform_distribution,
         boost::random::mt19937 &random_number_generator,
         ParticleDataLeaf &data_leaf, Bombardment &bombardment,
         size_t generation, double birth_time_seconds)
    : Particle(z_number, mass_number, simulation, uniform_distribution,
               random_number_generator, data_leaf, bombardment, generation,
               birth_time_seconds),
      custom_electronic_coefficient_(
          generation == 0
              ? simulation.get_settings(
                    &Simulation::Settings::custom_electronic_coefficient_)
              : simulation.get_settings(
                    &Simulation::Settings::
                        custom_recoil_electronic_coefficient_)),
      custom_electronic_exponent_(
          generation == 0
              ? simulation.get_settings(
                    &Simulation::Settings::custom_electronic_exponent_)
              : simulation.get_settings(
                    &Simulation::Settings::custom_recoil_electronic_exponent_)) {}

length_quantity Ion::screening_length(dimensionless_quantity z_1,
                    dimensionless_quantity z_2) const {
  return simulation_.nuclear_screening_length_function_(
      z_1, z_2,
      simulation_.get_settings(
          &Simulation::Settings::custom_nuclear_screening_scale_));
}

length_quantity Ion::screening_length_for_model(dimensionless_quantity z_1,
                                                dimensionless_quantity z_2,
                                                bool use_quickbasic) {
  if (use_quickbasic) {
    // Firsov screening length 0.8854 a_0 / (Z1^(1/2) + Z2^(1/2))^(2/3).
    auto numerator = 0.8854 * bohr_radius;
    auto denom = pow(pow(z_1, 0.5) + pow(z_2, 0.5), 2.0 / 3.0);
    return numerator / denom;
  }
  return 0.25 * cbrt(9.0 * pi * pi / 2.0) * bohr_radius /
         (pow(z_1, 0.23) + pow(z_2, 0.23));
}

// The "Thomas-Fermi" model is the four-exponential universal-type function
// with the 0.4028 exponent of the original code. Its last coefficient was
// 0.2817 (a transcription error of ZBL's 0.02817), which made phi(0) = 1.25
// instead of 1; fixed in v0.1.0 (docs/KNOWN_ISSUES.md, physics fix P4).
dimensionless_quantity Simulation::thomas_fermi_screening(
    dimensionless_quantity x) {
  return 0.1818 * exp(-3.2 * x) + 0.5099 * exp(-0.9423 * x) +
         0.2802 * exp(-0.4028 * x) + 0.02817 * exp(-0.2016 * x);
}

dimensionless_quantity Simulation::zbl_screening(dimensionless_quantity x) {
  return 0.1818 * exp(-3.2 * x) + 0.5099 * exp(-0.9423 * x) +
         0.2802 * exp(-0.4029 * x) + 0.02817 * exp(-0.2016 * x);
}

dimensionless_quantity Simulation::moliere_screening(dimensionless_quantity x) {
  return 0.35 * exp(-0.3 * x) + 0.55 * exp(-1.2 * x) + 0.10 * exp(-6.0 * x);
}

dimensionless_quantity Simulation::thomas_fermi_screening_derivative(
    dimensionless_quantity x) {
  return -3.2 * 0.1818 * exp(-3.2 * x) -
         0.9423 * 0.5099 * exp(-0.9423 * x) -
         0.4028 * 0.2802 * exp(-0.4028 * x) -
         0.2016 * 0.02817 * exp(-0.2016 * x);
}

dimensionless_quantity Simulation::zbl_screening_derivative(
    dimensionless_quantity x) {
  return -3.2 * 0.1818 * exp(-3.2 * x) -
         0.9423 * 0.5099 * exp(-0.9423 * x) -
         0.4029 * 0.2802 * exp(-0.4029 * x) -
         0.2016 * 0.02817 * exp(-0.2016 * x);
}

dimensionless_quantity Simulation::moliere_screening_derivative(
    dimensionless_quantity x) {
  return -0.3 * 0.35 * exp(-0.3 * x) -
         1.2 * 0.55 * exp(-1.2 * x) -
         6.0 * 0.10 * exp(-6.0 * x);
}

length_quantity Simulation::thomas_fermi_screening_length(
    dimensionless_quantity z_1, dimensionless_quantity z_2, double) {
  return 0.25 * cbrt(9.0 * pi * pi / 2.0) * bohr_radius /
         (pow(z_1, 0.23) + pow(z_2, 0.23));
}

length_quantity Simulation::firsov_screening_length(
    dimensionless_quantity z_1, dimensionless_quantity z_2, double) {
  return 0.8854 * bohr_radius /
         pow(pow(z_1, 0.5) + pow(z_2, 0.5), 2.0 / 3.0);
}

length_quantity Simulation::zbl_screening_length(
    dimensionless_quantity z_1, dimensionless_quantity z_2, double) {
  return 0.8854 * bohr_radius /
         (pow(z_1, 0.23) + pow(z_2, 0.23));
}

length_quantity Simulation::moliere_screening_length(
    dimensionless_quantity z_1, dimensionless_quantity z_2, double) {
  return 0.8853 * bohr_radius /
         sqrt(pow(z_1, 2.0 / 3.0) + pow(z_2, 2.0 / 3.0));
}

length_quantity Simulation::custom_screening_length(
    dimensionless_quantity z_1, dimensionless_quantity z_2, double scale) {
  return scale * thomas_fermi_screening_length(z_1, z_2, 1.0);
}

energy_quantity Ion::screening_potential(length_quantity radius,
                                         dimensionless_quantity z_1,
                                         dimensionless_quantity z_2) const {
  const auto length = screening_length(z_1, z_2);
  const auto reduced_length = radius / length;
  auto result = simulation_.nuclear_screening_function_(reduced_length) *
                z_1 * z_2 * coulomb_energy_length / radius;
  VLOG(2) << "Screening potential with radius " << radius << ", z_1 " << z_1
          << ", and z_2 " << z_2 << " is " << result;
  return result;
}

force_quantity Ion::screening_potential_derivative(length_quantity radius,
                                    dimensionless_quantity z_1,
                                    dimensionless_quantity z_2) const {
  auto s_length = screening_length(z_1, z_2);
  auto reduced_radius = radius / s_length;
  auto result =
      (simulation_.nuclear_screening_derivative_function_(reduced_radius) /
           s_length -
       simulation_.nuclear_screening_function_(reduced_radius) / radius) *
      z_1 * z_2 * coulomb_energy_length / radius;
  VLOG(2) << "Screening potential derivative with radius " << radius << ", z_1 " << z_1
          << ", and z_2 " << z_2 << " is " << result;
  return result;
}

dimensionless_quantity Ion::reduced_energy(length_quantity screening_length,
                                           energy_quantity cm_energy,
                                           dimensionless_quantity z_1,
                                           dimensionless_quantity z_2) const {
  auto result = screening_length * cm_energy /
                (z_1 * z_2 * coulomb_energy_length);
  VLOG(2) << "Reduced energy with screening length " << screening_length
          << ", center of mass energy " << cm_energy << ", z_1 " << z_1
          << ", z_2 " << z_2 << ", result is " << result;
  return result;
}

length_quantity Ion::free_flying_path_length(mass_quantity m_1, mass_quantity m_2,
                             dimensionless_quantity reduced_energy,
                             length_quantity screening_length,
                             number_density_quantity number_density) const {
  // One interatomic spacing, N^(-1/3).
  const length_quantity L = 1.0 / root<3>(number_density);
  VLOG(3) << "Low energy free flying path length with m_1 " << m_1 << ", m_2 "
          << m_2 << ", reduced_energy " << reduced_energy
          << ", screening_length " << screening_length
          << ", and number_density " << number_density << " with result " << L;
  return L;
}

length_quantity Ion::collision_diameter(dimensionless_quantity z_1,
                                        dimensionless_quantity z_2,
                                        mass_quantity rd_mass,
                                        velocity_quantity speed) const {
  auto result =
      4.0 * z_1 * z_2 * coulomb_energy_length / (rd_mass * speed * speed);
  VLOG(2) << "Collision diameter given z_1 " << z_1 << ", z_2 " << z_2
          << " reduced mass " << rd_mass << ", and speed " << speed << " is "
          << result;
  return result;
}

length_quantity Ion::closest_approach(dimensionless_quantity z_1,
                                      dimensionless_quantity z_2,
                                      energy_quantity cm_energy,
                                      length_quantity impact_param) const {
  auto F = [&](length_quantity closest_approach) {
    return 1 - screening_potential(closest_approach, z_1, z_2) / cm_energy -
           pow<2>(impact_param / closest_approach);
  };
  auto dF = [&](length_quantity closest_approach) {
    return -screening_potential_derivative(closest_approach, z_1, z_2) / cm_energy +
           2.0 * impact_param * impact_param /
               pow<static_rational<3>>(closest_approach);
  };
  // Newton iteration starting from r0 = b. If F changes sign between two
  // consecutive iterates, switch to bisection on that bracket.
  auto old_result = impact_param;
  VLOG(3) << "Closest approach initial guess " << old_result;
  auto old_pot = F(old_result);
  auto result = old_result - old_pot / dF(old_result);
  auto pot = F(result);
  size_t iterations = 0;
  constexpr size_t max_iterations = 100;
  constexpr double rel_tolerance = 1e-3;
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

    if ((old_pot < 0.0 && pot > 0.0) || (old_pot > 0.0 && pot < 0.0)) {
      // Bisection on [old_result, result]. Previously neither rel_diff nor
      // the iteration counter was updated inside this loop, so entering it
      // never terminated. The bracket width is now re-evaluated each step
      // and the loop is capped at max_iterations.
      size_t binary_search_iterations = 0;
      while (binary_search_iterations < max_iterations) {
        if (rel_diff < rel_tolerance) {
          VLOG(2) << "Closest approach found after " << iterations << " and "
                  << binary_search_iterations
                  << " binary search iterations. Result " << result;
          return result;
        }
        auto mid_result = old_result + (result - old_result) / 2.0;
        auto mid_pot = F(mid_result);
        if (copysign(old_pot, mid_pot) == old_pot) {
          old_result = mid_result;
          old_pot = mid_pot;
        } else {
          result = mid_result;
          pot = mid_pot;
        }
        binary_search_iterations++;
        rel_diff = abs(result - old_result) / result;
      }
      LOG(ERROR) << "Closest approach bisection did not converge. Result "
                 << result;
      return result;
    }

    old_result = result;
    result = old_result - pot / d_pot;
    iterations++;
  }

  LOG(ERROR) << "Convergence error. Please adjust closest approach parameters. Previous result "
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
  // G = 1 / (gamma ((1 + A^2)^(1/2) - A)) as in TRIM-85; the SRIM book's
  // Eq. 7-11 prints gamma in the numerator, which underestimates
  // sin^2(theta/2) by up to ~40 % (physics fix P6, docs/KNOWN_ISSUES.md).
  auto G = 1.0 / (gamma * (sqrt(1 + A * A) - A));
  auto R_o = closest_approach / screening_length;
  auto result = A * (R_o - B) / (1 + G);
  VLOG(2) << "Magic formula correction parameter with reduced_energy "
          << reduced_energy << ", impact_parameter " << impact_parameter
          << ", screening_length " << screening_length << " result is "
          << result;
  return result;
}

angle_quantity Ion::magic_formula_scattering_angle(length_quantity impact_parameter,
                                    length_quantity closest_approach,
                                    length_quantity radius_of_curvature,
                                    dimensionless_quantity correction_factor,
                                    length_quantity screening_distance) const {
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

energy_quantity Ion::nuclear_energy_loss(mass_quantity atom_mass, mass_quantity target_mass,
                         energy_quantity incident_energy,
                         angle_quantity scattering_angle) const {
  auto result = simulation_.nuclear_stopping_function_(
      atom_mass, target_mass, incident_energy, scattering_angle);
  VLOG(2) << "Nuclear scattering energy loss atom mass " << atom_mass
          << " target_mass " << target_mass << " incident_energy "
          << incident_energy << " scattering_angle " << scattering_angle
          << " has result " << result;
  return result;
}

energy_quantity Simulation::binary_collision_nuclear_stopping(
    mass_quantity atom_mass, mass_quantity target_mass,
    energy_quantity incident_energy, angle_quantity scattering_angle) {
  return 4.0 * atom_mass * target_mass /
         pow<2>(atom_mass + target_mass) * incident_energy *
         pow<2>(sin(scattering_angle));
}

angle_quantity Ion::laboratory_scattering_angle(angle_quantity cm_scattering_angle,
                                 mass_quantity atom_mass,
                                 mass_quantity target_mass) const {
  auto result = atan(sin(cm_scattering_angle) /
                     (cos(cm_scattering_angle) + (atom_mass / target_mass)));

  VLOG(2) << "Laboratory scattering angle with cm_scattering_angle "
          << cm_scattering_angle << " atom_mass " << atom_mass
          << " target_mass " << target_mass << " with result " << result;
  return result;
}

void Ion::fire() try {
  VLOG(1) << "Firing ion with initial energy " << velocity_.energy_;
  auto ion_stopping_energy = simulation_.get_settings(
      &Simulation::Settings::ion_stopping_energy_);
  VLOG(2) << "Ion stopping energy: " << ion_stopping_energy;

  VLOG(2) << "Particle generation: " << generation_;

  auto &energy = velocity_.energy_;
  VLOG(1) << "Initial energy: " << energy;
  auto max_depth = simulation_.get_volume().get_depth();
  VLOG(2) << "Max depth: " << max_depth;
  auto knock_on_threshold = simulation_.get_settings(&Simulation::Settings::ion_displacement_energy_);
  VLOG(2) << "Knock-on threshold: " << knock_on_threshold;
  auto enable_knock_on = simulation_.get_settings(&Simulation::Settings::enable_damage_cascade_);
  VLOG(2) << "Knock-on enabled: " << enable_knock_on;

  while (velocity_.energy_ > ion_stopping_energy) {
    const auto current_layer_pair = simulation_.get_volume().get_layer(coordinate_.z_);
    VLOG(3) << "Current layer: " << current_layer_pair.second;
    const Layer &current_layer = *(current_layer_pair.first);
    const size_t atom_index = current_layer.get_random_index(random_number_generator_);
    VLOG(3) << "Selected atom index: " << atom_index;
    const auto layer_properties = current_layer.get_property(atom_index);
    auto z_1 = properties_.charge_;
    VLOG(3) << "Ion charge: " << z_1;
    auto z_2 = layer_properties.charge_;
    VLOG(3) << "Target atom charge: " << z_2;
    auto m_1 = properties_.mass_;
    VLOG(3) << "Ion mass: " << m_1;
    auto m_2 = layer_properties.mass_;
    VLOG(3) << "Target atom mass: " << m_2;
    auto number_density = current_layer.get_number_density();
    VLOG(3) << "Layer number density: " << number_density;

    auto e_energy_loss = std::min(
        electronic_stopping_energy(z_1, m_1, energy, number_density,
                                   layer_properties),
        energy);
    VLOG(2) << "Electronic stopping energy loss: " << e_energy_loss;
    energy -= e_energy_loss;
    VLOG(2) << "Energy after electronic stopping: " << energy;
    bombardment_.record_energy_loss((long double)e_energy_loss.value(), coordinate_.z_.value(), Bombardment::EnergyLossType::electronic);
    if (energy <= ion_stopping_energy) {
      bombardment_.record_energy_loss(
          static_cast<long double>(e_energy_loss.value()),
          coordinate_.z_.value(), Bombardment::EnergyLossType::total);
      break;
    }

    auto speed = speed_from_energy(energy, m_1);
    VLOG(3) << "Speed from energy: " << speed;
    // Centre-of-mass energy (1/2) mu v^2 with the reduced mass mu.
    auto c_mass_energy = cm_energy(reduced_mass(m_1, m_2), speed);
    VLOG(3) << "Center of mass energy: " << c_mass_energy;
    auto screen_length = screening_length(z_1, z_2);
    VLOG(3) << "Screening length: " << screen_length;
    auto rd_energy = reduced_energy(screen_length, c_mass_energy, z_1, z_2);
    VLOG(3) << "Reduced energy: " << rd_energy;
    auto path_length = free_flying_path_length(m_1, m_2, rd_energy, screen_length, number_density);
    VLOG(3) << "Free flying path length: " << path_length;
    auto i_param = impact_parameter(number_density, path_length, rd_energy);
    VLOG(3) << "Impact parameter: " << i_param;
    auto closest_appr = closest_approach(z_1, z_2, c_mass_energy, i_param);
    VLOG(3) << "Closest approach: " << closest_appr;
    auto rad_curv = radius_of_curvature(closest_appr, c_mass_energy, z_1, z_2);
    VLOG(3) << "Radius of curvature: " << rad_curv;
    auto magic_correction = magic_formula_correction_parameter(rd_energy, i_param, screen_length, closest_appr);
    VLOG(3) << "Magic formula correction parameter: " << magic_correction;
    auto cm_scatter_angle = magic_formula_scattering_angle(i_param, closest_appr, rad_curv, magic_correction, screen_length);
    VLOG(3) << "Center of mass scattering angle: " << cm_scatter_angle;
    auto nuclear_e_loss = nuclear_energy_loss(m_1, m_2, energy, cm_scatter_angle);
    VLOG(2) << "Nuclear energy loss: " << nuclear_e_loss;

    energy -= nuclear_e_loss;
    bombardment_.record_energy_loss((long double)nuclear_e_loss.value(), coordinate_.z_.value(), Bombardment::EnergyLossType::nuclear);
    {
      // Record the collision with the post-collision speed (pre-rotation
      // direction). The time step is the upcoming flight length divided by
      // that speed.
      const double sp =
          speed_from_energy(velocity_.energy_, properties_.mass_).value();
      const auto &q = velocity_.quaternion_;
      const double vx = q.x_ * sp;
      const double vy = q.y_ * sp;
      const double vz = q.z_ * sp;
      const double path_m = (path_length / si::meter).value();
      const double delta_time_seconds = sp > 0.0 ? path_m / sp : 0.0;
      CollisionPoint cp{coordinate_, vx, vy, vz, generation_};
      cp.z_number = static_cast<size_t>(properties_.charge_.value());
      cp.atomic_mass_number = properties_.mass_number_.value();
      cp.energy_eV = (velocity_.energy_ / electron_volt).value();
      // Pass delta time to Bombardment which will accumulate and set cumulative time
      record_collision_point(cp, delta_time_seconds);
    }
    bombardment_.record_energy_loss((long double)(e_energy_loss + nuclear_e_loss).value(), coordinate_.z_.value(), Bombardment::EnergyLossType::total);

    // magic_formula_scattering_angle() returns theta_cm / 2.
    auto laboratory_scatter_angle =
        laboratory_scattering_angle(2.0 * cm_scatter_angle, m_1, m_2);
    // Per-element displacement energy of the struck atom, if the layer
    // defines one; otherwise the simulation-wide threshold.
    const energy_quantity threshold =
        current_layer.get_displacement_energy(atom_index).value_or(knock_on_threshold);
    if (nuclear_e_loss > threshold && enable_knock_on) {
      VLOG(2) << "Knock-on particle created with energy " << nuclear_e_loss;
      data_leaf_.add_child(ParticleDataLeaf());
      bombardment_.record_displacement(coordinate_.z_.value());
      auto knock_on_particle =
          Ion(static_cast<size_t>(std::lround(z_2)),
            static_cast<size_t>(std::round((m_2 / atomic_mass_unit).value())),
            simulation_, uniform_distribution_, random_number_generator_,
            data_leaf_.get_children().back(), bombardment_, generation_ + 1,
            current_time_seconds());
      // Recoil angle from transverse momentum conservation. NOTE: the asin
      // argument can exceed 1 through round-off (-> NaN), and the recoil's
      // azimuth is sampled independently of the projectile's (it should be
      // opposite). See docs/KNOWN_ISSUES.md.
      auto scattering_angle_2 = asin(sin(laboratory_scatter_angle) *
                                     sqrt(energy * m_1 / nuclear_e_loss / m_2));
      auto velocity = Velocity(velocity_.quaternion_, nuclear_e_loss);
      knock_on_particle.set_velocity(velocity);
      knock_on_particle.set_coordinate(coordinate_);
      knock_on_particle.rotate_inplace(scattering_angle_2);
      knock_on_particle.fire(); // Depth-first recursion through the cascade.
    }
    rotate_inplace(laboratory_scatter_angle);
    auto &velocity = velocity_.quaternion_;

    VLOG(2) << "Original coordinate (" << coordinate_.x_ << ", "
            << coordinate_.y_ << ", " << coordinate_.z_ << ")";
    VLOG(2) << "Original velocity (" << velocity.x_ << ", " << velocity.y_
            << ", " << velocity.z_ << ")";

    data_leaf_.add_coordinate(coordinate_);
    coordinates_.emplace_back(coordinate_);
    coordinate_.x_ += path_length * velocity.x_;
    coordinate_.y_ += path_length * velocity.y_;
    coordinate_.z_ += path_length * velocity.z_;

    if (coordinate_.z_ < 0.0 * angstrom) {
      data_leaf_.sputter();
      bombardment_.record_sputtered();
      break;
    }

    if (coordinate_.z_.value() > max_depth.value()) {
      data_leaf_.escape();
      break;
    }
  }
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

length_quantity Ion::impact_parameter(number_density_quantity number_density,
                      length_quantity length,
                      dimensionless_quantity reduced_energy) const {
  length_quantity result;
  auto random_number = uniform_distribution_(random_number_generator_);
  if (reduced_energy > 10.0) {
    result = sqrt(-log(random_number) / (pi * number_density * length));
    VLOG(2) << "Impact parameter at high energy for number density "
            << number_density << ", length " << length
            << ", and reduced energy " << reduced_energy << " is " << result;
  } else {
    result = sqrt(random_number /
                  (pi * pow<static_rational<2, 3>>(number_density)));
    VLOG(2) << "Impact parameter at low energy for number density "
            << number_density << ", length " << length
            << ", and reduced energy " << reduced_energy << " is " << result;
  }
  return result;
}

length_quantity Ion::interatomic_spacing(number_density_quantity number_density) const {
  return 1.0 / root<3>(number_density);
}

energy_quantity Ion::electronic_stopping_energy(dimensionless_quantity charge,
                                mass_quantity mass, energy_quantity energy,
                                number_density_quantity number_density,
                                const Properties &properties) const {
  if (simulation_.electronic_stopping_function_ == nullptr) {
    // Curve-based models (Bethe-Bloch, built-in, SRIM tables): the
    // stopping cross section of the struck element times N times one step.
    const auto *curve = simulation_.stopping_curve(
        static_cast<unsigned>(charge.value()), mass,
        static_cast<unsigned>(properties.charge_.value()));
    if (curve == nullptr) {
      throw std::logic_error("No electronic stopping curve for this projectile/element pair");
    }
    return curve->cross_section(energy) * number_density *
           interatomic_spacing(number_density);
  }
  return simulation_.electronic_stopping_function_(
      charge, mass, energy, number_density,
      interatomic_spacing(number_density), properties,
      custom_electronic_coefficient_, custom_electronic_exponent_);
}

// k_L sqrt(E) stopping (the ZBL low-energy / LSS-type form). NOTE: the
// origin of the c_2 = 1.59 A^3 factor is undocumented; see
// docs/KNOWN_ISSUES.md.
energy_quantity Simulation::zbl_electronic_stopping(
    dimensionless_quantity charge, mass_quantity mass, energy_quantity energy,
    number_density_quantity number_density, length_quantity travel_length,
    const Particle::Properties &target, double, double) {
  VLOG(2) << "Calculating ion electronic stopping energy with charge " << charge
          << ", mass " << mass << ", and energy " << energy;
  auto c_1 =
      1.212 * root<2>(constants::m_u) * root<2>(electron_volt) / angstrom;
  VLOG(3) << "Value of c_1 " << c_1;
  const auto k_l = c_1 * pow<static_rational<7, 6>>(charge) *
                   target.charge_ /
                   (pow<static_rational<3, 2>>(
                        pow<static_rational<2, 3>>(charge) +
                        pow<static_rational<2, 3>>(target.charge_)) *
                    root<2>(mass));
  VLOG(3) << "Value of k_l " << k_l;
  const auto stopping_energy = k_l * root<2>(energy);
  VLOG(3) << "Value of stopping_energy " << stopping_energy;
  const auto c_2 = 1.59 * pow<3>(angstrom);
  VLOG(3) << "Value of c_2 " << c_2;
  VLOG(3) << "Value of travel_length " << travel_length;
  const auto energy_loss =
      c_2 * travel_length * number_density * stopping_energy;
  VLOG(3) << "Value of energy_loss " << energy_loss;
  return energy_loss;
}

energy_quantity Simulation::lss_electronic_stopping(
    dimensionless_quantity charge, mass_quantity mass, energy_quantity energy,
    number_density_quantity number_density, length_quantity travel_length,
    const Particle::Properties &target, double, double) {
  const auto combined_charge =
      pow(pow(charge, 2.0 / 3.0) + pow(target.charge_, 2.0 / 3.0), 1.5);
  const auto screening_factor = pow(charge, 1.0 / 6.0);
  const auto speed = sqrt(2.0 * energy / mass);
  const auto bohr_speed = fine_structure_constant * si::constants::codata::c;
  const auto stopping_cross_section =
      8.0 * pi * coulomb_energy_length * bohr_radius * screening_factor *
      charge * target.charge_ / combined_charge * (speed / bohr_speed);
  return stopping_cross_section * travel_length * number_density;
}

// S_e/N = C * E[keV]^(n/100) in units of 1e-18 keV cm^2 (= 1e-19 eV m^2).
energy_quantity Simulation::custom_power_law_electronic_stopping(
    dimensionless_quantity, mass_quantity, energy_quantity energy,
    number_density_quantity number_density, length_quantity travel_length,
    const Particle::Properties &, double coefficient, double exponent) {
  const auto energy_keV = energy / kilo_electron_volt;
  const auto stopping_cross_section =
      coefficient * pow(energy_keV, exponent / 100.0) * 1e-19 *
      electron_volt * pow<2>(si::meter);
  return stopping_cross_section * travel_length * number_density;
}

} // namespace common
} // namespace amcset
