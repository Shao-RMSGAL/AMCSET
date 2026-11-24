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

/*! \file electron.cpp
 *  \brief Electron transport: Mott elastic scattering with tabulated
 *  parameters plus continuous ionization and bremsstrahlung losses.
 *
 *  All physical quantities use Boost.Units. The Mott and screening
 *  parameter tables are dimensionless fits in which the momentum transfer is
 *  expressed in atomic units (hbar / a_0).
 *
 *  Physics fixes in v0.1.0 (docs/KNOWN_ISSUES.md, fix P5): target rest
 *  energy (was A * 931, an MeV number mixed with keV), table index Z - 1
 *  (was Z), cross sections in m^2 with number densities in m^-3 (were cm^2
 *  and mol/m^3), collision stopping proportional to the electron density
 *  N Z (the Z/A factor was missing), bremsstrahlung proportional to the
 *  atomic density, and relativistic speeds.
 */

#include "amcset_common.h"

#include <absl/log/log.h>

#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace common {

namespace {
// Process-global parameter tables, filled once by Electron::initialize_tables.
//
// g_mott_scattering_params[i] holds the 5x6 Mott coefficient block of the
// (i+1)-th row group in the CSV, i.e. index 0 is Z = 1 (the file starts at
// hydrogen). g_electron_screening_params[i] likewise holds one row per element
// starting at Z = 1. Both are therefore indexed with Z - 1.
std::vector<std::vector<std::array<double, 6>>> g_mott_scattering_params;
std::vector<std::array<double, 6>> g_electron_screening_params;
std::once_flag g_electron_tables_initialized;

} // namespace

void Electron::initialize_tables(const std::string &screening_path,
                               const std::string &mott_path) {
  // An empty path selects the bundled table under <cwd>/common/assets.
  // Relative paths are resolved against the current working directory.
  const std::filesystem::path assets_path =
      std::filesystem::current_path() / "common" / "assets";
  const auto resolve = [](const std::string &path,
                          const std::filesystem::path &fallback) {
    if (path.empty()) {
      return fallback;
    }
    const std::filesystem::path p(path);
    return p.is_absolute() ? p : std::filesystem::current_path() / p;
  };
  const auto screening_table_path =
      resolve(screening_path, assets_path / "electron_screening_parameters.csv");
  const auto mott_table_path =
      resolve(mott_path, assets_path / "mott_scattering_parameters.csv");

  // If the lambda throws, call_once does not mark the flag as done, so the
  // next call retries.
  std::call_once(g_electron_tables_initialized, [screening_table_path,
                                                 mott_table_path]() {
    std::ifstream mott_stream(mott_table_path);
    if (!mott_stream.is_open()) {
      throw std::runtime_error("Unable to open Mott parameter table: " + mott_table_path.string());
    }
    std::ifstream screening_stream(screening_table_path);
    if (!screening_stream.is_open()) {
      throw std::runtime_error("Unable to open electron screening data: " + screening_table_path.string());
    }

    g_mott_scattering_params.clear();
    g_electron_screening_params.clear();
    std::string line;
    size_t element_index = 0;
    size_t row_index = 0;
    constexpr size_t mott_rows_per_element = 5;
    constexpr size_t mott_cols_per_row = 6;
    while (std::getline(mott_stream, line)) {
      if (line.empty()) continue;
      std::stringstream iss(line);
      std::array<double, mott_cols_per_row> vals{};
      std::string token;
      size_t col = 0;
      while (std::getline(iss, token, ',')) {
        if (col >= mott_cols_per_row) break;
        vals[col++] = std::stod(token);
      }
      if (col == 0) continue;
      if (row_index % mott_rows_per_element == 0 && row_index > 0) {
        ++element_index;
      }
      if (g_mott_scattering_params.size() <= element_index) {
        g_mott_scattering_params.resize(element_index + 1);
      }
      auto &element = g_mott_scattering_params[element_index];
      while (element.size() <= row_index % mott_rows_per_element) {
        element.push_back({});
      }
      for (size_t i = 0; i < col; ++i) {
        element[row_index % mott_rows_per_element][i] = vals[i];
      }
      ++row_index;
    }
    while (std::getline(screening_stream, line)) {
      if (line.empty()) continue;
      std::stringstream iss(line);
      std::array<double, 6> vals{};
      std::string token;
      size_t col = 0;
      while (std::getline(iss, token, ',')) {
        if (col >= 6) break;
        vals[col++] = std::stod(token);
      }
      if (col == 0) continue;
      g_electron_screening_params.push_back(vals);
    }
  });
}

double Electron::find_angle(size_t i, size_t divisor_count) {
  if (divisor_count == 0) return 0.0;
  return -pi / 9.0 * (1.0 - std::pow(10.0, static_cast<double>(i + 1) /
                                      static_cast<double>(divisor_count)));
}

namespace {

//! Classical electron radius r_e = e^2 / (4 pi epsilon_0 m_e c^2).
const length_quantity classical_electron_radius =
    coulomb_energy_length / electron_rest_energy;

//! Kinetic energy -> beta and gamma of an electron.
struct ElectronKinematics {
  double gamma;
  double beta;
  explicit ElectronKinematics(energy_quantity energy) {
    gamma = 1.0 + (energy / electron_rest_energy).value();
    beta = std::sqrt(std::max(0.0, 1.0 - 1.0 / (gamma * gamma)));
  }
};

} // namespace

area_quantity Electron::mott_differential_cross_section(double theta_rad,
                                                        double substrate_charge,
                                                        energy_quantity energy) const {
  const energy_quantity corrected_energy =
      std::max(energy, 1.0e-6 * kilo_electron_volt);
  const ElectronKinematics k(corrected_energy);
  const double beta_squared = k.beta * k.beta;
  constexpr double beta_ave = 0.7181287; // expansion point of the Mott fits
  const double one_minus_cos = 1.0 - std::cos(theta_rad);
  const double beta_quarted = beta_squared * beta_squared;
  // Relativistic Rutherford cross section.
  const area_quantity rutherford =
      substrate_charge * substrate_charge * classical_electron_radius *
      classical_electron_radius * (1.0 - beta_squared) /
      (beta_quarted * one_minus_cos * one_minus_cos);
  // Momentum transfer q = 2 p sin(theta/2) in atomic units (hbar/a_0):
  // p = beta gamma m_e c, and m_e c = 1/alpha in atomic units.
  const double moment = 2.0 * k.beta * k.gamma / fine_structure_constant.value() *
                        std::sin(theta_rad / 2.0);
  const double moment_squared = moment * moment;

  // Tables start at Z = 1 in slot 0.
  const size_t z_index =
      std::clamp<size_t>(static_cast<size_t>(std::lround(std::max(1.0, substrate_charge))),
                         1U, g_mott_scattering_params.size()) -
      1U;
  double rmott = 0.0;
  const auto &mott_coeffs = g_mott_scattering_params[z_index];
  for (size_t j = 0; j < mott_coeffs.size(); ++j) {
    double alpha1 = 0.0;
    for (size_t kk = 0; kk < 6; ++kk) {
      alpha1 += mott_coeffs[j][kk] * std::pow(k.beta - beta_ave, static_cast<int>(kk));
    }
    rmott += alpha1 * std::pow(std::sqrt(one_minus_cos), static_cast<int>(j));
  }

  double feq = 0.0;
  const auto &screening = g_electron_screening_params[std::min(
      z_index, g_electron_screening_params.size() - 1U)];
  for (size_t i = 0; i < 3; ++i) {
    feq += screening[i] * screening[i + 3] * screening[i + 3] /
           (screening[i + 3] * screening[i + 3] + moment_squared);
  }

  return rmott * rutherford * (1.0 - feq) * (1.0 - feq);
}

area_quantity Electron::mott_total_cross_section(energy_quantity energy,
                                                 double substrate_charge,
                                                 size_t divisor_count) const {
  area_quantity total = 0.0 * si::square_meter;
  if (divisor_count < 2) {
    return total;
  }
  std::vector<double> division_angles;
  division_angles.reserve(divisor_count);
  for (size_t i = 0; i < divisor_count; ++i) {
    division_angles.push_back(find_angle(i, divisor_count));
  }
  for (size_t i = 0; i + 1 < divisor_count; ++i) {
    const double theta_left = division_angles[i];
    const double theta_right = division_angles[i + 1];
    const auto ds_left = mott_differential_cross_section(theta_left, substrate_charge, energy);
    const auto ds_right = mott_differential_cross_section(theta_right, substrate_charge, energy);
    total += 2.0 * pi * 0.5 *
             (ds_left * std::sin(theta_left) + ds_right * std::sin(theta_right)) *
             (theta_right - theta_left);
  }
  return total;
}

length_quantity Electron::generate_flying_distance(area_quantity total_cross_section,
                                                   number_density_quantity atomic_density) const {
  if (total_cross_section <= 0.0 * si::square_meter ||
      atomic_density <= 0.0 * per_cubic_meter) {
    return 0.0 * si::meter;
  }
  const length_quantity mean_free_path = 1.0 / (atomic_density * total_cross_section);
  // The model adds one interatomic spacing to every flight segment.
  const length_quantity atomic_spacing = 1.0 / root<3>(atomic_density);
  const double sample =
      -std::log(std::max(uniform_distribution_(random_number_generator_), 1e-12));
  return mean_free_path * sample + atomic_spacing;
}

double Electron::sample_scattering_angle(area_quantity total_cross_section,
                                         double substrate_charge,
                                         energy_quantity energy,
                                         size_t divisor_count) const {
  if (total_cross_section <= 0.0 * si::square_meter || divisor_count < 2) {
    return 0.0;
  }
  std::vector<double> division_angles;
  division_angles.reserve(divisor_count);
  for (size_t i = 0; i < divisor_count; ++i) {
    division_angles.push_back(find_angle(i, divisor_count));
  }

  const double target = std::max(uniform_distribution_(random_number_generator_), 1e-12);
  double cumulative = 0.0;
  for (size_t i = 0; i + 1 < divisor_count; ++i) {
    const double theta_left = division_angles[i];
    const double theta_right = division_angles[i + 1];
    const auto ds_left = mott_differential_cross_section(theta_left, substrate_charge, energy);
    const auto ds_right = mott_differential_cross_section(theta_right, substrate_charge, energy);
    const area_quantity segment =
        2.0 * pi * 0.5 *
        (ds_left * std::sin(theta_left) + ds_right * std::sin(theta_right)) *
        (theta_right - theta_left);
    if (segment <= 0.0 * si::square_meter) continue;
    cumulative += (segment / total_cross_section).value();
    if (cumulative >= target) {
      return 0.5 * (theta_left + theta_right);
    }
  }
  return division_angles.back();
}

energy_quantity Electron::elastic_energy_loss(double scattering_angle_rad,
                                              energy_quantity energy,
                                              mass_quantity target_mass) const {
  // Exact relativistic energy transfer from an electron of kinetic energy E
  // scattered by theta to a nucleus of rest energy M c^2 initially at rest.
  const energy_quantity mc2 = electron_rest_energy;
  const energy_quantity Mc2 = target_mass * speed_of_light * speed_of_light;
  const double s = std::sin(scattering_angle_rad);
  const double c = std::cos(scattering_angle_rad);
  const auto numerator =
      ((energy + mc2) * s * s + Mc2 * (1.0 - c)) * energy * (energy + 2.0 * mc2);
  const auto denominator =
      (energy + Mc2) * (energy + Mc2) - energy * (energy + 2.0 * mc2) * c * c;
  return numerator / denominator;
}

stopping_power_quantity Electron::ionization_stopping_power(
    energy_quantity energy, double substrate_charge,
    number_density_quantity atomic_density) const {
  const ElectronKinematics k(energy);
  const double beta2 = k.beta * k.beta;
  if (beta2 <= 0.0) return 0.0 * electron_volt_per_angstrom;
  // Mean ionization energy J (keV) with the Joy-Luo (1989) correction
  // J / (1 + k J / E). The empirical fits take J and E in keV.
  const double energy_keV = in_units(energy, kilo_electron_volt);
  const double k_for_loss = 0.734 * std::pow(substrate_charge, 0.037);
  const double j_for_loss =
      (substrate_charge < 13.0)
          ? 0.0115 * substrate_charge
          : 0.00976 * substrate_charge + 0.0585 * std::pow(substrate_charge, -0.19);
  const double corrected_j = j_for_loss / (1.0 + k_for_loss * j_for_loss / energy_keV);
  const double mc2_keV = in_units(electron_rest_energy, kilo_electron_volt);
  const double bracket =
      std::log(mc2_keV * beta2 * energy_keV /
               (2.0 * corrected_j * corrected_j * (1.0 - beta2))) -
      std::log(2.0) * (2.0 * std::sqrt(1.0 - beta2) - 1.0 + beta2) + 1.0 - beta2 +
      0.125 * std::pow(1.0 - std::sqrt(1.0 - beta2), 2.0);
  // 2 pi r_e^2 m_e c^2 n_e / beta^2 with the electron density n_e = N Z.
  const number_density_quantity electron_density = atomic_density * substrate_charge;
  const stopping_power_quantity s = 2.0 * pi * classical_electron_radius *
                                    classical_electron_radius * electron_rest_energy *
                                    electron_density / beta2 * bracket;
  return std::max(s, 0.0 * electron_volt_per_angstrom);
}

stopping_power_quantity Electron::bremsstrahlung_stopping_power(
    energy_quantity energy, double substrate_charge,
    number_density_quantity atomic_density) const {
  // High-energy (Bethe-Heitler, no screening) radiative stopping power
  // 4 alpha r_e^2 N Z (Z + 1) E_tot (ln(2 E_tot / m_e c^2) - 1/3).
  const energy_quantity total_energy = energy + electron_rest_energy;
  const stopping_power_quantity s =
      4.0 * fine_structure_constant.value() * classical_electron_radius *
      classical_electron_radius * atomic_density * substrate_charge *
      (substrate_charge + 1.0) * total_energy *
      (std::log(2.0 * (total_energy / electron_rest_energy).value()) - 1.0 / 3.0);
  return std::max(s, 0.0 * electron_volt_per_angstrom);
}

void Electron::fire() try {
  initialize_tables(
      simulation_.get_settings(&Simulation::Settings::electron_screening_parameters_path_),
      simulation_.get_settings(&Simulation::Settings::mott_scattering_parameters_path_));
  const auto stopping_threshold = simulation_.get_settings(&Simulation::Settings::electron_stopping_energy_);

  while (velocity_.energy_ > stopping_threshold) {
    const auto layer_pair = simulation_.get_volume().get_layer(coordinate_.z_);
    const Layer &layer = *layer_pair.first;
    const auto atom_index = layer.get_random_index(random_number_generator_);
    const auto &target = layer.get_property(atom_index);
    const double substrate_charge = target.charge_.value();
    const number_density_quantity atomic_density = layer.get_number_density();
    const energy_quantity energy = velocity_.energy_;
    const size_t divisors =
        simulation_.get_settings(&Simulation::Settings::divisor_angle_number_);
    const area_quantity total_cross_section =
        mott_total_cross_section(energy, substrate_charge, divisors);
    const length_quantity step_length =
        generate_flying_distance(total_cross_section, atomic_density);
    const double scattering_angle =
        sample_scattering_angle(total_cross_section, substrate_charge, energy, divisors);

    const energy_quantity elastic_loss =
        elastic_energy_loss(scattering_angle, energy, target.mass_);
    const energy_quantity ionization_loss =
        ionization_stopping_power(energy, substrate_charge, atomic_density) * step_length;
    const energy_quantity brems_loss =
        bremsstrahlung_stopping_power(energy, substrate_charge, atomic_density) * step_length;
    const energy_quantity total_loss = elastic_loss + ionization_loss + brems_loss;
    const energy_quantity next_energy =
        std::max(energy - total_loss, 0.0 * electron_volt);

    bombardment_.record_energy_loss(static_cast<long double>(elastic_loss.value()),
                                    coordinate_.z_.value(),
                                    Bombardment::EnergyLossType::nuclear);
    bombardment_.record_energy_loss(
        static_cast<long double>((ionization_loss + brems_loss).value()),
        coordinate_.z_.value(), Bombardment::EnergyLossType::electronic);
    bombardment_.record_energy_loss(static_cast<long double>(total_loss.value()),
                                    coordinate_.z_.value(),
                                    Bombardment::EnergyLossType::total);

    const auto step = step_length;
    const auto dir_norm = std::sqrt(velocity_.quaternion_.x_ * velocity_.quaternion_.x_ +
                                   velocity_.quaternion_.y_ * velocity_.quaternion_.y_ +
                                   velocity_.quaternion_.z_ * velocity_.quaternion_.z_);
    const double ux = (dir_norm > 0.0) ? velocity_.quaternion_.x_ / dir_norm : 1.0;
    const double uy = (dir_norm > 0.0) ? velocity_.quaternion_.y_ / dir_norm : 0.0;
    const double uz = (dir_norm > 0.0) ? velocity_.quaternion_.z_ / dir_norm : 0.0;
    coordinate_.x_ += step * ux;
    coordinate_.y_ += step * uy;
    coordinate_.z_ += step * uz;

    velocity_.energy_ = next_energy;
    rotate_inplace(angle_quantity(scattering_angle * radian));

    // Record the post-step state. The time step uses the relativistic speed
    // after this step's energy loss.
    const auto &q = velocity_.quaternion_;
    const ElectronKinematics after(velocity_.energy_);
    const double sp = in_units(after.beta * speed_of_light, meters_per_second);
    const double vx = q.x_ * sp;
    const double vy = q.y_ * sp;
    const double vz = q.z_ * sp;
    const double step_in_meters = (step / si::meter).value();
    const double delta_time_seconds = sp > 0.0 ? step_in_meters / sp : 0.0;

    CollisionPoint cp{coordinate_, vx, vy, vz, generation_};
    cp.z_number = static_cast<size_t>(properties_.charge_.value());
    cp.atomic_mass_number = properties_.mass_number_.value();
    cp.energy_eV = (velocity_.energy_ / electron_volt).value();
    record_collision_point(cp, delta_time_seconds);

    data_leaf_.add_coordinate(coordinate_);
    coordinates_.emplace_back(coordinate_);
    data_leaf_.add_sample(coordinate_, velocity_, "electron_state");

    if (coordinate_.z_ < 0.0 * angstrom) {
      data_leaf_.sputter();
      bombardment_.record_sputtered();
      break;
    }

    if (coordinate_.z_.value() > simulation_.get_volume().get_depth().value()) {
      data_leaf_.escape();
      break;
    }

    if (velocity_.energy_ <= stopping_threshold) {
      break;
    }
  }
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

} // namespace common
} // namespace amcset
