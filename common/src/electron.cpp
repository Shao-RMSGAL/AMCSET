/*! \file electron.cpp
 *  \brief Electron transport implementation.
 */

#include "amcset_common.h"

#include <absl/log/log.h>

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

namespace {
std::vector<std::vector<std::array<double, 6>>> g_mott_scattering_params;
std::vector<std::array<double, 6>> g_electron_screening_params;
std::once_flag g_electron_tables_initialized;
std::mutex g_electron_tables_mutex;

} // namespace

void Electron::initialize_tables(){
                                initialize_tables("temp.csv","temp.csv");
                               }

void Electron::initialize_tables(const std::string &screening_path,
                               const std::string &mott_path) {
  const std::filesystem::path screening_file =
      std::filesystem::path(screening_path).is_absolute()
          ? std::filesystem::path(screening_path)
          : std::filesystem::current_path() / screening_path;
  const std::filesystem::path mott_file =
      std::filesystem::path(mott_path).is_absolute()
          ? std::filesystem::path(mott_path)
          : std::filesystem::current_path() / mott_path;

  std::call_once(g_electron_tables_initialized, [screening_file, mott_file]() {
    std::lock_guard<std::mutex> lock(g_electron_tables_mutex);
    const std::filesystem::path assets_path =
        std::filesystem::current_path() / "common" / "assets";
    const auto mott_table_path = mott_file.empty() ? assets_path / "mott_scattering_parameters.csv" : mott_file;
    const auto screening_table_path = screening_file.empty() ? assets_path / "electron_screening_parameters.csv" : screening_file;

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

const std::vector<std::array<double, 6>> &Electron::mott_table_for(size_t z_number) {
  initialize_tables();
  const size_t index = std::min(z_number, g_mott_scattering_params.size() - 1);
  return g_mott_scattering_params[index];
}

const std::vector<std::array<double, 6>> &Electron::screening_table_for(size_t z_number) {
  initialize_tables();
  static std::vector<std::array<double, 6>> fallback = {{}};
  if (g_electron_screening_params.empty()) {
    return fallback;
  }
  const size_t index = std::min(z_number, g_electron_screening_params.size() - 1U);
  static std::vector<std::array<double, 6>> singleton;
  singleton.clear();
  singleton.push_back(g_electron_screening_params[index]);
  return singleton;
}

double Electron::mott_differential_cross_section(double theta_rad, double substrate_charge,
                                                double energy_keV) const {
  const double corrected_energy = std::max(energy_keV, 1.0e-6);
  const double beta = std::sqrt(1.0 - 1.0 / ((corrected_energy / 511.0 + 1.0) *
                                            (corrected_energy / 511.0 + 1.0)));
  const double beta_squared = beta * beta;
  constexpr double beta_ave = 0.7181287;
  const double gamma = corrected_energy / 511.0 + 1.0;
  constexpr double elec_mass = 1.0;
  constexpr double vc = 137.0;
  constexpr double re1 = 2.817938E-13;
  const double one_minus_cos = 1.0 - std::cos(theta_rad);
  const double beta_quarted = beta_squared * beta_squared;
  const double dsc = std::pow(substrate_charge * re1, 2.0) * (1.0 - beta_squared) /
                     (beta_quarted * one_minus_cos * one_minus_cos);
  const double moment = 2.0 * beta * gamma * elec_mass * vc * std::sin(theta_rad / 2.0);
  const double moment_squared = moment * moment;

  double rmott = 0.0;
  const size_t z_index = std::clamp<size_t>(static_cast<size_t>(std::lround(std::max(1.0, substrate_charge))),
                                           1U, g_mott_scattering_params.size() - 1U);
  const auto &mott_coeffs = g_mott_scattering_params[z_index];
  for (size_t j = 0; j < mott_coeffs.size(); ++j) {
    double alpha1 = 0.0;
    for (size_t k = 0; k < 6; ++k) {
      alpha1 += mott_coeffs[j][k] * std::pow(beta - beta_ave, static_cast<int>(k));
    }
    rmott += alpha1 * std::pow(std::sqrt(one_minus_cos), static_cast<int>(j));
  }

  double feq = 0.0;
  const auto &screening = g_electron_screening_params[z_index];
  for (size_t i = 0; i < 3; ++i) {
    feq += screening[i] * screening[i + 3] * screening[i + 3] /
           (screening[i + 3] * screening[i + 3] + moment_squared);
  }

  return rmott * dsc * (1.0 - feq) * (1.0 - feq);
}

double Electron::mott_total_cross_section(double energy_keV, double substrate_charge,
                                         size_t divisor_count) const {
  if (divisor_count < 2) {
    return 0.0;
  }
  std::vector<double> division_angles;
  division_angles.reserve(divisor_count);
  for (size_t i = 0; i < divisor_count; ++i) {
    division_angles.push_back(find_angle(i, divisor_count));
  }

  double total_cross_section = 0.0;
  for (size_t i = 0; i + 1 < divisor_count; ++i) {
    const double theta_left = division_angles[i];
    const double theta_right = division_angles[i + 1];
    const double ds_left = mott_differential_cross_section(theta_left, substrate_charge, energy_keV);
    const double ds_right = mott_differential_cross_section(theta_right, substrate_charge, energy_keV);
    total_cross_section += 2.0 * pi * 0.5 * (ds_left * std::sin(theta_left) + ds_right * std::sin(theta_right)) *
                           (theta_right - theta_left);
  }
  return total_cross_section;
}

double Electron::generate_flying_distance(double total_cross_section,
                                          double number_density_m3_inv) const {
  if (total_cross_section <= 0.0 || number_density_m3_inv <= 0.0) {
    return 0.0;
  }

  // The layer stores number density in mol/m^3, while the Mott cross-section
  // formulas use atoms/m^3. The legacy implementation converted density to an
  // atom count before computing the free path and also added the interatomic
  // spacing offset to every flight segment.
  const double atoms_per_m3 = (number_density_m3_inv * constants::N_A).value();
  const double mean_free_path = 1.0 / (atoms_per_m3 * total_cross_section);
  const double atomic_spacing_m = std::pow(atoms_per_m3, -1.0 / 3.0);
  const double sample = -std::log(std::max(uniform_distribution_(random_number_generator_), 1e-12));
  return mean_free_path * sample + atomic_spacing_m;
}

double Electron::sample_scattering_angle(double total_cross_section, double substrate_charge,
                                         double energy_keV, size_t divisor_count) const {
  if (total_cross_section <= 0.0 || divisor_count < 2) {
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
    const double ds_left = mott_differential_cross_section(theta_left, substrate_charge, energy_keV);
    const double ds_right = mott_differential_cross_section(theta_right, substrate_charge, energy_keV);
    const double segment = 2.0 * pi * 0.5 * (ds_left * std::sin(theta_left) + ds_right * std::sin(theta_right)) *
                           (theta_right - theta_left);
    if (segment <= 0.0) continue;
    cumulative += segment / total_cross_section;
    if (cumulative >= target) {
      return 0.5 * (theta_left + theta_right);
    }
  }
  return division_angles.back();
}

double Electron::elastic_energy_loss(double scattering_angle_rad, double energy_keV,
                                    double substrate_mass_number) const {
  const double m_target = substrate_mass_number * 931.0;
  const double sin2 = std::sin(scattering_angle_rad) * std::sin(scattering_angle_rad);
  const double cos2 = std::cos(scattering_angle_rad) * std::cos(scattering_angle_rad);

  const double numerator = ((energy_keV + 511.0) * sin2 + m_target * (1.0 - std::cos(scattering_angle_rad))) *
                          energy_keV * (energy_keV + 2.0 * 511.0);
  const double denominator = (energy_keV + m_target) * (energy_keV + m_target) -
                            energy_keV * (energy_keV + 2.0 * 511.0) * cos2;
  return numerator / denominator;
}

double Electron::ionization_energy_loss(double energy_keV, double substrate_charge,
                                       double substrate_density_kg_m3) const {
  const double beta = std::sqrt(1.0 - 1.0 / ((energy_keV / 511.0 + 1.0) * (energy_keV / 511.0 + 1.0)));
  const double k_for_loss = 0.734 * std::pow(substrate_charge, 0.037);
  const double j_for_loss = (substrate_charge < 13.0)
                               ? 0.0115 * substrate_charge
                               : 0.00976 * substrate_charge + 0.0585 * std::pow(substrate_charge, -0.19);
  const double corrected_j = j_for_loss / (1.0 + k_for_loss * j_for_loss / energy_keV);
  const double ionization_energy = 153.55e24 * substrate_density_kg_m3 / (6.02214076e23) /
                                  (beta * beta) *
                                  (std::log(511.0 * beta * beta * energy_keV /
                                            (2.0 * corrected_j * corrected_j * (1.0 - beta * beta))) -
                                   std::log(2.0) * (2.0 * std::sqrt(1.0 - beta * beta) - 1.0 + beta * beta) +
                                   1.0 - beta * beta +
                                   0.125 * std::pow(1.0 - std::sqrt(1.0 - beta * beta), 2.0));
  return std::max(0.0, ionization_energy);
}

double Electron::bremsstrahlung_loss(double energy_keV, double substrate_charge,
                                    double substrate_density_kg_m3) const {
  const double brems = 1.4e-4 * 1.0e24 * substrate_density_kg_m3 / (6.02214076e23) *
                       substrate_charge * (substrate_charge + 1.0) *
                       (energy_keV + 511.0) *
                       (4.0 * std::log(2.0 * (energy_keV + 511.0) / 511.0) - 4.0 / 3.0);
  return std::max(0.0, brems);
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
    const double substrate_mass = target.mass_number_.value();
    const double number_density_mol_per_m3 = layer.get_number_density().value();
    const double energy_keV = (velocity_.energy_ / kilo_electron_volt).value();
    const double total_cross_section = mott_total_cross_section(
        energy_keV, substrate_charge,
        simulation_.get_settings(&Simulation::Settings::divisor_angle_number_));
    const double step_m = generate_flying_distance(total_cross_section,
                                                  number_density_mol_per_m3);
    const double scattering_angle = sample_scattering_angle(
        total_cross_section, substrate_charge, energy_keV,
        simulation_.get_settings(&Simulation::Settings::divisor_angle_number_));

    const double elastic_loss_keV = elastic_energy_loss(
        scattering_angle, energy_keV, substrate_mass);
    const double ionization_loss_keV = ionization_energy_loss(
        energy_keV, substrate_charge, layer.get_mass_density().value())*step_m;
    const double brems_loss_keV = bremsstrahlung_loss(
        energy_keV, substrate_charge, layer.get_mass_density().value())*step_m;
    const double total_loss_keV = elastic_loss_keV + ionization_loss_keV + brems_loss_keV;
    const double next_energy_keV = std::max(0.0, energy_keV - total_loss_keV);

    bombardment_.record_energy_loss(
        (long double)energy_quantity(elastic_loss_keV * kilo_electron_volt).value(),
        coordinate_.z_.value(), Bombardment::EnergyLossType::nuclear);
    bombardment_.record_energy_loss(
        (long double)energy_quantity((ionization_loss_keV + brems_loss_keV) * kilo_electron_volt).value(),
        coordinate_.z_.value(), Bombardment::EnergyLossType::electronic);
    bombardment_.record_energy_loss(
        (long double)energy_quantity(total_loss_keV * kilo_electron_volt).value(),
        coordinate_.z_.value(), Bombardment::EnergyLossType::total);

    const auto step = length_quantity(step_m * si::meter);
    const auto dir_norm = std::sqrt(velocity_.quaternion_.x_ * velocity_.quaternion_.x_ +
                                   velocity_.quaternion_.y_ * velocity_.quaternion_.y_ +
                                   velocity_.quaternion_.z_ * velocity_.quaternion_.z_);
    const double ux = (dir_norm > 0.0) ? velocity_.quaternion_.x_ / dir_norm : 1.0;
    const double uy = (dir_norm > 0.0) ? velocity_.quaternion_.y_ / dir_norm : 0.0;
    const double uz = (dir_norm > 0.0) ? velocity_.quaternion_.z_ / dir_norm : 0.0;
    coordinate_.x_ += step * ux;
    coordinate_.y_ += step * uy;
    coordinate_.z_ += step * uz;

    velocity_.energy_ = energy_quantity(next_energy_keV * kilo_electron_volt);
    rotate_inplace(angle_quantity(scattering_angle * radian));

    const auto &q = velocity_.quaternion_;
    const auto speed = speed_from_energy(velocity_.energy_, properties_.mass_);
    const double vx = q.x_ * speed.value();
    const double vy = q.y_ * speed.value();
    const double vz = q.z_ * speed.value();
    // Compute delta time for this step using `step` length and particle speed
    double delta_time_seconds = 0.0;
    try {
      const double step_in_meters = (step / si::meter).value();
      const auto speed_q = speed_from_energy(velocity_.energy_, properties_.mass_);
      const double sp = speed_q.value();
      if (sp > 0.0) delta_time_seconds = step_in_meters / sp;
    } catch (...) {
      delta_time_seconds = 0.0;
    }

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
