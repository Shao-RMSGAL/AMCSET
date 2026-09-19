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
                   size_t generation) try
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
      generation_(generation),
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

namespace {
std::vector<std::vector<std::array<double, 6>>> g_mott_scattering_params;
std::vector<std::array<double, 6>> g_electron_screening_params;
std::once_flag g_electron_tables_initialized;
std::mutex g_electron_tables_mutex;

void load_electron_tables() {
  const std::filesystem::path assets_path =
      std::filesystem::current_path() / "common" / "assets";
  const auto mott_path = assets_path / "mott_scattering_parameters.csv";
  const auto screening_path = assets_path / "electron_screening_parameters.csv";

  std::ifstream mott_stream(mott_path);
  if (!mott_stream.is_open()) {
    throw std::runtime_error("Unable to open Mott parameter table: " + mott_path.string());
  }
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
    std::vector<std::array<double, mott_cols_per_row>> &element =
        g_mott_scattering_params[element_index];
    while (element.size() <= row_index % mott_rows_per_element) {
      element.push_back({});
    }
    for (size_t i = 0; i < col; ++i) {
      element[row_index % mott_rows_per_element][i] = vals[i];
    }
    ++row_index;
  }

  std::ifstream screening_stream(screening_path);
  if (!screening_stream.is_open()) {
    throw std::runtime_error("Unable to open electron screening data: " + screening_path.string());
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
}

} // namespace

void Electron::initialize_tables() {
  std::call_once(g_electron_tables_initialized, []() {
    std::lock_guard<std::mutex> lock(g_electron_tables_mutex);
    load_electron_tables();
  });
}

double Electron::find_angle(size_t i, size_t divisor_count) {
  if (divisor_count == 0) return 0.0;
  return -pi / 9.0 * (1.0 - std::pow(10.0, static_cast<double>(i + 1) / static_cast<double>(divisor_count)));
}

const std::vector<std::array<double, 6>> &Electron::mott_table_for(size_t z_number) {
  initialize_tables();
  const size_t index = std::min(z_number, g_mott_scattering_params.size() - 1);
  return g_mott_scattering_params[index];
}

const std::vector<std::array<double, 6>> &Electron::screening_table_for(size_t z_number) {
  initialize_tables();
  static std::vector<std::array<double, 6>> fallback = {{} };
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
  const size_t z_index = std::clamp<size_t>(static_cast<size_t>(std::lround(std::max(1.0, substrate_charge))), 1U, g_mott_scattering_params.size() - 1U);
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
  const double mean_free_path = 1.0 / (number_density_m3_inv * total_cross_section);
  const double sample = -std::log(std::max(uniform_distribution_(random_number_generator_), 1e-12));
  return mean_free_path * sample;
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
  const double energy_mev = energy_keV * 1e-3;
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
  initialize_tables();
  const auto stopping_threshold = simulation_.get_settings(&Simulation::Settings::electron_stopping_energy_);

  while (velocity_.energy_ > stopping_threshold) {
    const auto layer_pair = simulation_.get_volume().get_layer(coordinate_.z_);
    const Layer &layer = *layer_pair.first;
    const auto atom_index = layer.get_random_index(random_number_generator_);
    const auto &target = layer.get_property(atom_index);
    const double substrate_charge = target.charge_.value();
    const double substrate_mass = target.mass_number_.value();
    const double number_density = layer.get_number_density().value();
    const double energy_keV = (velocity_.energy_ / kilo_electron_volt).value();
    const double total_cross_section = mott_total_cross_section(energy_keV, substrate_charge,
                                                                simulation_.get_settings(&Simulation::Settings::divisor_angle_number_));
    const double step_m = generate_flying_distance(total_cross_section, number_density);
    const double scattering_angle = sample_scattering_angle(total_cross_section, substrate_charge,
                                                           energy_keV,
                                                           simulation_.get_settings(&Simulation::Settings::divisor_angle_number_));
    const double elastic_loss = elastic_energy_loss(scattering_angle, energy_keV, substrate_mass);
    const double ionization_loss = ionization_energy_loss(energy_keV, substrate_charge,
                                                        layer.get_mass_density().value());
    const double brems_loss = bremsstrahlung_loss(energy_keV, substrate_charge,
                                                layer.get_mass_density().value());
    const double next_energy_keV = std::max(0.0, energy_keV - elastic_loss - ionization_loss - brems_loss);
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

    auto updated_velocity = velocity_;
    updated_velocity.energy_ = energy_quantity(next_energy_keV * kilo_electron_volt);
    updated_velocity.quaternion_.normalize_inplace();
    rotate_inplace(angle_quantity(scattering_angle * radian));
    velocity_.energy_ = updated_velocity.energy_;

    const auto &q = velocity_.quaternion_;
    const auto speed = speed_from_energy(velocity_.energy_, properties_.mass_);
    const double vx = q.x_ * speed.value();
    const double vy = q.y_ * speed.value();
    const double vz = q.z_ * speed.value();
    bombardment_.add_collision_point(CollisionPoint{coordinate_, vx, vy, vz, generation_});

    if (coordinate_.z_ < length_quantity(0.0 * angstrom)) {
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

dimensionless_quantity Ion::screening_function(dimensionless_quantity reduced_radius) const {
  const auto result = 0.1818 * exp(-3.2 * reduced_radius) +
                      0.5099 * exp(-0.9423 * reduced_radius) +
                      0.2802 * exp(-0.4028 * reduced_radius) +
                      0.2817 * exp(-0.2016 * reduced_radius);
  VLOG(2) << "Screening function with reduced radius " << reduced_radius
          << " is " << result;
  return result;
}

length_quantity Ion::screening_length(dimensionless_quantity z_1,
                    dimensionless_quantity z_2) const {
  return Ion::screening_length_for_model(
    z_1, z_2,
    simulation_.get_settings(&Simulation::Settings::screening_model_) ==
      Simulation::Settings::ScreeningModel::QuickBasic);
}

length_quantity Ion::screening_length_for_model(dimensionless_quantity z_1,
                                                dimensionless_quantity z_2,
                                                bool use_quickbasic) {
  if (use_quickbasic) {
    auto numerator = 0.8854 * 0.529 * angstrom;
    auto denom = pow(pow(z_1, 0.5) + pow(z_2, 0.5), 2.0 / 3.0);
    return numerator / denom;
  }
  return 0.25 * cbrt(9.0 * pi * pi / 2.0) * bohr_radius /
         (pow(z_1, 0.23) + pow(z_2, 0.23));
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

voltage_quantity Ion::screening_potential_derivative(length_quantity radius,
                                    dimensionless_quantity z_1,
                                    dimensionless_quantity z_2) const {
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
  VLOG(2) << "Screening potential derivative with radius " << radius << ", z_1 " << z_1
          << ", and z_2 " << z_2 << " is " << result;
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

length_quantity Ion::free_flying_path_length(mass_quantity m_1, mass_quantity m_2,
                             dimensionless_quantity reduced_energy,
                             length_quantity screening_length,
                             number_density_quantity number_density) const {
  length_quantity L;
  L = pow<static_rational<1, 3>>(1.0 / number_density / si::constants::codata::N_A);
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
      4 * z_1 * z_2 * e_statcoulomb * e_statcoulomb / (rd_mass * speed * speed);
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
        if (copysign(old_pot, mid_pot) == old_pot) {
          old_result = mid_result;
          old_pot = mid_pot;
        } else {
          result = mid_result;
          pot = mid_pot;
        }
        binary_search_iterations++;
      }
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
  auto G = gamma / (sqrt(1 + A * A) - A);
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
  auto result = 4.0 * atom_mass * target_mass /
                pow<2>(atom_mass + target_mass) * incident_energy *
                pow<2>(sin(scattering_angle));
  VLOG(2) << "Nuclear scattering energy loss atom mass " << atom_mass
          << " target_mass " << target_mass << " incident_energy "
          << incident_energy << " scattering_angle " << scattering_angle
          << " has result " << result;
  return result;
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

  size_t incident = (velocity_.energy_ == simulation_.get_settings(&Simulation::Settings::incident_energy_)) ? 0 : 1;
  if (incident == 1) {
    VLOG(2) << "Incident particle is a knock-on particle.";
  } else {
    VLOG(2) << "Incident particle is the primary particle.";
  }

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

    auto e_energy_loss = electronic_stopping_energy(
        z_1, m_1, energy, number_density, layer_properties);
    VLOG(2) << "Electronic stopping energy loss: " << e_energy_loss;
    energy -= e_energy_loss;
    VLOG(2) << "Energy after electronic stopping: " << energy;
    bombardment_.record_energy_loss((long double)e_energy_loss.value(), coordinate_.z_.value(), Bombardment::EnergyLossType::electronic);

    auto speed = speed_from_energy(energy, m_1);
    VLOG(3) << "Speed from energy: " << speed;
    auto c_mass_energy = cm_energy(m_1, speed);
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
      double vx = 0.0, vy = 0.0, vz = 0.0;
      try {
        auto speed_v = speed_from_energy(velocity_.energy_, properties_.mass_);
        double sp = speed_v.value();
        const auto &q = velocity_.quaternion_;
        vx = q.x_ * sp;
        vy = q.y_ * sp;
        vz = q.z_ * sp;
      } catch (...) {
        vx = vy = vz = 0.0;
      }
      CollisionPoint cp{coordinate_, vx, vy, vz, generation_};
      bombardment_.add_collision_point(cp);
    }
    bombardment_.record_energy_loss((long double)(e_energy_loss + nuclear_e_loss).value(), coordinate_.z_.value(), Bombardment::EnergyLossType::total);

    auto laboratory_scatter_angle = laboratory_scattering_angle(cm_scatter_angle, m_1, m_2);
    if (nuclear_e_loss > knock_on_threshold && enable_knock_on) {
      VLOG(2) << "Knock-on particle created with energy " << nuclear_e_loss;
      data_leaf_.add_child(ParticleDataLeaf());
      bombardment_.record_displacement(coordinate_.z_.value());
      auto knock_on_particle =
          Ion(z_2,
            static_cast<size_t>(std::round((m_2 / atomic_mass_unit).value())),
            simulation_, uniform_distribution_, random_number_generator_,
            data_leaf_.get_children().back(), bombardment_, generation_ + 1);
      auto scattering_angle_2 = asin(sin(laboratory_scatter_angle) *
                                     sqrt(energy * m_1 / nuclear_e_loss / m_2));
      auto velocity = Velocity(velocity_.quaternion_, nuclear_e_loss);
      knock_on_particle.set_velocity(velocity);
      knock_on_particle.set_coordinate(coordinate_);
      knock_on_particle.rotate_inplace(scattering_angle_2);
      knock_on_particle.fire();
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

length_quantity Ion::interatomic_spacing(number_density_quantity number_density) const {
  return root<-3>(number_density * constants::N_A);
}

energy_quantity Ion::electronic_stopping_energy(dimensionless_quantity charge,
                                mass_quantity mass, energy_quantity energy,
                                number_density_quantity number_density,
                                const Properties &properties) const try {
  VLOG(2) << "Calculating ion electronic stopping energy with charge " << charge
          << ", mass " << mass << ", and energy " << energy;
  auto c_1 =
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
  return energy_loss;
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

} // namespace common
} // namespace amcset
