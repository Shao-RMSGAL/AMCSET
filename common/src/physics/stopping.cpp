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

/*! \file stopping.cpp
 *  \brief Electronic stopping models, SRIM table loading and StoppingCurve.
 */

#include "physics/stopping.h"

#include <map>
#include <mutex>
#include <tuple>

#include "element_data.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <fstream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string_view>

namespace amcset {
namespace physics {

using namespace common; // units, constants

namespace {

constexpr double ln10 = 2.302585092994046;

double parse_double(std::string_view token) {
  // std::from_chars for floating point is available in GCC 11+ and MSVC.
  double value = 0.0;
  const auto *first = token.data();
  const auto *last = token.data() + token.size();
  if (!token.empty() && *first == '+') ++first;
  const auto result = std::from_chars(first, last, value);
  if (result.ec != std::errc() || result.ptr != last) {
    throw std::runtime_error("Not a number: '" + std::string(token) + "'");
  }
  return value;
}

std::vector<std::string_view> split_ws(std::string_view line) {
  std::vector<std::string_view> out;
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i])))
      ++i;
    std::size_t j = i;
    while (j < line.size() &&
           !std::isspace(static_cast<unsigned char>(line[j])))
      ++j;
    if (j > i) out.emplace_back(line.substr(i, j - i));
    i = j;
  }
  return out;
}

//! Interpolate y(x) on increasing x (log-log tables); clamps index range.
double interp(const std::vector<double> &x, const std::vector<double> &y,
              double xv) {
  const auto it = std::upper_bound(x.begin(), x.end(), xv);
  std::size_t i = static_cast<std::size_t>(it - x.begin());
  if (i == 0) i = 1;
  if (i >= x.size()) i = x.size() - 1;
  const double t = (xv - x[i - 1]) / (x[i] - x[i - 1]);
  return y[i - 1] + t * (y[i] - y[i - 1]);
}

//! Sternheimer-Peierls (1971) general density-effect correction delta.
double density_effect_delta(double beta_gamma, energy_quantity mean_excitation,
                            number_density_quantity electron_density) {
  if (beta_gamma <= 0.0 || electron_density <= 0.0 * per_cubic_meter) return 0.0;
  // Plasma energy hbar*omega_p = hbar sqrt(n_e e^2 / (eps0 m_e)).
  const auto omega_squared = electron_density * elementary_charge *
                             elementary_charge /
                             (constants::epsilon_0 * electron_mass);
  const energy_quantity plasma_energy =
      constants::hbar * boost::units::sqrt(omega_squared);
  const double c_bar = 2.0 * std::log((mean_excitation / plasma_energy).value()) + 1.0;
  double x0 = 0.0;
  double x1 = 0.0;
  if (mean_excitation < 100.0 * electron_volt) {
    x0 = c_bar < 3.681 ? 0.2 : 0.326 * c_bar - 1.0;
    x1 = 2.0;
  } else {
    x0 = c_bar < 5.215 ? 0.2 : 0.326 * c_bar - 1.5;
    x1 = 3.0;
  }
  constexpr double m = 3.0;
  const double a = (c_bar - 2.0 * ln10 * x0) / std::pow(x1 - x0, m);
  const double x = std::log10(beta_gamma);
  if (x < x0) return 0.0;
  if (x < x1) return 2.0 * ln10 * x - c_bar + a * std::pow(x1 - x, m);
  return 2.0 * ln10 * x - c_bar;
}

/*!
 * Bethe-Bloch cross section. The stopping number is
 * L = ln(x) - beta^2 - delta/2 with x = (2 m_e c^2 beta^2 gamma^2 W_max)^(1/2) / I.
 * With \p low_velocity_term c > 0, ln(x) is replaced by
 * ln(1 + x + c I/(m_e c^2 beta^2)): the extra terms keep the high-energy
 * branch finite and large at low velocities (like the A4/E term of the
 * Andersen-Ziegler fits), which is what builtin_stopping() needs.
 */
stopping_cross_section_quantity
bethe_bloch_impl(unsigned z1, mass_quantity m1, unsigned z2,
                 energy_quantity mean_excitation, energy_quantity energy,
                 number_density_quantity electron_density,
                 const BetheBlochOptions &options, double low_velocity_term) {
  if (energy <= 0.0 * electron_volt) return 0.0 * srim_stopping_unit;
  const double beta = relativistic_beta(energy, m1);
  const double gamma = lorentz_gamma(energy, m1);
  if (beta <= 0.0) return 0.0 * srim_stopping_unit;
  const double beta2 = beta * beta;
  double z = static_cast<double>(z1);
  if (options.effective_charge) {
    z *= 1.0 - std::exp(-125.0 * beta * std::pow(static_cast<double>(z1), -2.0 / 3.0));
  }
  const double mass_ratio = (electron_mass / m1).value();
  const energy_quantity w_max =
      2.0 * electron_rest_energy * beta2 * gamma * gamma /
      (1.0 + 2.0 * gamma * mass_ratio + mass_ratio * mass_ratio);
  const double arg = (2.0 * electron_rest_energy * beta2 * gamma * gamma *
                      w_max / (mean_excitation * mean_excitation))
                         .value();
  const double delta =
      options.density_effect
          ? density_effect_delta(beta * gamma, mean_excitation, electron_density)
          : 0.0;
  double logarithm = 0.5 * std::log(arg); // = ln(sqrt(arg))
  if (low_velocity_term > 0.0) {
    const double low = low_velocity_term *
                       (mean_excitation / (electron_rest_energy * beta2)).value();
    logarithm = std::log(1.0 + std::sqrt(arg) + low);
  }
  const double bracket = logarithm - beta2 - 0.5 * delta;
  if (!(bracket > 0.0)) return 0.0 * srim_stopping_unit;
  return 4.0 * pi_v * coulomb_energy_length * coulomb_energy_length * z * z *
         static_cast<double>(z2) / (electron_rest_energy * beta2) * bracket;
}

} // namespace

// ---------------------------------------------------------------------------
// Kinematics
// ---------------------------------------------------------------------------

double lorentz_gamma(energy_quantity energy, mass_quantity mass) {
  return 1.0 + (energy / (mass * speed_of_light * speed_of_light)).value();
}

double relativistic_beta(energy_quantity energy, mass_quantity mass) {
  const double gamma = lorentz_gamma(energy, mass);
  if (gamma <= 1.0) return 0.0;
  return std::sqrt(1.0 - 1.0 / (gamma * gamma));
}

velocity_quantity relativistic_speed(energy_quantity energy,
                                     mass_quantity mass) {
  return relativistic_beta(energy, mass) * speed_of_light;
}

// ---------------------------------------------------------------------------
// Analytic models
// ---------------------------------------------------------------------------

stopping_cross_section_quantity lindhard_scharff_stopping(unsigned z1,
                                                          mass_quantity m1,
                                                          unsigned z2,
                                                          energy_quantity energy) {
  if (energy <= 0.0 * electron_volt) return 0.0 * srim_stopping_unit;
  const double Z1 = static_cast<double>(z1);
  const double Z2 = static_cast<double>(z2);
  const double xi_e = std::pow(Z1, 1.0 / 6.0);
  const double z_factor =
      Z1 * Z2 / std::pow(std::pow(Z1, 2.0 / 3.0) + std::pow(Z2, 2.0 / 3.0), 1.5);
  const double v_over_v0 =
      (relativistic_speed(energy, m1) / bohr_velocity).value();
  return xi_e * 8.0 * pi_v * coulomb_energy_length * bohr_radius * z_factor *
         v_over_v0;
}

stopping_cross_section_quantity
bethe_bloch_stopping(unsigned z1, mass_quantity m1, unsigned z2,
                     energy_quantity mean_excitation_energy,
                     energy_quantity energy,
                     number_density_quantity electron_density,
                     const BetheBlochOptions &options) {
  return bethe_bloch_impl(z1, m1, z2, mean_excitation_energy, energy,
                          electron_density, options, 0.0);
}

stopping_cross_section_quantity
builtin_stopping(unsigned z1, mass_quantity m1, unsigned z2,
                 energy_quantity energy,
                 number_density_quantity electron_density,
                 const BetheBlochOptions &options) {
  // Constants fitted once to the SRIM-2013 tables of H, He, O and Fe ions in
  // ten elements (1 keV/u .. 100 MeV/u): rms deviation ~25 % below
  // 100 keV/u, < 6 % above 1 MeV/u.
  constexpr double low_velocity_term = 0.5;
  constexpr double p = 4.0;
  const auto s_low = lindhard_scharff_stopping(z1, m1, z2, energy);
  const auto s_high = bethe_bloch_impl(
      z1, m1, z2, common::Element::by_z(z2).mean_excitation_energy(), energy,
      electron_density, options, low_velocity_term);
  if (s_low <= 0.0 * srim_stopping_unit) return s_high;
  if (s_high <= 0.0 * srim_stopping_unit) return s_low;
  const double a = in_units(s_low, srim_stopping_unit);
  const double b = in_units(s_high, srim_stopping_unit);
  return std::pow(std::pow(a, -p) + std::pow(b, -p), -1.0 / p) * srim_stopping_unit;
}

energy_squared_per_length_quantity
trim_electronic_straggling(unsigned z1, mass_quantity m1, unsigned z2,
                           energy_quantity energy,
                           number_density_quantity atomic_density) {
  const double Z1 = static_cast<double>(z1);
  const double Z2 = static_cast<double>(z2);
  const double beta = relativistic_beta(energy, m1);
  // Relativistically corrected Bohr straggling; b = c^2 / v^2.
  double relativistic = 1.0;
  if (beta > 0.0) {
    const double b = 1.0 / (beta * beta);
    relativistic = (b - 0.5) / (b - 1.0);
  }
  const energy_squared_per_length_quantity bohr =
      4.0 * pi_v * Z1 * Z1 * coulomb_energy_length * coulomb_energy_length *
      Z2 * atomic_density * relativistic;
  // Low-energy form of Eq. 7-47; units as printed in the SRIM book
  // (E in keV, N in atoms/Angstrom^3, M1 in amu, result in eV^2/Angstrom).
  const double low_value = std::pow(Z1 + Z2, 8.0 / 3.0) *
                           in_units(energy, kilo_electron_volt) *
                           in_units(atomic_density, per_cubic_angstrom) /
                           (3.2 * (m1 / atomic_mass_unit).value());
  const energy_squared_per_length_quantity low =
      low_value * electron_volt * electron_volt / angstrom;
  return std::min(bohr, low);
}

// ---------------------------------------------------------------------------
// StoppingTable
// ---------------------------------------------------------------------------

energy_quantity StoppingTable::min_energy() const {
  return std::exp(log_energy.front()) * kilo_electron_volt;
}

energy_quantity StoppingTable::max_energy() const {
  return std::exp(log_energy.back()) * kilo_electron_volt;
}

stopping_cross_section_quantity
StoppingTable::electronic(energy_quantity energy, mass_quantity mass) const {
  if (energy <= 0.0 * electron_volt) return 0.0 * srim_stopping_unit;
  // Same velocity <=> same energy per mass.
  const energy_quantity scaled = energy * (m1 / mass).value();
  const double le = std::log(in_units(scaled, kilo_electron_volt));
  if (le < log_energy.front()) {
    // Velocity-proportional extrapolation: S ~ sqrt(E).
    return std::exp(log_electronic.front() + 0.5 * (le - log_energy.front())) *
           srim_stopping_unit;
  }
  if (le > log_energy.back()) {
    return std::exp(log_electronic.back()) * srim_stopping_unit;
  }
  return std::exp(interp(log_energy, log_electronic, le)) * srim_stopping_unit;
}

stopping_cross_section_quantity
StoppingTable::nuclear(energy_quantity energy) const {
  if (energy <= 0.0 * electron_volt) return 0.0 * srim_stopping_unit;
  const double le = std::log(in_units(energy, kilo_electron_volt));
  const double clamped = std::clamp(le, log_energy.front(), log_energy.back());
  return std::exp(interp(log_energy, log_nuclear, clamped)) * srim_stopping_unit;
}

// ---------------------------------------------------------------------------
// SrimStoppingDatabase
// ---------------------------------------------------------------------------

namespace {

/*!
 * What the interpolation relies on: at least two points, equal lengths,
 * strictly increasing finite energies (interp() divides by their
 * differences) and finite values. \returns an empty string or the problem.
 */
std::string table_problem(const StoppingTable &t) {
  const std::size_t n = t.log_energy.size();
  if (n < 2 || t.log_electronic.size() != n || t.log_nuclear.size() != n) {
    return "needs at least two energies and one stopping value per energy";
  }
  for (std::size_t i = 0; i < n; ++i) {
    if (!std::isfinite(t.log_energy[i]) || !std::isfinite(t.log_electronic[i]) ||
        !std::isfinite(t.log_nuclear[i])) {
      return "has a value that is not a finite positive number (entry " + std::to_string(i + 1) + ")";
    }
    if (i > 0 && !(t.log_energy[i] > t.log_energy[i - 1])) {
      return "energies must increase (entry " + std::to_string(i + 1) + ")";
    }
  }
  if (!(t.m1 > 0.0 * kilogram) || t.z1 < 1 || t.z1 > Element::max_z || t.z2 < 1 || t.z2 > Element::max_z) {
    return "has an invalid projectile or target (Z or mass)";
  }
  return {};
}

} // namespace

void SrimStoppingDatabase::add(StoppingTable table) {
  if (const auto problem = table_problem(table); !problem.empty()) {
    throw std::invalid_argument("Malformed stopping table for Z1=" + std::to_string(table.z1) +
                                " Z2=" + std::to_string(table.z2) + ": " + problem);
  }
  const auto key = std::make_tuple(table.z1, table.z2, table.gas);
  std::lock_guard<std::mutex> lock(mutex_);
  unparsed_.erase(key);
  tables_[key] = std::move(table);
}

std::size_t SrimStoppingDatabase::size() const {
  std::lock_guard<std::mutex> lock(mutex_);
  // Parsed and unparsed keys are disjoint (find_locked moves a key).
  return tables_.size() + unparsed_.size();
}

const StoppingTable *SrimStoppingDatabase::find(unsigned z1, unsigned z2,
                                                bool gas) const {
  std::lock_guard<std::mutex> lock(mutex_);
  const StoppingTable *t = find_locked(std::make_tuple(z1, z2, gas));
  if (t == nullptr && gas) t = find_locked(std::make_tuple(z1, z2, false));
  return t;
}

namespace {

StoppingTable make_table(unsigned z1, double m1_amu, unsigned z2, bool gas,
                         const std::vector<double> &energy_keV,
                         const std::vector<double> &se,
                         const std::vector<double> &sn) {
  StoppingTable t;
  t.z1 = z1;
  t.m1 = m1_amu * atomic_mass_unit;
  t.z2 = z2;
  t.gas = gas;
  // Floor at a tiny positive value so that the logarithms are defined.
  const auto safe_log = [](double v) { return std::log(std::max(v, 1e-30)); };
  for (std::size_t i = 0; i < energy_keV.size(); ++i) {
    t.log_energy.push_back(std::log(energy_keV[i]));
    t.log_electronic.push_back(safe_log(se[i]));
    t.log_nuclear.push_back(safe_log(sn[i]));
  }
  return t;
}

void load_sr_directory(const std::filesystem::path &dir, SrimStoppingDatabase &db) {
  static const std::regex name_rx(R"(sr_(\d{3})_(\d{3})_([sg])\.txt)",
                                  std::regex::icase);
  static const std::regex mass_rx(R"(Mass\s*=\s*([0-9.]+))");
  for (const auto &entry : std::filesystem::directory_iterator(dir)) {
    if (!entry.is_regular_file()) continue;
    const auto name = entry.path().filename().string();
    std::smatch m;
    if (!std::regex_match(name, m, name_rx)) continue;
    const unsigned z1 = static_cast<unsigned>(std::stoi(m[1]));
    const unsigned z2 = static_cast<unsigned>(std::stoi(m[2]));
    const bool gas = m[3] == "g" || m[3] == "G";
    std::ifstream in(entry.path(), std::ios::binary);
    std::string line;
    double mass = 0.0;
    bool units_ok = false;
    std::vector<double> e;
    std::vector<double> se;
    std::vector<double> sn;
    while (std::getline(in, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      std::smatch mm;
      if (mass == 0.0 && line.find("Ion") != std::string::npos &&
          std::regex_search(line, mm, mass_rx)) {
        mass = std::stod(mm[1]);
      }
      if (line.find("Stopping Units") != std::string::npos) {
        units_ok = line.find("eV/(1E15 atoms/cm2)") != std::string::npos;
      }
      auto tok = split_ws(line);
      if (tok.size() == 3) {
        try {
          const double a = parse_double(tok[0]);
          const double b = parse_double(tok[1]);
          const double c = parse_double(tok[2]);
          e.push_back(a);
          se.push_back(b);
          sn.push_back(c);
        } catch (const std::runtime_error &) {
        }
      }
    }
    if (!units_ok) {
      throw std::runtime_error(entry.path().string() +
                               ": SR output must use stopping units eV/(1E15 atoms/cm2)");
    }
    if (mass <= 0.0 || e.size() < 2) {
      throw std::runtime_error("Cannot parse SR output " + entry.path().string());
    }
    db.add(make_table(z1, mass, z2, gas, e, se, sn));
  }
}

} // namespace

namespace {

//! The line [pos, end of line) of \p text without its '\r', and the start of
//! the next line (text.size() at the end).
std::string_view line_at(const std::string &text, std::size_t pos, std::size_t &next) {
  const std::size_t eol = text.find('\n', pos);
  const std::size_t end = eol == std::string::npos ? text.size() : eol;
  next = eol == std::string::npos ? text.size() : eol + 1;
  std::string_view line(text.data() + pos, end - pos);
  if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
  return line;
}

} // namespace

void SrimStoppingDatabase::load_packed(const std::filesystem::path &path) {
  {
    // One read of the whole file (a stringstream copy goes character by
    // character through the stream buffer and takes ~0.1 s for 42 MB).
    std::ifstream in(path, std::ios::binary | std::ios::ate);
    if (!in) throw std::runtime_error("Cannot open stopping table " + path.string());
    const std::streamoff size = in.tellg();
    if (size < 0) throw std::runtime_error("Cannot read stopping table " + path.string());
    text_.resize(static_cast<std::size_t>(size));
    in.seekg(0);
    if (size > 0 && !in.read(text_.data(), size)) {
      throw std::runtime_error("Cannot read stopping table " + path.string());
    }
  }
  bool header_ok = false;
  std::size_t line_number = 0;
  std::size_t pos = 0;
  std::lock_guard<std::mutex> lock(mutex_);
  while (pos < text_.size()) {
    std::size_t next = 0;
    const std::string_view line = line_at(text_, pos, next);
    ++line_number;
    pos = next;
    if (line.empty() || line[0] == '#') continue;
    const auto tok = split_ws(line);
    if (tok.empty()) continue;
    if (tok[0] == "amcset-stopping-table") {
      if (tok.size() < 2 || tok[1] != "1") {
        throw std::runtime_error("Unsupported stopping table version in " + path.string());
      }
      header_ok = true;
    } else if (tok[0] == "source") {
      source_ = std::string(line.substr(std::min(line.size(), std::size_t{7})));
    } else if (tok[0] == "energy_unit") {
      if (tok.size() < 2 || tok[1] != "keV") {
        throw std::runtime_error("Stopping table energy unit must be keV");
      }
    } else if (tok[0] == "stopping_unit") {
      if (tok.size() < 2 || tok[1] != "eV/(1e15") {
        throw std::runtime_error("Stopping table unit must be eV/(1e15 atoms/cm2)");
      }
    } else if (tok[0] == "energies") {
      std::vector<double> energies;
      for (std::size_t i = 2; i < tok.size(); ++i) energies.push_back(parse_double(tok[i]));
      grids_.push_back(std::move(energies));
    } else if (tok[0] == "table") {
      if (!header_ok) throw std::runtime_error("Missing header in " + path.string());
      if (tok.size() < 5) throw std::runtime_error("Bad table line " + std::to_string(line_number));
      if (grids_.empty()) {
        throw std::runtime_error("Stopping table before its energies line (line " +
                                 std::to_string(line_number) + " of " + path.string() + ")");
      }
      Unparsed u;
      // Z must be checked before the integer conversion (NaN or 1e30 has no
      // unsigned value).
      const auto atomic_number = [&](std::string_view t) {
        const double z = parse_double(t);
        if (!(z >= 1.0 && z <= static_cast<double>(Element::max_z)) || z != std::floor(z)) {
          throw std::runtime_error("Bad table line " + std::to_string(line_number) + " of " + path.string() +
                                   ": Z must be a whole number from 1 to " + std::to_string(Element::max_z));
        }
        return static_cast<unsigned>(z);
      };
      const unsigned z1 = atomic_number(tok[1]);
      u.m1_amu = parse_double(tok[2]);
      const unsigned z2 = atomic_number(tok[3]);
      const bool gas = tok[4] == "g";
      u.data = pos;
      u.grid = grids_.size() - 1;
      u.line = line_number;
      // The Se and Sn lines follow; they are parsed when the table is needed.
      for (int k = 0; k < 2; ++k) {
        if (pos >= text_.size()) {
          throw std::runtime_error("Truncated stopping table file " + path.string());
        }
        line_at(text_, pos, next);
        pos = next;
        ++line_number;
      }
      const Key key{z1, z2, gas};
      tables_.erase(key); // a later table replaces an earlier one
      unparsed_[key] = u;
    }
  }
  if (!header_ok) {
    throw std::runtime_error(path.string() + " is not an AMCSET stopping table");
  }
}

const StoppingTable *SrimStoppingDatabase::find_locked(const Key &key) const {
  if (const auto it = tables_.find(key); it != tables_.end()) return &it->second;
  const auto u = unparsed_.find(key);
  if (u == unparsed_.end()) return nullptr;
  const std::vector<double> &energies = grids_[u->second.grid];
  std::size_t next = 0;
  const std::string_view se_line = line_at(text_, u->second.data, next);
  const std::string_view sn_line = line_at(text_, next, next);
  const auto se_tok = split_ws(se_line);
  const auto sn_tok = split_ws(sn_line);
  if (se_tok.empty() || se_tok[0] != "Se" || sn_tok.empty() || sn_tok[0] != "Sn" ||
      se_tok.size() != energies.size() + 1 || sn_tok.size() != energies.size() + 1) {
    throw std::runtime_error("Malformed table near line " + std::to_string(u->second.line + 2) +
                             " of " + path_.string());
  }
  std::vector<double> se;
  std::vector<double> sn;
  se.reserve(energies.size());
  sn.reserve(energies.size());
  for (std::size_t i = 1; i < se_tok.size(); ++i) {
    se.push_back(parse_double(se_tok[i]));
    sn.push_back(parse_double(sn_tok[i]));
  }
  StoppingTable t = make_table(std::get<0>(key), u->second.m1_amu, std::get<1>(key),
                               std::get<2>(key), energies, se, sn);
  if (const auto problem = table_problem(t); !problem.empty()) {
    throw std::runtime_error("Malformed stopping table Z1=" + std::to_string(std::get<0>(key)) +
                             " Z2=" + std::to_string(std::get<1>(key)) + " near line " +
                             std::to_string(u->second.line) + " of " + path_.string() + ": " + problem);
  }
  unparsed_.erase(u);
  return &(tables_[key] = std::move(t));
}

std::shared_ptr<const SrimStoppingDatabase>
SrimStoppingDatabase::load(const std::filesystem::path &path) {
  // Loaded databases are kept (the bundled file has 16928 tables): a file is
  // read again only when its size or time changes.
  static std::mutex cache_mutex;
  static std::map<std::string, std::tuple<std::uintmax_t, std::filesystem::file_time_type,
                                          std::shared_ptr<const SrimStoppingDatabase>>>
      cache;
  std::error_code ec;
  const bool is_file = std::filesystem::is_regular_file(path, ec);
  std::string key;
  std::uintmax_t size = 0;
  std::filesystem::file_time_type time{};
  if (is_file) {
    key = std::filesystem::weakly_canonical(path, ec).string();
    if (key.empty()) key = path.string();
    size = std::filesystem::file_size(path, ec);
    time = std::filesystem::last_write_time(path, ec);
    std::lock_guard<std::mutex> lock(cache_mutex);
    const auto it = cache.find(key);
    if (it != cache.end() && std::get<0>(it->second) == size && std::get<1>(it->second) == time) {
      return std::get<2>(it->second);
    }
  }
  auto db = load_uncached(path);
  if (is_file) {
    std::lock_guard<std::mutex> lock(cache_mutex);
    cache[key] = {size, time, db};
  }
  return db;
}

std::shared_ptr<const SrimStoppingDatabase>
SrimStoppingDatabase::load_uncached(const std::filesystem::path &path) {
  auto db = std::make_shared<SrimStoppingDatabase>();
  db->path_ = path;
  if (std::filesystem::is_directory(path)) {
    load_sr_directory(path, *db);
    db->source_ = "SRIM SR Module output directory " + path.string();
  } else if (std::filesystem::is_regular_file(path)) {
    db->load_packed(path);
  } else {
    throw std::runtime_error("Stopping table not found: " + path.string());
  }
  if (db->size() == 0) {
    throw std::runtime_error("No stopping tables found in " + path.string());
  }
  return db;
}

// ---------------------------------------------------------------------------
// Layer stopping
// ---------------------------------------------------------------------------

std::string to_string(ElectronicStoppingModel model) {
  switch (model) {
  case ElectronicStoppingModel::srim_tables: return "srim";
  case ElectronicStoppingModel::builtin: return "builtin";
  case ElectronicStoppingModel::bethe_bloch: return "bethe-bloch";
  case ElectronicStoppingModel::lindhard_scharff: return "lss";
  }
  return "builtin";
}

ElectronicStoppingModel parse_electronic_stopping_model(const std::string &s) {
  std::string t;
  for (char c : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  if (t == "srim" || t == "srim_tables" || t == "srim-tables" || t == "srim2013") {
    return ElectronicStoppingModel::srim_tables;
  }
  if (t == "builtin" || t == "built-in" || t == "default") return ElectronicStoppingModel::builtin;
  if (t == "bethe-bloch" || t == "bethe_bloch" || t == "bethe" || t == "bb") {
    return ElectronicStoppingModel::bethe_bloch;
  }
  if (t == "lss" || t == "lindhard-scharff" || t == "lindhard_scharff" || t == "ls") {
    return ElectronicStoppingModel::lindhard_scharff;
  }
  throw std::invalid_argument("Unknown electronic stopping model '" + s +
                              "' (expected srim, builtin, bethe-bloch or lss)");
}

namespace {

stopping_cross_section_quantity
component_stopping(ElectronicStoppingModel model,
                   const SrimStoppingDatabase *database, unsigned z1,
                   mass_quantity m1, unsigned z2, bool gas,
                   energy_quantity energy,
                   number_density_quantity electron_density,
                   const BetheBlochOptions &bethe) {
  switch (model) {
  case ElectronicStoppingModel::lindhard_scharff:
    return lindhard_scharff_stopping(z1, m1, z2, energy);
  case ElectronicStoppingModel::bethe_bloch:
    return bethe_bloch_stopping(z1, m1, z2,
                                Element::by_z(z2).mean_excitation_energy(),
                                energy, electron_density, bethe);
  case ElectronicStoppingModel::builtin:
    return builtin_stopping(z1, m1, z2, energy, electron_density, bethe);
  case ElectronicStoppingModel::srim_tables: {
    if (database == nullptr) {
      throw std::invalid_argument("SRIM stopping model selected but no stopping table loaded");
    }
    const StoppingTable *table = database->find(z1, z2, gas);
    if (table == nullptr) {
      throw std::runtime_error("No SRIM stopping table for Z1=" + std::to_string(z1) +
                               " in Z2=" + std::to_string(z2));
    }
    // Above the table: continue with Bethe-Bloch, scaled to be continuous.
    const energy_quantity table_top = table->max_energy() * (m1 / table->m1).value();
    if (energy > table_top) {
      const auto i_value = Element::by_z(z2).mean_excitation_energy();
      const auto s_top = table->electronic(table_top, m1);
      const auto bb_top =
          bethe_bloch_stopping(z1, m1, z2, i_value, table_top, electron_density, bethe);
      const auto bb = bethe_bloch_stopping(z1, m1, z2, i_value, energy, electron_density, bethe);
      if (bb_top > 0.0 * srim_stopping_unit && bb > 0.0 * srim_stopping_unit) {
        return s_top * (bb / bb_top).value();
      }
      return s_top;
    }
    return table->electronic(energy, m1);
  }
  }
  return 0.0 * srim_stopping_unit;
}

} // namespace

stopping_cross_section_quantity
layer_stopping(ElectronicStoppingModel model,
               const SrimStoppingDatabase *database, unsigned z1,
               mass_quantity m1, const std::vector<StoppingComponent> &layer,
               bool gas, double compound_correction, energy_quantity energy,
               number_density_quantity atomic_density,
               const BetheBlochOptions &bethe) {
  double mean_z = 0.0;
  for (const auto &c : layer) mean_z += c.atom_fraction * static_cast<double>(c.z);
  const number_density_quantity electron_density = atomic_density * mean_z;
  stopping_cross_section_quantity s = 0.0 * srim_stopping_unit;
  for (const auto &c : layer) {
    s += c.atom_fraction * component_stopping(model, database, z1, m1, c.z, gas,
                                              energy, electron_density, bethe);
  }
  return compound_correction * s;
}

// ---------------------------------------------------------------------------
// StoppingCurve
// ---------------------------------------------------------------------------

StoppingCurve::StoppingCurve(ElectronicStoppingModel model,
                             const SrimStoppingDatabase *database, unsigned z1,
                             mass_quantity m1,
                             const std::vector<StoppingComponent> &layer,
                             bool gas, double compound_correction,
                             number_density_quantity atomic_density,
                             const BetheBlochOptions &bethe,
                             energy_quantity min_energy,
                             energy_quantity max_energy,
                             unsigned points_per_decade)
    : density_(atomic_density) {
  if (!(min_energy > 0.0 * electron_volt) || !(max_energy > min_energy) ||
      points_per_decade == 0) {
    throw std::invalid_argument("Invalid StoppingCurve grid");
  }
  log_e0_ = std::log(in_units(min_energy, electron_volt));
  const double step = ln10 / static_cast<double>(points_per_decade);
  inv_step_ = 1.0 / step;
  const double log_e1 = std::log(in_units(max_energy, electron_volt));
  const auto n = static_cast<std::size_t>(std::ceil((log_e1 - log_e0_) * inv_step_)) + 1;
  log_s_.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const energy_quantity e =
        std::exp(log_e0_ + step * static_cast<double>(i)) * electron_volt;
    const stopping_power_quantity sp =
        atomic_density * layer_stopping(model, database, z1, m1, layer, gas,
                                        compound_correction, e, atomic_density, bethe);
    log_s_[i] = std::log(std::max(in_units(sp, electron_volt_per_angstrom), 1e-300));
  }
}

stopping_power_quantity StoppingCurve::stopping_power(energy_quantity energy) const {
  if (log_s_.empty() || energy <= 0.0 * electron_volt) {
    return 0.0 * electron_volt_per_angstrom;
  }
  const double x = (std::log(in_units(energy, electron_volt)) - log_e0_) * inv_step_;
  double value = 0.0;
  if (x <= 0.0) {
    // Below the grid: velocity-proportional, S ~ E^(1/2).
    value = log_s_.front() + 0.5 * x / inv_step_;
  } else {
    const auto last = static_cast<double>(log_s_.size() - 1);
    if (x >= last) {
      // Above the grid: continue the slope of the last interval.
      value = log_s_.back() + (x - last) * (log_s_.back() - log_s_[log_s_.size() - 2]);
    } else {
      const auto i = static_cast<std::size_t>(x);
      const double t = x - static_cast<double>(i);
      value = log_s_[i] + t * (log_s_[i + 1] - log_s_[i]);
    }
  }
  return std::exp(value) * electron_volt_per_angstrom;
}

void StoppingCurve::enable_fast_lookup() {
  if (log_s_.empty()) return;
  // 1 meV (the transport's lower limit of the energy) to 2^37 eV = 137 GeV
  // (above the highest ion energy), 128 bins per factor of two.
  fast_ = fast::OctaveTable(-10, 37, [this](double e_eV) {
    return in_units(stopping_power(e_eV * electron_volt), electron_volt_per_angstrom);
  });
}

stopping_cross_section_quantity
StoppingCurve::cross_section(energy_quantity energy) const {
  return stopping_power(energy) / density_;
}

} // namespace physics
} // namespace amcset
