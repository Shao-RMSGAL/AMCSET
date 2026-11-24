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

/*! \file run_config.h
 *  \brief Front-end-neutral description of a calculation (JSON scenario
 *  files), quantity parsing with units, and engine construction.
 *
 * A RunConfig fully describes one calculation in either engine:
 *  - Mode::srim   -- the TRIM-compatible engine (trim::TrimSimulation);
 *  - Mode::native -- the original AMCSET engine (common::Simulation).
 *
 * It is read from / written to JSON scenario files shared by the CLI and the
 * GUI. Physical values are written as strings with units ("1 MeV",
 * "1.2 um", "7.87 g/cm3", "40 eV"); plain numbers are accepted for keys with
 * a unit suffix ("energy_keV": 1000). See docs/CLI.md for the full schema.
 */

#pragma once

#include "amcset_units.h"
#include "simulation_engine.h"
#include "trim/trim.h"

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace config {

using common::angle_quantity;
using common::energy_quantity;
using common::length_quantity;
using common::mass_density_quantity;
using common::mass_quantity;

// ---------------------------------------------------------------------------
// Quantities with units
// ---------------------------------------------------------------------------

/*!
 * \brief Parse "<number> <unit>" (the space is optional).
 *
 * Energy: eV, keV, MeV, GeV, J. Length: A (Angstrom), nm, um, mm, cm, m.
 * Mass: amu (u, Da), kg, g. Mass density: g/cm3 (g/cm^3, g/cc), kg/m3.
 * Angle: deg, rad. A bare number uses \p default_unit (if given).
 * \throws std::invalid_argument naming \p what on malformed input or a unit
 *         of the wrong dimension.
 */
energy_quantity parse_energy(std::string_view text, std::string_view what = "energy");
/*!
 * \brief A list of energies: "1 keV", "0.5,1,2 keV" (a number without a
 * unit takes the unit of the next entry that has one) or "500 eV, 1 keV".
 * \throws std::invalid_argument naming \p what
 */
std::vector<energy_quantity> parse_energy_list(const std::string &text, const std::string &what);
length_quantity parse_length(std::string_view text, std::string_view what = "length");
mass_quantity parse_mass(std::string_view text, std::string_view what = "mass");
mass_density_quantity parse_density(std::string_view text, std::string_view what = "density");
angle_quantity parse_angle(std::string_view text, std::string_view what = "angle");
//! Temperature: K, or C (degrees Celsius; converted to kelvin).
common::temperature_quantity parse_temperature(std::string_view text, std::string_view what = "temperature");

//! Readable strings with a sensible unit ("1.5 MeV", "1200 A", "7.87 g/cm3").
std::string format_energy(energy_quantity e);
std::string format_length(length_quantity l);
std::string format_mass(mass_quantity m);
std::string format_density(mass_density_quantity d);
std::string format_angle(angle_quantity a);
std::string format_temperature(common::temperature_quantity t);

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

//! Which engine runs the calculation.
enum class Mode {
  srim,  //!< TRIM-compatible SRIM mode (default)
  native //!< Original AMCSET binary-collision engine (ions and electrons)
};

std::string to_string(Mode m);
Mode parse_mode(std::string_view s);

/*!
 * \brief Settings that only the native AMCSET engine uses. Names follow
 * common::Simulation::Settings; see its documentation.
 */
struct NativeOptions {
  bool electron = false;            //!< Electron bombardment instead of ions
  unsigned ion_mass_number = 0;     //!< Isotope of the ion (0 = from ion mass)
  bool damage_cascade = true;       //!< Follow recoils
  std::string screening_model = "thomas-fermi"; //!< thomas-fermi, quickbasic, zbl, moliere, custom
  std::string electronic_stopping_model = "zbl"; //!< zbl, lss, custom, bethe-bloch, builtin, srim
  double custom_electronic_coefficient = 1.0;
  double custom_electronic_exponent = 100.0;
  double custom_recoil_electronic_coefficient = 1.0;
  double custom_recoil_electronic_exponent = 100.0;
  double custom_nuclear_screening_scale = 1.0;
  energy_quantity ion_stop_energy = 40.0 * common::electron_volt;
  //! Default displacement energy (per-element Ed of the layers take precedence).
  energy_quantity displacement_energy = 40.0 * common::electron_volt;
  energy_quantity electron_stop_energy = 40.0 * common::electron_volt;
  std::size_t angle_divisors = 1000;
  std::size_t flying_distances = 1000;
  length_quantity range = 0.0 * common::angstrom; //!< Histogram extent (0 = target depth)
  std::size_t bins = 100;
  bool per_bombardment_output = false;
  std::string electron_screening_table = "common/assets/electron_screening_parameters.csv";
  std::string mott_table = "common/assets/mott_scattering_parameters.csv";
};

/*!
 * \brief One complete calculation.
 *
 * The ion, the target (with per-element Ed, Eb, Es), counts, seed, threads,
 * stopping tables and output options live in \c trim (a trim::TrimConfig)
 * for both modes; \c native holds what only the native engine needs.
 */
struct RunConfig {
  Mode mode = Mode::srim;
  //! SRIM mode defaults to TRIM's stopping: the SRIM-2013 tables that come
  //! with AMCSET (stopping model srim, empty table path).
  trim::TrimConfig trim = trim_mode_defaults();
  NativeOptions native;

  //! A TrimConfig with the defaults of a SRIM-mode scenario.
  static trim::TrimConfig trim_mode_defaults() {
    trim::TrimConfig t;
    t.stopping_model = physics::ElectronicStoppingModel::srim_tables;
    return t;
  }
};

/*!
 * \brief Parse a JSON scenario. Relative paths inside it (stopping table,
 * SRIM data, output directory) are resolved against \p base_dir.
 * \throws std::invalid_argument with a message naming the offending key.
 */
RunConfig parse_config(std::string_view json_text,
                       const std::filesystem::path &base_dir = {});
//! parse_config() on a file (base_dir = the file's directory).
RunConfig load_config(const std::filesystem::path &path);
//! Serialize to JSON (canonical form with unit strings).
std::string to_json(const RunConfig &config);
//! Write to_json() to a file.
void save_config(const RunConfig &config, const std::filesystem::path &path);

//! A layer's "crystal" object (docs/CLI.md), parsed.
struct CrystalDescription {
  crystal::Crystal crystal;
  /*!
   * A built-in structure given without a lattice constant: the cell is built
   * with a = 1 Angstrom and takes its size from the layer's density
   * (crystal::UnitCell::scale_to_density).
   */
  bool lattice_from_density = false;
};

/*!
 * \brief Parse a "crystal" object: a built-in structure, a cell file (CIF,
 * POSCAR, LAMMPS data, extended XYZ; relative to \p base_dir) or a lattice
 * with sites, plus orientation, temperature and thermal vibrations.
 * \throws std::invalid_argument naming the key (prefixed with \p where),
 *         std::runtime_error if a cell file cannot be read.
 */
CrystalDescription parse_crystal(std::string_view json_text, const std::filesystem::path &base_dir = {},
                                 const std::string &where = "crystal");
//! The canonical "crystal" object (lattice vectors and sites written out).
std::string crystal_to_json(const crystal::Crystal &crystal);

/*!
 * \brief Apply one "key=value" override, e.g. "ion.energy=2 MeV",
 * "ions=5000", "damage=full", "layers.0.elements.0.Ed=25 eV".
 * \throws std::invalid_argument for unknown keys or bad values.
 */
void apply_override(RunConfig &config, std::string_view assignment);

/*!
 * \brief A scenario template: Fe ions into Fe in SRIM mode, or the native
 * equivalent.
 */
RunConfig example_config(Mode mode);

/*!
 * \brief Build the engine described by \p config (validates it first).
 * \throws std::invalid_argument, std::runtime_error
 */
std::unique_ptr<common::SimulationEngine> make_engine(const RunConfig &config);

} // namespace config
} // namespace amcset
