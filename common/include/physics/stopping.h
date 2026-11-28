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

/*! \file stopping.h
 *  \brief Electronic stopping: analytic models (Lindhard-Scharff, relativistic
 *  Bethe-Bloch), tabulated SRIM stopping, Bragg's rule for compounds, and the
 *  fast per-(projectile, layer) StoppingCurve used in the transport loops.
 *
 * All stopping cross sections S = -(1/N) dE/dx are
 * stopping_cross_section_quantity (energy * area per atom); stopping powers
 * N*S are stopping_power_quantity (energy per length).
 */

#pragma once

#include "amcset_units.h"
#include "physics/fast_math.h"

#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <vector>

namespace amcset {
namespace physics {

using common::energy_quantity;
using common::energy_squared_per_length_quantity;
using common::length_quantity;
using common::mass_quantity;
using common::number_density_quantity;
using common::stopping_cross_section_quantity;
using common::stopping_power_quantity;
using common::velocity_quantity;

// ---------------------------------------------------------------------------
// Kinematics helpers
// ---------------------------------------------------------------------------

//! Relativistic speed of a particle of rest mass \p mass and kinetic energy \p energy.
velocity_quantity relativistic_speed(energy_quantity energy, mass_quantity mass);
//! beta = v / c for kinetic energy \p energy.
double relativistic_beta(energy_quantity energy, mass_quantity mass);
//! Lorentz factor gamma = 1 + E / (M c^2).
double lorentz_gamma(energy_quantity energy, mass_quantity mass);

// ---------------------------------------------------------------------------
// Analytic electronic stopping models
// ---------------------------------------------------------------------------

/*!
 * \brief Lindhard-Scharff low-velocity electronic stopping cross section.
 *
 * S_e = xi_e 8 pi e^2 a_0 Z1 Z2 / (Z1^(2/3) + Z2^(2/3))^(3/2) * v / v_0 with
 * xi_e = Z1^(1/6) (Lindhard & Scharff, Phys. Rev. 124, 128 (1961)).
 * Valid for v < v_0 Z1^(2/3).
 */
stopping_cross_section_quantity lindhard_scharff_stopping(unsigned z1,
                                                          mass_quantity m1,
                                                          unsigned z2,
                                                          energy_quantity energy);

//! Options of the Bethe-Bloch formula.
struct BetheBlochOptions {
  //! Use the Barkas-Pierce-Blann effective charge
  //! z* = Z1 (1 - exp(-125 beta Z1^(-2/3))) instead of the bare charge.
  bool effective_charge = true;
  //! Apply the Sternheimer-Peierls density-effect correction (needs the
  //! target electron density). Matters above ~1 GeV/u.
  bool density_effect = true;
};

/*!
 * \brief Corrections of the SRIM-2013 stopping tables (not TRIM: off by
 * default, so that SRIM mode reproduces TRIM.exe). docs/VALIDATION.md says
 * what each corrects and the evidence.
 */
struct StoppingCorrections {
  /*!
   * \brief Ni projectiles: the electronic stopping is the geometric mean of
   * Co's and Cu's at the same velocity, in the same target. SRIM-2013's Ni
   * stopping is up to 40 % below that of both neighbours between about 0.02
   * and 0.3 MeV/u, in every target (tools/validation/srim_table_scan.py), and
   * the measured ranges of Ni ions in Fe and Mo are 1.3-1.5 times shorter
   * than SRIM-2013's (Gurbich et al., Appl. Radiat. Isot. 215, 111590
   * (2025)). Ni recoils take it too.
   */
  bool nickel_ions = false;

  bool any() const noexcept { return nickel_ions; }
  bool operator==(const StoppingCorrections &) const = default;
};

/*!
 * \brief Relativistic Bethe-Bloch electronic stopping cross section of one
 * target atom.
 *
 * S = 4 pi (e^2/(4 pi eps0))^2 z*^2 Z2 / (m_e c^2 beta^2)
 *     [ 1/2 ln(2 m_e c^2 beta^2 gamma^2 W_max / I^2) - beta^2 - delta/2 ]
 * with W_max the maximum energy transfer to a free electron (PDG, Rev. Part.
 * Phys. "Passage of particles through matter"). Shell corrections are not
 * included, so the formula is only used where v >> Z2 v_0 (see
 * builtin_stopping and StoppingCurve). Returns zero where the logarithm
 * becomes negative.
 *
 * \param mean_excitation_energy The target's I value.
 * \param electron_density Target electron density (for the density effect;
 *        ignored when \p options.density_effect is false).
 */
stopping_cross_section_quantity
bethe_bloch_stopping(unsigned z1, mass_quantity m1, unsigned z2,
                     energy_quantity mean_excitation_energy,
                     energy_quantity energy,
                     number_density_quantity electron_density,
                     const BetheBlochOptions &options = {});

/*!
 * \brief Built-in stopping without data tables: Lindhard-Scharff and a
 * Bethe-Bloch form with a low-velocity term, joined by the generalized mean
 * S = (S_LS^-4 + S_BB^-4)^(-1/4).
 *
 * Compared with SRIM-2013 it deviates by ~25 % (rms) below 100 keV/u and
 * by < 6 % above 1 MeV/u. For TRIM-equivalent results use the SRIM
 * stopping tables (SrimStoppingDatabase).
 */
stopping_cross_section_quantity
builtin_stopping(unsigned z1, mass_quantity m1, unsigned z2,
                 energy_quantity energy,
                 number_density_quantity electron_density,
                 const BetheBlochOptions &options = {});

/*!
 * \brief TRIM's electronic-energy straggling per unit path length,
 * Omega^2 / L (Ziegler, Biersack & Ziegler, "SRIM - The Stopping and Range
 * of Ions in Matter", Eq. 7-47).
 *
 * The minimum of relativistic Bohr straggling,
 * 4 pi Z1^2 e^4 Z2 N (b - 1/2)/(b - 1) with b = c^2/v^2, and the low-energy
 * expression (Z1 + Z2)^(8/3) E N / (3.2 M1). The book gives the latter with
 * E in keV, N in atoms/Angstrom^3 and M1 in amu and the result in
 * eV^2/Angstrom; those units are attached here explicitly.
 */
energy_squared_per_length_quantity
trim_electronic_straggling(unsigned z1, mass_quantity m1, unsigned z2,
                           energy_quantity energy,
                           number_density_quantity atomic_density);

// ---------------------------------------------------------------------------
// Tabulated (SRIM) stopping
// ---------------------------------------------------------------------------

/*!
 * \brief A stopping table for one projectile species in one element.
 *
 * Electronic and nuclear stopping cross sections on an increasing energy
 * grid, interpolated log-log (exactly as TRIM interpolates its own
 * tables). Electronic stopping depends on the projectile velocity, so for a
 * projectile mass \p m different from the table mass the lookup energy is
 * scaled by m_table / m.
 */
struct StoppingTable {
  unsigned z1 = 0;                   //!< Projectile atomic number
  mass_quantity m1;                  //!< Projectile mass of the table
  unsigned z2 = 0;                   //!< Target atomic number
  bool gas = false;                  //!< Gas-phase target
  std::vector<double> log_energy;    //!< ln(E / keV), increasing
  std::vector<double> log_electronic; //!< ln(S_e / srim_stopping_unit)
  std::vector<double> log_nuclear;    //!< ln(S_n / srim_stopping_unit)

  //! Lowest and highest tabulated energies (for the table mass).
  energy_quantity min_energy() const;
  energy_quantity max_energy() const;

  /*!
   * \brief Electronic stopping cross section of a projectile of mass \p mass
   * at kinetic energy \p energy.
   *
   * Below the table, S_e is extrapolated proportional to the velocity; above
   * it, the last value is returned (callers that need high energies use
   * StoppingCurve, which continues with Bethe-Bloch).
   */
  stopping_cross_section_quantity electronic(energy_quantity energy,
                                             mass_quantity mass) const;
  //! Tabulated nuclear stopping cross section (table mass only).
  stopping_cross_section_quantity nuclear(energy_quantity energy) const;
};

/*!
 * \brief A collection of StoppingTables keyed by (Z1, Z2, phase).
 *
 * Load from either
 *  - an AMCSET stopping-table file ("amcset-stopping-table 1" text format,
 *    see tools/srim/pack_stopping.py), or
 *  - a directory of raw SRIM SR-Module output files named
 *    sr_<Z1>_<Z2>_<s|g>.txt (three-digit Z, s = solid, g = gas).
 *
 * AMCSET comes with the SRIM-2013 tables of every ion in every element
 * (data/srim2013, found with data::bundled_srim_stopping_table()); the
 * scripts in tools/srim/ regenerate them from a SRIM-2013 installation.
 *
 * A packed file is indexed when it is loaded and each table is parsed the
 * first time find() asks for it (a calculation needs a few of the 16928
 * tables of the bundled file), so a malformed table is reported by find().
 * The class is thread-safe.
 */
class SrimStoppingDatabase {
public:
  SrimStoppingDatabase() = default;
  SrimStoppingDatabase(const SrimStoppingDatabase &) = delete;
  SrimStoppingDatabase &operator=(const SrimStoppingDatabase &) = delete;

  //! Load a database file or directory (a file once per process while it
  //! is unchanged). \throws std::runtime_error
  static std::shared_ptr<const SrimStoppingDatabase>
  load(const std::filesystem::path &path);

  /*!
   * \brief Table for projectile \p z1 in element \p z2.
   *
   * For \p gas targets without a gas-phase table the solid table is used.
   * \returns nullptr if the pair is not tabulated. \throws
   * std::runtime_error if the table in the file is malformed.
   */
  const StoppingTable *find(unsigned z1, unsigned z2, bool gas) const;

  //! Source description recorded in the file.
  const std::string &source() const noexcept { return source_; }
  //! Path the database was loaded from.
  const std::filesystem::path &path() const noexcept { return path_; }
  //! Number of tables.
  std::size_t size() const;

  //! Add or replace a table (used by the loaders and by tests).
  void add(StoppingTable table);

  //! Read a packed table file (indexing only; tables parse on demand).
  //! \throws std::runtime_error
  void load_packed(const std::filesystem::path &path);

private:
  using Key = std::tuple<unsigned, unsigned, bool>;
  //! A table of the packed file that is not parsed yet.
  struct Unparsed {
    double m1_amu = 0.0;
    std::size_t data = 0; //!< Offset of its "Se" line in text_
    std::size_t grid = 0; //!< Index into grids_ (the energies line in force)
    std::size_t line = 0; //!< Line number of the "table" line (messages)
  };
  static std::shared_ptr<const SrimStoppingDatabase> load_uncached(const std::filesystem::path &path);
  //! The parsed table for \p key, parsing it if needed (caller holds mutex_).
  const StoppingTable *find_locked(const Key &key) const;

  mutable std::mutex mutex_;
  mutable std::map<Key, StoppingTable> tables_;
  mutable std::map<Key, Unparsed> unparsed_;
  std::string text_;                       //!< Packed file contents (lazy tables)
  std::vector<std::vector<double>> grids_; //!< Energy grids of the packed file (keV)
  std::string source_;
  std::filesystem::path path_;
};

// ---------------------------------------------------------------------------
// Layer stopping (Bragg's rule) and the fast evaluation curve
// ---------------------------------------------------------------------------

//! Electronic stopping model selection.
enum class ElectronicStoppingModel {
  srim_tables,    //!< SRIM-2013 tables (+ Bethe-Bloch above the tables)
  builtin,        //!< Lindhard-Scharff / Bethe-Bloch interpolation (no data)
  bethe_bloch,    //!< Bethe-Bloch only (high-energy particles)
  lindhard_scharff //!< Lindhard-Scharff only (low-energy heavy ions)
};

//! Name used in configuration files ("srim", "builtin", "bethe-bloch", "lss").
std::string to_string(ElectronicStoppingModel model);
//! Inverse of to_string(). \throws std::invalid_argument
ElectronicStoppingModel parse_electronic_stopping_model(const std::string &s);

//! One element of a layer, for the purpose of stopping calculations.
struct StoppingComponent {
  unsigned z = 0;             //!< Target atomic number
  double atom_fraction = 0.0; //!< Stoichiometric atom fraction (sum 1)
};

/*!
 * \brief Electronic stopping cross section of a layer (per target atom),
 * Bragg's additivity rule times the layer's compound correction factor.
 *
 * \param database Required for ElectronicStoppingModel::srim_tables.
 * \param gas Use gas-phase tables where available.
 * \param atomic_density Total atomic density of the layer (sets the electron
 *        density for the Bethe-Bloch density effect).
 * \param corrections Corrections of the SRIM tables (srim_tables only).
 */
stopping_cross_section_quantity
layer_stopping(ElectronicStoppingModel model,
               const SrimStoppingDatabase *database, unsigned z1,
               mass_quantity m1, const std::vector<StoppingComponent> &layer,
               bool gas, double compound_correction, energy_quantity energy,
               number_density_quantity atomic_density,
               const BetheBlochOptions &bethe = {},
               const StoppingCorrections &corrections = {});

/*!
 * \brief Fast electronic stopping power N*S_e(E) of one projectile species
 * in one layer, tabulated on a logarithmic energy grid and interpolated
 * log-log.
 *
 * Built once per (species, layer) before transport; evaluation is two table
 * lookups and an exp/log. Energies outside the grid are extrapolated
 * (velocity-proportional below, Bethe-Bloch shape above).
 *
 * In srim_tables mode, energies above the SRIM table are continued with
 * Bethe-Bloch scaled to be continuous at the table's upper end, which is how
 * AMCSET supports high-energy (relativistic) particles.
 */
class StoppingCurve {
public:
  StoppingCurve() = default;

  /*!
   * \brief Build the curve.
   * \param min_energy Lower grid limit (default 1 eV).
   * \param max_energy Upper grid limit (default 100 GeV).
   * \param points_per_decade Grid density (default 50).
   * \param corrections Corrections of the SRIM tables (layer_stopping()).
   */
  StoppingCurve(ElectronicStoppingModel model,
                const SrimStoppingDatabase *database, unsigned z1,
                mass_quantity m1, const std::vector<StoppingComponent> &layer,
                bool gas, double compound_correction,
                number_density_quantity atomic_density,
                const BetheBlochOptions &bethe = {},
                energy_quantity min_energy = 1.0 * common::electron_volt,
                energy_quantity max_energy = 100.0 * common::giga_electron_volt,
                unsigned points_per_decade = 50,
                const StoppingCorrections &corrections = {});

  //! Electronic stopping power N*S_e at \p energy.
  stopping_power_quantity stopping_power(energy_quantity energy) const;
  /*!
   * \brief Tabulate the curve for stopping_power_fast() (fast
   * approximations; 6016 nodes, 96 kB, about 0.3 ms).
   */
  void enable_fast_lookup();
  /*!
   * \brief stopping_power() from a table indexed by the bits of the energy
   * (fast::OctaveTable, 1 meV to 137 GeV): no log or exp; relative error
   * below 4e-5 (largest where a bin straddles a node of the curve). The exact interpolation if enable_fast_lookup() was not
   * called or the energy is outside the table.
   */
  stopping_power_quantity stopping_power_fast(energy_quantity energy) const {
    double y = 0.0;
    if (fast_.lookup(common::in_units(energy, common::electron_volt), y)) {
      return y * common::electron_volt_per_angstrom;
    }
    return stopping_power(energy);
  }
  //! Electronic stopping cross section S_e at \p energy.
  stopping_cross_section_quantity cross_section(energy_quantity energy) const;
  //! The layer's atomic density N.
  number_density_quantity atomic_density() const noexcept { return density_; }
  //! Whether the curve has been built.
  bool empty() const noexcept { return log_s_.empty(); }
  //! The interpolation table stopping_power() reads (the GPU backends copy
  //! it): ln(E_min / eV), 1 / step in ln E, and ln(N S_e / (eV/Angstrom)).
  double log_e0() const noexcept { return log_e0_; }
  double inv_step() const noexcept { return inv_step_; }
  const std::vector<double> &log_table() const noexcept { return log_s_; }

private:
  double log_e0_ = 0.0;    // ln(E_min / eV)
  double inv_step_ = 1.0;  // 1 / (grid step in ln E)
  std::vector<double> log_s_; // ln(N S_e / (eV/Angstrom))
  number_density_quantity density_;
  fast::OctaveTable fast_;    // N S_e in eV/Angstrom over 2^-10 .. 2^37 eV (fast approximations)
};

} // namespace physics
} // namespace amcset
