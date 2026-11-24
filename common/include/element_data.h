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

/*! \file element_data.h
 *  \brief Per-element reference data (Z = 1..92): masses, densities, surface
 *  binding energies and mean excitation energies.
 *
 * The raw table is generated from common/assets/element_data.csv (see that
 * file for sources). The numbers in ElementRecord are stored exactly as
 * tabulated, in the units named by each field; always read them through the
 * typed accessors of Element, which attach Boost.Units dimensions.
 */

#pragma once

#include "amcset_units.h"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace amcset {
namespace common {

/*!
 * \brief One row of the element table, as tabulated (plain numbers).
 *
 * The unit of each field is part of its name. Use Element for typed access.
 */
struct ElementRecord {
  unsigned z;                       //!< Atomic number
  const char *symbol;               //!< Chemical symbol, e.g. "Fe"
  const char *name;                 //!< English name, e.g. "Iron"
  unsigned mai_mass_number;         //!< Mass number of the most abundant isotope
  double mai_mass_amu;              //!< Mass of the most abundant isotope (amu)
  double natural_weight_amu;        //!< Standard atomic weight (amu)
  double density_g_cm3;             //!< Solid (or condensed) density (g/cm^3)
  double atomic_density_cm3;        //!< Solid atomic density (atoms/cm^3)
  double fermi_velocity_v0;         //!< Fermi velocity (units of the Bohr velocity)
  double sublimation_energy_eV;     //!< Heat of sublimation (eV)
  double gas_density_g_cm3;         //!< Gas density at STP (g/cm^3)
  double gas_atomic_density_cm3;    //!< Gas atomic density at STP (atoms/cm^3)
  double mean_excitation_energy_eV; //!< Mean excitation energy I (eV)
};

/*!
 * \brief Typed view of one element of the table.
 */
class Element {
public:
  static constexpr unsigned max_z = 92; //!< Highest tabulated atomic number

  //! Element with atomic number \p z. \throws std::out_of_range
  static const Element &by_z(unsigned z);
  //! Element by (case-insensitive) symbol or English name; nullopt if unknown.
  static std::optional<unsigned> find_z(std::string_view symbol_or_name);
  //! Element by symbol or name. \throws std::invalid_argument if unknown.
  static const Element &by_symbol(std::string_view symbol_or_name);
  //! All tabulated elements (index = Z - 1).
  static const std::array<Element, max_z> &all();

  unsigned z() const noexcept { return record_->z; }
  std::string_view symbol() const noexcept { return record_->symbol; }
  std::string_view name() const noexcept { return record_->name; }
  unsigned most_abundant_mass_number() const noexcept {
    return record_->mai_mass_number;
  }
  //! Mass of the most abundant isotope (the ion mass TRIM uses by default).
  mass_quantity most_abundant_isotope_mass() const noexcept {
    return record_->mai_mass_amu * atomic_mass_unit;
  }
  //! Standard atomic weight (the target-atom mass TRIM uses by default).
  mass_quantity natural_mass() const noexcept {
    return record_->natural_weight_amu * atomic_mass_unit;
  }
  //! Density of the element in its solid (condensed) form.
  mass_density_quantity solid_density() const noexcept {
    return record_->density_g_cm3 * gram_per_cubic_centimeter;
  }
  //! Atomic density of the element in its solid (condensed) form.
  number_density_quantity solid_atomic_density() const noexcept {
    return record_->atomic_density_cm3 * per_cubic_centimeter;
  }
  //! Density of the gaseous element at standard conditions.
  mass_density_quantity gas_density() const noexcept {
    return record_->gas_density_g_cm3 * gram_per_cubic_centimeter;
  }
  //! Whether the element is a gas at standard conditions.
  bool is_gas_at_stp() const noexcept {
    return record_->gas_density_g_cm3 < 0.01;
  }
  //! Fermi velocity of the solid.
  velocity_quantity fermi_velocity() const noexcept {
    return record_->fermi_velocity_v0 * bohr_velocity;
  }
  /*!
   * \brief Heat of sublimation, the customary default surface binding energy.
   *
   * Zero for elements that are gases at standard conditions.
   */
  energy_quantity sublimation_energy() const noexcept {
    return record_->sublimation_energy_eV * electron_volt;
  }
  //! Mean excitation energy I of the Bethe-Bloch formula.
  energy_quantity mean_excitation_energy() const noexcept {
    return record_->mean_excitation_energy_eV * electron_volt;
  }
  //! The raw table row.
  const ElementRecord &record() const noexcept { return *record_; }

  constexpr explicit Element(const ElementRecord *record) : record_(record) {}

private:
  const ElementRecord *record_;
};

/*!
 * \brief Default displacement energy used when none is given.
 *
 * 25 eV is TRIM's generic default for metals; ASTM E521 recommends 40 eV for
 * Fe and steels. Always set Ed explicitly for damage calculations.
 */
inline constexpr energy_quantity default_displacement_energy =
    25.0 * electron_volt;
//! Default lattice binding energy (TRIM's generic default).
inline constexpr energy_quantity default_lattice_binding_energy =
    3.0 * electron_volt;

} // namespace common
} // namespace amcset
