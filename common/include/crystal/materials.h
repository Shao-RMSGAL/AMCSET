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

/*! \file materials.h
 *  \brief Common crystals: structure, lattice constant and Debye temperature.
 *
 * The values are those of D. S. Gemmell, "Channeling and related effects in
 * the motion of charged particles through crystals", Rev. Mod. Phys. 46, 129
 * (1974), Appendix B: Table IV (the most commonly studied cubic elements:
 * Debye temperatures, lattice constants at room temperature) and Table V
 * (cubic diatomic compounds: lattice constants). Nothing else is built in.
 */

#pragma once

#include "amcset_units.h"
#include "crystal/crystal.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace crystal {

//! One row of Gemmell's Table IV or V.
struct ReferenceCrystal {
  std::string_view name;      //!< "Si", "GaAs", "UO2", ...
  std::string_view structure; //!< a builtin_structures() key
  double lattice_constant_A;  //!< at room temperature
  /*!
   * The elements of the structure's sublattices, in its order (one element
   * per sublattice: Z).
   */
  std::vector<unsigned> sublattices;
  //! Debye temperature (K) and rms vibration amplitude along one axis at
  //! 293 K (Angstrom) as tabulated; 0 for the compounds (Table V has none).
  double debye_temperature_K = 0.0;
  double u1_293K_A = 0.0;
  std::string_view note; //!< e.g. "above 120 C" for cubic BaTiO3
};

//! Gemmell's Tables IV and V (elements first, in the order of Z).
const std::vector<ReferenceCrystal> &reference_crystals();
//! The row called \p name (case-sensitive: "Si", "GaAs"), or null.
const ReferenceCrystal *find_reference_crystal(std::string_view name);
//! The tabulated Debye temperature of element \p z (Table IV), if any.
std::optional<common::temperature_quantity> reference_debye_temperature(unsigned z);

/*!
 * \brief The rms thermal displacement (one axis) the transport uses for
 * element \p z of \p crystal: the crystal's own (Crystal::vibrations); else,
 * in a crystal of one element that Table IV lists, the Debye model with the
 * tabulated Debye temperature (\p from_table is then set); else none (the
 * atoms stand still).
 */
std::optional<common::length_quantity> thermal_amplitude(const Crystal &crystal, unsigned z,
                                                         common::mass_quantity mass, bool *from_table = nullptr);

/*!
 * \brief Lindhard's characteristic angle of axial channeling,
 * psi_1 = sqrt(2 Z1 Z2 e^2 / (E d)), for atoms spaced \p spacing along the
 * row (Gemmell's Eq. A5: 0.307 sqrt(Z1 Z2 / (d[A] E[MeV])) degrees).
 * \throws std::invalid_argument for non-positive arguments.
 */
common::angle_quantity lindhard_angle(unsigned z1, unsigned z2, common::length_quantity spacing,
                                      common::energy_quantity energy);

/*!
 * \brief The crystal of a row: its cell at the tabulated lattice constant,
 * a (0 0 1) surface, 293 K, and for the elements the tabulated Debye
 * temperature.
 */
Crystal make_reference_crystal(const ReferenceCrystal &row);

} // namespace crystal
} // namespace amcset
