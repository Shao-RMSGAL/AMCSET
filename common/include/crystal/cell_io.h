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

/*! \file cell_io.h
 *  \brief Unit cells from and to CIF, VASP POSCAR, LAMMPS data and extended
 *  XYZ files.
 *
 * Readers:
 *  - CIF: the first data block; cell lengths and angles; the symmetry
 *    operations (`_symmetry_equiv_pos_as_xyz` or
 *    `_space_group_symop_operation_xyz`) applied to the atom sites
 *    (`_atom_site_fract_*` or `_atom_site_Cartn_*`, `_atom_site_type_symbol`
 *    or the label, `_atom_site_occupancy`). A CIF that names a space group
 *    other than P 1 without listing its operations is refused.
 *  - POSCAR/CONTCAR (VASP 4 and 5): scale (or volume), lattice, species,
 *    counts, Direct or Cartesian coordinates.
 *  - LAMMPS data: the box (with tilt factors), Masses, Atoms (atomic,
 *    charge, full, molecular, bond, angle styles); types become elements by
 *    a `# Symbol` comment in Masses, a given map, or the mass.
 *  - Extended XYZ: `Lattice="..."` and `Properties=...` (species and pos).
 * Positions closer than 0.01 Angstrom are one site (occupancies add up), so
 * the disordered sites of a CIF become one site with several occupants.
 *
 * Writers produce the same formats; only CIF can store partial occupancy
 * (and AMCSET's site kind, as `_amcset_site_kind`).
 */

#pragma once

#include "crystal/unit_cell.h"

#include <filesystem>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace crystal {

enum class CellFormat { cif, poscar, lammps_data, extxyz };

std::string to_string(CellFormat f);
//! "cif", "poscar"/"vasp", "lammps"/"lammps-data", "xyz"/"extxyz". \throws std::invalid_argument
CellFormat parse_cell_format(std::string_view s);
//! Format from the file name (.cif; POSCAR, CONTCAR, .vasp, .poscar; .data, .lmp, .lammps; .xyz, .extxyz).
std::optional<CellFormat> format_from_path(const std::filesystem::path &path);

//! What a file cannot say by itself.
struct CellReadOptions {
  std::vector<unsigned> elements;       //!< POSCAR without a species line: elements in the order of the counts
  std::map<int, unsigned> lammps_types; //!< LAMMPS atom type -> Z (wins over Masses)
  std::string lammps_atom_style;        //!< when the Atoms line has no "# style" comment
  double merge_distance_A = 0.01;       //!< positions closer than this are one site
};

/*!
 * \brief Read a cell file (format from \p format or the file name).
 * \throws std::runtime_error if the file cannot be read,
 *         std::invalid_argument for malformed content (with the line).
 */
UnitCell read_cell(const std::filesystem::path &path, std::optional<CellFormat> format = std::nullopt,
                   const CellReadOptions &options = {});
//! Parse cell text; \p source names it in messages.
UnitCell parse_cell(std::string_view text, CellFormat format, const CellReadOptions &options = {},
                    const std::string &source = "cell");
/*!
 * \brief The cell as text. \throws std::invalid_argument if the format
 * cannot store it (partial occupancy in anything but CIF; a left-handed or
 * otherwise unusable cell).
 */
std::string write_cell(const UnitCell &cell, CellFormat format);
//! write_cell() to a file. \throws std::runtime_error if it cannot be written.
void save_cell(const UnitCell &cell, const std::filesystem::path &path, std::optional<CellFormat> format = std::nullopt);

/*!
 * \brief Parse one CIF symmetry operation ("-x+1/2, y, z-y" style) into a
 * matrix (rows act on fractional x, y, z) and a translation.
 * \throws std::invalid_argument
 */
std::pair<Mat3, Vec3> parse_symmetry_operation(std::string_view text);

//! Element of a CIF type symbol or label ("Fe", "Fe3+", "O2-", "Si1"). nullopt if none.
std::optional<unsigned> element_of_symbol(std::string_view text);

} // namespace crystal
} // namespace amcset
