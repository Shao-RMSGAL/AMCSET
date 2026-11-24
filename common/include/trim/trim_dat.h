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

/*! \file trim_dat.h
 *  \brief Particle lists in TRIM's TRIM.DAT format.
 *
 * TRIM.DAT drives TRIM's cascade codes 4-7: every data line starts one
 * particle with its own energy, position and direction, either an ion (codes
 * 4 and 5) or a recoiling target atom (codes 6 and 7, damage from neutrons,
 * electrons or photons). The layout (SRIM manual, chapter 8): ten comment
 * lines, of which line 8 is the title copied into the output files, then one
 * line per particle:
 *
 *     Name  Z  Energy(eV)  X(A)  Y(A)  Z(A)  cos(X)  cos(Y)  cos(Z)
 *
 * with numbers separated by spaces or commas in any numeric format. TRIM's
 * BACKSCAT.txt, TRANSMIT.txt, SPUTTER.txt and TRIMOUT.txt use the same data
 * columns (with a one-letter kind and the ion number in place of the name),
 * so those files can be read back as particle sources too.
 */

#pragma once

#include "amcset_units.h"

#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

namespace amcset {
namespace trim {

//! One particle of a TRIM.DAT file.
struct SourceParticle {
  std::string name;                        //!< Event name (up to 5 characters)
  unsigned z = 0;                          //!< Atomic number
  common::energy_quantity energy;          //!< Kinetic energy
  common::length_quantity x;               //!< Depth of the starting point
  common::length_quantity y;               //!< Lateral position
  common::length_quantity z_position;      //!< Lateral position
  double cos_x = 1.0;                      //!< Direction cosines (x = depth)
  double cos_y = 0.0;
  double cos_z = 0.0;
};

//! A parsed TRIM.DAT file.
struct ParticleList {
  std::string title;                       //!< Line 8 of the file
  std::vector<SourceParticle> particles;
  //! Lines that could not be read (with their line numbers); TRIM skips them too.
  std::vector<std::string> warnings;
};

/*!
 * \brief Read a TRIM.DAT particle list. Malformed data lines are skipped and
 * reported in ParticleList::warnings; a trailing end-of-file character
 * (ASCII 26) is accepted. \throws std::runtime_error if the file has fewer
 * than ten header lines.
 */
ParticleList parse_trim_dat(std::istream &in);
//! parse_trim_dat() from a file. \throws std::runtime_error
ParticleList read_trim_dat(const std::filesystem::path &path);

/*!
 * \brief Read the exit records of BACKSCAT.txt, TRANSMIT.txt, SPUTTER.txt or
 * TRIMOUT.txt (lines "B|T|S ion Z E x y z cx cy cz") as particles, e.g. to
 * use sputtered atoms or transmitted ions as the source of a new
 * calculation. The particle names are "B0001"-style tags.
 */
ParticleList read_exit_file_as_particles(const std::filesystem::path &path);

//! Write \p list as a TRIM.DAT that TRIM accepts (CR LF line ends, EOF character).
void write_trim_dat(std::ostream &out, const ParticleList &list);

} // namespace trim
} // namespace amcset
