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

/*! \file grain_io.h
 *  \brief Grains from files: AMCSET's grain file, Neper tessellations and
 *  EBSD maps.
 *
 * Orientations in files are relative to the sample: x and y in the surface,
 * z the outward normal (as crystal::Orientation's Euler angles). A point
 * (xs, ys, zs) of the sample is at depth z_top - zs and at the lateral
 * position (xs, -ys) in the layer.
 *
 * **Grain file** (text; `#` starts a comment; lengths in `unit`):
 * \code
 *   unit nm                       # A (default), nm or um
 *   period 500 500                # lateral periods (y, z); or in the scenario
 *   columnar                      # optional: grains are columns
 *   grain x y z phi1 Phi phi2 [weight]   # one seed: depth, y, z; Bunge angles (degrees);
 *                                        # an optional Laguerre weight (a radius)
 * \endcode
 * or a grid of voxels:
 * \code
 *   unit nm
 *   period 500 500
 *   depth 200
 *   voxels nx ny nz
 *   orientation phi1 Phi phi2     # one per grain, numbered from 0
 *   ids                           # then nx*ny*nz grain numbers, the depth
 *   0 0 1 1 ...                   # index slowest and z fastest; -1: no grain
 * \endcode
 *
 * **Neper** `.tess`: the seeds with their weights (a Laguerre tessellation
 * by |r - p|^2 - w^2) and the cells' orientations (`*ori`: euler-bunge,
 * rodrigues, quaternion, axis-angle or rotmat, active unless marked
 * `:passive`); the box is that of the vertices and one unit of the file is
 * Polycrystal::scale long. The tessellation is repeated laterally, so cells
 * at the sides of a non-periodic tessellation differ from Neper's.
 *
 * **EBSD** maps, TSL/EDAX `.ang` (Euler angles in radians, x and y in
 * micrometres; square or hexagonal grids) and Oxford/HKL `.ctf` (degrees):
 * neighbouring points closer in orientation than
 * Polycrystal::misorientation, up to the lattice's rotations, are one
 * grain; every grain is a column into the depth with the orientation of
 * its first point; points that are not indexed (or of another phase than
 * Polycrystal::phase) take the grain of the nearest indexed point.
 */

#pragma once

#include "crystal/grain_map.h"

#include <filesystem>
#include <memory>
#include <string>
#include <string_view>

namespace amcset {
namespace crystal {

//! The format of a grain file by its name: .tess -> neper, .ang/.ctf -> ebsd, anything else the grain file.
Polycrystal::Source grain_source_from_path(const std::filesystem::path &path);

/*!
 * \brief Read the grains of \p p (a file source) from \p path for a crystal
 * with the unit cell \p cell.
 * \throws std::runtime_error if the file cannot be read,
 *         std::invalid_argument for malformed content (with the line).
 */
std::shared_ptr<const GrainMap> read_grain_map(const Polycrystal &p, const std::filesystem::path &path,
                                               const UnitCell &cell);

//! Parse the text of an AMCSET grain file; \p source names it in messages.
GrainMap parse_grain_file(std::string_view text, const Polycrystal &p, const UnitCell &cell,
                          const std::string &source = "grain file");
//! Parse a Neper tessellation.
GrainMap parse_neper_tess(std::string_view text, const Polycrystal &p, const UnitCell &cell,
                          const std::string &source = "tess");
//! Parse an EBSD map (\p ctf: Oxford's .ctf, otherwise TSL's .ang).
GrainMap parse_ebsd(std::string_view text, bool ctf, const Polycrystal &p, const UnitCell &cell,
                    const std::string &source = "EBSD map");

} // namespace crystal
} // namespace amcset
