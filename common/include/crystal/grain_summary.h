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

/*! \file grain_summary.h
 *  \brief Polycrystals: which grains touch, what the grains are like, and
 *  low-angle boundaries in registry.
 *
 * Not TRIM. Works on a GrainMap (grain_map.h) in the layer's frame, in
 * Angstrom (x = depth below the layer's top, y and z lateral and periodic).
 * Design: docs/CRYSTAL_DESIGN.md section 13.
 */

#pragma once

#include "crystal/grain_map.h"
#include "crystal/unit_cell.h"

#include <cstddef>
#include <string>
#include <vector>

namespace amcset {
namespace crystal {

//! Two grains that touch, and a point of their boundary.
struct GrainBoundary {
  int a = 0; //!< the grain with the smaller index
  int b = 0;
  //! On (or, from the grid, within half a grid step of) the boundary, inside
  //! the layer (0 <= x <= depth); laterally it may lie up to a period
  //! outside [0, period), next to the grain it was found from.
  Vec3 point_A;
  std::size_t contacts = 0; //!< seed pairs, grid neighbours or voxel faces found on this boundary
};

/*!
 * \brief The pairs of grains of \p map that touch inside a layer \p depth_A
 * thick, sorted by (a, b).
 *
 * Seeds: (1) for each pair of nearby seeds of different grains, the point
 * of equal power distance on the line between them, when no other seed
 * claims it; (2) a grid of at most about \p grid_points points over the
 * layer (over its surface for columns; 0: about 500 per grain, from
 * 5 x 10^4 to 2 x 10^6): neighbouring points in different grains, with the point
 * halfway between them. (1) is exact but misses the
 * faces that do not cross the line between their seeds (about a quarter of
 * a Poisson-Voronoi tessellation's); (2) misses faces smaller than the grid.
 * Voxels: the faces between voxels of different grains. Of several points
 * of one boundary the one nearest their middle is kept. Points of no grain
 * (grain -1) touch nothing.
 */
std::vector<GrainBoundary> grain_boundaries(const GrainMap &map, double depth_A, std::size_t grid_points = 0);

//! summarize_grains counts the boundaries of maps up to this many grains.
inline constexpr std::size_t max_grains_for_boundaries = 300'000;

//! What a layer's grains are like (summarize_grains).
struct GrainSummary {
  std::size_t grains = 0;          //!< in the map (generated grains include those just above and below the layer)
  std::size_t grains_in_layer = 0; //!< met by the sample points
  bool columnar = false;           //!< columns (sizes are of their sections)
  std::size_t samples = 0;         //!< points looked at: a grid over the layer (columns: over its surface)
  double sample_spacing_A = 0.0;
  //! Fewer than 32 sample points per grain: small grains are missed and
  //! the sizes are rough.
  bool coarse = false;
  double period_y_A = 0.0, period_z_A = 0.0;
  double ungrained_fraction = 0.0; //!< of the sample points: in no grain (amorphous)
  /*!
   * Equivalent diameter of the grains met (the sphere of their sampled volume,
   * the circle of their sampled section for columns): mean, median, 10th and
   * 90th percentile. Grains much smaller than the spacing are not resolved.
   */
  double diameter_mean_A = 0.0, diameter_median_A = 0.0, diameter_p10_A = 0.0, diameter_p90_A = 0.0;
  /*!
   * Texture, weighted by the sampled volume, against a reference orientation
   * (the crystal's own): the misorientation (smallest rotation, up to the
   * lattice's rotations), mean and rms; and the tilt of the crystal
   * direction along the surface normal from the reference's, mean.
   */
  double misorientation_mean_deg = 0.0, misorientation_rms_deg = 0.0, normal_tilt_mean_deg = 0.0;
  //! Boundaries (grain_boundaries): how many, their mean misorientation, the
  //! share below 15 degrees; not counted for maps of more than
  //! max_grains_for_boundaries grains.
  bool boundaries_counted = false;
  std::size_t boundaries = 0;
  double boundary_misorientation_mean_deg = 0.0;
  double low_angle_fraction = 0.0;
  //! Registry (Polycrystal::registry_below): the threshold (0: off), the groups of more than one grain and their grains.
  double registry_below_deg = 0.0;
  std::size_t registry_groups = 0;
  std::size_t registry_grains = 0;
};

/*!
 * \brief Summary of \p map in a layer \p depth_A thick, for a crystal of
 * \p lattice with the reference orientation \p reference (crystal ->
 * laboratory), from about 64 sample points per grain (at least 100 000),
 * at most about \p max_samples.
 */
GrainSummary summarize_grains(const GrainMap &map, const Lattice &lattice, const Mat3 &reference, double depth_A,
                              std::size_t max_samples = 1'000'000);

//! The summary as one JSON object (summary.json "crystals"[].grains, `amcset grains --json`).
std::string to_json(const GrainSummary &s);

//! What put_in_registry did.
struct RegistryResult {
  std::size_t groups = 0; //!< groups of more than one grain
  std::size_t grains = 0; //!< grains in those groups
};

/*!
 * \brief Low-angle boundaries in registry. Neighbouring grains
 * (grain_boundaries) misoriented by less than \p below_rad form groups; in
 * each group the grain with the smallest index keeps its lattice and the
 * others, in breadth-first order from it, take the orientation equivalent
 * to theirs (up to the lattice's rotations) nearest their parent's and pass
 * through the lattice point of their parent's lattice nearest their common
 * boundary point. The lattice then continues across that boundary and
 * drifts away from it by the misorientation times the distance, as at a
 * subgrain boundary; a group of one orientation is one lattice everywhere.
 * Boundaries not used to build a group (a group is a tree of them) and the
 * lateral repeats of the map are in registry only as far as that drift
 * allows. The map is otherwise unchanged.
 */
GrainMap put_in_registry(const GrainMap &map, const Lattice &lattice, double below_rad, double depth_A,
                         RegistryResult *result = nullptr);

} // namespace crystal
} // namespace amcset
