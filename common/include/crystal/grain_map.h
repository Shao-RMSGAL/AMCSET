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

/*! \file grain_map.h
 *  \brief Polycrystals: which grain, with which orientation, is at a point
 *  of a layer.
 *
 * Not TRIM. A layer's crystal can be a polycrystal: grains of the same unit
 * cell with different orientations. The grain map works in the layer's
 * frame, in Angstrom: x is the depth below the layer's top, y and z are
 * lateral and periodic (a representative patch of the surface). Grains are
 * the cells of a Voronoi or Laguerre tessellation of seed points (generated,
 * from a grain file, from Neper), columns under the points of an EBSD map,
 * or the voxels of a grid. Design: docs/CRYSTAL_DESIGN.md.
 */

#pragma once

#include "amcset_units.h"
#include "crystal/geometry.h"
#include "crystal/unit_cell.h"

#include <array>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace crystal {

//! One grain: how its crystal sits in the layer.
struct Grain {
  Mat3 rotation;  //!< crystal -> laboratory, as crystal::rotation()
  Vec3 origin_A;  //!< a lattice point of the grain (layer frame): grains are not in registry
};

/*!
 * \brief The grains of a layer.
 *
 * Built from seeds (the grain of a point is that of the seed with the
 * smallest |r - p|^2 - weight: a Voronoi tessellation for equal weights, a
 * Laguerre one otherwise; several seeds may belong to one grain, as the
 * points of an EBSD map do) or from a grid of voxels.
 */
class GrainMap {
public:
  struct Seed {
    Vec3 position_A;
    double weight_A2 = 0.0; //!< Laguerre weight (a squared radius)
    int grain = -1;         //!< index in grains(); -1: no crystal there (treated as amorphous)
  };
  //! What at() finds.
  struct At {
    int grain = -1; //!< -1: no grain (amorphous)
    //! Distance to the nearest boundary with another grain: exact where it
    //! is below the cap given to at(), at least the cap otherwise.
    double boundary_A = std::numeric_limits<double>::infinity();
  };

  GrainMap() = default;
  /*!
   * \brief From seeds. \p columnar: the grains are columns along the depth
   * (the seeds' x is not used).
   * \throws std::invalid_argument for no seeds, non-positive periods or a
   * seed whose grain does not exist.
   */
  static GrainMap from_seeds(std::vector<Grain> grains, std::vector<Seed> seeds, double period_y_A, double period_z_A,
                             bool columnar);
  /*!
   * \brief From voxels: \p ids holds n[0] x n[1] x n[2] grain indices (-1:
   * none), the depth index slowest and z fastest, over \p depth_A and the
   * periods. Below \p depth_A the last voxels continue.
   * \throws std::invalid_argument for sizes that do not match or ids out of range.
   */
  static GrainMap from_voxels(std::vector<Grain> grains, std::vector<int> ids, std::array<int, 3> n, double depth_A,
                              double period_y_A, double period_z_A);

  /*!
   * \brief The grain at \p r (layer frame, Angstrom) and its distance to the
   * nearest grain boundary, exact where below \p cap_A. For voxels the
   * boundaries are the faces between voxels of different grains.
   */
  At at(const Vec3 &r, double cap_A = 0.0) const;
  //! at() by looking at every seed and its lateral images (for tests).
  At at_brute_force(const Vec3 &r) const;

  /*!
   * \brief The same map with other grains (orientations and origins), as
   * many as before.
   * \throws std::invalid_argument for another number of grains.
   */
  GrainMap with_grains(std::vector<Grain> grains) const;

  const std::vector<Grain> &grains() const noexcept { return grains_; }
  const Grain &grain(int index) const { return grains_.at(static_cast<std::size_t>(index)); }
  std::size_t grain_count() const noexcept { return grains_.size(); }
  const std::vector<Seed> &seeds() const noexcept { return seeds_; }
  bool is_voxels() const noexcept { return !ids_.empty(); }
  bool is_columnar() const noexcept { return columnar_; }
  //! Voxels: the grid (depth, y, z), its depth and the grain of a voxel (-1: none).
  std::array<int, 3> voxel_grid() const noexcept { return n_; }
  double voxel_depth_A() const noexcept { return depth_A_; }
  int voxel(int i, int j, int k) const noexcept;
  double period_y_A() const noexcept { return period_[0]; }
  double period_z_A() const noexcept { return period_[1]; }
  //! Fraction of the map's points that belong to no grain (voxels and EBSD points).
  double ungrained_fraction() const noexcept { return ungrained_; }

private:
  At at_seeds(const Vec3 &r, double cap_A, bool brute_force) const;
  At at_voxels(const Vec3 &r) const;

  std::vector<Grain> grains_;
  std::vector<Seed> seeds_;
  bool columnar_ = false;
  double period_[2] = {1.0, 1.0};
  double ungrained_ = 0.0;
  // Seeds: buckets over (x, y, z) or (y, z).
  double x_lo_ = 0.0;                 // lowest seed depth
  double bucket_A_[3] = {1.0, 1.0, 1.0};
  int buckets_[3] = {1, 1, 1};
  std::vector<std::uint32_t> bucket_start_; // per bucket: its first entry in bucket_seed_ (+1 at the end)
  std::vector<std::uint32_t> bucket_seed_;
  double weight_max_ = 0.0;           // (weights are shifted so that the smallest is 0)
  // Voxels.
  std::vector<int> ids_;
  std::array<int, 3> n_{0, 0, 0};
  double depth_A_ = 0.0;
};

/*!
 * \brief A layer's crystal as a polycrystal: where its grains come from.
 *
 * `voronoi` generates the grains (a Poisson-Voronoi tessellation); the other
 * sources read them from Polycrystal::file, and #map then holds what was
 * read (config::parse_crystal and the GUI read it; see crystal/grain_io.h).
 */
struct Polycrystal {
  enum class Source { voronoi, grain_file, neper, ebsd };
  /*!
   * Orientations of generated grains: `random` (uniform), `fiber` (the
   * crystal's surface plane for every grain, turned at random about the
   * normal) or `mosaic` (the crystal's orientation), the last two with
   * #spread.
   */
  enum class Texture { random, fiber, mosaic };

  Source source = Source::voronoi;
  // --- generated grains ---
  common::length_quantity grain_size = 500.0 * common::angstrom; //!< mean diameter (sphere or, for columns, circle of the mean grain)
  bool columnar = false;
  Texture texture = Texture::random;
  common::angle_quantity spread = 0.0 * common::degree; //!< fiber: rms tilt of the axis; mosaic: rms rotation
  std::uint64_t seed = 1;                               //!< of the generator (not the run's seed)
  /*!
   * Lateral periods (y, z). 0: eight grain sizes for generated grains, the
   * extent of the file's map otherwise.
   */
  common::length_quantity period_y = 0.0 * common::angstrom;
  common::length_quantity period_z = 0.0 * common::angstrom;
  // --- every source ---
  //! Width of the disordered slab at grain boundaries, treated as amorphous (0: none).
  common::length_quantity boundary_width = 0.0 * common::angstrom;
  /*!
   * Neighbouring grains misoriented by less than this are put in registry
   * (put_in_registry): their lattices continue across the boundary, as at a
   * subgrain (low-angle) boundary. 0: every grain's lattice is its own.
   */
  common::angle_quantity registry_below = 0.0 * common::degree;
  // --- files ---
  std::string file;                                              //!< as written in the scenario
  common::length_quantity scale = 1.0 * common::micrometer;     //!< Neper: the length of one unit of the tessellation
  common::angle_quantity misorientation = 5.0 * common::degree; //!< EBSD: points closer in orientation than this are one grain
  int phase = 0;                                                 //!< EBSD: the phase that is this crystal (0: every indexed point)
  std::shared_ptr<const GrainMap> map;                           //!< file sources: the map read

  //! The description (the map itself is not compared).
  bool operator==(const Polycrystal &o) const {
    return source == o.source && grain_size == o.grain_size && columnar == o.columnar && texture == o.texture &&
           spread == o.spread && seed == o.seed && period_y == o.period_y && period_z == o.period_z &&
           boundary_width == o.boundary_width && registry_below == o.registry_below && file == o.file && scale == o.scale &&
           misorientation == o.misorientation && phase == o.phase;
  }
};

std::string to_string(Polycrystal::Source s);
//! "voronoi" ("generated"), "file" ("grains"), "neper" ("tess"), "ebsd". \throws std::invalid_argument
Polycrystal::Source parse_polycrystal_source(std::string_view s);
std::string to_string(Polycrystal::Texture t);
//! "random", "fiber" ("fibre"), "mosaic". \throws std::invalid_argument
Polycrystal::Texture parse_texture(std::string_view s);

//! Limits of a polycrystal (check_polycrystal).
inline constexpr std::size_t max_grains = 2'000'000;
inline constexpr double max_boundary_width_A = 100.0;
inline constexpr double max_registry_deg = 30.0;

/*!
 * \brief Check a polycrystal's description: sizes, periods, spread and
 * boundary width in range; a file for the file sources.
 * \throws std::invalid_argument naming the key (prefixed with \p where).
 */
void check_polycrystal(const Polycrystal &p, const std::string &where);

/*!
 * \brief Generated grains for a layer \p depth_A thick: Poisson-Voronoi
 * seeds (also a grain size above and below the layer, so that the grains at
 * its faces are like the others), orientations by the texture around
 * \p base (the crystal's own rotation), each grain's lattice through its
 * seed. The same description gives the same grains.
 * \throws std::invalid_argument if there would be more than max_grains.
 */
GrainMap generate_voronoi(const Polycrystal &p, const Mat3 &base, double depth_A);

/*!
 * \brief The proper rotations that map the lattice onto itself, in its
 * Cartesian frame, by its metric: cubic (24), hexagonal (12), tetragonal
 * (8), orthorhombic (4), otherwise the identity. (The lattice's symmetry,
 * which is at least that of the crystal: enough to tell grains apart.)
 */
std::vector<Mat3> lattice_rotations(const Lattice &lattice);

/*!
 * \brief The smallest rotation angle (radians) between two orientations
 * (crystal -> laboratory) of a crystal with the rotations \p symmetry.
 */
double misorientation_rad(const Mat3 &a, const Mat3 &b, const std::vector<Mat3> &symmetry);

//! A uniformly distributed rotation from three uniform numbers in [0, 1) (Shoemake's method).
Mat3 uniform_rotation(double u1, double u2, double u3) noexcept;

//! Sample frame (x, y in the surface, z the outward normal) -> laboratory: (depth, y, z) = (-zs, xs, -ys).
Mat3 sample_to_laboratory() noexcept;

} // namespace crystal
} // namespace amcset
