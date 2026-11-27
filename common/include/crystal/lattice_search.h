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

/*! \file lattice_search.h
 *  \brief The atoms of a crystal ahead of a moving particle: the collision
 *  partner search of a binary-collision transport in a lattice.
 *
 * After M. T. Robinson and I. M. Torrens, Phys. Rev. B 9, 5008 (1974)
 * (MARLOWE), Sec. IV: a target site T qualifies when its path coordinate
 * xi = lambda . (x_T - x_P) is positive (it lies ahead of the projectile P
 * moving along lambda) and its impact parameter s = sqrt(|x_T - x_P|^2 -
 * xi^2) is below a maximum s_max. Thermal vibrations are static, uncorrelated
 * Gaussian displacements along each Cartesian axis.
 *
 * The lattice is infinite and is not stored: the atom on a site (for sites
 * shared by several elements or partly vacant) and its thermal displacement
 * follow from a hash of the site's identity and a seed, so one seed is one
 * frozen configuration of the crystal (one ion history), the same however
 * often and from wherever a site is looked at.
 *
 * Lengths are Angstrom doubles in the frame the rotation maps to (the
 * laboratory frame of the transport, relative to the lattice origin).
 */

#pragma once

#include "crystal/geometry.h"
#include "crystal/unit_cell.h"

#include <cstdint>
#include <map>
#include <memory>
#include <vector>

namespace amcset {
namespace crystal {

//! A lattice site: its unit cell and its index in UnitCell::sites.
struct SiteId {
  std::int32_t i = 0;
  std::int32_t j = 0;
  std::int32_t k = 0;
  std::uint32_t site = 0;
  bool operator==(const SiteId &) const = default;
};

//! 64-bit hash of a site and a seed (well mixed: usable as a random stream's seed).
std::uint64_t site_hash(const SiteId &id, std::uint64_t seed) noexcept;

//! Hash for unordered containers of SiteId.
struct SiteIdHash {
  std::size_t operator()(const SiteId &id) const noexcept { return static_cast<std::size_t>(site_hash(id, 0)); }
};

//! An atom ahead of the projectile.
struct Candidate {
  SiteId id;
  unsigned z = 0;       //!< the element on the site
  Vec3 position;        //!< the atom's position, with its thermal displacement
  double xi = 0.0;      //!< path coordinate of the closest approach (> 0)
  double impact = 0.0;  //!< impact parameter
};

/*!
 * \brief Finds the atoms ahead of a projectile in an oriented crystal.
 */
class LatticeSearch {
public:
  /*!
   * \param cell          the unit cell (validated)
   * \param rotation      crystal -> search frame (crystal::rotation)
   * \param rms_A         rms thermal displacement along one axis per element
   *                      (Z -> Angstrom); elements without an entry are static
   * \param max_impact_A  largest impact parameter s_max
   * \param reach_A       how far ahead candidates() looks
   * \throws std::invalid_argument for non-positive distances.
   */
  LatticeSearch(const UnitCell &cell, const Mat3 &rotation, const std::map<unsigned, double> &rms_A,
                double max_impact_A, double reach_A);
  /*!
   * \brief The lattice of \p cell (the cell \p shape was built from) turned
   * by \p rotation instead, with \p shape's sites, occupants, thermal
   * amplitudes and distances. The site data is shared, not copied: the
   * grains of a polycrystal cost only their rotated vectors. The same as the
   * full constructor with \p rotation, to the last bit.
   */
  LatticeSearch(const LatticeSearch &shape, const UnitCell &cell, const Mat3 &rotation);

  /*!
   * \brief The atoms with 0 < xi <= reach and impact <= max_impact seen from
   * \p r along the unit vector \p dir, sorted by xi, into \p out (cleared
   * first). Only sites whose lattice position has x in [\p x_lo, \p x_hi)
   * exist (the layer's depth range). \p seed selects the configuration.
   */
  void candidates(const Vec3 &r, const Vec3 &dir, std::uint64_t seed, double x_lo, double x_hi,
                  std::vector<Candidate> &out) const;

  //! The lattice position of a site (no thermal displacement).
  Vec3 site_position(const SiteId &id) const noexcept;
  /*!
   * \brief The atom on a site in configuration \p seed: its element and
   * position (with thermal displacement). False if the site is vacant.
   */
  bool atom_at(const SiteId &id, std::uint64_t seed, unsigned &z, Vec3 &position) const noexcept;
  /*!
   * \brief The atom of element \p z (0: any) nearest to \p r in
   * configuration \p seed, among the sites with lattice x in [\p x_lo,
   * \p x_hi) within two cells of \p r: its site and position (with thermal
   * displacement). False if there is none. (Where a recoil that starts in
   * the crystal comes from.)
   */
  bool nearest_atom(const Vec3 &r, std::uint64_t seed, unsigned z, double x_lo, double x_hi, SiteId &id,
                    Vec3 &position) const noexcept;

  double max_impact_A() const noexcept { return max_impact_; }
  double reach_A() const noexcept { return reach_; }
  //! Largest thermal displacement an atom can have (the Gaussians are cut at 3.5 sigma per axis).
  double max_displacement_A() const noexcept { return margin_; }
  //! Largest number of sites one call can look at (its cost): the cells
  //! within reach of a point, times the occupied sites per cell.
  std::size_t sites_per_search() const noexcept;

  //! Thermal displacements are cut at this many standard deviations per axis.
  static constexpr double displacement_cut = 3.5;

private:
  struct SiteData {
    Vec3 fractional;
    std::vector<unsigned> z;         // occupants
    std::vector<double> cumulative;  // cumulative occupancy
    std::vector<double> sigma;       // rms displacement per occupant
  };
  bool occupant(const SiteData &sd, std::uint64_t &state, unsigned &z, double &sigma) const noexcept;
  //! The vectors, the search radius and the offsets for \p rotation (both constructors).
  void orient(const UnitCell &cell, const Mat3 &rotation);

  Mat3 vectors_;     // rows: lattice vectors in the search frame
  Mat3 reciprocal_;  // rows: reciprocal vectors in the search frame
  std::shared_ptr<const std::vector<SiteData>> sites_;         // the same in every orientation
  std::shared_ptr<const std::vector<std::uint32_t>> occupied_; // indices of the sites that can hold an atom
  std::vector<Vec3> offset_;            // their positions relative to the cell's origin (same order)
  double cell_reach_[3] = {0, 0, 0};    // the search radius in fractional coordinates, per axis
  double max_impact_ = 0.0;
  double reach_ = 0.0;
  double margin_ = 0.0;
};

/*!
 * \brief Radius of the largest sphere that fits between the atoms of the
 * cell: the largest distance from any point to its nearest lattice site
 * (interstitial sites do not count). A projectile is then never farther
 * than this from an atom, so a maximum impact parameter of at least this
 * radius leaves no part of the crystal free of collisions. Within about
 * 0.5 % (a grid search with refinement).
 */
double largest_hole_radius_A(const UnitCell &cell);

} // namespace crystal
} // namespace amcset
