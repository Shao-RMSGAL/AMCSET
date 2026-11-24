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

/*! \file unit_cell.h
 *  \brief Crystal lattices, unit cells with partially occupied lattice and
 *  interstitial sites, and the built-in crystal structures.
 *
 * Not TRIM: TRIM targets are amorphous. Design: docs/CRYSTAL_DESIGN.md.
 * Lengths are Angstrom doubles inside the library (see geometry.h).
 */

#pragma once

#include "amcset_units.h"
#include "crystal/geometry.h"

#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace crystal {

/*!
 * \brief A lattice: the vectors a, b, c (rows of a matrix, Angstrom).
 *
 * Built from parameters, a lies along x and b in the x-y plane (the usual
 * crystallographic Cartesian frame); read from a file, the vectors are kept
 * as the file gives them.
 */
class Lattice {
public:
  Lattice() = default;
  //! From the vectors (rows, Angstrom). \throws std::invalid_argument if degenerate or left-handed.
  static Lattice from_vectors_A(const Mat3 &rows);
  /*!
   * \brief From a, b, c and the angles alpha (b,c), beta (a,c), gamma (a,b).
   * \throws std::invalid_argument for non-positive lengths or impossible angles.
   */
  static Lattice from_parameters(common::length_quantity a, common::length_quantity b,
                                 common::length_quantity c, common::angle_quantity alpha,
                                 common::angle_quantity beta, common::angle_quantity gamma);

  const Mat3 &vectors_A() const noexcept { return m_; }
  common::length_quantity a() const noexcept { return norm(m_[0]) * common::angstrom; }
  common::length_quantity b() const noexcept { return norm(m_[1]) * common::angstrom; }
  common::length_quantity c() const noexcept { return norm(m_[2]) * common::angstrom; }
  common::angle_quantity alpha() const noexcept; //!< angle between b and c
  common::angle_quantity beta() const noexcept;  //!< angle between a and c
  common::angle_quantity gamma() const noexcept; //!< angle between a and b
  double volume_A3() const noexcept { return m_.determinant(); }
  common::volume_quantity volume() const noexcept {
    return volume_A3() * common::angstrom * common::angstrom * common::angstrom;
  }
  //! Cartesian position (Angstrom) of fractional coordinates \p f.
  Vec3 to_cartesian_A(const Vec3 &f) const noexcept { return m_.transposed() * f; }
  //! Fractional coordinates of a Cartesian position (Angstrom).
  Vec3 to_fractional(const Vec3 &r_A) const noexcept { return inv_t_ * r_A; }
  /*!
   * \brief Reciprocal vectors a*, b*, c* as rows (1/Angstrom, without 2 pi):
   * a* = (b x c)/V, so a_i . a*_j = delta_ij.
   */
  Mat3 reciprocal_per_A() const noexcept { return inv_t_; }
  //! A hexagonal lattice (a = b, alpha = beta = 90 deg, gamma = 120 deg), within 1e-6.
  bool is_hexagonal() const noexcept;
  //! The same lattice with every vector multiplied by \p s.
  Lattice scaled(double s) const;

private:
  explicit Lattice(const Mat3 &m);
  Mat3 m_ = Mat3::identity();
  Mat3 inv_t_ = Mat3::identity(); // (M^T)^-1 = rows a*, b*, c*
};

//! One element on a site with the probability that it is there.
struct Occupant {
  unsigned z = 0;
  double fraction = 1.0; //!< 0 < fraction <= 1; the fractions of a site add up to at most 1
  bool operator==(const Occupant &) const = default;
};

/*!
 * \brief Whether a site belongs to the host lattice (substitutional alloying
 * shares it) or is an interstitial position (e.g. an octahedral site that a
 * small solute may occupy).
 */
enum class SiteKind { lattice, interstitial };

std::string to_string(SiteKind k);
//! "lattice" / "substitutional" / "interstitial". \throws std::invalid_argument
SiteKind parse_site_kind(std::string_view s);

//! One site of the unit cell.
struct Site {
  std::string label;               //!< e.g. "Fe1" (files keep theirs; built-ins: element + number)
  Vec3 position;                   //!< fractional coordinates in [0, 1)
  SiteKind kind = SiteKind::lattice;
  std::vector<Occupant> occupants; //!< empty or below 1 in total: partly vacant

  //! Sum of the occupants' fractions.
  double occupancy() const noexcept;
  bool operator==(const Site &) const = default;
};

//! A distance shell around a site: \p count sites (weighted by occupancy) at \p distance_A.
struct NeighbourShell {
  double distance_A = 0.0;
  double count = 0.0;          //!< number of sites in the shell
  double occupied = 0.0;       //!< the same weighted by their occupancy
};

/*!
 * \brief A unit cell: a lattice and its sites.
 *
 * The cell does not know masses: composition() gives atoms per cell, and
 * density() takes the masses (the target element masses of the layer).
 */
class UnitCell {
public:
  std::string name;     //!< e.g. "bcc Fe" or the file's title
  Lattice lattice;
  std::vector<Site> sites;

  //! Atoms of each element (Z) per cell, lattice and interstitial sites together.
  std::map<unsigned, double> composition() const;
  //! Atoms per cell.
  double atoms_per_cell() const;
  //! Atoms per volume.
  common::number_density_quantity atomic_density() const;
  //! Mass density with the masses \p mass_of (Z -> mass).
  common::mass_density_quantity density(const std::map<unsigned, common::mass_quantity> &mass_of) const;
  //! Scale the lattice so that density() equals \p rho.
  void scale_to_density(common::mass_density_quantity rho,
                        const std::map<unsigned, common::mass_quantity> &mass_of);
  /*!
   * \brief Check the cell: a lattice with positive volume, at least one
   * occupied site, positions in [0,1), occupant fractions in (0,1] adding up
   * to at most 1 per site, known elements, and no two sites closer than
   * \p min_distance_A (minimum image).
   * \throws std::invalid_argument with \p where in the message.
   */
  void validate(const std::string &where, double min_distance_A = 0.3) const;
  //! Distance (Angstrom) between fractional positions, minimum image.
  double distance_A(const Vec3 &f1, const Vec3 &f2) const;
  //! Neighbour shells of position \p f (fractional) up to \p max_A, merged within \p tol_A.
  std::vector<NeighbourShell> neighbour_shells(const Vec3 &f, double max_A, double tol_A = 1e-3,
                                               std::optional<SiteKind> only = SiteKind::lattice) const;
};

//! Wrap a fractional coordinate into [0, 1) (values within 1e-9 of 1 become 0).
double wrap_fraction(double f) noexcept;
//! Wrap each coordinate of \p f into [0, 1).
Vec3 wrap_fractional(const Vec3 &f) noexcept;

// ---------------------------------------------------------------------------
// Built-in structures
// ---------------------------------------------------------------------------

//! Description of a built-in structure.
struct StructureInfo {
  std::string_view name;          //!< key, e.g. "bcc"
  std::string_view title;         //!< e.g. "body-centred cubic (A2)"
  std::string_view examples;      //!< e.g. "Fe, W, Mo, Cr"
  std::vector<std::string_view> sublattices; //!< names, in the order "sublattices" lists them
  bool hexagonal = false;         //!< hexagonal cell: a and c (c/a) are needed
  bool tetragonal = false;        //!< tetragonal cell: a and c
  bool has_u = false;             //!< internal parameter u (wurtzite)
  std::vector<std::string_view> interstitials; //!< built-in interstitial site sets
};

//! The built-in structures (sc, bcc, fcc, hcp, diamond, zincblende, wurtzite,
//! rocksalt, cesium_chloride, fluorite, perovskite, L1_2, L1_0).
const std::vector<StructureInfo> &builtin_structures();
//! \throws std::invalid_argument listing the known names.
const StructureInfo &structure_info(std::string_view name);

//! A set of interstitial sites of a built-in structure, with their occupants.
struct InterstitialSpec {
  std::string sites;               //!< "octahedral" or "tetrahedral"
  std::vector<Occupant> occupants; //!< usually one element with a small fraction
};

//! Everything a built-in structure needs.
struct StructureSpec {
  std::string name;                               //!< a builtin_structures() key
  double a_A = 0.0;                               //!< lattice constant a (Angstrom); 0: set later from the density
  double c_over_a = 0.0;                          //!< hexagonal / tetragonal c/a; 0: ideal (sqrt(8/3) / 1)
  double u = 0.375;                               //!< wurtzite internal parameter
  std::vector<std::vector<Occupant>> sublattices; //!< occupants of each sublattice
  std::vector<InterstitialSpec> interstitials;
};

/*!
 * \brief Build a built-in structure. With a_A = 0 the lattice constant is 1 Å
 * (scale it with UnitCell::scale_to_density).
 * \throws std::invalid_argument for unknown names, the wrong number of
 * sublattices, or unknown interstitial sets.
 */
UnitCell build_structure(const StructureSpec &spec);

} // namespace crystal
} // namespace amcset
