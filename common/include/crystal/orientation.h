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

/*! \file orientation.h
 *  \brief Orientation of a crystal in the laboratory frame, Miller indices,
 *  and the beam direction as tilt and twist.
 *
 * Laboratory frame (TRIM's): x is the depth (the inward surface normal), y
 * and z lie in the surface; the beam travels along (cos t, sin t, 0) for the
 * tilt t (`ion.angle`).
 *
 * A crystal orientation is a rotation R with v_lab = R v_crystal, where
 * v_crystal is in the lattice's Cartesian frame (Lattice::to_cartesian_A).
 * It is given as
 *  - Miller indices: the surface plane (hkl) (its normal becomes lab +x) and
 *    a reference direction [uvw] (its in-surface component becomes lab +y),
 *  - or Bunge Euler angles (phi1, Phi, phi2) of the crystal relative to a
 *    sample frame (xs, ys in the surface, zs the outward normal), with
 *    lab = (-zs, xs, -ys),
 * and then turned by the twist: a right-handed rotation of the crystal about
 * lab +x (the inward normal), as when a wafer is rotated on its holder.
 * Hexagonal cells accept four-index Miller-Bravais indices.
 */

#pragma once

#include "amcset_units.h"
#include "crystal/geometry.h"
#include "crystal/unit_cell.h"

#include <string>
#include <vector>

namespace amcset {
namespace crystal {

//! How an Orientation is given.
enum class OrientationKind { miller, euler };

//! A crystal orientation (see the file comment).
struct Orientation {
  OrientationKind kind = OrientationKind::miller;
  std::vector<double> surface{0, 0, 1};   //!< (hkl) or (hkil) of the surface plane
  std::vector<double> reference;          //!< [uvw] or [uvtw]; empty: automatic (see rotation())
  std::array<double, 3> euler_deg{0, 0, 0}; //!< Bunge phi1, Phi, phi2 (degrees)
  common::angle_quantity twist = 0.0 * common::degree;

  bool operator==(const Orientation &o) const {
    return kind == o.kind && surface == o.surface && reference == o.reference && euler_deg == o.euler_deg &&
           twist == o.twist;
  }
};

//! Three-index direction [UVW] from [uvw] or Miller-Bravais [uvtw] (hexagonal cells only).
//! \throws std::invalid_argument
Vec3 direction_indices(const Lattice &lattice, const std::vector<double> &uvw, const std::string &where);
//! Three-index plane (hkl) from (hkl) or (hkil) (hexagonal cells only, i = -(h+k)).
//! \throws std::invalid_argument
Vec3 plane_indices(const Lattice &lattice, const std::vector<double> &hkl, const std::string &where);
//! Cartesian vector (crystal frame, Angstrom) of the direction [UVW].
Vec3 direction_vector(const Lattice &lattice, const Vec3 &uvw);
//! Cartesian normal (crystal frame, 1/Angstrom) of the plane (hkl).
Vec3 plane_normal(const Lattice &lattice, const Vec3 &hkl);

/*!
 * \brief The rotation crystal -> laboratory.
 *
 * Without a reference direction the lattice vector a (or b, if a is within
 * 10 degrees of the surface normal) gives lab +y.
 * \throws std::invalid_argument (zero indices, reference parallel to the normal, ...).
 */
Mat3 rotation(const Lattice &lattice, const Orientation &o, const std::string &where = "orientation");

//! A beam direction as the tilt (`ion.angle`) and the crystal twist.
struct TiltTwist {
  common::angle_quantity tilt;
  common::angle_quantity twist;
};

/*!
 * \brief Tilt and twist that make the beam travel along the crystal
 * direction \p beam (crystal Cartesian frame, any length): the tilt is the
 * angle between the beam and the surface normal, the twist turns the crystal
 * so that the beam's in-surface part lies along lab +y. \p o's own twist is
 * replaced. With the beam along the normal the twist is kept.
 * \throws std::invalid_argument if the beam does not enter the surface
 * (90 degrees or more from the inward normal).
 */
TiltTwist tilt_twist_for_beam(const Lattice &lattice, const Orientation &o, const Vec3 &beam,
                              const std::string &where = "beam_direction");

//! Direction of the beam with tilt \p tilt in the crystal's Cartesian frame (unit vector).
Vec3 beam_in_crystal(const Lattice &lattice, const Orientation &o, common::angle_quantity tilt);

//! The low-index direction [UVW] (|U|,|V|,|W| <= \p max_index) closest to a
//! crystal-frame vector, and the angle between them.
struct NearestDirection {
  std::array<int, 3> uvw{0, 0, 1};
  common::angle_quantity off = 0.0 * common::degree;
};
NearestDirection nearest_direction(const Lattice &lattice, const Vec3 &v, int max_index = 12);
//! The low-index plane (hkl) (|h|,|k|,|l| <= \p max_index) whose normal is
//! closest to a crystal-frame vector (indices in NearestDirection::uvw).
NearestDirection nearest_plane(const Lattice &lattice, const Vec3 &v, int max_index = 12);

/*!
 * \brief The orientation that puts the crystal direction \p normal (crystal
 * Cartesian frame, any length) along the inward surface normal (lab +x) and
 * the part of \p lateral in the surface along lab +y.
 *
 * Given as Miller indices when a plane with indices up to \p max_index has
 * its normal within \p tolerance of \p normal: that plane is then the
 * surface exactly, and the reference is a direction with indices up to
 * \p max_index if one lies within the tolerance of \p lateral, otherwise the
 * automatic one with the twist that makes up the difference. Any other
 * \p normal gives Bunge Euler angles (no twist), which are exact.
 * \throws std::invalid_argument if \p normal is zero or \p lateral has no
 * part in the surface.
 */
Orientation orientation_for_frame(const Lattice &lattice, const Vec3 &normal, const Vec3 &lateral,
                                  int max_index = 6,
                                  common::angle_quantity tolerance = 0.05 * common::degree,
                                  const std::string &where = "orientation");

//! "[1 1 0]"-style text of indices (integers printed without decimals).
std::string format_indices(const std::vector<double> &v, char open = '[', char close = ']');

} // namespace crystal
} // namespace amcset
