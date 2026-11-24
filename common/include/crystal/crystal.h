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

/*! \file crystal.h
 *  \brief The crystal of a target layer: unit cell, orientation, temperature
 *  and thermal vibrations.
 *
 * Not TRIM (TRIM targets are amorphous); docs/CRYSTAL_DESIGN.md.
 */

#pragma once

#include "amcset_units.h"
#include "crystal/grain_map.h"
#include "crystal/orientation.h"
#include "crystal/unit_cell.h"

#include <optional>
#include <string>
#include <vector>

namespace amcset {
namespace crystal {

/*!
 * \brief Thermal vibrations of one element: a Debye temperature (the rms
 * amplitude then follows from the crystal temperature and the mass) or the
 * rms amplitude itself (one Cartesian component).
 */
struct Vibration {
  unsigned z = 0;
  std::optional<common::temperature_quantity> debye_temperature;
  std::optional<common::length_quantity> rms;
  bool operator==(const Vibration &) const = default;
};

//! The crystal of a layer.
struct Crystal {
  UnitCell cell;
  Orientation orientation;
  common::temperature_quantity temperature = 300.0 * common::kelvin;
  std::vector<Vibration> vibrations; //!< per element; elements without one have none given yet
  /*!
   * \brief Empty: a single crystal. Otherwise the layer is a polycrystal:
   * grains of this unit cell (grain_map.h). #orientation is then the
   * orientation the fiber and mosaic textures of generated grains are
   * built around; grains from files bring their own.
   */
  std::optional<Polycrystal> polycrystal;

  /*!
   * \brief rms thermal displacement (one Cartesian component) of element
   * \p z with mass \p mass; nullopt if the crystal gives none for it.
   */
  std::optional<common::length_quantity> rms_displacement(unsigned z, common::mass_quantity mass) const;
  //! The crystal -> laboratory rotation (orientation.h).
  Mat3 rotation() const { return crystal::rotation(cell.lattice, orientation); }

  bool operator==(const Crystal &o) const {
    return cell.name == o.cell.name && cell.lattice.vectors_A().row == o.cell.lattice.vectors_A().row &&
           cell.sites == o.cell.sites && orientation == o.orientation && temperature == o.temperature &&
           vibrations == o.vibrations && polycrystal == o.polycrystal;
  }
};

/*!
 * \brief The grains of a polycrystalline \p crystal in a layer \p depth_A
 * thick: generated (Polycrystal::Source::voronoi) or the map read from its
 * file.
 * \throws std::invalid_argument if the crystal is not a polycrystal or its
 * file has not been read (Polycrystal::map).
 */
std::shared_ptr<const GrainMap> make_grain_map(const Crystal &crystal, double depth_A);

/*!
 * \brief The Debye function Phi(x) = (1/x) Int_0^x t / (e^t - 1) dt
 * (Phi(0) = 1, Phi(x) -> pi^2 / (6 x) for large x).
 */
double debye_phi(double x);

/*!
 * \brief rms thermal displacement along one axis in the Debye model,
 * u = sqrt(3 hbar^2 / (M k_B Theta) * (Phi(Theta/T)/(Theta/T) + 1/4)):
 * 3 hbar^2 / (4 M k_B Theta) at T = 0 (zero-point motion), 3 hbar^2 T /
 * (M k_B Theta^2) for T >> Theta.
 * \throws std::invalid_argument for a non-positive mass or Debye temperature, or T < 0.
 */
common::length_quantity debye_rms_displacement(common::mass_quantity mass, common::temperature_quantity debye,
                                               common::temperature_quantity temperature);

} // namespace crystal
} // namespace amcset
