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

/*! \file layer_crystal.cpp
 *  \brief A target layer and its crystal: the layer's elements,
 *  stoichiometry and density follow from the unit cell (not TRIM).
 */

#include "trim/trim.h"

#include "element_data.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

namespace amcset {
namespace trim {

using namespace common;

namespace {

std::string num(double v) {
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

std::map<unsigned, mass_quantity> masses_of(const TargetLayer &layer) {
  std::map<unsigned, mass_quantity> m;
  for (const auto &e : layer.elements) m[e.z] = e.mass;
  return m;
}

} // namespace

void apply_crystal(TargetLayer &layer) {
  if (!layer.crystal) return;
  const auto composition = layer.crystal->cell.composition();
  if (composition.empty()) throw std::invalid_argument("crystal: the unit cell has no atoms");
  // Elements that left the crystal go; the others keep their energies.
  layer.elements.erase(std::remove_if(layer.elements.begin(), layer.elements.end(),
                                      [&](const TargetElement &e) { return composition.count(e.z) == 0; }),
                       layer.elements.end());
  for (const auto &[z, atoms] : composition) {
    auto it = std::find_if(layer.elements.begin(), layer.elements.end(),
                           [z = z](const TargetElement &e) { return e.z == z; });
    if (it == layer.elements.end()) {
      layer.elements.push_back(TargetElement::defaults(z, atoms));
    } else {
      it->stoichiometry = atoms;
    }
  }
  layer.density = layer.crystal->cell.density(masses_of(layer));
}

void check_layer_crystal(const TargetLayer &layer, const std::string &where) {
  if (!layer.crystal) return;
  const auto &c = *layer.crystal;
  const std::string cw = where + ": crystal";
  if (layer.gas) throw std::invalid_argument(cw + ": a gas layer cannot be a crystal");
  c.cell.validate(cw);
  // Elements and stoichiometry as in the unit cell.
  const auto composition = c.cell.composition();
  double layer_sum = 0.0;
  for (const auto &e : layer.elements) layer_sum += e.stoichiometry;
  const double cell_sum = c.cell.atoms_per_cell();
  for (const auto &e : layer.elements) {
    const auto it = composition.find(e.z);
    const std::string sym(Element::by_z(e.z).symbol());
    if (it == composition.end()) {
      if (e.stoichiometry > 0.0) throw std::invalid_argument(cw + ": the layer lists " + sym + ", which the unit cell does not contain");
      continue;
    }
    const double in_layer = e.stoichiometry / layer_sum;
    const double in_cell = it->second / cell_sum;
    if (std::abs(in_layer - in_cell) > 1e-6) {
      throw std::invalid_argument(cw + ": the layer's fraction of " + sym + " (" + num(in_layer) +
                                  ") is not the unit cell's (" + num(in_cell) + ")");
    }
  }
  for (const auto &[z, atoms] : composition) {
    if (std::none_of(layer.elements.begin(), layer.elements.end(), [z = z](const TargetElement &e) { return e.z == z; })) {
      throw std::invalid_argument(cw + ": the unit cell contains " + std::string(Element::by_z(z).symbol()) +
                                  ", which is not among the layer's elements");
    }
  }
  const double rho_cell = in_units(c.cell.density(masses_of(layer)), gram_per_cubic_centimeter);
  const double rho = in_units(layer.density, gram_per_cubic_centimeter);
  if (std::abs(rho - rho_cell) > 1e-6 * rho_cell) {
    throw std::invalid_argument(cw + ": the layer's density (" + num(rho) + " g/cm3) is not the unit cell's (" +
                                num(rho_cell) + " g/cm3)");
  }
  // Orientation (throws with the reason), temperature, vibrations.
  (void)crystal::rotation(c.cell.lattice, c.orientation, cw + ".orientation");
  if (c.polycrystal) crystal::check_polycrystal(*c.polycrystal, cw + ".polycrystal");
  const double t = in_units(c.temperature, kelvin);
  if (!(t >= 0.0) || !(t <= max_crystal_temperature_K)) {
    throw std::invalid_argument(cw + ".temperature: 0 to " + num(max_crystal_temperature_K) + " K, got " + num(t) + " K");
  }
  for (std::size_t i = 0; i < c.vibrations.size(); ++i) {
    const auto &v = c.vibrations[i];
    const std::string vw = cw + ".vibrations." + std::to_string(i);
    if (composition.count(v.z) == 0) {
      throw std::invalid_argument(vw + ": element Z = " + std::to_string(v.z) + " is not in the unit cell");
    }
    for (std::size_t k = 0; k < i; ++k) {
      if (c.vibrations[k].z == v.z) throw std::invalid_argument(vw + ": this element is listed twice");
    }
    if (!v.debye_temperature && !v.rms) throw std::invalid_argument(vw + ": give debye_temperature or rms");
    if (v.debye_temperature) {
      const double d = in_units(*v.debye_temperature, kelvin);
      if (!(d > 0.0) || !(d <= max_debye_temperature_K)) {
        throw std::invalid_argument(vw + ".debye_temperature: above 0 and at most " + num(max_debye_temperature_K) +
                                    " K (diamond: about 2200 K), got " + num(d) + " K");
      }
    }
    if (v.rms) {
      const double r = in_units(*v.rms, angstrom);
      if (!(r >= 0.0) || !(r <= max_vibration_rms_A)) {
        throw std::invalid_argument(vw + ".rms: 0 to " + num(max_vibration_rms_A) +
                                    " A (thermal amplitudes are a few hundredths of an Angstrom), got " + num(r) + " A");
      }
    }
  }
}

} // namespace trim
} // namespace amcset
