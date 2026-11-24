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

#include "crystal/materials.h"

#include <cmath>
#include <stdexcept>

namespace amcset {
namespace crystal {

using namespace common;

const std::vector<ReferenceCrystal> &reference_crystals() {
  // Gemmell, Rev. Mod. Phys. 46, 129 (1974), Appendix B.
  static const std::vector<ReferenceCrystal> rows{
      // Table IV: name, structure, lattice constant (A), element, Debye
      // temperature (K), u1 at 293 K (A).
      {"C", "diamond", 3.567, {6}, 2000.0, 0.04, "diamond"},
      {"Al", "fcc", 4.050, {13}, 390.0, 0.105, ""},
      {"Si", "diamond", 5.431, {14}, 543.0, 0.075, ""},
      {"V", "bcc", 3.024, {23}, 360.0, 0.082, ""},
      {"Cr", "bcc", 2.884, {24}, 485.0, 0.061, ""},
      {"Fe", "bcc", 2.867, {26}, 420.0, 0.068, ""},
      {"Ni", "fcc", 3.524, {28}, 425.0, 0.065, ""},
      {"Cu", "fcc", 3.615, {29}, 315.0, 0.084, ""},
      {"Ge", "diamond", 5.657, {32}, 290.0, 0.085, ""},
      {"Nb", "bcc", 3.300, {41}, 275.0, 0.079, ""},
      {"Mo", "bcc", 3.147, {42}, 380.0, 0.057, ""},
      {"Rh", "fcc", 3.803, {45}, 340.0, 0.061, ""},
      {"Pd", "fcc", 3.890, {46}, 275.0, 0.074, ""},
      {"Ag", "fcc", 4.086, {47}, 215.0, 0.093, ""},
      {"Ta", "bcc", 3.306, {73}, 245.0, 0.064, ""},
      {"W", "bcc", 3.165, {74}, 310.0, 0.050, ""},
      {"Pt", "fcc", 3.923, {78}, 225.0, 0.066, ""},
      {"Au", "fcc", 4.078, {79}, 170.0, 0.087, ""},
      {"Pb", "fcc", 4.951, {82}, 88.0, 0.164, ""},
      // Table V: name, structure, lattice constant (A), elements of the
      // sublattices (ZnS = zincblende: cation, anion; CaF2 = fluorite:
      // cation, anion; NaCl = rocksalt; CsCl; perovskite: A, B, X).
      {"AlSb", "zincblende", 6.135, {13, 51}, 0.0, 0.0, ""},
      {"BaF2", "fluorite", 6.200, {56, 9}, 0.0, 0.0, ""},
      {"BaTiO3", "perovskite", 3.996, {56, 22, 8}, 0.0, 0.0, "cubic above 120 C"},
      {"CaF2", "fluorite", 5.463, {20, 9}, 0.0, 0.0, ""},
      {"CsI", "cesium_chloride", 4.567, {55, 53}, 0.0, 0.0, ""},
      {"GaP", "zincblende", 5.451, {31, 15}, 0.0, 0.0, ""},
      {"GaAs", "zincblende", 5.654, {31, 33}, 0.0, 0.0, ""},
      {"GaSb", "zincblende", 6.118, {31, 51}, 0.0, 0.0, ""},
      {"InAs", "zincblende", 6.036, {49, 33}, 0.0, 0.0, ""},
      {"InSb", "zincblende", 6.478, {49, 51}, 0.0, 0.0, ""},
      {"KCl", "rocksalt", 6.293, {19, 17}, 0.0, 0.0, ""},
      {"KBr", "rocksalt", 6.600, {19, 35}, 0.0, 0.0, ""},
      {"KI", "rocksalt", 7.066, {19, 53}, 0.0, 0.0, ""},
      {"LiF", "rocksalt", 4.017, {3, 9}, 0.0, 0.0, ""},
      {"MgO", "rocksalt", 4.211, {12, 8}, 0.0, 0.0, ""},
      {"NaF", "rocksalt", 4.620, {11, 9}, 0.0, 0.0, ""},
      {"NaCl", "rocksalt", 5.640, {11, 17}, 0.0, 0.0, ""},
      {"NaI", "rocksalt", 6.473, {11, 53}, 0.0, 0.0, ""},
      {"PbS", "rocksalt", 5.936, {82, 16}, 0.0, 0.0, ""},
      {"RbBr", "rocksalt", 6.854, {37, 35}, 0.0, 0.0, ""},
      {"SrF2", "fluorite", 5.800, {38, 9}, 0.0, 0.0, ""},
      {"ThO2", "fluorite", 5.600, {90, 8}, 0.0, 0.0, ""},
      {"UC", "rocksalt", 4.959, {92, 6}, 0.0, 0.0, ""},
      {"UO2", "fluorite", 5.468, {92, 8}, 0.0, 0.0, ""},
  };
  return rows;
}

const ReferenceCrystal *find_reference_crystal(std::string_view name) {
  for (const auto &r : reference_crystals()) {
    if (r.name == name) return &r;
  }
  return nullptr;
}

std::optional<temperature_quantity> reference_debye_temperature(unsigned z) {
  for (const auto &r : reference_crystals()) {
    if (r.sublattices.size() == 1 && r.sublattices.front() == z && r.debye_temperature_K > 0.0) {
      return r.debye_temperature_K * kelvin;
    }
  }
  return std::nullopt;
}

std::optional<length_quantity> thermal_amplitude(const Crystal &crystal, unsigned z, mass_quantity mass,
                                                 bool *from_table) {
  if (from_table != nullptr) *from_table = false;
  if (const auto u = crystal.rms_displacement(z, mass)) return u;
  // A crystal of one element: the element's tabulated Debye temperature.
  const auto composition = crystal.cell.composition();
  if (composition.size() != 1 || composition.begin()->first != z) return std::nullopt;
  const auto debye = reference_debye_temperature(z);
  if (!debye) return std::nullopt;
  if (from_table != nullptr) *from_table = true;
  return debye_rms_displacement(mass, *debye, crystal.temperature);
}

angle_quantity lindhard_angle(unsigned z1, unsigned z2, length_quantity spacing, energy_quantity energy) {
  if (z1 == 0 || z2 == 0 || !(spacing > 0.0 * angstrom) || !(energy > 0.0 * electron_volt)) {
    throw std::invalid_argument("Lindhard angle: atomic numbers, spacing and energy must be positive");
  }
  const double ratio = (2.0 * static_cast<double>(z1) * static_cast<double>(z2) * coulomb_energy_length /
                        (energy * spacing))
                           .value();
  return std::sqrt(ratio) * radian;
}

Crystal make_reference_crystal(const ReferenceCrystal &row) {
  StructureSpec spec;
  spec.name = std::string(row.structure);
  spec.a_A = row.lattice_constant_A;
  for (unsigned z : row.sublattices) spec.sublattices.push_back({Occupant{z, 1.0}});
  Crystal c;
  c.cell = build_structure(spec);
  c.cell.name = std::string(row.name) + " (" + c.cell.name + ")";
  c.temperature = 293.0 * kelvin;
  if (row.debye_temperature_K > 0.0) {
    for (unsigned z : row.sublattices) c.vibrations.push_back({z, row.debye_temperature_K * kelvin, std::nullopt});
  }
  return c;
}

} // namespace crystal
} // namespace amcset
