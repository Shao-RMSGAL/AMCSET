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

#include "crystal/crystal_database.h"

#include "crystal/materials.h"
#include "crystal/unit_cell.h"
#include "element_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace amcset {
namespace crystal {

using namespace common;

const std::vector<DatabaseSource> &database_sources() {
  static const std::vector<DatabaseSource> sources{
      {"gemmell1974", "D. S. Gemmell, Rev. Mod. Phys. 46, 129 (1974), Appendix B, Tables IV and V"},
      {"wolfram", "Wolfram Research, Mathematica ElementData (lattice constants), as listed at "
                  "https://periodictable.com/Properties/A/LatticeConstants.html (retrieved 2026-11-26)"},
      {"kittel2004", "C. Kittel, Introduction to Solid State Physics, 8th ed. (Wiley, 2004), Debye temperatures as "
                     "tabulated in Wikipedia, \"Debye model\" (retrieved 2026-11-26)"},
      {"ioffe", "Ioffe Institute, \"Semiconductor NSM\" database (https://www.ioffe.ru/SVA/NSM/Semicond/), via "
                "Wikipedia, \"Lattice constant\" (retrieved 2026-11-26)"},
      {"aps", "Argonne National Laboratory, Advanced Photon Source, \"Lattice Constants\" "
              "(http://7id.xray.aps.anl.gov/calculators/crystal_lattice_parameters.html), via Wikipedia, "
              "\"Lattice constant\" (retrieved 2026-11-26)"},
      {"nist", "NIST, Fundamental physical constants (lattice parameter of silicon), via Wikipedia, \"Lattice "
               "constant\" (retrieved 2026-11-26)"},
      {"davey1925", "W. P. Davey, \"Precision Measurements of the Lattice Constants of Twelve Common Metals\", Phys. "
                    "Rev. 25, 753 (1925), via Wikipedia, \"Lattice constant\" (retrieved 2026-11-26)"},
      {"toth1967", "L. E. Toth, Transition Metal Carbides and Nitrides (Academic Press, 1967), via Wikipedia, "
                   "\"Lattice constant\" (retrieved 2026-11-26)"},
      {"morkoc2009", "H. Morkoc and U. Ozgur, Zinc Oxide: Fundamentals, Materials and Device Technology (Wiley-VCH, "
                     "2009), via Wikipedia, \"Lattice constant\" (retrieved 2026-11-26)"},
      {"spi", "SPI Supplies, \"Substrates\", via Wikipedia, \"Lattice constant\" (retrieved 2026-11-26)"},
      {"saha2010", "B. Saha et al., J. Appl. Phys. 107, 033715 (2010), via Wikipedia, \"Lattice constant\" "
                   "(retrieved 2026-11-26)"},
      {"landolt", "J. B. Goodenough and J. M. Longo, Landolt-Bornstein III/4a, 3.1.7 (perovskites), via Wikipedia, "
                  "\"Lattice constant\" (retrieved 2026-11-26)"},
      {"wikipedia", "Wikipedia, \"Lattice constant\", table of lattice constants at 300 K, row without a reference "
                    "(retrieved 2026-11-26)"},
  };
  return sources;
}

std::string_view database_citation(std::string_view key) {
  for (const auto &s : database_sources()) {
    if (s.key == key) return s.citation;
  }
  return {};
}

const std::vector<DatabaseCrystal> &crystal_database() {
  static const std::vector<DatabaseCrystal> rows = [] {
    std::vector<DatabaseCrystal> r;
    // Gemmell's Tables IV and V first, exactly as materials.h has them.
    for (const auto &g : reference_crystals()) {
      DatabaseCrystal d;
      d.name = g.name;
      d.structure = g.structure;
      d.a_A = g.lattice_constant_A;
      d.sublattices = g.sublattices;
      d.debye_temperature_K = g.debye_temperature_K;
      d.lattice_source = "gemmell1974";
      d.debye_source = g.debye_temperature_K > 0.0 ? "gemmell1974" : "";
      d.note = g.note;
      r.push_back(d);
    }
    // name, structure, a, c (A), sublattices, Debye (K), lattice source,
    // Debye source, note, preferred.
    const std::vector<DatabaseCrystal> more{
        // --- Elements: lattice constants from Mathematica's ElementData (pm
        // there, A here); structures as Wikipedia's "Periodic table (crystal
        // structure)" gives them at room temperature (W type = bcc, Cu type =
        // fcc, Mg type = hcp). Debye temperatures from Kittel where listed.
        // Left out: Rb (its constant is a low-temperature one: 6 % too dense
        // at room temperature) and Ac (a = 5.67 A gives 8.3 g/cm3, not ~10).
        {"Li", "bcc", 3.51, 0.0, {3}, 0.0, "wolfram", "", "", true},
        {"Na", "bcc", 4.2906, 0.0, {11}, 0.0, "wolfram", "", "", true},
        {"K", "bcc", 5.328, 0.0, {19}, 0.0, "wolfram", "", "", true},
        {"Cs", "bcc", 6.141, 0.0, {55}, 38.0, "wolfram", "kittel2004", "", true},
        {"Ba", "bcc", 5.028, 0.0, {56}, 0.0, "wolfram", "", "", true},
        {"Eu", "bcc", 4.581, 0.0, {63}, 0.0, "wolfram", "", "", true},
        {"Ra", "bcc", 5.148, 0.0, {88}, 0.0, "wolfram", "", "", true},
        {"Ca", "fcc", 5.5884, 0.0, {20}, 0.0, "wolfram", "", "", true},
        {"Sr", "fcc", 6.0849, 0.0, {38}, 0.0, "wolfram", "", "", true},
        {"Yb", "fcc", 5.4847, 0.0, {70}, 0.0, "wolfram", "", "", true},
        {"Ir", "fcc", 3.839, 0.0, {77}, 0.0, "wolfram", "", "", true},
        {"Th", "fcc", 5.0842, 0.0, {90}, 0.0, "wolfram", "", "", true},
        {"Be", "hcp", 2.2858, 3.5843, {4}, 1440.0, "wolfram", "kittel2004", "", true},
        {"Mg", "hcp", 3.2094, 5.2108, {12}, 0.0, "wolfram", "", "", true},
        {"Sc", "hcp", 3.309, 5.2733, {21}, 0.0, "wolfram", "", "", true},
        {"Ti", "hcp", 2.9508, 4.6855, {22}, 420.0, "wolfram", "kittel2004", "alpha-Ti", true},
        {"Co", "hcp", 2.5071, 4.0695, {27}, 0.0, "wolfram", "", "alpha-Co", true},
        {"Zn", "hcp", 2.6649, 4.9468, {30}, 327.0, "wolfram", "kittel2004", "", true},
        {"Y", "hcp", 3.6474, 5.7306, {39}, 0.0, "wolfram", "", "", true},
        {"Zr", "hcp", 3.232, 5.147, {40}, 0.0, "wolfram", "", "alpha-Zr", true},
        {"Tc", "hcp", 2.735, 4.388, {43}, 0.0, "wolfram", "", "", true},
        {"Ru", "hcp", 2.7059, 4.2815, {44}, 0.0, "wolfram", "", "", true},
        {"Cd", "hcp", 2.9794, 5.6186, {48}, 209.0, "wolfram", "kittel2004", "", true},
        {"Gd", "hcp", 3.636, 5.7826, {64}, 0.0, "wolfram", "", "", true},
        {"Tb", "hcp", 3.601, 5.6936, {65}, 0.0, "wolfram", "", "", true},
        {"Dy", "hcp", 3.593, 5.6537, {66}, 0.0, "wolfram", "", "", true},
        {"Ho", "hcp", 3.5773, 5.6158, {67}, 0.0, "wolfram", "", "", true},
        {"Er", "hcp", 3.5588, 5.5874, {68}, 0.0, "wolfram", "", "", true},
        {"Tm", "hcp", 3.5375, 5.5546, {69}, 0.0, "wolfram", "", "", true},
        {"Lu", "hcp", 3.5031, 5.5509, {71}, 0.0, "wolfram", "", "", true},
        {"Hf", "hcp", 3.1964, 5.0511, {72}, 0.0, "wolfram", "", "alpha-Hf", true},
        {"Re", "hcp", 2.761, 4.456, {75}, 0.0, "wolfram", "", "", true},
        {"Os", "hcp", 2.7344, 4.3173, {76}, 0.0, "wolfram", "", "", true},
        {"Tl", "hcp", 3.4566, 5.5248, {81}, 0.0, "wolfram", "", "alpha-Tl", true},
        {"Po", "sc", 3.359, 0.0, {84}, 0.0, "wolfram", "", "alpha-Po", true},

        // --- Compounds: Wikipedia, "Lattice constant", lattice constants at
        // 300 K, with the reference each row gives.
        // Zinc blende (cation, anion).
        {"AlAs", "zincblende", 5.6605, 0.0, {13, 33}, 0.0, "ioffe", "", "", true},
        {"AlP", "zincblende", 5.4510, 0.0, {13, 15}, 0.0, "ioffe", "", "", true},
        {"InP", "zincblende", 5.869, 0.0, {49, 15}, 0.0, "ioffe", "", "", true},
        {"CdS", "zincblende", 5.8320, 0.0, {48, 16}, 0.0, "aps", "", "zinc blende (hawleyite) form", false},
        {"CdSe", "zincblende", 6.050, 0.0, {48, 34}, 0.0, "aps", "", "zinc blende form", true},
        {"CdTe", "zincblende", 6.482, 0.0, {48, 52}, 0.0, "aps", "", "", true},
        {"ZnS", "zincblende", 5.420, 0.0, {30, 16}, 0.0, "aps", "", "sphalerite", true},
        {"BN", "zincblende", 3.6150, 0.0, {5, 7}, 0.0, "aps", "", "cubic boron nitride (c-BN)", true},
        {"BP", "zincblende", 4.5380, 0.0, {5, 15}, 0.0, "aps", "", "", true},
        // Wurtzite (A, B); u = 3/8 (ideal), c/a from the table.
        {"ZnO", "wurtzite", 3.25, 5.2, {30, 8}, 0.0, "morkoc2009", "", "zincite", true},
        {"CdS", "wurtzite", 4.160, 6.756, {48, 16}, 0.0, "aps", "", "greenockite", true},
        {"ZnS", "wurtzite", 3.82, 6.26, {30, 16}, 0.0, "aps", "", "wurtzite (high-temperature) form", false},
        {"AlN", "wurtzite", 3.112, 4.982, {13, 7}, 0.0, "ioffe", "", "", true},
        {"GaN", "wurtzite", 3.189, 5.185, {31, 7}, 0.0, "ioffe", "", "", true},
        {"InN", "wurtzite", 3.533, 5.693, {49, 7}, 0.0, "ioffe", "", "", true},
        // Rock salt (cation, anion). (Not HfN: the table gives it NbN's
        // 4.392 A, which makes it 9 % denser than HfN is.)
        {"ZnO", "rocksalt", 4.580, 0.0, {30, 8}, 0.0, "aps", "", "high-pressure phase", false},
        {"PbTe", "rocksalt", 6.4620, 0.0, {82, 52}, 0.0, "aps", "", "", true},
        {"LiCl", "rocksalt", 5.14, 0.0, {3, 17}, 0.0, "wikipedia", "", "", true},
        {"LiBr", "rocksalt", 5.50, 0.0, {3, 35}, 0.0, "wikipedia", "", "", true},
        {"LiI", "rocksalt", 6.01, 0.0, {3, 53}, 0.0, "wikipedia", "", "", true},
        {"NaBr", "rocksalt", 5.97, 0.0, {11, 35}, 0.0, "wikipedia", "", "", true},
        {"KF", "rocksalt", 5.34, 0.0, {19, 9}, 0.0, "wikipedia", "", "", true},
        {"RbF", "rocksalt", 5.65, 0.0, {37, 9}, 0.0, "wikipedia", "", "", true},
        {"RbCl", "rocksalt", 6.59, 0.0, {37, 17}, 0.0, "wikipedia", "", "", true},
        {"RbI", "rocksalt", 7.35, 0.0, {37, 53}, 0.0, "wikipedia", "", "", true},
        {"CsF", "rocksalt", 6.02, 0.0, {55, 9}, 0.0, "wikipedia", "", "", true},
        {"TiN", "rocksalt", 4.249, 0.0, {22, 7}, 0.0, "wikipedia", "", "nitride", true},
        {"ZrN", "rocksalt", 4.577, 0.0, {40, 7}, 0.0, "wikipedia", "", "nitride", true},
        {"VN", "rocksalt", 4.136, 0.0, {23, 7}, 0.0, "wikipedia", "", "nitride", true},
        {"CrN", "rocksalt", 4.149, 0.0, {24, 7}, 0.0, "wikipedia", "", "nitride", true},
        {"NbN", "rocksalt", 4.392, 0.0, {41, 7}, 0.0, "wikipedia", "", "nitride", true},
        {"ScN", "rocksalt", 4.52, 0.0, {21, 7}, 0.0, "saha2010", "", "nitride", true},
        {"TiC", "rocksalt", 4.328, 0.0, {22, 6}, 0.0, "toth1967", "", "carbide", true},
        {"ZrC", "rocksalt", 4.698, 0.0, {40, 6}, 0.0, "toth1967", "", "carbide; measured on ZrC0.97", true},
        {"HfC", "rocksalt", 4.640, 0.0, {72, 6}, 0.0, "toth1967", "", "carbide; measured on HfC0.99", true},
        {"VC", "rocksalt", 4.166, 0.0, {23, 6}, 0.0, "toth1967", "", "carbide; measured on VC0.97", true},
        {"NbC", "rocksalt", 4.470, 0.0, {41, 6}, 0.0, "toth1967", "", "carbide; measured on NbC0.99", true},
        {"TaC", "rocksalt", 4.456, 0.0, {73, 6}, 0.0, "toth1967", "", "carbide; measured on TaC0.99", true},
        // Cesium chloride (A, B).
        {"CsCl", "cesium_chloride", 4.123, 0.0, {55, 17}, 0.0, "wikipedia", "", "", true},
        {"CsBr", "cesium_chloride", 4.291, 0.0, {55, 35}, 0.0, "wikipedia", "", "", true},
        // Cubic perovskites (A, B, X).
        {"SrTiO3", "perovskite", 3.98805, 0.0, {38, 22, 8}, 0.0, "landolt", "", "", true},
        {"KTaO3", "perovskite", 3.9885, 0.0, {19, 73, 8}, 0.0, "landolt", "", "", true},
        {"SrVO3", "perovskite", 3.838, 0.0, {38, 23, 8}, 0.0, "landolt", "", "", true},
        {"CaVO3", "perovskite", 3.767, 0.0, {20, 23, 8}, 0.0, "landolt", "", "", true},
    };
    r.insert(r.end(), more.begin(), more.end());
    return r;
  }();
  return rows;
}

namespace {

std::string lower(std::string_view s) {
  std::string out(s);
  for (char &c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return out;
}

//! Everything a query word may match, lower case, separated by spaces.
std::string haystack(const DatabaseCrystal &row) {
  std::string h = lower(row.name) + " " + lower(row.structure) + " " + lower(row.note);
  for (const auto &s : builtin_structures()) {
    if (s.name == row.structure) h += " " + lower(s.title);
  }
  for (unsigned z : row.sublattices) {
    const auto &e = Element::by_z(z);
    h += " " + lower(e.symbol()) + " " + lower(e.name());
  }
  // (Common spellings: "zinc blende" for zincblende, "rock salt", "cesium/caesium chloride".)
  if (row.structure == "zincblende") h += " zinc blende sphalerite";
  if (row.structure == "rocksalt") h += " rock salt halite";
  if (row.structure == "cesium_chloride") h += " cesium chloride caesium";
  return h;
}

} // namespace

std::vector<const DatabaseCrystal *> search_crystal_database(std::string_view query) {
  std::vector<std::string> words;
  std::string w;
  for (const char c : lower(query) + " ") {
    if (std::isspace(static_cast<unsigned char>(c)) || c == ',') {
      if (!w.empty()) words.push_back(w);
      w.clear();
    } else {
      w += c;
    }
  }
  std::vector<const DatabaseCrystal *> out;
  for (const auto &row : crystal_database()) {
    const std::string h = haystack(row);
    // A word that is a whole element symbol or row name matches only as a
    // word ("n" must not find every "nitride"... but "Ni" finds nickel).
    const bool all = std::all_of(words.begin(), words.end(), [&](const std::string &q) {
      if (q.size() <= 2) {
        std::string token;
        for (const char c : h + " ") {
          if (c == ' ') {
            if (token == q) return true;
            token.clear();
          } else {
            token += c;
          }
        }
        return false;
      }
      return h.find(q) != std::string::npos;
    });
    if (all) out.push_back(&row);
  }
  return out;
}

const DatabaseCrystal &find_database_crystal(std::string_view name) {
  const std::string wanted = lower(name);
  const DatabaseCrystal *named = nullptr;
  for (const auto &row : crystal_database()) {
    if (lower(row.name) != wanted) continue;
    if (row.preferred) return row;
    if (named == nullptr) named = &row;
  }
  if (named != nullptr) return *named;
  const auto found = search_crystal_database(name);
  if (found.size() == 1) return *found.front();
  std::string list;
  for (std::size_t i = 0; i < found.size() && i < 12; ++i) {
    list += (i ? ", " : "") + std::string(found[i]->name) + " (" + std::string(found[i]->structure) + ")";
  }
  if (found.empty()) {
    throw std::invalid_argument("'" + std::string(name) + "' is not in the crystal database (amcset crystals lists it)");
  }
  throw std::invalid_argument("'" + std::string(name) + "' matches " + std::to_string(found.size()) +
                              " crystals of the database: " + list + (found.size() > 12 ? ", ..." : "") +
                              "; name one");
}

std::vector<std::pair<unsigned, double>> database_composition(const DatabaseCrystal &row) {
  const Crystal c = make_database_crystal(row);
  const auto m = c.cell.composition();
  return {m.begin(), m.end()};
}

const DatabaseCrystal *match_composition(const std::vector<std::pair<unsigned, double>> &composition,
                                         double tolerance) {
  std::map<unsigned, double> want;
  double total = 0.0;
  for (const auto &[z, amount] : composition) {
    if (amount > 0.0) {
      want[z] += amount;
      total += amount;
    }
  }
  if (want.empty() || !(total > 0.0)) return nullptr;
  const DatabaseCrystal *second = nullptr;
  for (const auto &row : crystal_database()) {
    const auto have = database_composition(row);
    if (have.size() != want.size()) continue;
    double cell_total = 0.0;
    for (const auto &[z, atoms] : have) cell_total += atoms;
    bool same = true;
    for (const auto &[z, atoms] : have) {
      const auto it = want.find(z);
      if (it == want.end() || std::abs(it->second / total - atoms / cell_total) > tolerance) {
        same = false;
        break;
      }
    }
    if (!same) continue;
    if (row.preferred) return &row;
    if (second == nullptr) second = &row;
  }
  return second;
}

Crystal make_database_crystal(const DatabaseCrystal &row) {
  StructureSpec spec;
  spec.name = std::string(row.structure);
  spec.a_A = row.a_A;
  if (row.c_A > 0.0) spec.c_over_a = row.c_A / row.a_A;
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
