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

#include "crystal/unit_cell.h"

#include "element_data.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace amcset {
namespace crystal {

using namespace common;

namespace {

std::string num(double v) {
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

double angle_between(const Vec3 &u, const Vec3 &v) {
  const double c = dot(u, v) / (norm(u) * norm(v));
  return std::acos(std::clamp(c, -1.0, 1.0));
}

} // namespace

// ---------------------------------------------------------------------------
// Lattice
// ---------------------------------------------------------------------------

Lattice::Lattice(const Mat3 &m) : m_(m), inv_t_(m.transposed().inverse()) {}

Lattice Lattice::from_vectors_A(const Mat3 &rows) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (!std::isfinite(rows[i][j])) throw std::invalid_argument("lattice vectors: not finite numbers");
    }
  }
  const double v = rows.determinant();
  const double scale = norm(rows[0]) * norm(rows[1]) * norm(rows[2]);
  if (!(scale > 0.0) || std::abs(v) < 1e-9 * scale) {
    throw std::invalid_argument("lattice vectors: degenerate (zero volume)");
  }
  if (v < 0.0) {
    throw std::invalid_argument("lattice vectors: left-handed (a . (b x c) < 0); swap two vectors");
  }
  return Lattice(rows);
}

Lattice Lattice::from_parameters(length_quantity a, length_quantity b, length_quantity c,
                                 angle_quantity alpha, angle_quantity beta, angle_quantity gamma) {
  const double aa = in_units(a, angstrom);
  const double bb = in_units(b, angstrom);
  const double cc = in_units(c, angstrom);
  if (!(aa > 0.0) || !(bb > 0.0) || !(cc > 0.0) || !std::isfinite(aa * bb * cc)) {
    throw std::invalid_argument("lattice: a, b and c must be positive, got " + num(aa) + ", " + num(bb) +
                                ", " + num(cc) + " A");
  }
  const double al = in_units(alpha, radian);
  const double be = in_units(beta, radian);
  const double ga = in_units(gamma, radian);
  for (double x : {al, be, ga}) {
    if (!(x > 0.0) || !(x < pi_v)) throw std::invalid_argument("lattice: angles must lie between 0 and 180 deg");
  }
  const double ca = std::cos(al);
  const double cb = std::cos(be);
  const double cg = std::cos(ga);
  const double sg = std::sin(ga);
  const double cy = (ca - cb * cg) / sg;
  const double rad = 1.0 - cb * cb - cy * cy;
  if (!(rad > 1e-12)) {
    throw std::invalid_argument("lattice: the angles alpha, beta, gamma do not form a cell (zero volume)");
  }
  // Cosines of exact right angles are ~6e-17: make them 0 so that cubic
  // cells have exactly axis-aligned vectors.
  const auto tidy = [](double x) { return std::abs(x) < 1e-15 ? 0.0 : x; };
  return Lattice(Mat3::from_rows({aa, 0.0, 0.0}, {tidy(bb * cg), bb * sg, 0.0},
                                 {tidy(cc * cb), tidy(cc * cy), cc * std::sqrt(rad)}));
}

angle_quantity Lattice::alpha() const noexcept { return angle_between(m_[1], m_[2]) * radian; }
angle_quantity Lattice::beta() const noexcept { return angle_between(m_[0], m_[2]) * radian; }
angle_quantity Lattice::gamma() const noexcept { return angle_between(m_[0], m_[1]) * radian; }

bool Lattice::is_hexagonal() const noexcept {
  const double a = norm(m_[0]);
  const double b = norm(m_[1]);
  const double tol = 1e-6;
  return std::abs(a - b) <= tol * a && std::abs(in_units(alpha(), radian) - pi_v / 2) < tol &&
         std::abs(in_units(beta(), radian) - pi_v / 2) < tol &&
         std::abs(in_units(gamma(), radian) - 2.0 * pi_v / 3.0) < tol;
}

Lattice Lattice::scaled(double s) const {
  if (!(s > 0.0) || !std::isfinite(s)) throw std::invalid_argument("lattice: scale factor must be positive");
  return Lattice(Mat3::from_rows(m_[0] * s, m_[1] * s, m_[2] * s));
}

// ---------------------------------------------------------------------------
// Sites and cells
// ---------------------------------------------------------------------------

std::string to_string(SiteKind k) { return k == SiteKind::lattice ? "lattice" : "interstitial"; }

SiteKind parse_site_kind(std::string_view s) {
  if (s == "lattice" || s == "substitutional") return SiteKind::lattice;
  if (s == "interstitial") return SiteKind::interstitial;
  throw std::invalid_argument("expected lattice (substitutional) or interstitial, got '" + std::string(s) + "'");
}

double Site::occupancy() const noexcept {
  double s = 0.0;
  for (const auto &o : occupants) s += o.fraction;
  return s;
}

double wrap_fraction(double f) noexcept {
  double w = f - std::floor(f);
  if (w >= 1.0 - 1e-9 || w < 1e-12) w = 0.0; // 0.9999999999 and -1e-17 are 0
  return w;
}

Vec3 wrap_fractional(const Vec3 &f) noexcept {
  return {wrap_fraction(f.x), wrap_fraction(f.y), wrap_fraction(f.z)};
}

std::map<unsigned, double> UnitCell::composition() const {
  std::map<unsigned, double> c;
  for (const auto &s : sites) {
    for (const auto &o : s.occupants) c[o.z] += o.fraction;
  }
  return c;
}

double UnitCell::atoms_per_cell() const {
  double n = 0.0;
  for (const auto &s : sites) n += s.occupancy();
  return n;
}

number_density_quantity UnitCell::atomic_density() const {
  return atoms_per_cell() / lattice.volume_A3() * per_cubic_angstrom;
}

mass_density_quantity UnitCell::density(const std::map<unsigned, mass_quantity> &mass_of) const {
  mass_quantity m = 0.0 * kilogram;
  for (const auto &[z, n] : composition()) {
    const auto it = mass_of.find(z);
    const mass_quantity mz = it != mass_of.end() ? it->second : Element::by_z(z).natural_mass();
    m += n * mz;
  }
  return m / lattice.volume();
}

void UnitCell::scale_to_density(mass_density_quantity rho, const std::map<unsigned, mass_quantity> &mass_of) {
  if (!(rho > 0.0 * gram_per_cubic_centimeter)) throw std::invalid_argument("crystal: density must be positive");
  const auto now = density(mass_of);
  if (!(now > 0.0 * gram_per_cubic_centimeter)) {
    throw std::invalid_argument("crystal: a cell without atoms has no density to scale to");
  }
  lattice = lattice.scaled(std::cbrt((now / rho).value()));
}

double UnitCell::distance_A(const Vec3 &f1, const Vec3 &f2) const {
  Vec3 d = f2 - f1;
  d = {d.x - std::round(d.x), d.y - std::round(d.y), d.z - std::round(d.z)};
  // Rounding gives the minimum image for orthogonal cells; check the
  // neighbouring images for skewed ones.
  double best = std::numeric_limits<double>::max();
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        const Vec3 r = lattice.to_cartesian_A(d + Vec3{double(i), double(j), double(k)});
        best = std::min(best, norm(r));
      }
    }
  }
  return best;
}

std::vector<NeighbourShell> UnitCell::neighbour_shells(const Vec3 &f, double max_A, double tol_A,
                                                       std::optional<SiteKind> only) const {
  // Images needed in each direction: the spacing of the lattice planes
  // normal to a*_i is 1/|a*_i|.
  const Mat3 rec = lattice.reciprocal_per_A();
  int n[3];
  for (int i = 0; i < 3; ++i) n[i] = static_cast<int>(std::ceil(max_A * norm(rec[i]))) + 1;
  std::vector<std::pair<double, double>> d; // (distance, occupancy)
  for (const auto &s : sites) {
    if (only && s.kind != *only) continue;
    for (int i = -n[0]; i <= n[0]; ++i) {
      for (int j = -n[1]; j <= n[1]; ++j) {
        for (int k = -n[2]; k <= n[2]; ++k) {
          const Vec3 r = lattice.to_cartesian_A(s.position + Vec3{double(i), double(j), double(k)} - f);
          const double dist = norm(r);
          if (dist > tol_A && dist <= max_A) d.emplace_back(dist, s.occupancy());
        }
      }
    }
  }
  std::sort(d.begin(), d.end());
  std::vector<NeighbourShell> shells;
  for (const auto &[dist, occ] : d) {
    if (shells.empty() || dist - shells.back().distance_A > tol_A) {
      shells.push_back({dist, 0.0, 0.0});
    }
    shells.back().count += 1.0;
    shells.back().occupied += occ;
  }
  return shells;
}

void UnitCell::validate(const std::string &where, double min_distance_A) const {
  if (!(lattice.volume_A3() > 0.0)) throw std::invalid_argument(where + ": the lattice has no volume");
  if (sites.empty()) throw std::invalid_argument(where + ": the cell has no sites");
  double atoms = 0.0;
  for (std::size_t i = 0; i < sites.size(); ++i) {
    const auto &s = sites[i];
    const std::string sw = where + ".sites." + std::to_string(i) + (s.label.empty() ? "" : " (" + s.label + ")");
    for (int c = 0; c < 3; ++c) {
      if (!(s.position[c] >= 0.0 && s.position[c] < 1.0)) {
        throw std::invalid_argument(sw + ": fractional coordinates must lie in [0, 1)");
      }
    }
    double sum = 0.0;
    for (const auto &o : s.occupants) {
      if (o.z < 1 || o.z > Element::max_z) throw std::invalid_argument(sw + ": unknown element Z = " + std::to_string(o.z));
      if (!(o.fraction > 0.0 && o.fraction <= 1.0)) {
        throw std::invalid_argument(sw + ": occupancy fractions must lie in (0, 1], got " + num(o.fraction));
      }
      sum += o.fraction;
    }
    if (sum > 1.0 + 1e-6) {
      throw std::invalid_argument(sw + ": occupancies add up to " + num(sum) + ", more than 1");
    }
    atoms += sum;
  }
  if (!(atoms > 0.0)) throw std::invalid_argument(where + ": no site is occupied");
  for (std::size_t i = 0; i < sites.size(); ++i) {
    for (std::size_t j = i + 1; j < sites.size(); ++j) {
      const double d = distance_A(sites[i].position, sites[j].position);
      if (d < min_distance_A) {
        throw std::invalid_argument(where + ": sites " + std::to_string(i) + " (" + sites[i].label + ") and " +
                                    std::to_string(j) + " (" + sites[j].label + ") are " + num(d) +
                                    " A apart (closer than " + num(min_distance_A) +
                                    " A); merge them into one site with several occupants");
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Built-in structures
// ---------------------------------------------------------------------------

namespace {

using Positions = std::vector<Vec3>;

const Positions fcc_translations{Vec3{0.0, 0.0, 0.0}, Vec3{0.5, 0.5, 0.0}, Vec3{0.5, 0.0, 0.5},
                                 Vec3{0.0, 0.5, 0.5}};

Positions shifted(const Positions &p, const Vec3 &t) {
  Positions out;
  for (const auto &x : p) out.push_back(wrap_fractional(x + t));
  return out;
}

Positions corners_of_quarter_cube() {
  // The eight (1/4|3/4)^3 positions: the tetrahedral holes of fcc.
  Positions out;
  for (double x : {0.25, 0.75}) {
    for (double y : {0.25, 0.75}) {
      for (double z : {0.25, 0.75}) out.push_back(Vec3{x, y, z});
    }
  }
  return out;
}

struct Layout {
  std::vector<Positions> sublattices;
  std::map<std::string, Positions, std::less<>> interstitials;
};

Positions joined(Positions a, const Positions &b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

Layout layout(std::string_view name, double u) {
  const Vec3 origin{0.0, 0.0, 0.0};
  const Vec3 centre{0.5, 0.5, 0.5};
  const Vec3 quarter{0.25, 0.25, 0.25};
  const Positions face_centres{Vec3{0.5, 0.5, 0.0}, Vec3{0.5, 0.0, 0.5}, Vec3{0.0, 0.5, 0.5}};
  const Positions edge_centres{Vec3{0.5, 0.0, 0.0}, Vec3{0.0, 0.5, 0.0}, Vec3{0.0, 0.0, 0.5}};
  // Octahedral holes of fcc: the cube centre and the edge centres.
  const Positions octa_fcc = joined(Positions{centre}, edge_centres);
  // Tetrahedral holes of diamond and zinc blende: an fcc lattice at
  // (1/2,1/2,1/2) and one at (3/4,3/4,3/4).
  const Positions tetra_diamond =
      joined(shifted(fcc_translations, centre), shifted(fcc_translations, Vec3{0.75, 0.75, 0.75}));
  Layout l;
  if (name == "sc") {
    l.sublattices.push_back(Positions{origin});
  } else if (name == "bcc") {
    l.sublattices.push_back(Positions{origin, centre});
    // Octahedral: face and edge centres (6 per cell). Tetrahedral:
    // (1/2, 1/4, 0) and its relatives on the cube faces (12 per cell).
    l.interstitials["octahedral"] = joined(face_centres, edge_centres);
    Positions tet;
    for (double a : {0.25, 0.75}) {
      tet.push_back(Vec3{0.5, a, 0.0});
      tet.push_back(Vec3{a, 0.5, 0.0});
      tet.push_back(Vec3{0.0, 0.5, a});
      tet.push_back(Vec3{0.0, a, 0.5});
      tet.push_back(Vec3{0.5, 0.0, a});
      tet.push_back(Vec3{a, 0.0, 0.5});
    }
    l.interstitials["tetrahedral"] = tet;
  } else if (name == "fcc") {
    l.sublattices.push_back(fcc_translations);
    l.interstitials["octahedral"] = octa_fcc;
    l.interstitials["tetrahedral"] = corners_of_quarter_cube();
  } else if (name == "hcp") {
    // Wyckoff 2c of P6_3/mmc; octahedral holes 2a; tetrahedral holes 4f
    // with z = 5/8 (exact for the ideal c/a).
    l.sublattices.push_back(Positions{Vec3{1.0 / 3, 2.0 / 3, 0.25}, Vec3{2.0 / 3, 1.0 / 3, 0.75}});
    l.interstitials["octahedral"] = Positions{origin, Vec3{0.0, 0.0, 0.5}};
    l.interstitials["tetrahedral"] = Positions{Vec3{1.0 / 3, 2.0 / 3, 0.625}, Vec3{1.0 / 3, 2.0 / 3, 0.875},
                                               Vec3{2.0 / 3, 1.0 / 3, 0.125}, Vec3{2.0 / 3, 1.0 / 3, 0.375}};
  } else if (name == "diamond") {
    l.sublattices.push_back(joined(fcc_translations, shifted(fcc_translations, quarter)));
    l.interstitials["tetrahedral"] = tetra_diamond;
  } else if (name == "zincblende") {
    l.sublattices.push_back(fcc_translations);
    l.sublattices.push_back(shifted(fcc_translations, quarter));
    l.interstitials["tetrahedral"] = tetra_diamond;
  } else if (name == "wurtzite") {
    l.sublattices.push_back(Positions{Vec3{1.0 / 3, 2.0 / 3, 0.0}, Vec3{2.0 / 3, 1.0 / 3, 0.5}});
    l.sublattices.push_back(Positions{wrap_fractional(Vec3{1.0 / 3, 2.0 / 3, u}),
                                      wrap_fractional(Vec3{2.0 / 3, 1.0 / 3, 0.5 + u})});
  } else if (name == "rocksalt") {
    l.sublattices.push_back(fcc_translations);
    l.sublattices.push_back(shifted(fcc_translations, Vec3{0.5, 0.0, 0.0}));
    l.interstitials["tetrahedral"] = corners_of_quarter_cube();
  } else if (name == "cesium_chloride") {
    l.sublattices.push_back(Positions{origin});
    l.sublattices.push_back(Positions{centre});
  } else if (name == "fluorite") {
    l.sublattices.push_back(fcc_translations);
    l.sublattices.push_back(corners_of_quarter_cube());
    l.interstitials["octahedral"] = octa_fcc;
  } else if (name == "perovskite") {
    l.sublattices.push_back(Positions{origin});
    l.sublattices.push_back(Positions{centre});
    l.sublattices.push_back(face_centres);
  } else if (name == "L1_2") {
    l.sublattices.push_back(Positions{origin});
    l.sublattices.push_back(face_centres);
    l.interstitials["octahedral"] = octa_fcc;
  } else if (name == "L1_0") {
    l.sublattices.push_back(Positions{origin, Vec3{0.5, 0.5, 0.0}});
    l.sublattices.push_back(Positions{Vec3{0.5, 0.0, 0.5}, Vec3{0.0, 0.5, 0.5}});
  } else {
    throw std::logic_error("crystal: no layout for structure " + std::string(name));
  }
  return l;
}

} // namespace

const std::vector<StructureInfo> &builtin_structures() {
  static const std::vector<StructureInfo> list{
      {"sc", "simple cubic", "Po", {"atoms"}, false, false, false, {}},
      {"bcc", "body-centred cubic (A2)", "Fe, W, Mo, Cr, V, Nb, Ta", {"atoms"}, false, false, false,
       {"octahedral", "tetrahedral"}},
      {"fcc", "face-centred cubic (A1)", "Cu, Ni, Al, Au, Ag, Pt, Pd, austenite", {"atoms"}, false, false, false,
       {"octahedral", "tetrahedral"}},
      {"hcp", "hexagonal close-packed (A3)", "Ti, Zr, Mg, Be, Zn, Co", {"atoms"}, true, false, false,
       {"octahedral", "tetrahedral"}},
      {"diamond", "diamond (A4)", "Si, Ge, C", {"atoms"}, false, false, false, {"tetrahedral"}},
      {"zincblende", "zinc blende (B3)", "GaAs, InP, 3C-SiC, ZnS", {"A (0,0,0)", "B (1/4,1/4,1/4)"}, false, false,
       false, {"tetrahedral"}},
      {"wurtzite", "wurtzite (B4)", "GaN, ZnO, AlN, 2H-SiC", {"A (1/3,2/3,0)", "B (1/3,2/3,u)"}, true, false, true,
       {}},
      {"rocksalt", "rock salt (B1)", "NaCl, MgO, NiO, UC, TiC", {"A (0,0,0)", "B (1/2,0,0)"}, false, false, false,
       {"tetrahedral"}},
      {"cesium_chloride", "cesium chloride (B2)", "CsCl, NiAl, FeAl, CuZn", {"A (0,0,0)", "B (1/2,1/2,1/2)"}, false,
       false, false, {}},
      {"fluorite", "fluorite (C1)", "CaF2, UO2, ThO2, CeO2", {"A (0,0,0) cation", "B (1/4,1/4,1/4) anion"}, false,
       false, false, {"octahedral"}},
      {"perovskite", "cubic perovskite (E2_1)", "SrTiO3, BaTiO3 (cubic)",
       {"A (0,0,0)", "B (1/2,1/2,1/2)", "X (1/2,1/2,0)"}, false, false, false, {}},
      {"L1_2", "L1_2 (Cu3Au)", "Ni3Al, Cu3Au, Pt3Fe", {"A (0,0,0)", "B (1/2,1/2,0)"}, false, false, false,
       {"octahedral"}},
      {"L1_0", "L1_0 (CuAu)", "FePt, TiAl, CuAu", {"A (0,0,0)", "B (1/2,0,1/2)"}, false, true, false, {}},
  };
  return list;
}

const StructureInfo &structure_info(std::string_view name) {
  for (const auto &s : builtin_structures()) {
    if (s.name == name) return s;
  }
  std::string known;
  for (const auto &s : builtin_structures()) known += (known.empty() ? "" : ", ") + std::string(s.name);
  throw std::invalid_argument("unknown crystal structure '" + std::string(name) + "' (built in: " + known + ")");
}

UnitCell build_structure(const StructureSpec &spec) {
  const StructureInfo &info = structure_info(spec.name);
  if (spec.sublattices.size() != info.sublattices.size()) {
    std::string names;
    for (auto n : info.sublattices) names += (names.empty() ? "" : ", ") + std::string(n);
    throw std::invalid_argument("structure " + spec.name + " has " + std::to_string(info.sublattices.size()) +
                                " sublattice(s) (" + names + "), got " + std::to_string(spec.sublattices.size()));
  }
  if (info.has_u && !(spec.u > 0.0 && spec.u < 1.0)) {
    throw std::invalid_argument("structure " + spec.name + ": u must lie between 0 and 1");
  }
  double c_over_a = spec.c_over_a;
  if (c_over_a == 0.0) c_over_a = info.hexagonal ? std::sqrt(8.0 / 3.0) : 1.0;
  if (!(c_over_a > 0.0) || !std::isfinite(c_over_a)) throw std::invalid_argument("structure " + spec.name + ": c/a must be positive");
  if (!(spec.a_A >= 0.0) || !std::isfinite(spec.a_A)) throw std::invalid_argument("structure " + spec.name + ": a must be positive");
  const double a = spec.a_A > 0.0 ? spec.a_A : 1.0;
  const bool uses_c = info.hexagonal || info.tetragonal;
  UnitCell cell;
  cell.name = std::string(info.title);
  cell.lattice = Lattice::from_parameters(a * angstrom, a * angstrom, (uses_c ? c_over_a : 1.0) * a * angstrom,
                                          90.0 * degree, 90.0 * degree,
                                          (info.hexagonal ? 120.0 : 90.0) * degree);
  const Layout l = layout(info.name, spec.u);
  // Labels: the element and a running number (Fe1, Fe2; Ga1, As1); a
  // sublattice shared by several elements (or empty) gets its letter (A1).
  std::map<std::string, int> numbers;
  for (std::size_t s = 0; s < l.sublattices.size(); ++s) {
    const auto &occupants = spec.sublattices[s];
    const std::string prefix = occupants.size() == 1 && occupants.front().z >= 1 &&
                                       occupants.front().z <= Element::max_z
                                   ? std::string(Element::by_z(occupants.front().z).symbol())
                                   : std::string(1, static_cast<char>('A' + s));
    int &n = numbers[prefix];
    for (const auto &p : l.sublattices[s]) {
      Site site;
      site.label = prefix + std::to_string(++n);
      site.position = wrap_fractional(p);
      site.kind = SiteKind::lattice;
      site.occupants = spec.sublattices[s];
      cell.sites.push_back(std::move(site));
    }
  }
  for (const auto &is : spec.interstitials) {
    const auto it = l.interstitials.find(is.sites);
    if (it == l.interstitials.end()) {
      std::string names;
      for (auto n : info.interstitials) names += (names.empty() ? "" : ", ") + std::string(n);
      throw std::invalid_argument("structure " + spec.name + " has no built-in '" + is.sites + "' interstitial sites" +
                                  (names.empty() ? std::string(" (none built in; add sites explicitly)")
                                                 : " (built in: " + names + ")"));
    }
    const std::string prefix = is.sites.substr(0, 3);
    int n = 0;
    for (const auto &p : it->second) {
      Site site;
      site.label = prefix + std::to_string(++n);
      site.position = wrap_fractional(p);
      site.kind = SiteKind::interstitial;
      site.occupants = is.occupants;
      cell.sites.push_back(std::move(site));
    }
  }
  return cell;
}

} // namespace crystal
} // namespace amcset
