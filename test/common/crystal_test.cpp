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

/*! \file crystal_test.cpp
 *  \brief The crystal library: lattices, built-in structures, orientations,
 *  thermal vibrations and cell files. References are textbook geometry
 *  (coordination numbers and distances) and hand calculations.
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>

#include "crystal/cell_io.h"
#include "crystal/crystal.h"
#include "crystal/crystal_database.h"
#include "crystal/grain_io.h"
#include "crystal/grain_map.h"
#include "crystal/lattice_search.h"
#include "crystal/materials.h"
#include "crystal/orientation.h"
#include "crystal/unit_cell.h"
#include "element_data.h"

using namespace amcset;
using namespace amcset::common;
using namespace amcset::crystal;

namespace {

const double deg = pi_v / 180.0;

UnitCell simple(const std::string &structure, std::vector<unsigned> z, double a_A = 4.0, double c_over_a = 0.0) {
  StructureSpec s;
  s.name = structure;
  s.a_A = a_A;
  s.c_over_a = c_over_a;
  for (unsigned zi : z) s.sublattices.push_back({Occupant{zi, 1.0}});
  return build_structure(s);
}

double angle_deg(const Vec3 &a, const Vec3 &b) {
  return std::acos(std::clamp(dot(a, b) / (norm(a) * norm(b)), -1.0, 1.0)) / deg;
}

void expect_vec(const Vec3 &a, const Vec3 &b, double tol = 1e-12) {
  EXPECT_NEAR(a.x, b.x, tol);
  EXPECT_NEAR(a.y, b.y, tol);
  EXPECT_NEAR(a.z, b.z, tol);
}

//! Same lattice parameters, same sites (position, kind, occupants) in any order.
void expect_same_cell(const UnitCell &a, const UnitCell &b, double tol_A = 1e-6) {
  EXPECT_NEAR(in_units(a.lattice.a(), angstrom), in_units(b.lattice.a(), angstrom), tol_A);
  EXPECT_NEAR(in_units(a.lattice.b(), angstrom), in_units(b.lattice.b(), angstrom), tol_A);
  EXPECT_NEAR(in_units(a.lattice.c(), angstrom), in_units(b.lattice.c(), angstrom), tol_A);
  EXPECT_NEAR(in_units(a.lattice.alpha(), degree), in_units(b.lattice.alpha(), degree), 1e-5);
  EXPECT_NEAR(in_units(a.lattice.beta(), degree), in_units(b.lattice.beta(), degree), 1e-5);
  EXPECT_NEAR(in_units(a.lattice.gamma(), degree), in_units(b.lattice.gamma(), degree), 1e-5);
  ASSERT_EQ(a.sites.size(), b.sites.size());
  for (const auto &sa : a.sites) {
    bool found = false;
    for (const auto &sb : b.sites) {
      if (b.distance_A(sa.position, sb.position) > tol_A) continue;
      found = true;
      EXPECT_EQ(sa.kind, sb.kind) << sa.label;
      ASSERT_EQ(sa.occupants.size(), sb.occupants.size()) << sa.label;
      for (const auto &oa : sa.occupants) {
        bool has = false;
        for (const auto &ob : sb.occupants) has = has || (oa.z == ob.z && std::abs(oa.fraction - ob.fraction) < 1e-7);
        EXPECT_TRUE(has) << sa.label << " Z=" << oa.z;
      }
    }
    EXPECT_TRUE(found) << sa.label;
  }
}

} // namespace

// ---------------------------------------------------------------------------
// Geometry and lattices
// ---------------------------------------------------------------------------

TEST(CrystalGeometry, VectorsAndMatrices) {
  const Vec3 a{1, 2, 3}, b{-2, 0.5, 4};
  EXPECT_DOUBLE_EQ(dot(a, b), -2 + 1 + 12);
  expect_vec(cross(a, b), {2 * 4 - 3 * 0.5, 3 * -2 - 1 * 4, 1 * 0.5 - 2 * -2});
  EXPECT_NEAR(dot(cross(a, b), a), 0.0, 1e-12);
  const Mat3 m = Mat3::from_rows({2, 1, 0}, {0, 3, 1}, {1, 0, 4});
  EXPECT_DOUBLE_EQ(m.determinant(), 2 * 12 - 1 * (0 - 1) + 0); // 25
  const Mat3 id = m * m.inverse();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) EXPECT_NEAR(id[i][j], i == j ? 1.0 : 0.0, 1e-14);
  }
  expect_vec(m * Vec3{1, 1, 1}, {3, 4, 5});
  expect_vec(m.transposed() * Vec3{1, 0, 0}, {2, 1, 0});
  // A rotation by 90 degrees about z takes x to y; about x it takes y to z.
  expect_vec(rotation_about({0, 0, 1}, 90 * deg) * Vec3{1, 0, 0}, {0, 1, 0}, 1e-15);
  expect_vec(rotation_about({1, 0, 0}, 90 * deg) * Vec3{0, 1, 0}, {0, 0, 1}, 1e-15);
  EXPECT_LT(rotation_error(rotation_about({1, 2, 3}, 0.7)), 1e-14);
  EXPECT_GT(rotation_error(m), 0.1);
}

TEST(CrystalLattice, ParametersVectorsAndVolume) {
  const auto cubic = Lattice::from_parameters(4.0 * angstrom, 4.0 * angstrom, 4.0 * angstrom, 90.0 * degree, 90.0 * degree,
                                              90.0 * degree);
  EXPECT_DOUBLE_EQ(cubic.volume_A3(), 64.0);
  expect_vec(cubic.vectors_A()[1], {0, 4, 0}, 0.0); // exactly axis-aligned
  expect_vec(cubic.vectors_A()[2], {0, 0, 4}, 0.0);
  EXPECT_FALSE(cubic.is_hexagonal());

  const auto hex = Lattice::from_parameters(3.0 * angstrom, 3.0 * angstrom, 5.0 * angstrom, 90.0 * degree, 90.0 * degree,
                                            120.0 * degree);
  EXPECT_NEAR(hex.volume_A3(), 9.0 * 5.0 * std::sqrt(3.0) / 2.0, 1e-12); // a^2 c sin 120
  EXPECT_TRUE(hex.is_hexagonal());
  EXPECT_NEAR(in_units(hex.gamma(), degree), 120.0, 1e-12);

  // Triclinic: V = abc sqrt(1 - cos^2 a - cos^2 b - cos^2 g + 2 cos a cos b cos g).
  const double al = 80 * deg, be = 95 * deg, ga = 105 * deg;
  const auto tri = Lattice::from_parameters(3.0 * angstrom, 4.0 * angstrom, 5.0 * angstrom, al * radian, be * radian,
                                            ga * radian);
  const double v = 60.0 * std::sqrt(1 - std::cos(al) * std::cos(al) - std::cos(be) * std::cos(be) -
                                    std::cos(ga) * std::cos(ga) + 2 * std::cos(al) * std::cos(be) * std::cos(ga));
  EXPECT_NEAR(tri.volume_A3(), v, 1e-12);
  EXPECT_NEAR(in_units(tri.a(), angstrom), 3.0, 1e-13);
  EXPECT_NEAR(in_units(tri.b(), angstrom), 4.0, 1e-13);
  EXPECT_NEAR(in_units(tri.c(), angstrom), 5.0, 1e-13);
  EXPECT_NEAR(in_units(tri.alpha(), radian), al, 1e-13);
  EXPECT_NEAR(in_units(tri.beta(), radian), be, 1e-13);
  EXPECT_NEAR(in_units(tri.gamma(), radian), ga, 1e-13);
  // Fractional <-> Cartesian, and a_i . a*_j = delta_ij.
  const Vec3 f{0.1, 0.7, 0.35};
  expect_vec(tri.to_fractional(tri.to_cartesian_A(f)), f, 1e-14);
  const Mat3 rec = tri.reciprocal_per_A();
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) EXPECT_NEAR(dot(tri.vectors_A()[i], rec[j]), i == j ? 1.0 : 0.0, 1e-14);
  }
  EXPECT_NEAR(tri.scaled(2.0).volume_A3(), 8.0 * v, 1e-10);
}

TEST(CrystalLattice, RejectsImpossibleCells) {
  const auto a = 3.0 * angstrom;
  EXPECT_THROW(Lattice::from_parameters(0.0 * angstrom, a, a, 90.0 * degree, 90.0 * degree, 90.0 * degree),
               std::invalid_argument);
  EXPECT_THROW(Lattice::from_parameters(a, a, a, 180.0 * degree, 90.0 * degree, 90.0 * degree), std::invalid_argument);
  // 30 + 30 < 120: the three angles cannot meet.
  EXPECT_THROW(Lattice::from_parameters(a, a, a, 30.0 * degree, 30.0 * degree, 120.0 * degree), std::invalid_argument);
  EXPECT_THROW(Lattice::from_vectors_A(Mat3::from_rows({1, 0, 0}, {0, 1, 0}, {1, 1, 0})), std::invalid_argument);
  EXPECT_THROW(Lattice::from_vectors_A(Mat3::from_rows({1, 0, 0}, {0, 0, 1}, {0, 1, 0})), std::invalid_argument); // left-handed
  EXPECT_NO_THROW(Lattice::from_vectors_A(Mat3::from_rows({0, 2, 2}, {2, 0, 2}, {2, 2, 0}))); // fcc primitive
  EXPECT_DOUBLE_EQ(wrap_fraction(1.25), 0.25);
  EXPECT_DOUBLE_EQ(wrap_fraction(-0.25), 0.75);
  EXPECT_DOUBLE_EQ(wrap_fraction(1.0 - 1e-12), 0.0);
  EXPECT_DOUBLE_EQ(wrap_fraction(-1e-17), 0.0);
}

// ---------------------------------------------------------------------------
// Built-in structures
// ---------------------------------------------------------------------------

TEST(CrystalStructures, AtomsPerCellAndNeighbours) {
  // Atoms per conventional cell; nearest-neighbour count and distance (in a)
  // seen from the first site; all textbook values.
  struct Row {
    const char *name;
    std::vector<unsigned> z;
    double atoms;
    double nn_count;
    double nn_distance_a;
  };
  const double ideal = std::sqrt(8.0 / 3.0);
  const std::vector<Row> rows{
      {"sc", {84}, 1, 6, 1.0},
      {"bcc", {26}, 2, 8, std::sqrt(3.0) / 2},
      {"fcc", {29}, 4, 12, 1 / std::sqrt(2.0)},
      {"hcp", {22}, 2, 12, 1.0},
      {"diamond", {14}, 8, 4, std::sqrt(3.0) / 4},
      {"zincblende", {31, 33}, 8, 4, std::sqrt(3.0) / 4},
      {"wurtzite", {31, 7}, 4, 4, 0.375 * ideal}, // u c = 3c/8 = a sqrt(3/8)
      {"rocksalt", {11, 17}, 8, 6, 0.5},
      {"cesium_chloride", {55, 17}, 2, 8, std::sqrt(3.0) / 2},
      {"fluorite", {92, 8}, 12, 8, std::sqrt(3.0) / 4},
      {"perovskite", {38, 22, 8}, 5, 12, 1 / std::sqrt(2.0)}, // A site: 12 X
      {"L1_2", {79, 29}, 4, 12, 1 / std::sqrt(2.0)},
      {"L1_0", {29, 79}, 4, 12, 1 / std::sqrt(2.0)}, // c/a = 1
  };
  ASSERT_EQ(builtin_structures().size(), rows.size());
  const double a = 4.0;
  for (const auto &r : rows) {
    const auto cell = simple(r.name, r.z, a);
    EXPECT_NO_THROW(cell.validate(r.name)) << r.name;
    EXPECT_DOUBLE_EQ(cell.atoms_per_cell(), r.atoms) << r.name;
    const auto shells = cell.neighbour_shells(cell.sites.front().position, 1.2 * a);
    ASSERT_FALSE(shells.empty()) << r.name;
    EXPECT_NEAR(shells.front().distance_A, r.nn_distance_a * a, 1e-9) << r.name;
    EXPECT_DOUBLE_EQ(shells.front().count, r.nn_count) << r.name;
  }
  // Second shells: bcc 6 at a, fcc 6 at a, diamond 12 at a/sqrt 2.
  EXPECT_NEAR(simple("bcc", {26}, a).neighbour_shells({0, 0, 0}, 1.2 * a)[1].distance_A, a, 1e-9);
  EXPECT_DOUBLE_EQ(simple("bcc", {26}, a).neighbour_shells({0, 0, 0}, 1.2 * a)[1].count, 6);
  EXPECT_DOUBLE_EQ(simple("fcc", {29}, a).neighbour_shells({0, 0, 0}, 1.2 * a)[1].count, 6);
  EXPECT_NEAR(simple("diamond", {14}, a).neighbour_shells({0, 0, 0}, 1.2 * a)[1].distance_A, a / std::sqrt(2.0), 1e-9);
  EXPECT_DOUBLE_EQ(simple("diamond", {14}, a).neighbour_shells({0, 0, 0}, 1.2 * a)[1].count, 12);
  // Perovskite B site: an octahedron of 6 X at a/2.
  const auto pv = simple("perovskite", {38, 22, 8}, a);
  EXPECT_DOUBLE_EQ(pv.neighbour_shells({0.5, 0.5, 0.5}, a)[0].count, 6);
  EXPECT_NEAR(pv.neighbour_shells({0.5, 0.5, 0.5}, a)[0].distance_A, a / 2, 1e-9);
  // Labels: element and number; a sublattice shared by several elements gets its letter.
  const auto gaas = simple("zincblende", {31, 33}, a);
  EXPECT_EQ(gaas.sites[0].label, "Ga1");
  EXPECT_EQ(gaas.sites[3].label, "Ga4");
  EXPECT_EQ(gaas.sites[4].label, "As1");
  StructureSpec alloy;
  alloy.name = "cesium_chloride";
  alloy.sublattices = {{Occupant{26, 0.5}, Occupant{27, 0.5}}, {Occupant{13, 1.0}}};
  EXPECT_EQ(build_structure(alloy).sites[0].label, "A1");
  EXPECT_EQ(build_structure(alloy).sites[1].label, "Al1");
  // hcp with the ideal c/a by default; another c/a splits the 12 into 6 + 6.
  EXPECT_NEAR(in_units(simple("hcp", {22}, a).lattice.c(), angstrom), ideal * a, 1e-12);
  EXPECT_DOUBLE_EQ(simple("hcp", {22}, a, 1.587).neighbour_shells({1.0 / 3, 2.0 / 3, 0.25}, 1.2 * a)[0].count, 6);
}

TEST(CrystalStructures, InterstitialSites) {
  // Number of sites per cell and their lattice-atom neighbours (textbook
  // hole geometry): count and distance in units of a.
  struct Row {
    const char *name;
    const char *sites;
    std::size_t per_cell;
    double count;
    double distance_a;
  };
  const std::vector<Row> rows{
      {"fcc", "octahedral", 4, 6, 0.5},
      {"fcc", "tetrahedral", 8, 4, std::sqrt(3.0) / 4},
      {"bcc", "octahedral", 6, 2, 0.5},                 // 2 at a/2 (then 4 at a/sqrt 2)
      {"bcc", "tetrahedral", 12, 4, std::sqrt(5.0) / 4},
      {"hcp", "octahedral", 2, 6, 1 / std::sqrt(2.0)},
      {"hcp", "tetrahedral", 4, 4, std::sqrt(3.0 / 8.0)},
      {"diamond", "tetrahedral", 8, 4, std::sqrt(3.0) / 4},
  };
  const double a = 4.0;
  for (const auto &r : rows) {
    StructureSpec s;
    s.name = r.name;
    s.a_A = a;
    s.sublattices = {{Occupant{26, 1.0}}};
    s.interstitials = {{r.sites, {Occupant{6, 0.05}}}};
    const auto cell = build_structure(s);
    EXPECT_NO_THROW(cell.validate(r.name)) << r.name << " " << r.sites;
    std::size_t n = 0;
    for (const auto &site : cell.sites) {
      if (site.kind != SiteKind::interstitial) continue;
      ++n;
      // Every site of the set has the same surroundings.
      const auto shells = cell.neighbour_shells(site.position, a);
      ASSERT_FALSE(shells.empty());
      EXPECT_NEAR(shells.front().distance_A, r.distance_a * a, 1e-9) << r.name << " " << r.sites << " " << site.label;
      EXPECT_DOUBLE_EQ(shells.front().count, r.count) << r.name << " " << r.sites << " " << site.label;
    }
    EXPECT_EQ(n, r.per_cell) << r.name << " " << r.sites;
  }
  // bcc octahedral: the second shell is 4 atoms at a/sqrt 2.
  StructureSpec s;
  s.name = "bcc";
  s.a_A = a;
  s.sublattices = {{Occupant{26, 0.9}, Occupant{24, 0.1}}};
  s.interstitials = {{"octahedral", {Occupant{6, 0.01}}}};
  const auto cell = build_structure(s);
  const auto shells = cell.neighbour_shells({0.5, 0.5, 0.0}, a);
  ASSERT_GE(shells.size(), 2u);
  EXPECT_DOUBLE_EQ(shells[1].count, 4);
  EXPECT_NEAR(shells[1].distance_A, a / std::sqrt(2.0), 1e-9);
  // Composition: 2 x (0.9 Fe + 0.1 Cr) + 6 x 0.01 C.
  const auto comp = cell.composition();
  EXPECT_NEAR(comp.at(26), 1.8, 1e-12);
  EXPECT_NEAR(comp.at(24), 0.2, 1e-12);
  EXPECT_NEAR(comp.at(6), 0.06, 1e-12);
  EXPECT_NEAR(cell.atoms_per_cell(), 2.06, 1e-12);
  // Fluorite's empty cube centres: 8 anions at a sqrt(3)/4, then 6 cations at a/2.
  StructureSpec f;
  f.name = "fluorite";
  f.a_A = a;
  f.sublattices = {{Occupant{92, 1.0}}, {Occupant{8, 1.0}}};
  f.interstitials = {{"octahedral", {Occupant{54, 0.01}}}};
  const auto uo2 = build_structure(f);
  const auto fs = uo2.neighbour_shells({0.5, 0.5, 0.5}, a);
  ASSERT_GE(fs.size(), 2u);
  EXPECT_DOUBLE_EQ(fs[0].count, 8);
  EXPECT_NEAR(fs[0].distance_A, a * std::sqrt(3.0) / 4, 1e-9);
  EXPECT_DOUBLE_EQ(fs[1].count, 6);
  EXPECT_NEAR(fs[1].distance_A, a / 2, 1e-9);
}

TEST(CrystalStructures, DensityAndLatticeConstant) {
  // rho = n M / a^3 with 1 amu = 1.66053907e-24 g: Si (diamond, a = 5.431 A).
  auto si = simple("diamond", {14}, 5.431);
  const double m_si = in_units(Element::by_z(14).natural_mass(), atomic_mass_unit);
  const double rho = 8.0 * m_si * 1.66053907 / (5.431 * 5.431 * 5.431); // g/cm3
  EXPECT_NEAR(in_units(si.density({}), gram_per_cubic_centimeter), rho, 1e-6);
  EXPECT_NEAR(rho, 2.329, 0.005); // the density of silicon
  EXPECT_NEAR(in_units(si.atomic_density(), per_cubic_angstrom), 8.0 / std::pow(5.431, 3), 1e-12);
  // An isotope mass instead of the standard atomic weight.
  EXPECT_NEAR(in_units(si.density({{14u, 28.0 * atomic_mass_unit}}), gram_per_cubic_centimeter),
              rho * 28.0 / m_si, 1e-6);
  // Lattice constant from the density: fcc Cu at 8.96 g/cm3.
  auto cu = simple("fcc", {29}, 0.0);
  EXPECT_NEAR(in_units(cu.lattice.a(), angstrom), 1.0, 1e-12);
  cu.scale_to_density(8.96 * gram_per_cubic_centimeter, {});
  EXPECT_NEAR(in_units(cu.density({}), gram_per_cubic_centimeter), 8.96, 1e-9);
  const double m_cu = in_units(Element::by_z(29).natural_mass(), atomic_mass_unit);
  EXPECT_NEAR(in_units(cu.lattice.a(), angstrom), std::cbrt(4.0 * m_cu * 1.66053907 / 8.96), 1e-6);
  EXPECT_NEAR(in_units(cu.lattice.a(), angstrom), 3.61, 0.01);
  EXPECT_THROW(cu.scale_to_density(0.0 * gram_per_cubic_centimeter, {}), std::invalid_argument);
}

TEST(CrystalStructures, SpecsAndCellsAreValidated) {
  StructureSpec s;
  s.name = "bcd";
  s.sublattices = {{Occupant{26, 1.0}}};
  EXPECT_THROW(build_structure(s), std::invalid_argument);
  s.name = "zincblende"; // two sublattices
  EXPECT_THROW(build_structure(s), std::invalid_argument);
  s.name = "sc";
  s.interstitials = {{"octahedral", {Occupant{1, 0.1}}}}; // none built in
  EXPECT_THROW(build_structure(s), std::invalid_argument);
  s.interstitials.clear();
  s.c_over_a = -1.0;
  EXPECT_THROW(build_structure(s), std::invalid_argument);
  EXPECT_THROW(structure_info("nope"), std::invalid_argument);
  EXPECT_EQ(structure_info("hcp").hexagonal, true);
  EXPECT_EQ(parse_site_kind("substitutional"), SiteKind::lattice);
  EXPECT_EQ(parse_site_kind(to_string(SiteKind::interstitial)), SiteKind::interstitial);
  EXPECT_THROW(parse_site_kind("vacancy"), std::invalid_argument);

  auto cell = simple("bcc", {26}, 2.87);
  EXPECT_NO_THROW(cell.validate("c"));
  auto bad = cell;
  bad.sites[0].occupants = {Occupant{26, 0.7}, Occupant{24, 0.4}}; // 1.1
  EXPECT_THROW(bad.validate("c"), std::invalid_argument);
  bad = cell;
  bad.sites[0].occupants = {Occupant{0, 1.0}};
  EXPECT_THROW(bad.validate("c"), std::invalid_argument);
  bad = cell;
  bad.sites[1].position = {0.02, 0.0, 0.0}; // 0.06 A from site 0
  EXPECT_THROW(bad.validate("c"), std::invalid_argument);
  bad = cell;
  bad.sites[1].position = {1.0, 0.5, 0.5};
  EXPECT_THROW(bad.validate("c"), std::invalid_argument);
  bad = cell;
  for (auto &site : bad.sites) site.occupants.clear();
  EXPECT_THROW(bad.validate("c"), std::invalid_argument);
  // Minimum image in a skewed cell: (0,0,0) and (0.9,0.9,0.9) are 0.1 of the diagonal apart.
  EXPECT_NEAR(cell.distance_A({0, 0, 0}, {0.9, 0.9, 0.9}), 0.1 * std::sqrt(3.0) * 2.87, 1e-12);
}

// ---------------------------------------------------------------------------
// Orientation
// ---------------------------------------------------------------------------

TEST(CrystalOrientation, MillerIndicesGiveTheLaboratoryFrame) {
  const auto cubic = simple("fcc", {29}, 3.615).lattice;
  Orientation o;
  o.surface = {0, 0, 1};
  o.reference = {1, 0, 0};
  Mat3 r = rotation(cubic, o);
  EXPECT_LT(rotation_error(r), 1e-14);
  expect_vec(r * Vec3{0, 0, 1}, {1, 0, 0}); // the surface normal is the depth axis
  expect_vec(r * Vec3{1, 0, 0}, {0, 1, 0});
  expect_vec(r * Vec3{0, 1, 0}, {0, 0, 1});
  // (111) surface, [1 -1 0] in the surface along lab y.
  o.surface = {1, 1, 1};
  o.reference = {1, -1, 0};
  r = rotation(cubic, o);
  EXPECT_LT(rotation_error(r), 1e-14);
  expect_vec(r * (Vec3{1, 1, 1} / std::sqrt(3.0)), {1, 0, 0}, 1e-14);
  expect_vec(r * (Vec3{1, -1, 0} / std::sqrt(2.0)), {0, 1, 0}, 1e-14);
  // A reference out of the surface is projected into it: [1 0 0] on (1 1 0) gives [1 -1 0].
  o.surface = {1, 1, 0};
  o.reference = {1, 0, 0};
  r = rotation(cubic, o);
  expect_vec(r * (Vec3{1, -1, 0} / std::sqrt(2.0)), {0, 1, 0}, 1e-14);
  // Automatic reference: a, or b when a is the normal.
  o.reference.clear();
  o.surface = {0, 0, 1};
  expect_vec(rotation(cubic, o) * Vec3{1, 0, 0}, {0, 1, 0});
  o.surface = {1, 0, 0};
  expect_vec(rotation(cubic, o) * Vec3{0, 1, 0}, {0, 1, 0});
  // The twist turns the crystal about the depth axis: lab y -> lab z at 90 degrees.
  o.surface = {0, 0, 1};
  o.reference = {1, 0, 0};
  o.twist = 90.0 * degree;
  r = rotation(cubic, o);
  expect_vec(r * Vec3{0, 0, 1}, {1, 0, 0}, 1e-14);
  expect_vec(r * Vec3{1, 0, 0}, {0, 0, 1}, 1e-14);
  // Errors.
  o.twist = 0.0 * degree;
  o.reference = {0, 0, 2};
  EXPECT_THROW(rotation(cubic, o), std::invalid_argument); // parallel to the normal
  o.reference = {1, 0, 0};
  o.surface = {0, 0, 0};
  EXPECT_THROW(rotation(cubic, o), std::invalid_argument);
  o.surface = {0, 0, 0, 1};
  EXPECT_THROW(rotation(cubic, o), std::invalid_argument); // four indices need a hexagonal cell
  o.surface = {0, 1};
  EXPECT_THROW(rotation(cubic, o), std::invalid_argument);
}

TEST(CrystalOrientation, NonCubicPlanesUseTheReciprocalLattice) {
  // In a tetragonal cell (c = 2a) the normal of (1 0 1) is a* + c* = (1/a, 0, 1/(2a)),
  // not the direction [1 0 1] = (a, 0, 2a).
  const auto tet = Lattice::from_parameters(3.0 * angstrom, 3.0 * angstrom, 6.0 * angstrom, 90.0 * degree, 90.0 * degree,
                                            90.0 * degree);
  const Vec3 n = plane_normal(tet, {1, 0, 1});
  EXPECT_NEAR(angle_deg(n, {2, 0, 1}), 0.0, 1e-6);
  EXPECT_NEAR(angle_deg(direction_vector(tet, {1, 0, 1}), {1, 0, 2}), 0.0, 1e-6);
  // The plane contains the directions [0 1 0] and [1 0 -1] (h u + k v + l w = 0).
  EXPECT_NEAR(dot(n, direction_vector(tet, {0, 1, 0})), 0.0, 1e-14);
  EXPECT_NEAR(dot(n, direction_vector(tet, {1, 0, -1})), 0.0, 1e-14);
}

TEST(CrystalOrientation, HexagonalFourIndexNotation) {
  const auto hex = simple("hcp", {22}, 2.95, 1.587).lattice;
  // [2 -1 -1 0] is the a axis; [1 0 -1 0] lies 30 degrees from it; [0 0 0 1] is c.
  const Vec3 a = direction_vector(hex, direction_indices(hex, {2, -1, -1, 0}, "d"));
  EXPECT_NEAR(angle_deg(a, hex.vectors_A()[0]), 0.0, 1e-6);
  const Vec3 m = direction_vector(hex, direction_indices(hex, {1, 0, -1, 0}, "d"));
  EXPECT_NEAR(angle_deg(m, a), 30.0, 1e-9);
  EXPECT_NEAR(angle_deg(direction_vector(hex, direction_indices(hex, {0, 0, 0, 1}, "d")), hex.vectors_A()[2]), 0.0,
              1e-6);
  // Prism and basal planes: the plane normal is parallel to the direction with the same indices.
  for (const auto &idx : {std::vector<double>{1, 0, -1, 0}, std::vector<double>{1, 1, -2, 0},
                          std::vector<double>{0, 0, 0, 1}}) {
    const Vec3 n = plane_normal(hex, plane_indices(hex, idx, "p"));
    const Vec3 d = direction_vector(hex, direction_indices(hex, idx, "d"));
    EXPECT_NEAR(angle_deg(n, d), 0.0, 1e-6) << format_indices(idx);
  }
  EXPECT_THROW(direction_indices(hex, {1, 0, 0, 0}, "d"), std::invalid_argument); // t != -(u+v)
  EXPECT_THROW(plane_indices(hex, {1, 0, 0, 0}, "p"), std::invalid_argument);
  // A basal-plane surface with a prism direction in it.
  Orientation o;
  o.surface = {0, 0, 0, 1};
  o.reference = {1, 0, -1, 0};
  const Mat3 r = rotation(hex, o);
  EXPECT_LT(rotation_error(r), 1e-14);
  expect_vec(r * normalized(hex.vectors_A()[2]), {1, 0, 0}, 1e-14);
  expect_vec(r * normalized(m), {0, 1, 0}, 1e-14);
  EXPECT_EQ(format_indices({1, 0, -1, 0}), "[1 0 -1 0]");
  EXPECT_EQ(format_indices({1, 1, 1}, '(', ')'), "(1 1 1)");
  EXPECT_EQ(format_indices({0.5, 1, 0}), "[0.5 1 0]");
}

TEST(CrystalOrientation, BungeEulerAngles) {
  const auto cubic = simple("bcc", {26}, 2.8665).lattice;
  Orientation o;
  o.kind = OrientationKind::euler;
  // Cube texture (0,0,0): the crystal axes are the sample axes, so the
  // outward normal is [0 0 1] and the depth axis is [0 0 -1].
  Mat3 r = rotation(cubic, o);
  EXPECT_LT(rotation_error(r), 1e-14);
  expect_vec(r * Vec3{0, 0, -1}, {1, 0, 0});
  expect_vec(r * Vec3{1, 0, 0}, {0, 1, 0});
  // Goss texture {110}<001> is (0, 45, 0): normal direction [0 1 1], rolling direction [1 0 0].
  o.euler_deg = {0, 45, 0};
  r = rotation(cubic, o);
  expect_vec(r * (Vec3{0, 1, 1} / std::sqrt(2.0)), {-1, 0, 0}, 1e-14); // the outward normal
  expect_vec(r * Vec3{1, 0, 0}, {0, 1, 0}, 1e-14);                      // xs = lab y
  // phi1 turns the crystal about the normal: its x axis is then 30 degrees from xs towards ys.
  o.euler_deg = {30, 0, 0};
  r = rotation(cubic, o);
  expect_vec(r * Vec3{1, 0, 0}, {0, std::cos(30 * deg), -std::sin(30 * deg)}, 1e-14); // lab = (-zs, xs, -ys)
  EXPECT_LT(rotation_error(r), 1e-14);
  o.euler_deg = {12.5, 77.0, 301.0};
  EXPECT_LT(rotation_error(rotation(cubic, o)), 1e-14);
}

TEST(CrystalOrientation, TiltAndTwistForABeamDirection) {
  const auto cubic = simple("diamond", {14}, 5.431).lattice;
  Orientation o;
  o.surface = {0, 0, 1};
  o.reference = {1, 1, 0};
  // Along the normal: no tilt, the twist stays.
  o.twist = 22.0 * degree;
  auto tt = tilt_twist_for_beam(cubic, o, {0, 0, 3});
  EXPECT_NEAR(in_units(tt.tilt, degree), 0.0, 1e-12);
  EXPECT_NEAR(in_units(tt.twist, degree), 22.0, 1e-12);
  // [0 1 1] is 45 degrees from [0 0 1]; [1 1 1] is acos(1/sqrt 3) = 54.7356 degrees.
  EXPECT_NEAR(in_units(tilt_twist_for_beam(cubic, o, {0, 1, 1}).tilt, degree), 45.0, 1e-12);
  EXPECT_NEAR(in_units(tilt_twist_for_beam(cubic, o, {1, 1, 1}).tilt, degree), 54.735610317245346, 1e-10);
  // The beam's in-surface part is already along the reference: no twist.
  EXPECT_NEAR(in_units(tilt_twist_for_beam(cubic, o, {1, 1, 1}).twist, degree), 0.0, 1e-10);
  // With the tilt and twist applied, the beam travels along the requested direction.
  const auto tri = Lattice::from_parameters(3.0 * angstrom, 4.0 * angstrom, 5.0 * angstrom, 80.0 * degree, 95.0 * degree,
                                            105.0 * degree);
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  int checked = 0;
  for (const auto *lat : {&cubic, &tri}) {
    Orientation base;
    base.surface = {1, 2, 3};
    base.reference = {1, 0, 0};
    const Vec3 inward = rotation(*lat, base).transposed() * Vec3{1, 0, 0};
    for (int i = 0; i < 200; ++i) {
      const Vec3 beam = normalized({u(rng), u(rng), u(rng)});
      if (dot(beam, inward) < 0.05) {
        continue;
      }
      const auto t = tilt_twist_for_beam(*lat, base, beam * 2.5);
      Orientation with = base;
      with.twist = t.twist;
      expect_vec(beam_in_crystal(*lat, with, t.tilt), beam, 1e-12);
      EXPECT_NEAR(in_units(t.tilt, degree), angle_deg(beam, inward), 1e-9);
      ++checked;
    }
  }
  EXPECT_GT(checked, 100);
  EXPECT_THROW(tilt_twist_for_beam(cubic, o, {0, 0, -1}), std::invalid_argument); // leaves the target
  EXPECT_THROW(tilt_twist_for_beam(cubic, o, {1, 0, 0}), std::invalid_argument);  // grazing: 90 degrees
  EXPECT_THROW(tilt_twist_for_beam(cubic, o, {0, 0, 0}), std::invalid_argument);
  // Naming a direction.
  auto near = nearest_direction(cubic, {1.0, 1.0, 0.02});
  EXPECT_EQ(near.uvw, (std::array<int, 3>{1, 1, 0}));
  EXPECT_NEAR(in_units(near.off, degree), std::atan(0.02 / std::sqrt(2.0)) / deg, 1e-9);
  near = nearest_direction(cubic, {-2.0, 4.0, 6.0});
  EXPECT_EQ(near.uvw, (std::array<int, 3>{-1, 2, 3}));
  EXPECT_NEAR(in_units(near.off, degree), 0.0, 1e-6);
}

// An orientation from a frame (the GUI's "Turn the crystal to this view"): the
// given normal becomes lab +x and the lateral direction's in-surface part
// lab +y, as Miller indices where a low-index plane fits, else as Euler angles.
TEST(CrystalOrientation, OrientationForAFrame) {
  const auto cubic = simple("diamond", {14}, 5.431).lattice;
  const auto tet = Lattice::from_parameters(3.0 * angstrom, 3.0 * angstrom, 6.0 * angstrom, 90.0 * degree, 90.0 * degree,
                                            90.0 * degree);
  const auto tri = Lattice::from_parameters(3.0 * angstrom, 4.0 * angstrom, 5.0 * angstrom, 80.0 * degree, 95.0 * degree,
                                            105.0 * degree);
  // Looking along [1 1 0] with [0 0 1] up: the (1 1 0) surface with the reference [0 0 1].
  auto o = orientation_for_frame(cubic, {2, 2, 0}, {0, 0, 5});
  EXPECT_EQ(o.kind, OrientationKind::miller);
  EXPECT_EQ(o.surface, (std::vector<double>{1, 1, 0}));
  EXPECT_EQ(o.reference, (std::vector<double>{0, 0, 1}));
  EXPECT_EQ(in_units(o.twist, degree), 0.0);
  // Within the tolerance of the plane's normal: the plane itself.
  o = orientation_for_frame(cubic, {1.0, 1.0, 0.0005}, {0, 0, 1});
  EXPECT_EQ(o.surface, (std::vector<double>{1, 1, 0}));
  expect_vec(rotation(cubic, o) * (Vec3{1, 1, 0} / std::sqrt(2.0)), {1, 0, 0}, 1e-14);
  // The view rolled by 17 degrees about [0 0 1]: the same surface, the rest in the twist.
  const Vec3 up = rotation_about({0, 0, 1}, 17.0 * deg) * Vec3{1, 0, 0};
  o = orientation_for_frame(cubic, {0, 0, 1}, up);
  EXPECT_EQ(o.kind, OrientationKind::miller);
  EXPECT_EQ(o.surface, (std::vector<double>{0, 0, 1}));
  EXPECT_TRUE(o.reference.empty());
  EXPECT_NEAR(std::abs(in_units(o.twist, degree)), 17.0, 1e-9);
  expect_vec(rotation(cubic, o) * up, {0, 1, 0}, 1e-12);
  // A tetragonal cell: the plane whose normal is the view, not the direction of the same indices.
  o = orientation_for_frame(tet, {2, 0, 1}, {0, 1, 0});
  EXPECT_EQ(o.surface, (std::vector<double>{1, 0, 1}));
  EXPECT_EQ(o.reference, (std::vector<double>{0, 1, 0}));
  // Any frame, in any lattice: the rotation is the frame's.
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  int euler = 0;
  for (const auto *lat : {&cubic, &tet, &tri}) {
    for (int i = 0; i < 100; ++i) {
      const Vec3 x = normalized({u(rng), u(rng), u(rng)});
      const Vec3 l{u(rng), u(rng), u(rng)};
      if (norm(cross(x, l)) < 0.1) continue;
      const Vec3 y = normalized(l - x * dot(l, x));
      const auto f = orientation_for_frame(*lat, x * 3.0, l);
      const Mat3 r = rotation(*lat, f);
      EXPECT_LT(rotation_error(r), 1e-12);
      // (A direction that fits a low-index plane is moved onto it: up to the tolerance.)
      const double tol = f.kind == OrientationKind::euler ? 1e-9 : 2e-3;
      expect_vec(r * x, {1, 0, 0}, tol);
      expect_vec(r * y, {0, 1, 0}, tol);
      euler += f.kind == OrientationKind::euler ? 1 : 0;
    }
  }
  EXPECT_GT(euler, 200);
  // Phi = 0 and 180 degrees (the view along the sample normal): still the frame.
  for (const Vec3 &x : {Vec3{0, 0, -1}, Vec3{0, 0, 1}}) {
    const Vec3 l = rotation_about({0, 0, 1}, 33.0 * deg) * Vec3{1, 0, 0};
    const auto f = orientation_for_frame(tri, x, l, 0); // no Miller indices: Euler angles
    EXPECT_EQ(f.kind, OrientationKind::euler);
    const Mat3 r = rotation(tri, f);
    expect_vec(r * x, {1, 0, 0}, 1e-9);
    expect_vec(r * l, {0, 1, 0}, 1e-9);
  }
  EXPECT_THROW(orientation_for_frame(cubic, {0, 0, 0}, {1, 0, 0}), std::invalid_argument);
  EXPECT_THROW(orientation_for_frame(cubic, {0, 0, 2}, {0, 0, 1}), std::invalid_argument); // no part in the surface
  EXPECT_THROW(orientation_for_frame(cubic, {0, 0, 2}, {0, 0, 0}), std::invalid_argument);
  // Naming a plane.
  const auto p = nearest_plane(tet, {2.0, 0.0, 1.01});
  EXPECT_EQ(p.uvw, (std::array<int, 3>{1, 0, 1}));
}

// ---------------------------------------------------------------------------
// Thermal vibrations
// ---------------------------------------------------------------------------

TEST(CrystalVibrations, DebyeFunction) {
  EXPECT_DOUBLE_EQ(debye_phi(0.0), 1.0);
  // Series 1 - x/4 + x^2/36 - x^4/3600 for small x.
  for (double x : {1e-7, 1e-3, 0.1, 0.5}) {
    EXPECT_NEAR(debye_phi(x), 1 - x / 4 + x * x / 36 - x * x * x * x / 3600, 2e-6) << x;
  }
  // Abramowitz and Stegun, Table 27.1: (1/x) Int_0^x t/(e^t - 1) dt = 0.777505 at x = 1.
  EXPECT_NEAR(debye_phi(1.0), 0.777505, 1e-6);
  // Large x: the integral tends to pi^2/6.
  EXPECT_NEAR(debye_phi(40.0) * 40.0, pi_v * pi_v / 6.0, 1e-9);
  EXPECT_NEAR(debye_phi(100.0) * 100.0, pi_v * pi_v / 6.0, 1e-12);
  EXPECT_THROW(debye_phi(-1.0), std::invalid_argument);
}

TEST(CrystalVibrations, DebyeAmplitudeIsTheOscillatorAverage) {
  // Independent of the closed form: a harmonic oscillator of frequency w has
  // <x^2> = hbar/(2 M w) coth(hbar w / 2 k T); the Debye spectrum is
  // g(w) = 3 w^2 / w_D^3 up to w_D = k Theta / hbar. Integrate numerically.
  const double hbar = 1.054571817e-34, kb = 1.380649e-23, amu = 1.66053907e-27;
  struct Case {
    double mass_amu, theta, t;
  };
  for (const Case &c : {Case{28.0855, 645.0, 300.0}, Case{183.84, 310.0, 300.0}, Case{55.845, 470.0, 77.0},
                        Case{63.546, 343.0, 1000.0}, Case{196.97, 165.0, 4.0}}) {
    const double wd = kb * c.theta / hbar;
    const int n = 200000;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      const double w = (i + 0.5) * wd / n; // midpoint rule
      sum += 3.0 * w * w / (wd * wd * wd) * hbar / (2.0 * c.mass_amu * amu * w) / std::tanh(hbar * w / (2.0 * kb * c.t));
    }
    const double u_ref_A = std::sqrt(sum * wd / n) * 1e10;
    const double u = in_units(debye_rms_displacement(c.mass_amu * atomic_mass_unit, c.theta * kelvin, c.t * kelvin),
                              angstrom);
    EXPECT_NEAR(u, u_ref_A, 2e-6 * u_ref_A) << c.mass_amu << " " << c.theta << " " << c.t;
    EXPECT_GT(u, 0.01);
    EXPECT_LT(u, 0.2); // thermal amplitudes are a few hundredths of an Angstrom
  }
  // Limits: zero-point motion 3 hbar^2/(4 M k Theta); classical 3 hbar^2 T/(M k Theta^2).
  const auto m = 28.0855 * atomic_mass_unit;
  const double zp = std::sqrt(3.0 * hbar * hbar / (4.0 * 28.0855 * amu * kb * 645.0)) * 1e10;
  // (2e-6: Boost.Units carries the CODATA 2006 constants, the reference the 2018 ones.)
  EXPECT_NEAR(in_units(debye_rms_displacement(m, 645.0 * kelvin, 0.0 * kelvin), angstrom), zp, 2e-6 * zp);
  const double hot = std::sqrt(3.0 * hbar * hbar * 20000.0 / (28.0855 * amu * kb * 645.0 * 645.0)) * 1e10;
  EXPECT_NEAR(in_units(debye_rms_displacement(m, 645.0 * kelvin, 20000.0 * kelvin), angstrom), hot, 1e-4 * hot);
  // Hotter and lighter atoms vibrate more.
  EXPECT_GT(debye_rms_displacement(m, 645.0 * kelvin, 600.0 * kelvin),
            debye_rms_displacement(m, 645.0 * kelvin, 300.0 * kelvin));
  EXPECT_GT(debye_rms_displacement(m, 645.0 * kelvin, 300.0 * kelvin),
            debye_rms_displacement(2.0 * m, 645.0 * kelvin, 300.0 * kelvin));
  EXPECT_THROW(debye_rms_displacement(m, 0.0 * kelvin, 300.0 * kelvin), std::invalid_argument);
  EXPECT_THROW(debye_rms_displacement(m, 645.0 * kelvin, -1.0 * kelvin), std::invalid_argument);
  EXPECT_THROW(debye_rms_displacement(0.0 * kilogram, 645.0 * kelvin, 300.0 * kelvin), std::invalid_argument);

  // A crystal gives the amplitude of the elements it has data for.
  Crystal c;
  c.cell = simple("diamond", {14}, 5.431);
  c.temperature = 300.0 * kelvin;
  EXPECT_FALSE(c.rms_displacement(14, m).has_value());
  c.vibrations.push_back({14, 645.0 * kelvin, std::nullopt});
  ASSERT_TRUE(c.rms_displacement(14, m).has_value());
  EXPECT_EQ(*c.rms_displacement(14, m), debye_rms_displacement(m, 645.0 * kelvin, 300.0 * kelvin));
  c.vibrations[0].rms = 0.05 * angstrom; // a given amplitude wins
  EXPECT_NEAR(in_units(*c.rms_displacement(14, m), angstrom), 0.05, 1e-15);
  EXPECT_FALSE(c.rms_displacement(32, m).has_value());
}

// ---------------------------------------------------------------------------
// Cell files
// ---------------------------------------------------------------------------

namespace {

// Rutile TiO2, space group P 4_2/m n m (136): Ti on 2a, O on 4f (u = 0.305);
// the sixteen operations of International Tables A.
const char *rutile_cif = R"(#------------------------------------------------
# rutile, for the reader test
data_rutile
_chemical_name_mineral 'Rutile'
_chemical_formula_sum 'Ti O2'
_publ_section_comment
;
 A text field with a line that starts with data_ inside:
 data_not_a_block loop_ _not_a_tag
;
_cell_length_a 4.594(1)
_cell_length_b 4.594(1)
_cell_length_c 2.959
_cell_angle_alpha 90
_cell_angle_beta 90.
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'P 42/m n m'
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
 1 'x, y, z'
 2 '-x, -y, z'
 3 '-y+1/2, x+1/2, z+1/2'
 4 'y+1/2, -x+1/2, z+1/2'
 5 '-x+1/2, y+1/2, -z+1/2'
 6 'x+1/2, -y+1/2, -z+1/2'
 7 'y, x, -z'
 8 '-y, -x, -z'
 9 '-x, -y, -z'
 10 'x, y, -z'
 11 'y+1/2, -x+1/2, -z+1/2'
 12 '-y+1/2, x+1/2, -z+1/2'
 13 'x+1/2, -y+1/2, z+1/2'
 14 '-x+1/2, y+1/2, z+1/2'
 15 '-y, -x, z'
 16 'y, x, z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Ti1 Ti4+ 0 0 0 1.0
O1 O2- 0.305(2) 0.305(2) 0 1
)";

// A disordered P 1 cell: Fe and Cr share one site; a carbon interstitial.
const char *alloy_cif = R"(data_alloy
_cell_length_a 2.87
_cell_length_b 2.87
_cell_length_c 2.87
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_space_group_name_H-M_alt "P 1"
loop_
_atom_site_label _atom_site_fract_x _atom_site_fract_y _atom_site_fract_z _atom_site_occupancy _amcset_site_kind
Fe1 0 0 0 0.9 lattice
Cr1 0 0 0 0.10(1) lattice
Fe2 0.5 0.5 0.5 ? .
C1 0.5 0.5 0 0.02 interstitial
)";

} // namespace

TEST(CrystalFiles, SymmetryOperations) {
  auto [m, t] = parse_symmetry_operation("-x+1/2, y+1/2, -z");
  expect_vec(m[0], {-1, 0, 0});
  expect_vec(m[1], {0, 1, 0});
  expect_vec(m[2], {0, 0, -1});
  expect_vec(t, {0.5, 0.5, 0});
  std::tie(m, t) = parse_symmetry_operation("X-Y, x, z+1/6");
  expect_vec(m[0], {1, -1, 0});
  expect_vec(m[1], {1, 0, 0});
  expect_vec(t, {0, 0, 1.0 / 6});
  std::tie(m, t) = parse_symmetry_operation("1/2+x,-y,0.25-z");
  expect_vec(m[0], {1, 0, 0});
  expect_vec(t, {0.5, 0, 0.25});
  expect_vec(m[2], {0, 0, -1});
  std::tie(m, t) = parse_symmetry_operation("+x,+y,+z");
  expect_vec(m * Vec3{0.1, 0.2, 0.3} + t, {0.1, 0.2, 0.3});
  EXPECT_THROW(parse_symmetry_operation("x,y"), std::invalid_argument);
  EXPECT_THROW(parse_symmetry_operation("x,q,z"), std::invalid_argument);
  EXPECT_THROW(parse_symmetry_operation("x,,z"), std::invalid_argument);
  EXPECT_THROW(parse_symmetry_operation("x,y,z+"), std::invalid_argument);

  EXPECT_EQ(element_of_symbol("Fe3+"), 26u);
  EXPECT_EQ(element_of_symbol("O2-"), 8u);
  EXPECT_EQ(element_of_symbol("Si1"), 14u);
  EXPECT_EQ(element_of_symbol("CA2"), 20u);
  EXPECT_EQ(element_of_symbol("OW"), 8u); // not an element: the first letter
  EXPECT_EQ(element_of_symbol("U"), 92u);
  EXPECT_FALSE(element_of_symbol("X1").has_value());
  EXPECT_FALSE(element_of_symbol("12").has_value());
  EXPECT_FALSE(element_of_symbol("").has_value());
}

TEST(CrystalFiles, CifWithSymmetry) {
  const auto cell = parse_cell(rutile_cif, CellFormat::cif);
  EXPECT_EQ(cell.name, "Ti O2 (rutile)");
  EXPECT_NEAR(in_units(cell.lattice.a(), angstrom), 4.594, 1e-12);
  EXPECT_NEAR(in_units(cell.lattice.c(), angstrom), 2.959, 1e-12);
  // 2 Ti + 4 O from one Ti and one O of the asymmetric unit.
  ASSERT_EQ(cell.sites.size(), 6u);
  const auto comp = cell.composition();
  EXPECT_DOUBLE_EQ(comp.at(22), 2.0);
  EXPECT_DOUBLE_EQ(comp.at(8), 4.0);
  // Ti at (0,0,0) and (1/2,1/2,1/2); O at +-(u,u,0) and (1/2 -+ u, 1/2 +- u, 1/2).
  const double u = 0.305;
  for (const Vec3 &p : {Vec3{0, 0, 0}, Vec3{0.5, 0.5, 0.5}, Vec3{u, u, 0}, Vec3{1 - u, 1 - u, 0},
                        Vec3{0.5 - u, 0.5 + u, 0.5}, Vec3{0.5 + u, 0.5 - u, 0.5}}) {
    bool found = false;
    for (const auto &s : cell.sites) found = found || cell.distance_A(s.position, p) < 1e-9;
    EXPECT_TRUE(found) << p.x << " " << p.y << " " << p.z;
  }
  // The TiO6 octahedron of rutile: two apical O at sqrt(2) u a = 1.98 A and
  // four equatorial at 1.95 A (the textbook bond lengths).
  std::vector<double> d;
  for (const auto &s : cell.sites) {
    if (s.occupants.front().z != 8) continue;
    for (int i = -1; i <= 1; ++i) {
      for (int j = -1; j <= 1; ++j) {
        for (int k = -1; k <= 1; ++k) {
          const double r = norm(cell.lattice.to_cartesian_A(s.position + Vec3{double(i), double(j), double(k)}));
          if (r < 2.5) d.push_back(r);
        }
      }
    }
  }
  std::sort(d.begin(), d.end());
  ASSERT_EQ(d.size(), 6u);
  EXPECT_NEAR(d[0], 1.948, 0.002);
  EXPECT_NEAR(d[3], 1.948, 0.002);
  EXPECT_NEAR(d[4], std::sqrt(2.0) * u * 4.594, 1e-9);
  EXPECT_NEAR(d[5], 1.9815, 0.001);
}

TEST(CrystalFiles, CifOccupancyAndErrors) {
  const auto cell = parse_cell(alloy_cif, CellFormat::cif);
  ASSERT_EQ(cell.sites.size(), 3u);
  // Fe1 and Cr1 are one site with two occupants.
  ASSERT_EQ(cell.sites[0].occupants.size(), 2u);
  EXPECT_EQ(cell.sites[0].occupants[0], (Occupant{26, 0.9}));
  EXPECT_EQ(cell.sites[0].occupants[1], (Occupant{24, 0.1}));
  EXPECT_EQ(cell.sites[1].occupants[0], (Occupant{26, 1.0})); // "?" = fully occupied
  EXPECT_EQ(cell.sites[1].kind, SiteKind::lattice);
  EXPECT_EQ(cell.sites[2].kind, SiteKind::interstitial);
  EXPECT_NEAR(cell.sites[2].occupants[0].fraction, 0.02, 1e-15);
  EXPECT_NEAR(cell.atoms_per_cell(), 2.02, 1e-12);
  // Round trip through AMCSET's own CIF (occupancies and site kinds are kept).
  expect_same_cell(parse_cell(write_cell(cell, CellFormat::cif), CellFormat::cif), cell);

  const std::string text(alloy_cif);
  const auto with = [&](const std::string &from, const std::string &to) {
    std::string s = text;
    const auto p = s.find(from);
    EXPECT_NE(p, std::string::npos) << from;
    return s.replace(p, from.size(), to);
  };
  // A space group without its operations cannot be expanded.
  EXPECT_THROW(parse_cell(with("\"P 1\"", "'I m -3 m'"), CellFormat::cif), std::invalid_argument);
  EXPECT_THROW(parse_cell(with("_space_group_name_H-M_alt \"P 1\"", "_space_group_IT_number 229"), CellFormat::cif),
               std::invalid_argument);
  EXPECT_NO_THROW(parse_cell(with("_space_group_name_H-M_alt \"P 1\"", "_space_group_IT_number 1"), CellFormat::cif));
  EXPECT_THROW(parse_cell(with("_cell_length_a 2.87", ""), CellFormat::cif), std::invalid_argument);
  EXPECT_THROW(parse_cell(with("0.10(1)", "0.3"), CellFormat::cif), std::invalid_argument);      // 1.2 on one site
  EXPECT_THROW(parse_cell(with("Fe2 0.5 0.5 0.5 ? .", "Fe2 0.5 0.5 ? ."), CellFormat::cif), std::invalid_argument); // short row
  EXPECT_THROW(parse_cell(with("Fe2 0.5", "Xx2 0.5"), CellFormat::cif), std::invalid_argument);  // no element
  EXPECT_THROW(parse_cell(with("C1 0.5 0.5 0 0.02", "C1 0.5 abc 0 0.02"), CellFormat::cif), std::invalid_argument);
  EXPECT_THROW(parse_cell("no cif here", CellFormat::cif), std::invalid_argument);
  EXPECT_THROW(parse_cell("data_x\n_cell_length_a 'unterminated\n", CellFormat::cif), std::invalid_argument);
}

TEST(CrystalFiles, Poscar) {
  // VASP 5, Direct, with selective dynamics and a POTCAR-style species name.
  const char *direct = "bcc Fe\n1.0\n 2.8665 0 0\n 0 2.8665 0\n 0 0 2.8665\n Fe_pv\n 2\nSelective dynamics\nDirect\n"
                       " 0 0 0 T T T\n 0.5 0.5 0.5 F F F\n";
  auto cell = parse_cell(direct, CellFormat::poscar);
  EXPECT_EQ(cell.name, "bcc Fe");
  expect_same_cell(cell, simple("bcc", {26}, 2.8665));
  // Cartesian coordinates and lattice share the scale factor.
  const char *cart = "NaCl\n5.64\n 1 0 0\n 0 1 0\n 0 0 1\n Na Cl\n 4 4\nCartesian\n"
                     "0 0 0\n0.5 0.5 0\n0.5 0 0.5\n0 0.5 0.5\n0.5 0 0\n0 0.5 0\n0 0 0.5\n0.5 0.5 0.5\n";
  expect_same_cell(parse_cell(cart, CellFormat::poscar), simple("rocksalt", {11, 17}, 5.64));
  // A negative scale is the cell volume.
  const char *volume = "Cu\n-47.24\n 1 0 0\n 0 1 0\n 0 0 1\n Cu\n 4\nD\n0 0 0\n0.5 0.5 0\n0.5 0 0.5\n0 0.5 0.5\n";
  EXPECT_NEAR(parse_cell(volume, CellFormat::poscar).lattice.volume_A3(), 47.24, 1e-9);
  // VASP 4: no species line.
  const char *v4 = "Si\n5.431\n 1 0 0\n 0 1 0\n 0 0 1\n 2\nDirect\n0 0 0\n0.25 0.25 0.25\n";
  EXPECT_THROW(parse_cell(v4, CellFormat::poscar), std::invalid_argument);
  CellReadOptions opt;
  opt.elements = {14};
  EXPECT_EQ(parse_cell(v4, CellFormat::poscar, opt).composition().at(14), 2.0);
  // Errors: too few atoms, unknown element, bad numbers, wrong mode.
  EXPECT_THROW(parse_cell("x\n1.0\n 1 0 0\n 0 1 0\n 0 0 1\n Fe\n 2\nDirect\n 0 0 0\n", CellFormat::poscar),
               std::invalid_argument);
  EXPECT_THROW(parse_cell("x\n1.0\n 1 0 0\n 0 1 0\n 0 0 1\n Qq\n 1\nDirect\n 0 0 0\n", CellFormat::poscar),
               std::invalid_argument);
  EXPECT_THROW(parse_cell("x\n1.0\n 1 0 0\n 0 1 0\n 0 0\n Fe\n 1\nDirect\n 0 0 0\n", CellFormat::poscar),
               std::invalid_argument);
  EXPECT_THROW(parse_cell("x\n\n 1 0 0\n 0 1 0\n 0 0 1\n Fe\n 1\nDirect\n 0 0 0\n", CellFormat::poscar),
               std::invalid_argument);
  EXPECT_THROW(parse_cell("x\n1.0\n 1 0 0\n 0 1 0\n 0 0 1\n Fe\n 1\nReciprocal\n 0 0 0\n", CellFormat::poscar),
               std::invalid_argument);
  EXPECT_THROW(parse_cell("x\n1.0\n 1 0 0\n 0 1 0\n 0 0 1\n Fe Cr\n 1\nDirect\n 0 0 0\n", CellFormat::poscar),
               std::invalid_argument);
}

TEST(CrystalFiles, LammpsData) {
  // Atomic style, an origin that is not zero, elements from the Masses comments.
  const char *atomic = "bcc Fe-Cr supercell\n\n2 atoms\n2 atom types\n\n-1.0 1.8665 xlo xhi\n0.0 2.8665 ylo yhi\n"
                       "0.0 2.8665 zlo zhi\n\nMasses\n\n1 55.845 # Fe\n2 51.9961 # Cr\n\nAtoms # atomic\n\n"
                       "1 1 -1.0 0.0 0.0\n2 2 0.43325 1.43325 1.43325 0 0 0\n\nVelocities\n\n1 0 0 0\n2 0 0 0\n";
  auto cell = parse_cell(atomic, CellFormat::lammps_data);
  expect_same_cell(cell, simple("cesium_chloride", {26, 24}, 2.8665));
  // Triclinic box, full style, elements from the masses alone.
  const char *full = "hcp Ti\n\n2 atoms\n1 atom types\n\n0 2.95 xlo xhi\n0 2.554774941 ylo yhi\n0 4.68 zlo zhi\n"
                     "-1.475 0 0 xy xz yz\n\nMasses\n\n1 47.867\n\nAtoms # full\n\n"
                     "1 1 1 0.0 0.0 1.703183294 1.17\n2 1 1 0.0 1.475 0.851591647 3.51\n";
  cell = parse_cell(full, CellFormat::lammps_data);
  EXPECT_TRUE(cell.lattice.is_hexagonal());
  expect_same_cell(cell, simple("hcp", {22}, 2.95, 4.68 / 2.95), 1e-5);
  // A given type map wins; an unknown style or too few columns is an error.
  CellReadOptions opt;
  opt.lammps_types = {{1, 29u}, {2, 79u}};
  EXPECT_EQ(parse_cell(atomic, CellFormat::lammps_data, opt).composition().count(79), 1u);
  std::string s(atomic);
  const auto swap = [&](const std::string &from, const std::string &to) {
    std::string t = s;
    const auto p = t.find(from);
    EXPECT_NE(p, std::string::npos);
    return t.replace(p, from.size(), to);
  };
  EXPECT_THROW(parse_cell(swap("Atoms # atomic", "Atoms # sphere"), CellFormat::lammps_data), std::invalid_argument);
  EXPECT_THROW(parse_cell(swap("2 atoms", "3 atoms"), CellFormat::lammps_data), std::invalid_argument);
  EXPECT_THROW(parse_cell(swap("0.0 2.8665 zlo zhi\n", ""), CellFormat::lammps_data), std::invalid_argument);
  EXPECT_THROW(parse_cell(swap("1 55.845 # Fe", "1 12.5"), CellFormat::lammps_data), std::invalid_argument); // no such element
  EXPECT_THROW(parse_cell(swap("2 2 0.43325", "2 3 0.43325"), CellFormat::lammps_data), std::invalid_argument); // type 3 of 2
  // Six columns could be charge or molecular: the style must be given.
  const char *six = "x\n\n1 atoms\n1 atom types\n\n0 3 xlo xhi\n0 3 ylo yhi\n0 3 zlo zhi\n\nMasses\n\n1 55.845\n\n"
                    "Atoms\n\n1 1 0.0 0 0 0\n";
  EXPECT_THROW(parse_cell(six, CellFormat::lammps_data), std::invalid_argument);
  opt = {};
  opt.lammps_atom_style = "charge";
  EXPECT_EQ(parse_cell(six, CellFormat::lammps_data, opt).sites.size(), 1u);
}

TEST(CrystalFiles, ExtendedXyz) {
  const char *xyz = "2\nLattice=\"2.8665 0 0 0 2.8665 0 0 0 2.8665\" Properties=species:S:1:pos:R:3:forces:R:3 "
                    "energy=-8.2 pbc=\"T T T\"\nFe 0 0 0 0.1 0.2 0.3\nFe 1.43325 1.43325 1.43325 0 0 0\n";
  expect_same_cell(parse_cell(xyz, CellFormat::extxyz), simple("bcc", {26}, 2.8665));
  // Atomic numbers instead of symbols, columns in another order.
  const char *by_z = "1\nProperties=pos:R:3:Z:I:1 Lattice='3 0 0 0 3 0 0 0 3'\n1.5 1.5 1.5 74\n";
  const auto cell = parse_cell(by_z, CellFormat::extxyz);
  EXPECT_EQ(cell.composition().count(74), 1u);
  expect_vec(cell.sites[0].position, {0.5, 0.5, 0.5}, 1e-12);
  EXPECT_THROW(parse_cell("1\ncomment only\nFe 0 0 0\n", CellFormat::extxyz), std::invalid_argument); // no cell
  EXPECT_THROW(parse_cell("2\nLattice=\"3 0 0 0 3 0 0 0 3\"\nFe 0 0 0\n", CellFormat::extxyz), std::invalid_argument);
  EXPECT_THROW(parse_cell("1\nLattice=\"3 0 0 0 3 0 0 0\"\nFe 0 0 0\n", CellFormat::extxyz), std::invalid_argument);
  EXPECT_THROW(parse_cell("1\nLattice=\"3 0 0 0 3 0 0 0 3\"\nQq 0 0 0\n", CellFormat::extxyz), std::invalid_argument);
  EXPECT_THROW(parse_cell("1\nLattice=\"3 0 0 0 3 0 0 0 3\"\nFe 0 0\n", CellFormat::extxyz), std::invalid_argument);
  EXPECT_THROW(parse_cell("", CellFormat::extxyz), std::invalid_argument);
}

TEST(CrystalFiles, WritersRoundTrip) {
  // A triclinic cell with three atoms of two elements.
  UnitCell cell;
  cell.name = "test cell";
  cell.lattice = Lattice::from_parameters(3.1 * angstrom, 4.2 * angstrom, 5.3 * angstrom, 82.0 * degree, 97.0 * degree,
                                          104.0 * degree);
  cell.sites = {{"W1", {0.1, 0.2, 0.3}, SiteKind::lattice, {Occupant{74, 1.0}}},
                {"W2", {0.6, 0.7, 0.05}, SiteKind::lattice, {Occupant{74, 1.0}}},
                {"Re1", {0.35, 0.9, 0.55}, SiteKind::lattice, {Occupant{75, 1.0}}}};
  for (CellFormat f : {CellFormat::cif, CellFormat::poscar, CellFormat::lammps_data, CellFormat::extxyz}) {
    const std::string text = write_cell(cell, f);
    const auto back = parse_cell(text, f);
    expect_same_cell(back, cell, 2e-6);
    // Writing what was read gives the same text (stable output).
    EXPECT_EQ(write_cell(back, f), write_cell(parse_cell(write_cell(back, f), f), f)) << to_string(f);
  }
  // Only CIF stores partial occupancy.
  cell.sites[2].occupants = {Occupant{75, 0.5}, Occupant{76, 0.25}};
  EXPECT_NO_THROW(write_cell(cell, CellFormat::cif));
  expect_same_cell(parse_cell(write_cell(cell, CellFormat::cif), CellFormat::cif), cell, 2e-6);
  EXPECT_THROW(write_cell(cell, CellFormat::poscar), std::invalid_argument);
  EXPECT_THROW(write_cell(cell, CellFormat::lammps_data), std::invalid_argument);
  EXPECT_THROW(write_cell(cell, CellFormat::extxyz), std::invalid_argument);

  // Files: format from the name, read back.
  const auto dir = std::filesystem::temp_directory_path() / "amcset_test_crystal_files";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  const auto si = simple("diamond", {14}, 5.431);
  for (const char *name : {"si.cif", "POSCAR", "si.vasp", "si.data", "data.si", "si.xyz"}) {
    save_cell(si, dir / name);
    expect_same_cell(read_cell(dir / name), si);
  }
  EXPECT_EQ(format_from_path("CONTCAR"), CellFormat::poscar);
  EXPECT_EQ(format_from_path("a/b/Fe.CIF"), CellFormat::cif);
  EXPECT_EQ(format_from_path("x.extxyz"), CellFormat::extxyz);
  EXPECT_EQ(format_from_path("x.lmp"), CellFormat::lammps_data);
  EXPECT_FALSE(format_from_path("x.txt").has_value());
  EXPECT_THROW(read_cell(dir / "x.txt"), std::invalid_argument);
  EXPECT_THROW(read_cell(dir / "missing.cif"), std::runtime_error);
  EXPECT_NO_THROW(read_cell(dir / "si.vasp", CellFormat::poscar));
  EXPECT_EQ(parse_cell_format("VASP"), CellFormat::poscar);
  EXPECT_EQ(parse_cell_format(to_string(CellFormat::lammps_data)), CellFormat::lammps_data);
  EXPECT_THROW(parse_cell_format("pdb"), std::invalid_argument);
  std::filesystem::remove_all(dir);
}

// ---------------------------------------------------------------------------
// The collision-partner search of the crystal transport
// ---------------------------------------------------------------------------

TEST(CrystalSearch, FindsTheAtomsAheadLikeABruteForceSearch) {
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> u(-1.0, 1.0);
  struct Case {
    UnitCell cell;
    std::map<unsigned, double> rms;
  };
  StructureSpec alloy;
  alloy.name = "bcc";
  alloy.a_A = 2.9;
  alloy.sublattices = {{Occupant{26, 0.6}, Occupant{24, 0.3}}}; // a tenth of the sites vacant
  const std::vector<Case> cases{
      {simple("bcc", {26}, 2.8665), {{26u, 0.06}}},
      {simple("diamond", {14}, 5.431), {{14u, 0.075}}},
      {simple("hcp", {22}, 2.95, 1.587), {{22u, 0.08}}},
      {simple("fcc", {29}, 3.615), {}}, // static
      {build_structure(alloy), {{26u, 0.06}}}, // Cr static, Fe vibrating
  };
  const double pmax = 1.9, reach = 2.3;
  std::size_t total = 0;
  for (const auto &c : cases) {
    Orientation o;
    o.kind = OrientationKind::euler;
    o.euler_deg = {180.0 * (u(rng) + 1.0), 90.0 * (u(rng) + 1.0), 180.0 * (u(rng) + 1.0)};
    const Mat3 rot = rotation(c.cell.lattice, o);
    const LatticeSearch search(c.cell, rot, c.rms, pmax, reach);
    EXPECT_DOUBLE_EQ(search.max_impact_A(), pmax);
    EXPECT_DOUBLE_EQ(search.reach_A(), reach);
    std::vector<Candidate> found;
    for (int trial = 0; trial < 60; ++trial) {
      const Vec3 r{40.0 * u(rng), 40.0 * u(rng), 40.0 * u(rng)};
      const Vec3 dir = normalized({u(rng), u(rng), u(rng)});
      const std::uint64_t seed = 7919u * static_cast<unsigned>(trial) + 3u;
      const double x_lo = -15.0, x_hi = 22.0; // the layer's depth range
      search.candidates(r, dir, seed, x_lo, x_hi, found);
      // Every site of the cells around r, one by one.
      const Vec3 f = c.cell.lattice.to_fractional(rot.transposed() * r);
      const int c0[3] = {static_cast<int>(std::floor(f.x)), static_cast<int>(std::floor(f.y)),
                         static_cast<int>(std::floor(f.z))};
      std::vector<Candidate> expected;
      for (int i = c0[0] - 4; i <= c0[0] + 4; ++i) {
        for (int j = c0[1] - 4; j <= c0[1] + 4; ++j) {
          for (int k = c0[2] - 4; k <= c0[2] + 4; ++k) {
            for (std::uint32_t s = 0; s < c.cell.sites.size(); ++s) {
              const SiteId id{i, j, k, s};
              const Vec3 site = search.site_position(id);
              if (site.x < x_lo - 1e-7 || site.x >= x_hi - 1e-7) continue;
              unsigned z = 0;
              Vec3 p;
              if (!search.atom_at(id, seed, z, p)) continue;
              const Vec3 d = p - r;
              const double xi = dot(d, dir);
              if (!(xi > 0.0) || xi > reach) continue;
              const double s2 = dot(d, d) - xi * xi;
              if (s2 > pmax * pmax) continue;
              expected.push_back({id, z, p, xi, std::sqrt(std::max(0.0, s2))});
            }
          }
        }
      }
      std::sort(expected.begin(), expected.end(), [](const Candidate &a, const Candidate &b) { return a.xi < b.xi; });
      ASSERT_EQ(found.size(), expected.size()) << c.cell.name << " trial " << trial;
      for (std::size_t n = 0; n < found.size(); ++n) {
        EXPECT_TRUE(found[n].id == expected[n].id) << c.cell.name << " trial " << trial;
        EXPECT_EQ(found[n].z, expected[n].z);
        EXPECT_NEAR(found[n].xi, expected[n].xi, 1e-9);
        EXPECT_NEAR(found[n].impact, expected[n].impact, 1e-9);
        EXPECT_LT(norm(found[n].position - expected[n].position), 1e-9);
        if (n > 0) {
          EXPECT_GE(found[n].xi, found[n - 1].xi);
        }
      }
      total += found.size();
    }
  }
  EXPECT_GT(total, 150u); // the comparison was not of empty lists
  EXPECT_THROW(LatticeSearch(cases[0].cell, Mat3::identity(), {}, 0.0, 1.0), std::invalid_argument);
  EXPECT_THROW(LatticeSearch(cases[0].cell, Mat3::identity(), {}, 1.0, 0.0), std::invalid_argument);
  EXPECT_THROW(LatticeSearch(cases[0].cell, Mat3::identity(), {{26u, -0.1}}, 1.0, 1.0), std::invalid_argument);
}

TEST(CrystalSearch, ThermalDisplacementsAndOccupancy) {
  // One seed is one frozen crystal: 70 % Fe, 20 % Cr, 10 % vacancies, with
  // Gaussian displacements of 0.06 A (Fe) and 0.10 A (Cr) along each axis.
  StructureSpec spec;
  spec.name = "sc";
  spec.a_A = 2.5;
  spec.sublattices = {{Occupant{26, 0.7}, Occupant{24, 0.2}}};
  const UnitCell cell = build_structure(spec);
  const LatticeSearch search(cell, Mat3::identity(), {{26u, 0.06}, {24u, 0.10}}, 1.5, 1.5);
  const int n = 48; // 110 592 sites
  double count[3] = {0, 0, 0}; // Fe, Cr, vacant
  double sum[2][3] = {{0, 0, 0}, {0, 0, 0}}, sum2[2][3] = {{0, 0, 0}, {0, 0, 0}};
  double largest[2] = {0, 0};
  double covariance = 0.0, pairs = 0.0; // of the x displacements of neighbouring Fe atoms
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k) {
        const SiteId id{i, j, k, 0};
        unsigned z = 0;
        Vec3 p;
        if (!search.atom_at(id, 5, z, p)) {
          count[2] += 1;
          continue;
        }
        const int e = z == 26 ? 0 : 1;
        ASSERT_TRUE(z == 26 || z == 24);
        count[e] += 1;
        const Vec3 d = p - search.site_position(id);
        for (int c = 0; c < 3; ++c) {
          sum[e][c] += d[c];
          sum2[e][c] += d[c] * d[c];
          largest[e] = std::max(largest[e], std::abs(d[c]));
        }
        unsigned z2 = 0;
        Vec3 p2;
        const SiteId next{i + 1, j, k, 0};
        if (e == 0 && search.atom_at(next, 5, z2, p2) && z2 == 26) {
          covariance += d.x * (p2 - search.site_position(next)).x;
          pairs += 1;
        }
      }
    }
  }
  const double sites = static_cast<double>(n) * n * n;
  EXPECT_NEAR(count[0] / sites, 0.7, 0.005);
  EXPECT_NEAR(count[1] / sites, 0.2, 0.005);
  EXPECT_NEAR(count[2] / sites, 0.1, 0.004);
  const double sigma[2] = {0.06, 0.10};
  for (int e = 0; e < 2; ++e) {
    for (int c = 0; c < 3; ++c) {
      EXPECT_NEAR(sum[e][c] / count[e], 0.0, 5.0 * sigma[e] / std::sqrt(count[e])) << e << " " << c;
      EXPECT_NEAR(std::sqrt(sum2[e][c] / count[e]), sigma[e], 0.015 * sigma[e]) << e << " " << c;
    }
    // Cut at 3.5 standard deviations per axis, and reached.
    EXPECT_LE(largest[e], LatticeSearch::displacement_cut * sigma[e] + 1e-12);
    EXPECT_GT(largest[e], 3.0 * sigma[e]);
  }
  EXPECT_NEAR(search.max_displacement_A(), 3.5 * std::sqrt(3.0) * 0.10, 1e-12);
  // Neighbours are uncorrelated (five standard errors of the correlation).
  EXPECT_LT(std::abs(covariance / pairs) / (0.06 * 0.06), 5.0 / std::sqrt(pairs));
  // The same seed is the same crystal; another seed is another one.
  unsigned za = 0, zb = 0;
  Vec3 pa, pb;
  int differ = 0;
  for (int i = 0; i < 200; ++i) {
    const SiteId id{i, -3 * i, 7, 0};
    const bool a = search.atom_at(id, 5, za, pa);
    EXPECT_EQ(search.atom_at(id, 5, zb, pb), a);
    if (a) {
      EXPECT_EQ(za, zb);
      EXPECT_TRUE(pa == pb);
    }
    const bool other = search.atom_at(id, 6, zb, pb);
    if (a != other || (a && !(pa == pb))) ++differ;
  }
  EXPECT_GT(differ, 150);
  // Far from the origin the lattice is the same lattice.
  const SiteId far{1000000, -2000000, 3000000, 0};
  expect_vec(search.site_position(far), {2.5e6, -5e6, 7.5e6}, 1e-6);
  EXPECT_NE(site_hash({1, 2, 3, 0}, 9), site_hash({2, 1, 3, 0}, 9));
  EXPECT_NE(site_hash({1, 2, 3, 0}, 9), site_hash({1, 2, 3, 1}, 9));
  EXPECT_NE(site_hash({1, 2, 3, 0}, 9), site_hash({1, 2, 3, 0}, 10));
  EXPECT_NE(site_hash({-1, 0, 0, 0}, 9), site_hash({1, 0, 0, 0}, 9));
}

TEST(CrystalSearch, LargestHoleRadius) {
  // The largest interstitial hole of each structure (textbook geometry), in
  // units of a: the default maximum impact parameter of the transport.
  const double a = 4.0;
  struct Row {
    const char *name;
    std::vector<unsigned> z;
    double radius_a;
  };
  for (const Row &r : {Row{"sc", {84}, std::sqrt(3.0) / 2},            // the cube centre
                       Row{"bcc", {26}, std::sqrt(5.0) / 4},           // tetrahedral site
                       Row{"fcc", {29}, 0.5},                          // octahedral site
                       Row{"diamond", {14}, std::sqrt(3.0) / 4},       // tetrahedral site
                       Row{"hcp", {22}, 1 / std::sqrt(2.0)},           // octahedral site (ideal c/a)
                       Row{"rocksalt", {11, 17}, std::sqrt(3.0) / 4},  // tetrahedral site
                       Row{"cesium_chloride", {55, 17}, std::sqrt(5.0) / 4}}) {
    EXPECT_NEAR(largest_hole_radius_A(simple(r.name, r.z, a)), r.radius_a * a, 0.005 * r.radius_a * a) << r.name;
  }
  // Atoms on interstitial sites do not make the holes smaller (the hole is the lattice's).
  StructureSpec s;
  s.name = "fcc";
  s.a_A = a;
  s.sublattices = {{Occupant{26, 1.0}}};
  s.interstitials = {{"octahedral", {Occupant{6, 0.05}}}};
  EXPECT_NEAR(largest_hole_radius_A(build_structure(s)), 0.5 * a, 0.0025 * a);
  UnitCell empty = simple("sc", {84}, a);
  empty.sites[0].occupants.clear();
  EXPECT_THROW(largest_hole_radius_A(empty), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Common crystals (Gemmell, Rev. Mod. Phys. 46, 129 (1974), Tables IV and V)
// ---------------------------------------------------------------------------

TEST(CrystalMaterials, DebyeAmplitudesMatchGemmellsTable) {
  // Table IV lists u1(293 K) for each element next to its Debye
  // temperature: the Debye-model function reproduces every one of them
  // (the table rounds to 0.001 A, to 0.01 A for diamond).
  int elements = 0;
  for (const auto &row : reference_crystals()) {
    if (row.sublattices.size() != 1) continue;
    ++elements;
    const unsigned z = row.sublattices.front();
    const double u = in_units(debye_rms_displacement(Element::by_z(z).natural_mass(), row.debye_temperature_K * kelvin,
                                                     293.0 * kelvin),
                              angstrom);
    EXPECT_NEAR(u, row.u1_293K_A, row.name == "C" ? 0.006 : 0.0012) << row.name;
    EXPECT_EQ(in_units(*reference_debye_temperature(z), kelvin), row.debye_temperature_K);
  }
  EXPECT_EQ(elements, 19);
  // The worked examples of Gemmell's Appendix A: u1 = 0.085 A for Ge at
  // 20 C (Debye temperature 290 K), 0.088 A for Au at 25 C (170 K).
  EXPECT_NEAR(in_units(debye_rms_displacement(72.59 * atomic_mass_unit, 290.0 * kelvin, 293.15 * kelvin), angstrom),
              0.085, 0.001);
  EXPECT_NEAR(in_units(debye_rms_displacement(196.97 * atomic_mass_unit, 170.0 * kelvin, 298.15 * kelvin), angstrom),
              0.088, 0.001);
  EXPECT_FALSE(reference_debye_temperature(22).has_value()); // Ti is not in the table
}

TEST(CrystalMaterials, RowsBuildTheirStructures) {
  EXPECT_EQ(reference_crystals().size(), 19u + 24u);
  for (const auto &row : reference_crystals()) {
    const Crystal c = make_reference_crystal(row);
    EXPECT_NO_THROW(c.cell.validate(std::string(row.name))) << row.name;
    EXPECT_NEAR(in_units(c.cell.lattice.a(), angstrom), row.lattice_constant_A, 1e-12) << row.name;
    EXPECT_EQ(c.cell.composition().size(), row.sublattices.size() == 1 ? 1u : c.cell.composition().size());
    EXPECT_EQ(c.vibrations.size(), row.debye_temperature_K > 0.0 ? 1u : 0u) << row.name;
    EXPECT_LT(rotation_error(c.rotation()), 1e-12);
  }
  // Densities from the tabulated lattice constants (g/cm3): Si 2.33, W 19.3,
  // GaAs 5.32, UO2 10.96 (handbook values, to 1 %).
  const auto density = [](const char *name) {
    return in_units(make_reference_crystal(*find_reference_crystal(name)).cell.density({}), gram_per_cubic_centimeter);
  };
  EXPECT_NEAR(density("Si"), 2.33, 0.02);
  EXPECT_NEAR(density("W"), 19.3, 0.2);
  EXPECT_NEAR(density("GaAs"), 5.32, 0.05);
  EXPECT_NEAR(density("UO2"), 10.96, 0.1);
  EXPECT_EQ(find_reference_crystal("GaAs")->structure, "zincblende");
  EXPECT_EQ(find_reference_crystal("Unobtainium"), nullptr);
  // Stoichiometry of the compounds: UO2 has two O per U, BaTiO3 three.
  const auto uo2 = make_reference_crystal(*find_reference_crystal("UO2")).cell.composition();
  EXPECT_DOUBLE_EQ(uo2.at(8) / uo2.at(92), 2.0);
  const auto bto = make_reference_crystal(*find_reference_crystal("BaTiO3")).cell.composition();
  EXPECT_DOUBLE_EQ(bto.at(8) / bto.at(22), 3.0);
}

TEST(CrystalMaterials, ThermalAmplitudeFallsBackToTheTable) {
  // What the crystal says wins; a pure element in the table gets its
  // tabulated Debye temperature; anything else stands still.
  Crystal si;
  si.cell = simple("diamond", {14}, 5.431);
  si.temperature = 293.0 * kelvin;
  const auto m = Element::by_z(14).natural_mass();
  bool from_table = false;
  auto u = thermal_amplitude(si, 14, m, &from_table);
  ASSERT_TRUE(u.has_value());
  EXPECT_TRUE(from_table);
  EXPECT_NEAR(in_units(*u, angstrom), 0.075, 0.001); // Table IV
  si.vibrations.push_back({14, std::nullopt, 0.05 * angstrom});
  u = thermal_amplitude(si, 14, m, &from_table);
  EXPECT_FALSE(from_table);
  EXPECT_NEAR(in_units(*u, angstrom), 0.05, 1e-15);
  Crystal gaas;
  gaas.cell = simple("zincblende", {31, 33}, 5.654);
  EXPECT_FALSE(thermal_amplitude(gaas, 31, Element::by_z(31).natural_mass(), &from_table).has_value()); // a compound
  Crystal ti;
  ti.cell = simple("hcp", {22}, 2.95);
  EXPECT_FALSE(thermal_amplitude(ti, 22, Element::by_z(22).natural_mass()).has_value()); // not in the table
  // Lindhard's angle: 1 MeV He along <110> of Ge (d = a/sqrt 2 = 4.00 A) is
  // 1.23 degrees (Gemmell's example); 0.307 sqrt(Z1 Z2/(d E)) in general.
  EXPECT_NEAR(in_units(lindhard_angle(2, 32, 4.00 * angstrom, 1.0 * mega_electron_volt), degree), 1.23, 0.005);
  EXPECT_NEAR(in_units(lindhard_angle(2, 74, 3.165 * angstrom, 1.0 * mega_electron_volt), degree),
              0.307 * std::sqrt(2.0 * 74.0 / 3.165), 0.005);
  EXPECT_THROW(lindhard_angle(2, 74, 0.0 * angstrom, 1.0 * mega_electron_volt), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Polycrystals: grain maps, generated grains, grain files, Neper, EBSD
// ---------------------------------------------------------------------------

// The crystal database: every row builds its structure, names known
// sources, and (for the elements) gives the density of AMCSET's element
// table within 4 %: an independent check of the transcribed lattice
// constants and structures.
TEST(CrystalDatabase, EveryRowBuildsWithItsSources) {
  const auto &rows = crystal_database();
  EXPECT_GE(rows.size(), 122u);
  std::size_t elements = 0;
  for (const auto &r : rows) {
    SCOPED_TRACE(std::string(r.name) + " " + std::string(r.structure));
    EXPECT_FALSE(database_citation(r.lattice_source).empty());
    if (r.debye_temperature_K > 0.0) {
      EXPECT_FALSE(database_citation(r.debye_source).empty());
    }
    EXPECT_EQ(r.sublattices.size(), structure_info(r.structure).sublattices.size());
    EXPECT_EQ(structure_info(r.structure).hexagonal, r.c_A > 0.0);
    const Crystal c = make_database_crystal(r);
    EXPECT_NO_THROW(c.cell.validate("database"));
    // (Not for C, whose table density is graphite's, nor Ra, whose table
    // density, 5.0 g/cm3, is an estimate.)
    if (r.sublattices.size() == 1 && r.sublattices.front() != 6 && r.sublattices.front() != 88) {
      ++elements;
      const auto &e = Element::by_z(r.sublattices.front());
      const double rho = in_units(c.cell.density({{e.z(), e.natural_mass()}}), gram_per_cubic_centimeter);
      const double table = in_units(e.solid_density(), gram_per_cubic_centimeter);
      EXPECT_NEAR(rho / table, 1.0, 0.04) << rho << " against " << table << " g/cm3";
    }
  }
  EXPECT_GT(elements, 50u); // (52 of the 54 elements)
  // Gemmell's rows come first and unchanged.
  ASSERT_GE(rows.size(), reference_crystals().size());
  for (std::size_t i = 0; i < reference_crystals().size(); ++i) {
    EXPECT_EQ(rows[i].name, reference_crystals()[i].name);
    EXPECT_EQ(rows[i].a_A, reference_crystals()[i].lattice_constant_A);
  }
}

TEST(CrystalDatabase, SearchMatchAndFind) {
  const auto names = [](const std::vector<const DatabaseCrystal *> &v) {
    std::string s;
    for (const auto *r : v) s += std::string(r->name) + " ";
    return s;
  };
  EXPECT_EQ(names(search_crystal_database("GaAs")), "GaAs ");
  EXPECT_EQ(names(search_crystal_database("ga as")), "GaAs ");
  EXPECT_EQ(names(search_crystal_database("iron")), "Fe ");
  EXPECT_EQ(search_crystal_database("nitride rocksalt").size(), 6u);
  EXPECT_GE(search_crystal_database("wurtzite").size(), 6u);
  EXPECT_EQ(search_crystal_database("").size(), crystal_database().size());
  EXPECT_TRUE(search_crystal_database("unobtainium").empty());
  // A composition picks the form most common at room temperature.
  const auto match = [](std::vector<std::pair<unsigned, double>> c) -> std::string {
    const auto *r = match_composition(c);
    return r ? std::string(r->name) + " " + std::string(r->structure) : "none";
  };
  EXPECT_EQ(match({{26, 1.0}}), "Fe bcc");
  EXPECT_EQ(match({{22, 1.0}}), "Ti hcp");
  EXPECT_EQ(match({{31, 1.0}, {33, 1.0}}), "GaAs zincblende");
  EXPECT_EQ(match({{33, 0.5}, {31, 0.5}}), "GaAs zincblende");
  EXPECT_EQ(match({{30, 1.0}, {16, 1.0}}), "ZnS zincblende");
  EXPECT_EQ(match({{48, 1.0}, {16, 1.0}}), "CdS wurtzite");
  EXPECT_EQ(match({{30, 1.0}, {8, 1.0}}), "ZnO wurtzite");
  EXPECT_EQ(match({{92, 1.0}, {8, 2.0}}), "UO2 fluorite");
  EXPECT_EQ(match({{38, 1.0}, {22, 1.0}, {8, 3.0}}), "SrTiO3 perovskite");
  EXPECT_EQ(match({{26, 0.9}, {24, 0.1}}), "none"); // an alloy
  EXPECT_EQ(match({{92, 1.0}, {8, 1.9}}), "none");  // off stoichiometry
  // Names in scenarios.
  EXPECT_EQ(find_database_crystal("zns").structure, "zincblende");
  EXPECT_EQ(find_database_crystal("ZnS wurtzite").structure, "wurtzite");
  EXPECT_EQ(find_database_crystal("silicon").name, "Si");
  EXPECT_THROW(find_database_crystal("Xx"), std::invalid_argument);
  EXPECT_THROW(find_database_crystal("nitride"), std::invalid_argument);
  // The crystal: the row's constants, 293 K, Debye temperatures if known.
  const Crystal ti = make_database_crystal(find_database_crystal("Ti"));
  EXPECT_NEAR(in_units(ti.cell.lattice.c(), angstrom), 4.6855, 1e-9);
  ASSERT_EQ(ti.vibrations.size(), 1u);
  EXPECT_EQ(in_units(*ti.vibrations.front().debye_temperature, kelvin), 420.0);
  EXPECT_TRUE(make_database_crystal(find_database_crystal("GaN")).vibrations.empty());
}

TEST(CrystalGrains, LookupMatchesBruteForce) {
  // The grain of a point and its distance to the nearest grain boundary,
  // through the buckets and by looking at every seed and its lateral
  // images: Voronoi cells, columns, Laguerre cells, and a map of three
  // seeds (where the images matter).
  std::mt19937_64 rng(7);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  for (int variant = 0; variant < 4; ++variant) {
    const bool columnar = variant == 1;
    const bool weighted = variant == 2;
    const std::size_t n = variant == 3 ? 3 : 400;
    const double py = 900.0, pz = 700.0, depth = 500.0;
    std::vector<Grain> grains(n);
    std::vector<GrainMap::Seed> seeds(n);
    for (std::size_t i = 0; i < n; ++i) {
      seeds[i].position_A = {depth * u(rng), py * u(rng), pz * u(rng)};
      seeds[i].weight_A2 = weighted ? 900.0 * u(rng) : 0.0;
      seeds[i].grain = static_cast<int>(i);
    }
    const auto map = GrainMap::from_seeds(grains, seeds, py, pz, columnar);
    EXPECT_EQ(map.grain_count(), n);
    EXPECT_EQ(map.is_columnar(), columnar);
    const double cap = 30.0;
    int close = 0;
    for (int k = 0; k < 3000; ++k) {
      // (Also above and below the seeds, and outside one period.)
      const Vec3 r{depth * (1.4 * u(rng) - 0.2), py * (3.0 * u(rng) - 1.0), pz * (3.0 * u(rng) - 1.0)};
      const auto fast = map.at(r, cap);
      const auto slow = map.at_brute_force(r);
      ASSERT_EQ(fast.grain, slow.grain) << "variant " << variant << ", point " << k;
      if (slow.boundary_A < cap) {
        ++close;
        ASSERT_NEAR(fast.boundary_A, slow.boundary_A, 1e-7) << "variant " << variant << ", point " << k;
      } else {
        ASSERT_GE(fast.boundary_A, cap - 1e-9) << "variant " << variant << ", point " << k;
      }
      // Periodic: the same grain one period away.
      EXPECT_EQ(map.at({r.x, r.y + py, r.z - 2.0 * pz}).grain, fast.grain);
    }
    if (n > 3) {
      EXPECT_GT(close, 300) << variant; // boundaries were met
    }
  }
  // Across a boundary the grain changes: a step of the boundary distance
  // (and a little more) towards the other seed leaves the cell.
  std::vector<Grain> two(2);
  const auto pair = GrainMap::from_seeds(two, {{{0, 100, 500}, 0.0, 0}, {{0, 300, 500}, 0.0, 1}}, 1000.0, 1000.0, true);
  const auto here = pair.at({50.0, 180.0, 500.0}, 100.0);
  EXPECT_EQ(here.grain, 0);
  EXPECT_NEAR(here.boundary_A, 20.0, 1e-9); // the bisector is at y = 200
  EXPECT_EQ(pair.at({50.0, 180.0 + here.boundary_A + 1e-6, 500.0}).grain, 1);
  // Mistakes are refused.
  EXPECT_THROW(GrainMap::from_seeds(two, {}, 10.0, 10.0, false), std::invalid_argument);
  EXPECT_THROW(GrainMap::from_seeds(two, {{{0, 0, 0}, 0.0, 2}}, 10.0, 10.0, false), std::invalid_argument);
  EXPECT_THROW(GrainMap::from_seeds(two, {{{0, 0, 0}, 0.0, 0}}, 0.0, 10.0, false), std::invalid_argument);
}

TEST(CrystalGrains, VoxelsAndTheirFaces) {
  // 2 x 2 x 2 voxels of 50 x 100 x 100 A: grain 0 in the upper half, 1 and
  // 2 side by side below; the last voxels continue below the grid.
  std::vector<Grain> grains(3);
  const auto map = GrainMap::from_voxels(grains, {0, 0, 0, 0, 1, 1, 2, 2}, {2, 2, 2}, 100.0, 200.0, 200.0);
  EXPECT_TRUE(map.is_voxels());
  EXPECT_EQ(map.at({10.0, 30.0, 30.0}).grain, 0);
  EXPECT_EQ(map.at({60.0, 30.0, 30.0}).grain, 1);
  EXPECT_EQ(map.at({60.0, 130.0, 30.0}).grain, 2);
  EXPECT_EQ(map.at({500.0, 130.0, 30.0}).grain, 2);
  EXPECT_EQ(map.at({60.0, 130.0 - 200.0, 30.0 + 400.0}).grain, 2); // periodic
  EXPECT_NEAR(map.at({40.0, 30.0, 30.0}).boundary_A, 10.0, 1e-12);  // the face at depth 50
  EXPECT_NEAR(map.at({80.0, 95.0, 30.0}).boundary_A, 5.0, 1e-12);   // the face at y = 100
  EXPECT_NEAR(map.at({80.0, 3.0, 30.0}).boundary_A, 3.0, 1e-12);    // the periodic face at y = 0
  EXPECT_EQ(map.ungrained_fraction(), 0.0);
  const auto holes = GrainMap::from_voxels(grains, {0, -1, -1, 1}, {1, 2, 2}, 10.0, 20.0, 20.0);
  EXPECT_EQ(holes.at({5.0, 5.0, 15.0}).grain, -1);
  EXPECT_DOUBLE_EQ(holes.ungrained_fraction(), 0.5);
  EXPECT_THROW(GrainMap::from_voxels(grains, {0, 1, 2}, {2, 2, 2}, 100.0, 200.0, 200.0), std::invalid_argument);
  EXPECT_THROW(GrainMap::from_voxels(grains, {0, 0, 0, 0, 1, 1, 2, 3}, {2, 2, 2}, 100.0, 200.0, 200.0),
               std::invalid_argument);
}

TEST(CrystalGrains, BoundariesBetweenSeedsAndVoxels) {
  // Two columns side by side in a periodic strip: they touch at y = 200 and,
  // through the repeat, at y = 0; one boundary, with a point on one of
  // those lines (within half a grid step).
  std::vector<Grain> two(2);
  const auto pair = GrainMap::from_seeds(two, {{{0, 100, 500}, 0.0, 0}, {{0, 300, 500}, 0.0, 1}}, 400.0, 1000.0, true);
  const auto b = grain_boundaries(pair, 800.0, 40'000);
  ASSERT_EQ(b.size(), 1u);
  EXPECT_EQ(b[0].a, 0);
  EXPECT_EQ(b[0].b, 1);
  EXPECT_GT(b[0].contacts, 100u);
  const double y = std::fmod(b[0].point_A.y + 400.0, 400.0);
  EXPECT_TRUE(std::abs(y - 200.0) < 2.0 || y < 2.0 || y > 398.0) << b[0].point_A.y;

  // Random Voronoi cells: the boundary points are about as far from the
  // two seeds as from each other, and nearly all pairs of touching cells
  // are found (compared with points sampled on both sides of a boundary).
  std::mt19937_64 rng(11);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  const std::size_t n = 300;
  const double py = 900.0, pz = 700.0, depth = 600.0;
  std::vector<Grain> grains(n);
  std::vector<GrainMap::Seed> seeds(n);
  for (std::size_t i = 0; i < n; ++i) {
    seeds[i].position_A = {depth * u(rng), py * u(rng), pz * u(rng)};
    seeds[i].grain = static_cast<int>(i);
  }
  const auto map = GrainMap::from_seeds(grains, seeds, py, pz, false);
  const auto found = grain_boundaries(map, depth);
  std::set<std::pair<int, int>> pairs;
  const double step = std::cbrt(depth * py * pz / 1.5e5); // the grid's (500 points per grain)
  std::size_t on_one_of_them = 0;
  const auto distance = [&](const Vec3 &p, int grain) {
    const Vec3 s = map.seeds()[static_cast<std::size_t>(grain)].position_A;
    const Vec3 d{p.x - s.x, p.y - s.y - py * std::nearbyint((p.y - s.y) / py), p.z - s.z - pz * std::nearbyint((p.z - s.z) / pz)};
    return norm(d);
  };
  for (const auto &g : found) {
    ASSERT_LT(g.a, g.b);
    pairs.insert({g.a, g.b});
    EXPECT_GE(g.point_A.x, 0.0);
    EXPECT_LE(g.point_A.x, depth);
    EXPECT_NEAR(distance(g.point_A, g.a), distance(g.point_A, g.b), 2.0 * step);
    const int here = map.at(g.point_A).grain;
    on_one_of_them += here == g.a || here == g.b ? 1 : 0;
  }
  EXPECT_GT(on_one_of_them, found.size() * 95 / 100);
  std::set<std::pair<int, int>> touching;
  for (int k = 0; k < 20000; ++k) {
    const Vec3 r{depth * u(rng), py * u(rng), pz * u(rng)};
    const auto at = map.at(r, 3.0);
    if (at.boundary_A >= 3.0) continue;
    for (const Vec3 &s : {Vec3{3, 0, 0}, Vec3{-3, 0, 0}, Vec3{0, 3, 0}, Vec3{0, -3, 0}, Vec3{0, 0, 3}, Vec3{0, 0, -3}}) {
      const Vec3 q = r + s;
      if (q.x < 0.0 || q.x > depth) continue;
      const int other = map.at(q).grain;
      if (other != at.grain) touching.insert({std::min(at.grain, other), std::max(at.grain, other)});
    }
  }
  std::size_t missed = 0;
  for (const auto &t : touching) missed += pairs.count(t) == 0 ? 1 : 0;
  EXPECT_GT(touching.size(), 500u);
  EXPECT_LT(static_cast<double>(missed), 0.02 * static_cast<double>(touching.size())) << missed;

  // Voxels: grain 0 above 1 and 2 (1 | 2 side by side, periodic in y).
  std::vector<Grain> three(3);
  const auto vox = GrainMap::from_voxels(three, {0, 0, 0, 0, 1, 1, 2, 2}, {2, 2, 2}, 100.0, 200.0, 200.0);
  const auto vb = grain_boundaries(vox, 100.0);
  ASSERT_EQ(vb.size(), 3u); // 0-1, 0-2, 1-2
  EXPECT_EQ(vb[0].a, 0);
  EXPECT_EQ(vb[0].b, 1);
  EXPECT_EQ(vb[0].contacts, 2u);
  EXPECT_DOUBLE_EQ(vb[0].point_A.x, 50.0);
  EXPECT_EQ(vb[2].a, 1);
  EXPECT_EQ(vb[2].b, 2);
  EXPECT_EQ(vb[2].contacts, 4u); // at y = 100 and, through the repeat, at y = 200
}

TEST(CrystalGrains, SummaryOfGeneratedGrains) {
  Crystal c;
  c.cell = simple("fcc", {29}, 3.615);
  c.orientation.surface = {0, 0, 1};
  Polycrystal p;
  p.grain_size = 100.0 * angstrom;
  p.period_y = p.period_z = 800.0 * angstrom;
  c.polycrystal = p;
  const double depth = 800.0;
  const auto random = summarize_grains(*make_grain_map(c, depth), c.cell.lattice, c.rotation(), depth, 200'000);
  // Spheres of 100 A in 800^3 A^3 (and the grains just above and below).
  EXPECT_EQ(random.grains, static_cast<std::size_t>(std::llround(1000.0 * 800.0 * 800.0 / (pi_v / 6.0 * 1e6))));
  EXPECT_GT(random.grains_in_layer, random.grains * 7 / 10);
  EXPECT_LT(random.grains_in_layer, random.grains);
  EXPECT_EQ(random.ungrained_fraction, 0.0);
  EXPECT_FALSE(random.coarse);
  EXPECT_GT(random.diameter_mean_A, 70.0); // (grains cut by the layer's faces are smaller)
  EXPECT_LT(random.diameter_mean_A, 105.0);
  EXPECT_LE(random.diameter_p10_A, random.diameter_median_A);
  EXPECT_LE(random.diameter_median_A, random.diameter_p90_A);
  EXPECT_NEAR(random.misorientation_mean_deg, 40.7, 3.0); // Mackenzie's mean for cubic crystals
  EXPECT_TRUE(random.boundaries_counted);
  EXPECT_GT(random.boundaries, random.grains_in_layer * 4); // (about 7 per grain inside, fewer at the faces)
  EXPECT_LT(random.low_angle_fraction, 0.1);
  // A mosaic without spread: the crystal's orientation everywhere.
  c.polycrystal->texture = Polycrystal::Texture::mosaic;
  const auto same = summarize_grains(*make_grain_map(c, depth), c.cell.lattice, c.rotation(), depth, 100'000);
  EXPECT_NEAR(same.misorientation_mean_deg, 0.0, 1e-5);
  EXPECT_NEAR(same.normal_tilt_mean_deg, 0.0, 1e-5);
  EXPECT_DOUBLE_EQ(same.low_angle_fraction, 1.0);
  // A fiber texture: the normal direction stays, the rest turns.
  c.polycrystal->texture = Polycrystal::Texture::fiber;
  const auto fiber = summarize_grains(*make_grain_map(c, depth), c.cell.lattice, c.rotation(), depth, 100'000);
  EXPECT_NEAR(fiber.normal_tilt_mean_deg, 0.0, 1e-5);
  EXPECT_GT(fiber.misorientation_mean_deg, 10.0);
  // The JSON form names the fields.
  const auto j = to_json(fiber);
  EXPECT_NE(j.find("\"grains_in_layer\": "), std::string::npos);
  EXPECT_NE(j.find("\"low_angle_fraction\": "), std::string::npos);
  // Columns are measured over the surface.
  c.polycrystal->columnar = true;
  c.polycrystal->grain_size = 200.0 * angstrom;
  c.polycrystal->period_y = c.polycrystal->period_z = 1600.0 * angstrom;
  const auto columns = summarize_grains(*make_grain_map(c, 5000.0), c.cell.lattice, c.rotation(), 5000.0, 100'000);
  EXPECT_TRUE(columns.columnar);
  EXPECT_EQ(columns.grains_in_layer, columns.grains);
  EXPECT_NEAR(columns.diameter_mean_A, 200.0, 25.0);
}

TEST(CrystalGrains, LowAngleBoundariesInRegistry) {
  // Grains of one orientation (a mosaic without spread), put in registry:
  // every grain's lattice is one lattice, so the polycrystal is the single
  // crystal. Points in the layer, also beyond one period.
  Crystal c;
  c.cell = simple("diamond", {14}, 5.431);
  c.orientation.surface = {0, 0, 1};
  Polycrystal p;
  p.grain_size = 150.0 * angstrom;
  p.texture = Polycrystal::Texture::mosaic;
  p.registry_below = 2.0 * degree;
  c.polycrystal = p;
  const double depth = 600.0;
  RegistryResult rr;
  const auto map = make_grain_map(c, depth, &rr);
  EXPECT_EQ(rr.groups, 1u);
  const Mat3 r = c.rotation();
  const auto fractional = [&](const Vec3 &v) { return c.cell.lattice.to_fractional(r.transposed() * v); };
  std::mt19937_64 rng(5);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  const double py = map->period_y_A(), pz = map->period_z_A();
  const auto origin = [&](const Vec3 &at) { return map->grain(map->at(at).grain).origin_A; };
  const Vec3 reference = origin({depth / 2, py / 2, pz / 2});
  std::set<int> met;
  for (int k = 0; k < 2000; ++k) {
    const Vec3 q{depth * u(rng), py * (3.0 * u(rng) - 1.0), pz * (3.0 * u(rng) - 1.0)};
    met.insert(map->at(q).grain);
    const Vec3 f = fractional(origin(q) - reference);
    for (int i = 0; i < 3; ++i) ASSERT_NEAR(f[i], std::round(f[i]), 1e-6) << k;
  }
  EXPECT_GT(met.size(), 50u);
  EXPECT_GE(rr.grains, met.size());
  // Without the registry the grains' lattices go through their seeds.
  c.polycrystal->registry_below = 0.0 * degree;
  const auto free = make_grain_map(c, depth, &rr);
  EXPECT_EQ(rr.groups, 0u);
  EXPECT_EQ(free->grains()[3].origin_A, free->seeds()[3].position_A);

  // Random orientations: misoriented grains are not joined; a large
  // threshold joins some.
  c.polycrystal->texture = Polycrystal::Texture::random;
  c.polycrystal->registry_below = 1.0 * degree;
  make_grain_map(c, depth, &rr);
  EXPECT_EQ(rr.groups, 0u);
  c.polycrystal->registry_below = 25.0 * degree;
  make_grain_map(c, depth, &rr);
  EXPECT_GT(rr.groups, 0u);

  // A small spread: at the point of each boundary in a group the child's
  // lattice point found is the parent's: lattice points of both grains lie
  // within the drift of the misorientation over a cell there.
  c.polycrystal->texture = Polycrystal::Texture::mosaic;
  c.polycrystal->spread = 0.5 * degree;
  c.polycrystal->registry_below = 2.0 * degree;
  const auto spread = make_grain_map(c, depth, &rr);
  EXPECT_EQ(rr.groups, 1u);
  std::size_t boundaries = 0;
  std::size_t checked = 0;
  for (const auto &b : grain_boundaries(*spread, depth)) {
    const auto &ga = spread->grain(b.a);
    const auto &gb = spread->grain(b.b);
    ++boundaries;
    // The equivalent orientations chosen are those nearest each other.
    EXPECT_LT(misorientation_rad(ga.rotation, gb.rotation, {}), 4.0 * pi_v / 180.0);
    // b's lattice point nearest the boundary point, seen in a's lattice.
    const Vec3 fb = c.cell.lattice.to_fractional(gb.rotation.transposed() * (b.point_A - gb.origin_A));
    const Vec3 nb{std::round(fb.x), std::round(fb.y), std::round(fb.z)};
    const Vec3 point = gb.origin_A + gb.rotation * c.cell.lattice.to_cartesian_A(nb);
    const Vec3 fa = c.cell.lattice.to_fractional(ga.rotation.transposed() * (point - ga.origin_A));
    const Vec3 off = c.cell.lattice.to_cartesian_A(fa - Vec3{std::round(fa.x), std::round(fa.y), std::round(fa.z)});
    if (norm(off) < 0.5) ++checked; // in registry there
  }
  EXPECT_GT(boundaries, 50u);
  // The boundaries through which the group was built (grains - groups of
  // them) are exactly in registry; the others are where the lattice drifts.
  EXPECT_GE(checked, rr.grains - rr.groups);
  std::printf("[ registry ] %zu of %zu boundaries in registry (tree edges %zu)\n", checked, boundaries,
              rr.grains - rr.groups);
}

TEST(CrystalSearch, SharedSitesGiveTheSameLattice) {
  // A lattice built from another (shared sites) finds the same atoms as one
  // built in full with the same rotation, number for number.
  const UnitCell cell = simple("diamond", {14}, 5.431);
  const Mat3 first = rotation_about({1, 0, 0}, 0.3);
  const Mat3 second = rotation_about({1, 2, 3}, 1.1);
  const std::map<unsigned, double> rms{{14u, 0.08}};
  const LatticeSearch shape(cell, first, rms, 1.5, 1.5);
  const LatticeSearch shared(shape, cell, second);
  const LatticeSearch full(cell, second, rms, 1.5, 1.5);
  EXPECT_EQ(shared.sites_per_search(), full.sites_per_search());
  EXPECT_EQ(shared.max_displacement_A(), full.max_displacement_A());
  std::vector<Candidate> a, b;
  std::mt19937_64 rng(2);
  std::uniform_real_distribution<double> u(-20.0, 20.0);
  std::size_t seen = 0;
  for (int k = 0; k < 300; ++k) {
    const Vec3 r{u(rng), u(rng), u(rng)};
    const Vec3 d = normalized({u(rng), u(rng), u(rng)});
    shared.candidates(r, d, 77, -1e9, 1e9, a);
    full.candidates(r, d, 77, -1e9, 1e9, b);
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
      EXPECT_EQ(a[i].id, b[i].id);
      EXPECT_EQ(a[i].position, b[i].position);
      EXPECT_EQ(a[i].xi, b[i].xi);
    }
    seen += a.size();
  }
  EXPECT_GT(seen, 100u);
  EXPECT_THROW(LatticeSearch(shape, simple("fcc", {29}, 3.6), second), std::invalid_argument);
}

TEST(CrystalGrains, LatticeRotationsAndMisorientation) {
  const auto ninety = 90.0 * degree;
  const auto lattice = [&](double a, double b, double c) {
    return Lattice::from_parameters(a * angstrom, b * angstrom, c * angstrom, ninety, ninety, ninety);
  };
  const auto triclinic = Lattice::from_parameters(3.0 * angstrom, 4.0 * angstrom, 5.0 * angstrom, 80.0 * degree,
                                                  85.0 * degree, 95.0 * degree);
  const std::vector<std::pair<Lattice, std::size_t>> cases{{lattice(4, 4, 4), 24},
                                                           {simple("hcp", {22}, 2.95).lattice, 12},
                                                           {lattice(3, 3, 4), 8},
                                                           {lattice(3, 4, 5), 4},
                                                           {triclinic, 1}};
  for (const auto &[l, count] : cases) {
    const auto group = lattice_rotations(l);
    EXPECT_EQ(group.size(), count);
    for (const auto &s : group) {
      EXPECT_LT(rotation_error(s), 1e-9);
      // Each lattice vector is turned into a lattice vector.
      for (int i = 0; i < 3; ++i) {
        const Vec3 f = l.to_fractional(s * l.vectors_A()[i]);
        for (int k = 0; k < 3; ++k) EXPECT_NEAR(f[k], std::round(f[k]), 1e-9);
      }
    }
  }
  // A cubic crystal turned by 90 degrees about a cube axis is the same
  // crystal; by 45 degrees it is 45 degrees away. Without symmetry the
  // angle is the rotation's own.
  const auto cubic = lattice_rotations(lattice(4, 4, 4));
  const double d = pi_v / 180.0;
  const Mat3 a = rotation_about({1, 2, 3}, 0.7);
  EXPECT_NEAR(misorientation_rad(a, a * rotation_about({0, 0, 1}, 90.0 * d), cubic), 0.0, 1e-7);
  EXPECT_NEAR(misorientation_rad(a, a * rotation_about({0, 0, 1}, 45.0 * d), cubic), 45.0 * d, 1e-9);
  EXPECT_NEAR(misorientation_rad(a, a * rotation_about({0, 0, 1}, 90.0 * d), {}), 90.0 * d, 1e-9);
  EXPECT_NEAR(misorientation_rad(a, a * rotation_about({1, 1, 1}, 120.0 * d), cubic), 0.0, 1e-7);
  // No two cubic orientations are farther apart than 62.8 degrees (Mackenzie).
  std::mt19937_64 rng(3);
  std::uniform_real_distribution<double> u(0.0, 1.0);
  double largest = 0.0;
  double mean = 0.0;
  const int n = 4000;
  for (int i = 0; i < n; ++i) {
    const Mat3 r = uniform_rotation(u(rng), u(rng), u(rng));
    EXPECT_LT(rotation_error(r), 1e-12);
    const double m = misorientation_rad(a, r, cubic);
    largest = std::max(largest, m);
    mean += m / n;
  }
  EXPECT_LE(largest, 62.8 * d);
  EXPECT_GT(largest, 55.0 * d);
  EXPECT_NEAR(mean, 40.7 * d, 1.0 * d); // the mean of Mackenzie's distribution
}

TEST(CrystalGrains, GeneratedGrains) {
  Crystal c;
  c.cell = simple("fcc", {29}, 3.615);
  c.orientation.surface = {1, 1, 1};
  c.orientation.reference = {1, -1, 0};
  Polycrystal p;
  p.grain_size = 100.0 * angstrom;
  p.period_y = p.period_z = 1000.0 * angstrom;
  c.polycrystal = p;
  // Seeds also one grain size above and below the 800 A: 1000 A x 1000^2 A^2
  // of spheres of 100 A.
  const auto map = make_grain_map(c, 800.0);
  const std::size_t expected = static_cast<std::size_t>(std::llround(1000.0 * 1e6 / (pi_v / 6.0 * 1e6)));
  EXPECT_EQ(map->grain_count(), expected);
  EXPECT_EQ(map->seeds().size(), expected);
  EXPECT_DOUBLE_EQ(map->period_y_A(), 1000.0);
  // The same description, the same grains; another seed, others.
  const auto again = make_grain_map(c, 800.0);
  EXPECT_EQ(again->seeds()[17].position_A, map->seeds()[17].position_A);
  EXPECT_EQ(again->grains()[17].rotation.row, map->grains()[17].rotation.row);
  c.polycrystal->seed = 2;
  EXPECT_FALSE(make_grain_map(c, 800.0)->seeds()[17].position_A == map->seeds()[17].position_A);
  c.polycrystal->seed = 1;
  // Random orientations: rotations, with no preferred direction along the depth.
  double mean = 0.0;
  for (const auto &g : map->grains()) {
    EXPECT_LT(rotation_error(g.rotation), 1e-12);
    mean += g.rotation[0][0] / static_cast<double>(expected);
  }
  EXPECT_LT(std::abs(mean), 4.0 / std::sqrt(3.0 * static_cast<double>(expected)));
  // Every point of the layer is in a grain, and each grain's lattice goes
  // through its seed.
  EXPECT_GE(map->at({400.0, 123.0, 987.0}).grain, 0);
  EXPECT_EQ(map->grains()[5].origin_A, map->seeds()[5].position_A);

  // A fiber texture: (1 1 1) along the depth in every grain, turned about it.
  const Mat3 base = c.rotation();
  c.polycrystal->texture = Polycrystal::Texture::fiber;
  const auto fiber = make_grain_map(c, 800.0);
  int turned = 0;
  for (const auto &g : fiber->grains()) {
    EXPECT_LT(norm(g.rotation.transposed() * Vec3{1, 0, 0} - base.transposed() * Vec3{1, 0, 0}), 1e-12);
    turned += misorientation_rad(base, g.rotation, {}) > 0.01 ? 1 : 0;
  }
  EXPECT_GT(turned, static_cast<int>(expected) * 9 / 10);
  // A mosaic crystal: the crystal's orientation within the spread.
  c.polycrystal->texture = Polycrystal::Texture::mosaic;
  c.polycrystal->spread = 2.0 * degree;
  const auto mosaic = make_grain_map(c, 800.0);
  double square = 0.0;
  for (const auto &g : mosaic->grains()) {
    const double m = misorientation_rad(base, g.rotation, {});
    square += m * m / static_cast<double>(expected);
  }
  EXPECT_NEAR(std::sqrt(square) * 180.0 / pi_v, 2.0, 0.2);
  c.polycrystal->spread = 0.0 * degree;
  const auto exact = make_grain_map(c, 800.0); // (kept: the loop below refers into it)
  for (const auto &g : exact->grains()) EXPECT_EQ(g.rotation.row, base.row);

  // Columns: 4/pi of the area per grain; the default period is eight grains.
  Polycrystal columns;
  columns.grain_size = 200.0 * angstrom;
  columns.columnar = true;
  const auto col = generate_voronoi(columns, Mat3::identity(), 5000.0);
  EXPECT_DOUBLE_EQ(col.period_z_A(), 1600.0);
  EXPECT_EQ(col.grain_count(), static_cast<std::size_t>(std::llround(1600.0 * 1600.0 / (0.25 * pi_v * 200.0 * 200.0))));
  EXPECT_EQ(col.at({10.0, 300.0, 400.0}).grain, col.at({4990.0, 300.0, 400.0}).grain);

  // Limits and names.
  Polycrystal bad;
  bad.grain_size = 10.0 * angstrom;
  bad.period_y = bad.period_z = 1.0e6 * angstrom;
  EXPECT_THROW(generate_voronoi(bad, Mat3::identity(), 1.0e6), std::invalid_argument); // too many grains
  bad = Polycrystal{};
  bad.grain_size = 1.0 * angstrom;
  EXPECT_THROW(check_polycrystal(bad, "p"), std::invalid_argument);
  bad = Polycrystal{};
  bad.boundary_width = 500.0 * angstrom;
  EXPECT_THROW(check_polycrystal(bad, "p"), std::invalid_argument);
  bad = Polycrystal{};
  bad.source = Polycrystal::Source::ebsd;
  EXPECT_THROW(check_polycrystal(bad, "p"), std::invalid_argument); // no file
  EXPECT_NO_THROW(check_polycrystal(Polycrystal{}, "p"));
  for (auto s : {Polycrystal::Source::voronoi, Polycrystal::Source::grain_file, Polycrystal::Source::neper,
                 Polycrystal::Source::ebsd}) {
    EXPECT_EQ(parse_polycrystal_source(to_string(s)), s);
  }
  for (auto t : {Polycrystal::Texture::random, Polycrystal::Texture::fiber, Polycrystal::Texture::mosaic}) {
    EXPECT_EQ(parse_texture(to_string(t)), t);
  }
  EXPECT_THROW(parse_texture("cube"), std::invalid_argument);
  EXPECT_EQ(grain_source_from_path("map.CTF"), Polycrystal::Source::ebsd);
  EXPECT_EQ(grain_source_from_path("n8.tess"), Polycrystal::Source::neper);
  EXPECT_EQ(grain_source_from_path("grains.txt"), Polycrystal::Source::grain_file);
  Crystal single;
  single.cell = c.cell;
  EXPECT_THROW(make_grain_map(single, 100.0), std::invalid_argument);
}

TEST(CrystalGrains, GrainFiles) {
  const UnitCell cell = simple("bcc", {26}, 2.87);
  Polycrystal p;
  p.source = Polycrystal::Source::grain_file;
  // Seeds: two grains side by side in y, the second turned by 30 degrees
  // about the surface normal.
  const auto seeds = parse_grain_file("# two grains\n"
                                      "unit nm\n"
                                      "period 40 20\n"
                                      "grain 5 10 10   0 0 0\n"
                                      "grain 5, 30, 10, 30, 0, 0   # commas are fine\n",
                                      p, cell);
  EXPECT_EQ(seeds.grain_count(), 2u);
  EXPECT_DOUBLE_EQ(seeds.period_y_A(), 400.0);
  EXPECT_DOUBLE_EQ(seeds.period_z_A(), 200.0);
  EXPECT_EQ(seeds.at({50.0, 150.0, 100.0}).grain, 0);
  EXPECT_EQ(seeds.at({50.0, 250.0, 100.0}).grain, 1);
  EXPECT_NEAR(seeds.at({50.0, 150.0, 100.0}, 100.0).boundary_A, 50.0, 1e-9);
  // Bunge angles as in crystal::rotation(): (0, 0, 0) is the cube
  // orientation of the sample, and phi1 turns about the surface normal.
  Orientation cube;
  cube.kind = OrientationKind::euler;
  EXPECT_EQ(seeds.grains()[0].rotation.row, rotation(cell.lattice, cube).row);
  cube.euler_deg = {30.0, 0.0, 0.0};
  for (int i = 0; i < 3; ++i) {
    EXPECT_LT(norm(seeds.grains()[1].rotation[i] - rotation(cell.lattice, cube)[i]), 1e-12);
  }
  EXPECT_NEAR(misorientation_rad(seeds.grains()[0].rotation, seeds.grains()[1].rotation, {}), 30.0 * pi_v / 180.0, 1e-12);
  // The scenario's period wins over the file's.
  Polycrystal wider = p;
  wider.period_y = 800.0 * angstrom;
  EXPECT_DOUBLE_EQ(parse_grain_file("period 400 200\ngrain 0 0 0 0 0 0\n", wider, cell).period_y_A(), 800.0);

  // Voxels.
  const auto voxels = parse_grain_file("unit A\nperiod 200 200\ndepth 100\nvoxels 2 2 2\n"
                                       "orientation 0 0 0\norientation 45 0 0\norientation 0 45 0\n"
                                       "ids\n0 0 0 0\n1 1 2 -1\n",
                                       p, cell);
  EXPECT_TRUE(voxels.is_voxels());
  EXPECT_EQ(voxels.grain_count(), 3u);
  EXPECT_EQ(voxels.at({60.0, 30.0, 30.0}).grain, 1);
  EXPECT_EQ(voxels.at({60.0, 130.0, 130.0}).grain, -1);
  EXPECT_DOUBLE_EQ(voxels.ungrained_fraction(), 0.125);

  // Mistakes name the line.
  const auto message = [&](const char *text) {
    try {
      parse_grain_file(text, p, cell, "g.txt");
    } catch (const std::invalid_argument &e) {
      return std::string(e.what());
    }
    return std::string();
  };
  EXPECT_NE(message("period 10 10\ngrain 1 2 3 0 0\n").find("g.txt, line 2"), std::string::npos);
  EXPECT_NE(message("period 10 10\ngrain 1 2 x 0 0 0\n").find("'x' is not a number"), std::string::npos);
  EXPECT_NE(message("grain 1 2 3 0 0 0\n").find("periods are missing"), std::string::npos);
  EXPECT_NE(message("period 10 10\ncolour red\n").find("unknown entry 'colour'"), std::string::npos);
  EXPECT_NE(message("period 10 10\nunit furlong\n").find("unknown unit"), std::string::npos);
  EXPECT_NE(message("period 10 10\n").find("no grains"), std::string::npos);
  EXPECT_NE(message("period 10 10\ndepth 5\nvoxels 1 1 2\norientation 0 0 0\nids\n0\n").find("do not fill"),
            std::string::npos);
  EXPECT_NE(message("period 10 10\nvoxels 1 1 1\norientation 0 0 0\nids\n0\n").find("depth"), std::string::npos);
}

TEST(CrystalGrains, NeperTessellation) {
  // Two cells in a unit cube, cut at x = 0.5 (written after the format of
  // Neper's .tess files; not produced by Neper).
  const UnitCell cell = simple("fcc", {28}, 3.524);
  const auto tess = [](const std::string &ori) {
    return "***tess\n **format\n   3.4\n **general\n   3 standard\n **cell\n   2\n  *id\n   1 2\n  *seed\n"
           "   1 0.25 0.5 0.5 0.0\n   2 0.75 0.5 0.5 0.0\n  *ori\n" +
           ori +
           " **vertex\n   8\n   1 0 0 0 0\n   2 1 0 0 0\n   3 0 1 0 0\n   4 1 1 0 0\n   5 0 0 1 0\n   6 1 0 1 0\n"
           "   7 0 1 1 0\n   8 1 1 1 0\n***end\n";
  };
  Polycrystal p;
  p.source = Polycrystal::Source::neper;
  p.scale = 1000.0 * angstrom;
  const auto map = parse_neper_tess(tess("   euler-bunge\n   0 0 0\n   30 0 0\n"), p, cell);
  EXPECT_EQ(map.grain_count(), 2u);
  EXPECT_DOUBLE_EQ(map.period_y_A(), 1000.0);
  EXPECT_DOUBLE_EQ(map.period_z_A(), 1000.0);
  // Sample x is the layer's y: the cells meet at y = 500 A.
  EXPECT_EQ(map.at({500.0, 100.0, 500.0}).grain, 0);
  EXPECT_EQ(map.at({500.0, 900.0, 500.0}).grain, 1);
  EXPECT_NEAR(map.at({500.0, 400.0, 500.0}, 500.0).boundary_A, 100.0, 1e-9);
  // Sample z is the outward normal: the top of the box is the surface.
  EXPECT_NEAR(map.seeds()[0].position_A.x, 500.0, 1e-9);
  Orientation o;
  o.kind = OrientationKind::euler;
  o.euler_deg = {30.0, 0.0, 0.0};
  const Mat3 turned = rotation(cell.lattice, o);
  for (int i = 0; i < 3; ++i) EXPECT_LT(norm(map.grains()[1].rotation[i] - turned[i]), 1e-12);
  // The same orientations by the other descriptors, as Neper's documentation
  // gives them for its example (30 degrees about a sample axis, passive
  // unless marked; here about z): Rodrigues tan 15 deg, quaternion (cos 15,
  // sin 15), the matrix g = Bunge's, Kocks (Psi, Theta, phi) = (phi1 - 90,
  // Phi, 90 - phi2); the active values of the opposite rotation.
  const std::vector<std::string> same{
      "   euler-bunge:passive\n   0 0 0\n   30 0 0\n",
      "   euler-bunge:active\n   0 0 0\n   -30 0 0\n",
      "   euler-kocks\n   270 0 90\n   300 0 90\n",
      "   euler-roe\n   270 0 90\n   300 0 90\n",
      "   rodrigues\n   0 0 0\n   0 0 0.267949192431\n",
      "   rodrigues:passive\n   0 0 0\n   0 0 0.267949192431\n",
      "   rodrigues:active\n   0 0 0\n   0 0 -0.267949192431\n",
      "   quaternion:passive\n   1 0 0 0\n   0.965925826289 0 0 0.258819045103\n",
      "   axis-angle\n   0 0 1 0\n   0 0 1 30\n",
      "   rotmat\n   1 0 0 0 1 0 0 0 1\n   0.866025403784 0.5 0 -0.5 0.866025403784 0 0 0 1\n",
      "   rotmat:active\n   1 0 0 0 1 0 0 0 1\n   0.866025403784 -0.5 0 0.5 0.866025403784 0 0 0 1\n"};
  for (const auto &text : same) {
    const auto other = parse_neper_tess(tess(text), p, cell);
    for (int i = 0; i < 3; ++i) {
      EXPECT_LT(norm(other.grains()[1].rotation[i] - turned[i]), 1e-9) << text;
      EXPECT_LT(norm(other.grains()[0].rotation[i] - map.grains()[0].rotation[i]), 1e-9) << text;
    }
  }
  // Laguerre weights (radii): a larger first cell moves the boundary.
  std::string weighted = tess("   euler-bunge\n   0 0 0\n   30 0 0\n");
  weighted.replace(weighted.find("0.25 0.5 0.5 0.0"), 16, "0.25 0.5 0.5 0.3");
  const auto laguerre = parse_neper_tess(weighted, p, cell);
  // |y - 250|^2 - 300^2 = |y - 750|^2 at y = 590 A.
  EXPECT_NEAR(laguerre.at({500.0, 500.0, 500.0}, 500.0).boundary_A, 90.0, 1e-6);
  EXPECT_EQ(laguerre.at({500.0, 580.0, 500.0}).grain, 0);
  // Mistakes.
  EXPECT_THROW(parse_neper_tess("***tesr\n", p, cell), std::invalid_argument);
  EXPECT_THROW(parse_neper_tess(tess("   label\n   Cube\n   Goss\n"), p, cell), std::invalid_argument);
  EXPECT_THROW(parse_neper_tess(tess("   rodrigues:sideways\n   0 0 0\n   0 0 0\n"), p, cell), std::invalid_argument);
  EXPECT_THROW(parse_neper_tess(tess("   rotmat\n   1 0 0 0 1 0 0 0 1\n   2 0 0 0 1 0 0 0 1\n"), p, cell),
               std::invalid_argument);
  std::string no_ori = tess("");
  no_ori.erase(no_ori.find("  *ori\n"), 7);
  EXPECT_THROW(parse_neper_tess(no_ori, p, cell), std::invalid_argument);
}

namespace {

std::string read_text(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

} // namespace

TEST(CrystalGrains, PublicSampleFiles) {
  // Files written by Neper, by MTEX's EBSD test-data generator and from a
  // measured map (orix), in the folder AMCSET_PUBLIC_GRAIN_DATA names
  // (docs/TESTING.md lists the files and where they come from).
  const char *root_env = std::getenv("AMCSET_PUBLIC_GRAIN_DATA");
  if (!root_env) GTEST_SKIP() << "AMCSET_PUBLIC_GRAIN_DATA not set";
  const std::filesystem::path root(root_env);
  const UnitCell cubic = simple("fcc", {29}, 3.615);
  const double d = pi_v / 180.0;

  // Neper, "-ori fibre(x,1,1,1)": a <111> axis of every cell along sample x
  // (the layer's y; rows of crystal -> laboratory are the laboratory axes in
  // the crystal frame).
  Polycrystal neper;
  neper.source = Polycrystal::Source::neper;
  const auto fibre = parse_neper_tess(read_text(root / "neper/ori/ref.tess"), neper, cubic);
  ASSERT_EQ(fibre.grain_count(), 2u);
  for (const auto &g : fibre.grains()) {
    const Vec3 along = normalized(g.rotation[1]);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(std::abs(along[k]), 1.0 / std::sqrt(3.0), 1e-6);
  }
  // "-ori file(ori)": the cells keep the Rodrigues vectors of the file
  // (passive: R(r) takes the crystal's axes to the sample's).
  const auto given = parse_neper_tess(read_text(root / "neper/ori_file/ref.tess"), neper, cubic);
  ASSERT_EQ(given.grain_count(), 3u);
  const Vec3 r[3] = {{0.1, 0.2, 0.3}, {0.4, 0.5, 0.6}, {0.7, 0.8, 0.9}};
  for (int i = 0; i < 3; ++i) {
    const Mat3 expected = sample_to_laboratory() * rotation_about(r[i], 2.0 * std::atan(norm(r[i])));
    for (int k = 0; k < 3; ++k) EXPECT_LT(norm(given.grain(i).rotation[k] - expected[k]), 1e-9) << i;
  }
  // Quaternions and the 2D format (columns).
  EXPECT_EQ(parse_neper_tess(read_text(root / "neper/oridescriptor/ref.tess"), neper, cubic).grain_count(), 2u);
  const auto flat = parse_neper_tess(read_text(root / "neper/format_2d/ref.tess"), neper, cubic);
  EXPECT_TRUE(flat.is_columnar());
  EXPECT_EQ(flat.grain_count(), 2u);

  // MTEX's generated maps: 2 x 2 grains of 15 x 15 points, square and
  // hexagonal grids; each point's grain column and row are in its MAD and
  // BC columns, its Euler angles within 1 degree of its grain's.
  UnitCell monoclinic = cubic; // (glaucophane's lattice: no rotations of its own)
  monoclinic.lattice = Lattice::from_parameters(9.541 * angstrom, 17.74 * angstrom, 5.295 * angstrom, 90.0 * degree,
                                                103.67 * degree, 90.0 * degree);
  Polycrystal ebsd;
  ebsd.source = Polycrystal::Source::ebsd;
  for (const char *name : {"mtex/testdata_sqr.ctf", "mtex/testdata_hex.ctf"}) {
    const std::string text = read_text(root / name);
    const auto map = parse_ebsd(text, true, ebsd, monoclinic, name);
    ASSERT_EQ(map.grain_count(), 4u) << name;
    ASSERT_EQ(map.seeds().size(), 900u) << name;
    // The points in file order with their true grain (column + 2 row).
    std::vector<int> truth;
    std::vector<std::array<double, 3>> euler;
    std::istringstream in(text);
    std::string line;
    bool data = false;
    while (std::getline(in, line)) {
      if (!data) {
        data = line.rfind("Phase\tX", 0) == 0;
        continue;
      }
      std::istringstream f(line);
      double v[10];
      int count = 0;
      while (count < 10 && f >> v[count]) ++count;
      if (count < 10) continue;
      truth.push_back(static_cast<int>(v[8]) + 2 * static_cast<int>(v[9]));
      euler.push_back({v[5], v[6], v[7]});
    }
    ASSERT_EQ(truth.size(), 900u) << name;
    std::map<int, int> to_truth;
    for (std::size_t i = 0; i < truth.size(); ++i) {
      const int g = map.seeds()[i].grain;
      const auto it = to_truth.try_emplace(g, truth[i]).first;
      ASSERT_EQ(it->second, truth[i]) << name << ": point " << i << " is in the wrong grain";
    }
    EXPECT_EQ(to_truth.size(), 4u) << name;
    // A grain's orientation is one of its points': within 3.5 degrees of every other.
    for (std::size_t i = 0; i < truth.size(); i += 37) {
      Orientation o;
      o.kind = OrientationKind::euler;
      o.euler_deg = euler[i];
      const Mat3 point = rotation(monoclinic.lattice, o);
      EXPECT_LT(misorientation_rad(point, map.grain(map.seeds()[i].grain).rotation, {}), 3.5 * d) << name << " " << i;
    }
  }

  // A measured map (orix's super-duplex steel, dictionary indexed): two
  // phases on a square grid of 117 x 100 points 1.5 um apart.
  const std::string sdss = read_text(root / "orix/sdss_ferrite_austenite.ang");
  const auto both = parse_ebsd(sdss, false, ebsd, cubic, "sdss");
  EXPECT_EQ(both.seeds().size(), 11700u);
  EXPECT_NEAR(both.period_y_A(), 117 * 1.5e4, 1e-6);
  EXPECT_NEAR(both.period_z_A(), 100 * 1.5e4, 1e-6);
  Polycrystal austenite = ebsd;
  austenite.phase = 1;
  Polycrystal ferrite = ebsd;
  ferrite.phase = 2;
  const auto a = parse_ebsd(sdss, false, austenite, cubic, "sdss");
  const auto f = parse_ebsd(sdss, false, ferrite, cubic, "sdss");
  EXPECT_GT(both.grain_count(), 10u);
  EXPECT_LT(a.grain_count(), both.grain_count());
  EXPECT_LT(f.grain_count(), both.grain_count());
  std::printf("[ sdss ] grains: both phases %zu, austenite %zu, ferrite %zu\n", both.grain_count(), a.grain_count(),
              f.grain_count());
  // A map without a step in its header (MTEX's ACOM.ang): the step is the grid's.
  const auto acom = parse_ebsd(read_text(root / "mtex/ACOM.ang"), false, ebsd, cubic, "ACOM");
  EXPECT_EQ(acom.seeds().size(), 225u);
  EXPECT_NEAR(acom.period_y_A(), 15 * 2.0e4, 1e-6);
}

TEST(CrystalGrains, EbsdMaps) {
  // A map of 12 x 6 points, 0.5 um apart: the left half one orientation
  // (written in three ways that a cubic crystal cannot tell apart), the
  // right half turned by 20 degrees; two points are not indexed.
  const UnitCell cell = simple("fcc", {29}, 3.615);
  const double d = pi_v / 180.0;
  std::string ang = "# TEM_PIXperUM          1.000000\n# GRID: SqrGrid\n# XSTEP: 0.500000\n# YSTEP: 0.500000\n";
  std::string ctf = "Channel Text File\nPrj\ttest\nJobMode\tGrid\nXCells\t12\nYCells\t6\nXStep\t0.5\nYStep\t0.5\n"
                    "Phases\t1\n3.615;3.615;3.615\t90;90;90\tCopper\t11\t225\n"
                    "Phase\tX\tY\tBands\tError\tEuler1\tEuler2\tEuler3\tMAD\tBC\tBS\n";
  char line[200];
  for (int iy = 0; iy < 6; ++iy) {
    for (int ix = 0; ix < 12; ++ix) {
      // Left: phi1 = 10, 100 or 190 degrees (90 degrees apart about the
      // normal: the same cubic crystal). Right: 30 degrees.
      double phi1 = ix < 6 ? 10.0 + 90.0 * ((ix + iy) % 3) : 30.0;
      const bool missing = (ix == 2 && iy == 3) || (ix == 9 && iy == 1);
      std::snprintf(line, sizeof line, "%9.5f %9.5f %9.5f %12.5f %12.5f %.1f %6.3f %2d %6d %.3f\n",
                    missing ? 4.0 * pi_v : phi1 * d, missing ? 4.0 * pi_v : 0.0, missing ? 4.0 * pi_v : 0.0, 0.5 * ix,
                    0.5 * iy, 100.0, missing ? -1.0 : 0.9, 1, 0, 1.0);
      ang += line;
      std::snprintf(line, sizeof line, "%d\t%.4f\t%.4f\t8\t0\t%.4f\t0.0000\t0.0000\t0.5\t100\t200\n", missing ? 0 : 1,
                    0.5 * ix, 0.5 * iy, missing ? 0.0 : phi1);
      ctf += line;
    }
  }
  Polycrystal p;
  p.source = Polycrystal::Source::ebsd;
  for (const bool is_ctf : {false, true}) {
    const auto map = parse_ebsd(is_ctf ? ctf : ang, is_ctf, p, cell);
    EXPECT_EQ(map.grain_count(), 2u) << is_ctf;
    EXPECT_TRUE(map.is_columnar());
    EXPECT_EQ(map.seeds().size(), 72u);
    EXPECT_EQ(map.ungrained_fraction(), 0.0); // the two points took their neighbours' grain
    EXPECT_DOUBLE_EQ(map.period_y_A(), 6.0e4); // 12 points of 0.5 um
    EXPECT_DOUBLE_EQ(map.period_z_A(), 3.0e4);
    // Map x is the layer's y; the grains are columns.
    const int left = map.at({10.0, 1.0e4, 1.0e4}).grain;
    const int right = map.at({10.0, 4.5e4, 1.0e4}).grain;
    EXPECT_NE(left, right);
    EXPECT_EQ(map.at({9.0e5, 1.0e4, 1.0e4}).grain, left);
    // (The point at x = 1 um, y = 1.5 um, which is the one asked for above,
    // was not indexed; nor was the one at x = 4.5 um, y = 0.5 um.)
    EXPECT_EQ(map.at({10.0, 1.0e4, 1.5e4}).grain, left);
    EXPECT_EQ(map.at({10.0, 4.5e4, 2.0e4}).grain, right);
    // The boundary between the columns of points at 2.5 and 3.0 um.
    EXPECT_NEAR(map.at({10.0, 2.6e4, 1.0e4}, 1.0e4).boundary_A, 0.15e4, 1.0);
    EXPECT_NEAR(misorientation_rad(map.grain(left).rotation, map.grain(right).rotation, lattice_rotations(cell.lattice)),
                20.0 * d, 1e-4);
  }
  // A triclinic cell has no rotations of its own: the three ways of
  // writing the left half are three orientations, in many grains.
  UnitCell low = cell;
  low.lattice = Lattice::from_parameters(3.0 * angstrom, 4.0 * angstrom, 5.0 * angstrom, 80.0 * degree, 85.0 * degree,
                                         95.0 * degree);
  EXPECT_GT(parse_ebsd(ang, false, p, low).grain_count(), 10u);
  // A larger threshold joins the halves.
  Polycrystal coarse = p;
  coarse.misorientation = 25.0 * degree;
  EXPECT_EQ(parse_ebsd(ang, false, coarse, cell).grain_count(), 1u);
  // Another phase than the map's: nothing to use.
  Polycrystal other = p;
  other.phase = 2;
  EXPECT_THROW(parse_ebsd(ctf, true, other, cell), std::invalid_argument);
  EXPECT_THROW(parse_ebsd("# header only\n", false, p, cell), std::invalid_argument);
  EXPECT_THROW(parse_ebsd("Channel Text File\nXCells\t2\n", true, p, cell), std::invalid_argument);
}
