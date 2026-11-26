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

/*!
 * \file trim_test.cpp
 *
 * \brief Tests of the SRIM-mode (TRIM-compatible) engine: random numbers,
 * the Kinchin-Pease model, TRIM.IN parsing/writing, determinism, the
 * relation between the two damage modes, energy bookkeeping and output
 * files, plus regression checks against SRIM-2013 TRIM results.
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>

#include "run_config.h"
#include "data_files.h"
#include "amcset_common.h"
#include "crystal/crystal_database.h"
#include "trim/random.h"
#include "trim/single_precision.h"
#include "trim/trim.h"
#include "trim/trim_dat.h"
#include "trim/trim_in.h"
#include "trim/trim_tables.h"

#include <chrono>
#include <thread>

using namespace amcset;
using namespace amcset::common;
using namespace amcset::trim;

namespace {

double rel(double a, double b) { return std::abs(a - b) / std::abs(b); }

//! Fe -> Fe with TRIM's defaults (Ed 40 eV, Eb 0, Es 4.34 eV), built-in
//! stopping so that no data files are needed.
TrimConfig fe_on_fe(double energy_keV, DamageCalculation damage, std::size_t ions,
                    double width_A = 0.0) {
  TrimConfig c;
  c.ion = IonBeam::defaults(26, energy_keV * kilo_electron_volt);
  TargetLayer layer;
  layer.name = "Fe";
  layer.density = 7.8658 * gram_per_cubic_centimeter;
  layer.width = (width_A > 0.0 ? width_A : 20.0 * energy_keV + 500.0) * angstrom;
  auto fe = TargetElement::defaults(26);
  fe.displacement_energy = 40.0 * electron_volt;
  fe.lattice_binding_energy = 0.0 * electron_volt;
  layer.elements.push_back(fe);
  c.layers.push_back(layer);
  c.damage = damage;
  c.ion_count = ions;
  c.threads = 2;
  return c;
}

std::filesystem::path temp_dir(const std::string &name) {
  auto d = std::filesystem::temp_directory_path() / ("amcset_test_" + name);
  std::filesystem::remove_all(d);
  std::filesystem::create_directories(d);
  return d;
}

} // namespace

// ---------------------------------------------------------------------------
// Random numbers
// ---------------------------------------------------------------------------

TEST(TrimRandom, LecuyerGeneratorIsUniformAndReproducible) {
  LecuyerRandom a(1, 2), b(1, 2), c(1, 3);
  double mean = 0.0, var = 0.0;
  const int n = 200000;
  bool differs = false;
  for (int i = 0; i < n; ++i) {
    const double x = a();
    ASSERT_GT(x, 0.0);
    ASSERT_LT(x, 1.0);
    EXPECT_EQ(x, b());
    if (x != c()) differs = true;
    mean += x;
    var += (x - 0.5) * (x - 0.5);
  }
  EXPECT_TRUE(differs);
  EXPECT_NEAR(mean / n, 0.5, 0.005);
  EXPECT_NEAR(var / n, 1.0 / 12.0, 0.002);
}

TEST(TrimRandom, LecuyerRecurrence) {
  // One step of each component by hand.
  LecuyerRandom r(12345, 67890);
  const std::int64_t s1 = (40014LL * 12346) % 2147483563LL;
  const std::int64_t s2 = (40692LL * 67891) % 2147483399LL;
  std::int64_t z = s1 - s2;
  if (z < 1) z += 2147483562LL;
  EXPECT_DOUBLE_EQ(r(), static_cast<double>(z) / 2147483563.0);
}

TEST(TrimRandom, ModularReductionMatchesTheRemainder) {
  using R = LecuyerRandom;
  std::mt19937_64 rng(42);
  for (const auto [a, m] : {std::pair{R::a1, R::m1}, std::pair{R::a2, R::m2}}) {
    for (const std::int64_t s : {std::int64_t{0}, std::int64_t{1}, std::int64_t{2}, m - 2, m - 1, m, m + 1,
                                 (std::int64_t{1} << 31) - 1}) {
      ASSERT_EQ(R::mul_mod(a, s, m), (a * s) % m) << s;
    }
    for (int i = 0; i < 2000000; ++i) {
      const auto s = static_cast<std::int64_t>(rng() >> 33); // 0 .. 2^31 - 1
      ASSERT_EQ(R::mul_mod(a, s, m), (a * s) % m) << s;
    }
  }
  // A long stream is the textbook recurrence.
  LecuyerRandom g(716381, 12345);
  std::int64_t s1 = g.state1();
  std::int64_t s2 = g.state2();
  for (int i = 0; i < 1000000; ++i) {
    s1 = (R::a1 * s1) % R::m1;
    s2 = (R::a2 * s2) % R::m2;
    std::int64_t z = s1 - s2;
    if (z < 1) z += R::m1 - 1;
    ASSERT_EQ(g(), static_cast<double>(z) * (1.0 / static_cast<double>(R::m1)));
  }
}

TEST(TrimRandom, GaussianMoments) {
  LecuyerRandom r(7, 11);
  double m1 = 0.0, m2 = 0.0;
  const int n = 100000;
  for (int i = 0; i < n; ++i) {
    const double g = r.gaussian();
    m1 += g;
    m2 += g * g;
  }
  EXPECT_NEAR(m1 / n, 0.0, 0.02);
  EXPECT_NEAR(m2 / n, 1.0, 0.02);
}

TEST(TrimRandom, StreamsAreIndependentPerIon) {
  std::set<double> first;
  for (std::uint64_t ion = 0; ion < 1000; ++ion) {
    for (unsigned stream = 0; stream < 2; ++stream) {
      auto r = make_stream(716381, ion, stream);
      first.insert(r());
    }
  }
  EXPECT_EQ(first.size(), 2000u);
  auto x = make_stream(1, 5, 0);
  auto y = make_stream(1, 5, 0);
  EXPECT_EQ(x(), y());
}

// ---------------------------------------------------------------------------
// Kinchin-Pease (quick calculation)
// ---------------------------------------------------------------------------

TEST(TrimKinchinPease, MatchesTrimCollisionTable) {
  // Recoil energies and vacancy numbers from TRIM's COLLISON.txt, 100 keV Fe
  // into Fe (Ed = 40 eV), quick calculation.
  const mass_quantity m = 55.847 * atomic_mass_unit;
  const energy_quantity ed = 40.0 * electron_volt;
  const std::vector<std::pair<double, double>> trim{
      {472.24, 3.951}, {234.82, 2.001}, {3672.8, 28.704}, {23279.0, 165.947},
      {58.591, 1.0},   {13393.0, 98.620}, {1920.5, 15.377}, {317.43, 2.684}};
  for (const auto &[t, nu] : trim) {
    EXPECT_LT(rel(kinchin_pease_vacancies(26, m, ed, t * electron_volt), nu), 2e-3) << t;
  }
}

TEST(TrimKinchinPease, ThresholdIsOnTheDamageEnergy) {
  // TRIM's smallest tabulated recoil in that run is 45.269 eV, whose damage
  // energy is Ed = 40 eV: recoils with 40 < T < 45.27 eV are not displaced.
  const mass_quantity m = 55.847 * atomic_mass_unit;
  const energy_quantity ed = 40.0 * electron_volt;
  EXPECT_NEAR(in_units(robinson_damage_energy(26, m, 45.27 * electron_volt), electron_volt), 40.0, 0.01);
  EXPECT_DOUBLE_EQ(kinchin_pease_vacancies(26, m, ed, 45.3 * electron_volt), 1.0);
  EXPECT_DOUBLE_EQ(kinchin_pease_vacancies(26, m, ed, 45.2 * electron_volt), 0.0);
  EXPECT_DOUBLE_EQ(kinchin_pease_vacancies(26, m, ed, 41.0 * electron_volt), 0.0);
  // Damage energy is below T and tends to T at low energy.
  const auto d = robinson_damage_energy(26, m, 1.0 * mega_electron_volt);
  // 1 MeV Fe in Fe: eps_d = 5.05, g = 10.9, k_d = 0.157 -> E_d = 0.37 MeV.
  EXPECT_NEAR(in_units(d, mega_electron_volt), 0.370, 0.005);
}

// ---------------------------------------------------------------------------
// Configuration helpers
// ---------------------------------------------------------------------------

TEST(TrimConfig, DamageNamesAndValidation) {
  EXPECT_EQ(parse_damage_calculation("kp"), DamageCalculation::kinchin_pease);
  EXPECT_EQ(parse_damage_calculation("full"), DamageCalculation::full_cascade);
  EXPECT_EQ(parse_damage_calculation(to_string(DamageCalculation::monolayer)),
            DamageCalculation::monolayer);
  EXPECT_THROW(parse_damage_calculation("nope"), std::invalid_argument);

  auto c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 10);
  EXPECT_NO_THROW(c.validate());
  EXPECT_NEAR(in_units(c.target_thickness(), angstrom), 700.0, 1e-9);
  auto bad = c;
  bad.layers.clear();
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.ion.energy = 0.0 * electron_volt;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.layers[0].elements[0].stoichiometry = 0.0;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.stopping_model = physics::ElectronicStoppingModel::srim_tables;
  EXPECT_NO_THROW(bad.validate()); // no table: the bundled SRIM-2013 tables
  EXPECT_FALSE(c.describe().empty());
}

TEST(DataFiles, ElectronTablesAreFoundOutsideTheSourceTree) {
  // The native engine's electron tables: installed as share/amcset/electrons
  // (found from any working directory), and in common/assets of the source
  // tree for a build.
  for (const char *name : {"electron_screening_parameters.csv", "mott_scattering_parameters.csv"}) {
    const auto p = data::find_electron_table(name);
    ASSERT_FALSE(p.empty()) << name;
    EXPECT_EQ(p.filename().string(), name);
  }
  EXPECT_TRUE(data::find_electron_table("no_such_table.csv").empty());
}

TEST(TrimConfig, BundledSrimTablesAreTheDefault) {
  // TRIM.IN and scenario files default to the SRIM-2013 tables of data/srim2013.
  const auto path = data::bundled_srim_stopping_table();
  if (path.empty()) GTEST_SKIP() << "data/srim2013 not found";
  auto c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 20);
  c.stopping_model = physics::ElectronicStoppingModel::srim_tables;
  TrimSimulation s(c);
  EXPECT_EQ(s.config().stopping_table, path);
  std::ostringstream out;
  write_trim_in(out, c, extras_for(c));
  std::istringstream in(out.str());
  EXPECT_EQ(parse_trim_in(in).stopping_model, physics::ElectronicStoppingModel::srim_tables);
  EXPECT_EQ(config::RunConfig{}.trim.stopping_model, physics::ElectronicStoppingModel::srim_tables);
}

TEST(TrimBugs, EachFixChangesWhatItShould) {
  TrimBugs b;
  EXPECT_TRUE(b.all_reproduced());
  b.set_all(false);
  for (const auto &m : TrimBugs::members()) EXPECT_FALSE(b.*(m.second)) << m.first;
  EXPECT_EQ(TrimBugs::members().size(), 11u);

  // Replacements: with the fix the interstitials no longer exceed the
  // vacancies by the replacements.
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 200);
  TrimSimulation trim_like(c);
  trim_like.run();
  const auto &a = trim_like.tallies();
  ASSERT_GT(a.total_replacements, 0.0);
  c.bugs.replacement_interstitial = false;
  TrimSimulation fixed(c);
  fixed.run();
  const auto &f = fixed.tallies();
  EXPECT_DOUBLE_EQ(f.total_vacancies, a.total_vacancies); // same random numbers
  EXPECT_DOUBLE_EQ(f.total_interstitials, a.total_interstitials - a.total_replacements);

  // Phonons counted twice: the fix removes the second booking only.
  c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 200);
  TrimSimulation twice(c);
  twice.run();
  c.bugs.ion_phonons_counted_twice = false;
  TrimSimulation once(c);
  once.run();
  EXPECT_LT(once.tallies().total_phonons_recoils, twice.tallies().total_phonons_recoils);
  EXPECT_DOUBLE_EQ(once.tallies().total_phonons_ion, twice.tallies().total_phonons_ion);

  // Energy creation in recoil collisions: with the fix the energy balance of
  // a slow cascade closes better.
  c = fe_on_fe(1.0, DamageCalculation::full_cascade, 300);
  c.bugs.ion_phonons_counted_twice = false;
  auto balance = [](const TrimTallies &t, double e0) {
    const double out = t.total_ionization_ion + t.total_ionization_recoils + t.total_phonons_ion +
                       t.total_phonons_recoils + t.total_binding_energy_loss;
    return out / (e0 * static_cast<double>(t.ions));
  };
  TrimSimulation created(c);
  created.run();
  c.bugs.recoil_energy_created = false;
  TrimSimulation conserved(c);
  conserved.run();
  EXPECT_LT(balance(conserved.tallies(), 1000.0), balance(created.tallies(), 1000.0));

  // The configuration keeps the flags.
  config::RunConfig rc;
  rc.trim = fe_on_fe(10.0, DamageCalculation::full_cascade, 10);
  rc.trim.bugs.struck_atom_direction = false;
  const auto back = config::parse_config(config::to_json(rc));
  EXPECT_FALSE(back.trim.bugs.struck_atom_direction);
  EXPECT_TRUE(back.trim.bugs.stale_primary_energy);
  rc.trim.bugs.struck_atom_direction = true;
  std::string text = config::to_json(rc);
  const auto key = text.find("\"trim_bugs\"");
  ASSERT_NE(key, std::string::npos);
  const auto value = text.find("\"trim\"", key + 11);
  ASSERT_NE(value, std::string::npos);
  text.replace(value, 6, "\"fixed\"");
  const auto all = config::parse_config(text);
  for (const auto &m : TrimBugs::members()) EXPECT_FALSE(all.trim.bugs.*(m.second)) << m.first;
}

TEST(TrimConfig, ElementDefaults) {
  const auto fe = TargetElement::defaults(26);
  EXPECT_NEAR(in_units(fe.displacement_energy, electron_volt), 25.0, 1e-9);
  EXPECT_NEAR(in_units(fe.lattice_binding_energy, electron_volt), 3.0, 1e-9);
  EXPECT_NEAR(in_units(fe.surface_binding_energy, electron_volt), 4.34, 1e-9);
  const auto ion = IonBeam::defaults(2, 1.0 * mega_electron_volt);
  EXPECT_NEAR(in_units(ion.mass, atomic_mass_unit), 4.003, 1e-3);
}

// ---------------------------------------------------------------------------
// TRIM.IN
// ---------------------------------------------------------------------------

namespace {

// A TRIM.IN as written by SRIM-2013 (two elements, two layers).
const char *trim_in_text =
    "==> SRIM-2013.00 This file controls TRIM Calculations.\r\n"
    "Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.\r\n"
    "     2   4.003   2000   7   500   1   10000\r\n"
    "Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders\r\n"
    "                      2                                   123       0\r\n"
    "Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file\r\n"
    "                          1       0       1       1       0       0\r\n"
    "Target material : Number of Elements & Layers\r\n"
    "\"FeCr on Ni\"       3               2\r\n"
    "PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]\r\n"
    "       5                         0            0\r\n"
    "Target Elements:    Z   Mass(amu)\r\n"
    "Atom 1 = Fe =      26  55.847\r\n"
    "Atom 2 = Cr =      24  51.996\r\n"
    "Atom 3 = Ni =      28  58.69\r\n"
    "Layer   Layer Name /               Width Density    Fe(26)  Cr(24)  Ni(28)\r\n"
    "Numb.   Description                (Ang) (g/cm3)    Stoich  Stoich  Stoich\r\n"
    " 1      \"FeCr\"           5000  7.8  .91  .09  0\r\n"
    " 2      \"Ni\"           20000  8.9  0  0  1\r\n"
    "0  Target layer phases (0=Solid, 1=Gas)\r\n"
    "0 0 \r\n"
    "Target Compound Corrections (Bragg)\r\n"
    " 1  .95  \r\n"
    "Individual target atom displacement energies (eV)\r\n"
    "   40 35 30\r\n"
    "Individual target atom lattice binding energies (eV)\r\n"
    "   3 3 3\r\n"
    "Individual target atom surface binding energies (eV)\r\n"
    "   4.34 4.12 4.46\r\n"
    "Stopping Power Version (1=2011, 0=2011)\r\n"
    " 0 \r\n";

} // namespace

TEST(TrimIn, ParsesSrim2013Input) {
  std::istringstream in(trim_in_text);
  TrimInExtras extras;
  const auto c = parse_trim_in(in, &extras);
  EXPECT_EQ(c.ion.z, 2u);
  EXPECT_NEAR(in_units(c.ion.mass, atomic_mass_unit), 4.003, 1e-9);
  EXPECT_NEAR(in_units(c.ion.energy, kilo_electron_volt), 2000.0, 1e-9);
  EXPECT_NEAR(in_units(c.ion.angle, amcset::common::degree), 7.0, 1e-9);
  EXPECT_EQ(c.ion_count, 500u);
  EXPECT_EQ(c.damage, DamageCalculation::full_cascade);
  EXPECT_EQ(c.seed, 123u);
  EXPECT_EQ(extras.cascade_code, 2);
  EXPECT_FALSE(extras.backscatter_file);
  EXPECT_EQ(c.title, "FeCr on Ni");
  ASSERT_EQ(c.layers.size(), 2u);
  // Elements with zero stoichiometry in a layer are dropped from it.
  ASSERT_EQ(c.layers[0].elements.size(), 2u);
  ASSERT_EQ(c.layers[1].elements.size(), 1u);
  EXPECT_EQ(c.layers[0].name, "FeCr");
  EXPECT_NEAR(in_units(c.layers[0].width, angstrom), 5000.0, 1e-9);
  EXPECT_NEAR(in_units(c.layers[1].density, gram_per_cubic_centimeter), 8.9, 1e-12);
  EXPECT_NEAR(c.layers[1].compound_correction, 0.95, 1e-12);
  EXPECT_NEAR(c.layers[0].elements[1].stoichiometry, 0.09, 1e-12);
  EXPECT_NEAR(in_units(c.layers[0].elements[1].displacement_energy, electron_volt), 35.0, 1e-12);
  EXPECT_NEAR(in_units(c.layers[1].elements[0].surface_binding_energy, electron_volt), 4.46, 1e-12);
  EXPECT_EQ(c.layers[1].elements[0].z, 28u);
}

TEST(TrimIn, RoundTrip) {
  std::istringstream in(trim_in_text);
  TrimInExtras extras;
  const auto c = parse_trim_in(in, &extras);
  std::ostringstream out;
  write_trim_in(out, c, extras);
  std::istringstream again(out.str());
  TrimInExtras extras2;
  const auto d = parse_trim_in(again, &extras2);
  EXPECT_EQ(d.ion.z, c.ion.z);
  EXPECT_EQ(d.ion_count, c.ion_count);
  EXPECT_EQ(d.damage, c.damage);
  EXPECT_EQ(d.seed, c.seed);
  ASSERT_EQ(d.layers.size(), c.layers.size());
  for (std::size_t i = 0; i < c.layers.size(); ++i) {
    EXPECT_EQ(d.layers[i].name, c.layers[i].name);
    EXPECT_NEAR(in_units(d.layers[i].width, angstrom), in_units(c.layers[i].width, angstrom), 1e-6);
    ASSERT_EQ(d.layers[i].elements.size(), c.layers[i].elements.size());
    for (std::size_t j = 0; j < c.layers[i].elements.size(); ++j) {
      EXPECT_EQ(d.layers[i].elements[j].z, c.layers[i].elements[j].z);
      EXPECT_NEAR(d.layers[i].elements[j].stoichiometry, c.layers[i].elements[j].stoichiometry, 1e-9);
      EXPECT_NEAR(in_units(d.layers[i].elements[j].displacement_energy, electron_volt),
                  in_units(c.layers[i].elements[j].displacement_energy, electron_volt), 1e-9);
    }
  }
  EXPECT_EQ(extras2.backscatter_file, extras.backscatter_file);
  // TRIM needs DOS line endings.
  EXPECT_NE(out.str().find("\r\n"), std::string::npos);
}

TEST(TrimIn, KeepsTrimsAtomOrder) {
  // Atom 1 (Fe) appears only in the second layer; TRIM's order is kept.
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 10);
  TargetLayer cr;
  cr.name = "Cr";
  cr.width = 300.0 * angstrom;
  cr.density = 7.19 * gram_per_cubic_centimeter;
  cr.elements.push_back(TargetElement::defaults(24));
  c.layers.insert(c.layers.begin(), cr);
  c.layers[1].elements[0].trim_atom = 1;
  c.layers[0].elements[0].trim_atom = 2;
  std::ostringstream out;
  write_trim_in(out, c);
  EXPECT_LT(out.str().find("Atom 1 = Fe"), out.str().find("Atom 2 = Cr"));
  std::istringstream in(out.str());
  const auto d = parse_trim_in(in);
  EXPECT_EQ(d.layers[0].elements[0].trim_atom, 2u);
  EXPECT_EQ(d.layers[1].elements[0].trim_atom, 1u);
  TrimSimulation s(d);
  ASSERT_EQ(s.element_names().size(), 2u);
  EXPECT_EQ(s.element_names()[0], "Fe");
  EXPECT_EQ(s.element_names()[1], "Cr");
}

TEST(TrimIn, ReadsTrimDatModesAndRejectsUnknownCodes) {
  std::string text = trim_in_text;
  const auto pos = text.find("                      2      ");
  ASSERT_NE(pos, std::string::npos);
  for (char code : {'4', '5', '6', '7'}) {
    text[pos + 22] = code; // TRIM.DAT modes
    std::istringstream in(text);
    const auto c = parse_trim_in(in);
    EXPECT_EQ(c.trim_code(), code - '0');
    EXPECT_NE(c.source, TrimConfig::Source::beam);
    EXPECT_EQ(c.ion_count, 0u); // every TRIM.DAT line
  }
  text[pos + 22] = '8';
  std::istringstream in(text);
  EXPECT_THROW(parse_trim_in(in), std::exception);
  // read_trim_in needs TRIM.DAT next to TRIM.IN.
  text[pos + 22] = '5';
  const auto dir = temp_dir("trim_in_dat");
  {
    std::ofstream o(dir / "TRIM.IN", std::ios::binary);
    o << text;
  }
  EXPECT_THROW(read_trim_in(dir / "TRIM.IN"), std::runtime_error);
  {
    std::ofstream o(dir / "trim.dat", std::ios::binary);
    o << "1\n2\n3\n4\n5\n6\n7\nT\n9\n10\nA 5 1000 0 0 0 1 0 0\n";
  }
  const auto c = read_trim_in(dir / "TRIM.IN");
  ASSERT_EQ(c.particles.size(), 1u);
  EXPECT_EQ(c.particles_title, "T");
}

// ---------------------------------------------------------------------------
// Engine behaviour
// ---------------------------------------------------------------------------

TEST(TrimEngine, ResultsDoNotDependOnThreadCount) {
  auto c = fe_on_fe(20.0, DamageCalculation::full_cascade, 40);
  c.threads = 1;
  TrimSimulation one(c);
  one.run();
  c.threads = 4;
  TrimSimulation four(c);
  four.run();
  const auto &a = one.tallies();
  const auto &b = four.tallies();
  EXPECT_EQ(a.ions, 40u);
  EXPECT_DOUBLE_EQ(a.range.s1, b.range.s1);
  EXPECT_DOUBLE_EQ(a.range.s2, b.range.s2);
  EXPECT_DOUBLE_EQ(a.total_vacancies, b.total_vacancies);
  EXPECT_DOUBLE_EQ(a.total_replacements, b.total_replacements);
  EXPECT_DOUBLE_EQ(a.total_ionization_recoils, b.total_ionization_recoils);
  EXPECT_EQ(a.backscattered, b.backscattered);
}

TEST(TrimEngine, SeedChangesTheResult) {
  auto c = fe_on_fe(20.0, DamageCalculation::kinchin_pease, 40);
  TrimSimulation a(c);
  a.run();
  c.seed += 1;
  TrimSimulation b(c);
  b.run();
  EXPECT_NE(a.tallies().range.s1, b.tallies().range.s1);
}

// "seed": "random": every engine draws its own seed, config() and the
// outputs give it, and that seed reruns the same calculation.
TEST(TrimEngine, RandomSeedIsDrawnRecordedAndRepeatable) {
  auto c = fe_on_fe(20.0, DamageCalculation::kinchin_pease, 20);
  c.random_seed = true;
  c.output_directory = temp_dir("random_seed").string();
  c.outputs.formats = {output::Format::trim};
  TrimSimulation a(c);
  TrimSimulation b(c);
  EXPECT_NE(a.config().seed, b.config().seed); // (31 random bits: equal once in 2^31)
  EXPECT_LT(a.config().seed, 1u << 31);
  EXPECT_TRUE(a.config().random_seed);
  a.run();
  a.write_outputs(c.output_directory);
  std::ifstream in(std::filesystem::path(c.output_directory) / "summary.json");
  const std::string summary((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  EXPECT_NE(summary.find("\"seed\": " + std::to_string(a.config().seed) + ","), std::string::npos);
  EXPECT_NE(summary.find("\"seed_random\": true"), std::string::npos);
  auto again = c;
  again.random_seed = false;
  again.seed = a.config().seed;
  again.output_directory.clear();
  TrimSimulation r(again);
  EXPECT_EQ(r.config().seed, a.config().seed);
  r.run();
  EXPECT_EQ(r.tallies().range.s1, a.tallies().range.s1);
  EXPECT_EQ(r.tallies().total_vacancies, a.tallies().total_vacancies);
}

TEST(TrimEngine, DamageModesShareIonTrajectories) {
  // Refined transport: recoils use their own random stream, so the ion's path
  // (and thus the range distribution and the primary knock-ons) is the same
  // in both modes.
  auto kp = fe_on_fe(50.0, DamageCalculation::kinchin_pease, 60);
  kp.transport = TrimConfig::Transport::refined;
  auto fc = kp;
  fc.damage = DamageCalculation::full_cascade;
  TrimSimulation a(kp);
  a.run();
  TrimSimulation b(fc);
  b.run();
  // Only the last step of an ion that ends in a replacement collision (full
  // cascades) differs, and E2RECOIL uses each mode's displacement test.
  EXPECT_LT(rel(a.tallies().range.s1, b.tallies().range.s1), 1e-3);
  EXPECT_LT(rel(a.tallies().lateral_abs_y.s1, b.tallies().lateral_abs_y.s1), 2e-3);
  EXPECT_LT(rel(a.tallies().total_energy_to_recoils, b.tallies().total_energy_to_recoils), 0.01);
  // Full cascades give roughly twice the Kinchin-Pease vacancies for Fe
  // (TRIM: 1313 vs 658 per ion at 100 keV).
  const double ratio = b.tallies().total_vacancies / a.tallies().total_vacancies;
  EXPECT_GT(ratio, 1.6);
  EXPECT_LT(ratio, 2.6);
  EXPECT_GT(b.tallies().total_replacements, 0.0);
  EXPECT_DOUBLE_EQ(a.tallies().total_replacements, 0.0);
}

// (deposited + carried out) / incident energy.
double energy_balance(const TrimSimulation &s, double energy_eV) {
  const auto &t = s.tallies();
  const double e0 = energy_eV * static_cast<double>(t.ions);
  double sputtered = 0.0;
  for (double e : t.sputtered_energy) sputtered += e;
  const double deposited = t.total_ionization_ion + t.total_ionization_recoils + t.total_phonons_ion +
                           t.total_phonons_recoils + t.total_binding_energy_loss;
  const double out = t.backscattered_energy + t.transmitted_energy + sputtered;
  return (deposited + out) / e0;
}

TEST(TrimEngine, EnergyIsConserved) {
  for (auto damage : {DamageCalculation::kinchin_pease, DamageCalculation::full_cascade}) {
    auto c = fe_on_fe(30.0, damage, 50);
    c.transport = TrimConfig::Transport::refined;
    TrimSimulation s(c);
    s.run();
    EXPECT_NEAR(energy_balance(s, 30000.0), 1.0, 0.01) << to_string(damage);
  }
}

TEST(TrimEngine, Trim2013KeepsTrimsEnergyBookkeeping) {
  // TRIM.exe books some energy twice (transfers below 1 eV, the energy of a
  // primary slot that was not displaced, T + DEE above the mover's energy),
  // so its IONIZ/PHONON tables add up to slightly more than the incident
  // energy (Fe -> Fe: 100.8 % at 100 keV, 102.8 % at 10 keV, quick mode).
  for (auto damage : {DamageCalculation::kinchin_pease, DamageCalculation::full_cascade}) {
    auto c = fe_on_fe(30.0, damage, 50);
    ASSERT_EQ(c.transport, TrimConfig::Transport::trim2013);
    TrimSimulation s(c);
    s.run();
    const double balance = energy_balance(s, 30000.0);
    EXPECT_GT(balance, 0.995) << to_string(damage);
    EXPECT_LT(balance, 1.06) << to_string(damage);
  }
}

TEST(TrimEngine, Trim2013FlightRulesDependOnTheDamageMode) {
  // TRIM.exe randomises every ion flight in full-cascade mode but not in
  // quick mode, so the ion paths differ while the ranges agree on average
  // (TRIM: 346.5 A quick, 350.7 A full cascades for Fe -> Fe at 100 keV).
  auto kp = fe_on_fe(50.0, DamageCalculation::kinchin_pease, 400);
  auto fc = kp;
  fc.damage = DamageCalculation::full_cascade;
  TrimSimulation a(kp);
  a.run();
  TrimSimulation b(fc);
  b.run();
  EXPECT_NE(a.tallies().range.s1, b.tallies().range.s1);
  EXPECT_LT(rel(a.tallies().range.mean(), b.tallies().range.mean()), 0.12);
}

TEST(TrimEngine, Trim2013SputteringFollowsTrim) {
  // Fe -> Fe at 1 keV, full cascades: TRIM gives 1.98 sputtered atoms per ion
  // (AMCSET with SRIM's stopping tables: 2.01; the built-in stopping used
  // here gives about 1.5). Its SPUTTER.txt lists the energy inside the
  // surface (the binding energy is not subtracted) and the cosines before the
  // surface barrier.
  auto c = fe_on_fe(1.0, DamageCalculation::full_cascade, 400, 200.0);
  c.layers[0].elements[0].surface_binding_energy = 4.34 * electron_volt;
  TrimSimulation s(c);
  s.run();
  const auto &t = s.tallies();
  const double yield = t.sputtered[0] / static_cast<double>(t.ions);
  EXPECT_GT(yield, 1.2);
  EXPECT_LT(yield, 2.5);
  std::size_t listed = 0;
  for (const auto &r : t.exits) {
    if (r.kind != 'S') continue;
    ++listed;
    EXPECT_GE(r.energy_eV * r.cx * r.cx, 4.34 * (1.0 - 1e-6));
    EXPECT_LT(r.cx, 0.0);
    EXPECT_LE(std::abs(r.cx), 0.9999999);
  }
  EXPECT_GT(listed, 0u);
}

TEST(TrimEngine, ThinTargetTransmitsLightIons) {
  TrimConfig c;
  c.ion = IonBeam::defaults(1, 1.0 * mega_electron_volt);
  TargetLayer layer;
  layer.name = "Fe";
  layer.density = 7.8658 * gram_per_cubic_centimeter;
  layer.width = 1000.0 * angstrom;
  layer.elements.push_back(TargetElement::defaults(26));
  c.layers.push_back(layer);
  c.ion_count = 100;
  TrimSimulation s(c);
  s.run();
  EXPECT_EQ(s.tallies().transmitted, 100u);
  // Energy lost in 0.1 um of Fe at 1 MeV: ~ S x = 10.4 eV/A * 1000 A.
  const double lost = 1.0e6 - s.tallies().transmitted_energy / 100.0;
  EXPECT_GT(lost, 8000.0);
  EXPECT_LT(lost, 14000.0);
}

TEST(TrimEngine, PerElementDisplacementEnergyMatters) {
  auto low = fe_on_fe(20.0, DamageCalculation::kinchin_pease, 40);
  auto high = low;
  high.layers[0].elements[0].displacement_energy = 80.0 * electron_volt;
  TrimSimulation a(low);
  a.run();
  TrimSimulation b(high);
  b.run();
  // NRT: vacancies ~ 1/Ed for energetic recoils.
  const double ratio = a.tallies().total_vacancies / b.tallies().total_vacancies;
  EXPECT_GT(ratio, 1.6);
  EXPECT_LT(ratio, 2.4);
}

TEST(TrimEngine, ZeroEnergiesOfLaterAtomsAreReplacedAsInTrim) {
  // TRIM.exe sets Ed <= 0 to 20 eV, Es <= 0 to 4 eV and Eb <= 0 to 2 eV for
  // its target atoms 2..n, not for atom 1 (TDATA.txt shows the result).
  auto c = fe_on_fe(20.0, DamageCalculation::kinchin_pease, 1);
  auto cr = TargetElement::defaults(24, 0.5);
  cr.displacement_energy = 0.0 * electron_volt;
  cr.lattice_binding_energy = 0.0 * electron_volt;
  cr.surface_binding_energy = 0.0 * electron_volt;
  c.layers[0].elements.push_back(cr);
  c.layers[0].elements[0].surface_binding_energy = 0.0 * electron_volt;
  TrimSimulation s(c);
  const auto &fe_used = s.effective_layers()[0].elements[0];
  const auto &cr_used = s.effective_layers()[0].elements[1];
  EXPECT_DOUBLE_EQ(in_units(fe_used.lattice_binding_energy, electron_volt), 0.0);
  EXPECT_DOUBLE_EQ(in_units(fe_used.surface_binding_energy, electron_volt), 0.0);
  EXPECT_DOUBLE_EQ(in_units(cr_used.displacement_energy, electron_volt), 20.0);
  EXPECT_DOUBLE_EQ(in_units(cr_used.surface_binding_energy, electron_volt), 4.0);
  EXPECT_DOUBLE_EQ(in_units(cr_used.lattice_binding_energy, electron_volt), 2.0);
  ASSERT_EQ(s.warnings().size(), 1u);
  EXPECT_NE(s.warnings()[0].find("Target atom 2 (Cr)"), std::string::npos);
  // The configuration itself is unchanged.
  EXPECT_DOUBLE_EQ(in_units(s.config().layers[0].elements[1].lattice_binding_energy, electron_volt), 0.0);

  // TRIM's atom order (from a TRIM.IN) decides which atom is the first.
  auto d = c;
  d.layers[0].elements[0].trim_atom = 2;
  d.layers[0].elements[1].trim_atom = 1;
  TrimSimulation t(d);
  EXPECT_EQ(t.element_names().front(), "Cr");
  EXPECT_DOUBLE_EQ(in_units(t.effective_layers()[0].elements[1].lattice_binding_energy, electron_volt), 0.0);
  EXPECT_DOUBLE_EQ(in_units(t.effective_layers()[0].elements[0].lattice_binding_energy, electron_volt), 2.0);
}

TEST(TrimEngine, QuickCalculationBooksTheBindingEnergyOncePerDisplacedAtom) {
  // TRIM's E2RECOIL "from IONS" loses Eb once per displaced atom (not nu Eb);
  // the "absorbed" column keeps the whole transfer.
  auto c = fe_on_fe(50.0, DamageCalculation::kinchin_pease, 40);
  c.layers[0].elements[0].lattice_binding_energy = 3.0 * electron_volt;
  TrimSimulation s(c);
  s.run();
  const auto &t = s.tallies();
  double from_ions = 0.0, absorbed = 0.0, displaced = 0.0;
  for (std::size_t i = 0; i < t.energy_to_recoils.values.size(); ++i) {
    from_ions += t.energy_to_recoils.values[i];
    absorbed += t.energy_absorbed[0].values[i];
    displaced += t.vacancies_by_ion.values[i];
  }
  ASSERT_GT(displaced, 0.0);
  EXPECT_NEAR(absorbed - from_ions, 3.0 * displaced, 1e-6 * absorbed);
}

TEST(TrimTables, SinglePrecisionSumsMatchStepByStepAddition) {
  // The fast single-precision sum equals adding the steps one by one.
  std::mt19937_64 g(5);
  for (int t = 0; t < 12; ++t) {
    const std::uint64_t up = 1000 + g() % 400000;
    const std::uint64_t down = t % 3 == 0 ? 0 : g() % (up / 10);
    const float step = 1.0f / static_cast<float>(20 + g() % 300);
    const std::uint64_t total = up + down;
    const double ratio = static_cast<double>(down) / static_cast<double>(total);
    auto before = [&](std::uint64_t k) {
      return k >= total ? down : std::min<std::uint64_t>(down, static_cast<std::uint64_t>(std::floor(static_cast<double>(k) * ratio)));
    };
    float s = 0.0f;
    for (std::uint64_t k = 0; k < total; ++k) {
      const bool dec = before(k + 1) > before(k);
      s = static_cast<float>(static_cast<double>(s) + static_cast<double>(dec ? -step : step));
    }
    EXPECT_EQ(single_precision_sum(up, down, step), static_cast<double>(s)) << up << " " << down << " " << step;
  }
  // TRIM's VACANCY.txt for 2000 Fe ions at 2 MeV: a bin of 480000 vacancies
  // (180 A bins) comes out 0.65 % high.
  EXPECT_NEAR(trim_table_count(480000.0, 0.0, 180.0) / 480000.0, 1.0065, 0.0002);
  // Small counts stay (nearly) exact.
  EXPECT_NEAR(trim_table_count(37.0, 0.0, 180.0), 37.0, 1e-4);
}

TEST(TrimEngine, CollisionDetailsListCascadeAtoms) {
  // TRIM's Collisions = 2: every displacement and replacement of the
  // cascades, and a damage summary per ion.
  auto c = fe_on_fe(30.0, DamageCalculation::full_cascade, 4);
  c.collision_details_ions = 3;
  c.collision_recoil_details = true;
  TrimSimulation s(c);
  s.run();
  const auto &t = s.tallies();
  ASSERT_FALSE(t.cascades.empty());
  ASSERT_EQ(t.ion_summaries.size(), 3u);
  double disp = 0.0, repl = 0.0;
  for (const auto &r : t.cascades) {
    ASSERT_FALSE(r.atoms.empty());
    EXPECT_TRUE(r.atoms.front().vacancy);
    EXPECT_DOUBLE_EQ(r.atoms.front().energy_eV, r.recoil_energy_eV);
    double v = 0.0, p = 0.0;
    for (const auto &a : r.atoms) {
      v += a.vacancy ? 1.0 : 0.0;
      p += a.replacement ? 1.0 : 0.0;
    }
    // TRIM also counts sputtered atoms set in motion below Ed as
    // displacements (and vacancies), without listing them.
    EXPECT_LE(v, r.displacements);
    EXPECT_DOUBLE_EQ(p, r.replacements);
    EXPECT_DOUBLE_EQ(r.vacancies, r.displacements - r.replacements);
    disp += r.displacements;
    repl += r.replacements;
  }
  double sd = 0.0, sr = 0.0;
  for (const auto &i : t.ion_summaries) {
    sd += i.displacements;
    sr += i.replacements;
    EXPECT_DOUBLE_EQ(i.interstitials, i.displacements - i.sputtered - i.transmitted);
    EXPECT_DOUBLE_EQ(i.vacancies, i.displacements - i.replacements);
  }
  // A cascade whose primary is not displaced is not listed, but what it
  // does is in the ion's totals.
  EXPECT_GE(sd, disp);
  EXPECT_LE(sd, disp * 1.05);
  EXPECT_DOUBLE_EQ(sr, repl);
  // TRIM.IN round trip of the setting.
  std::ostringstream out;
  write_trim_in(out, c, extras_for(c));
  std::istringstream in(out.str());
  const auto d = parse_trim_in(in);
  EXPECT_TRUE(d.collision_recoil_details);
}

TEST(TrimEngine, CancellationStopsEarly) {
  auto c = fe_on_fe(100.0, DamageCalculation::full_cascade, 100000);
  TrimSimulation s(c);
  s.request_cancel();
  s.run();
  EXPECT_TRUE(s.is_cancelled());
  EXPECT_LT(s.completed_count(), 100000u);
}

TEST(TrimEngine, WritesTrimFormatOutputs) {
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 20);
  c.collision_details_ions = 2;
  TrimSimulation s(c);
  s.run();
  const auto dir = temp_dir("trim_outputs");
  s.write_outputs(dir);
  for (const char *f : {"RANGE.txt", "VACANCY.txt", "NOVAC.txt", "IONIZ.txt", "PHONON.txt",
                        "E2RECOIL.txt", "LATERAL.txt", "TDATA.txt", "BACKSCAT.txt",
                        "TRANSMIT.txt", "SPUTTER.txt", "TRIMOUT.txt", "RANGE_3D.txt",
                        "COLLISON.txt", "PKA_SPECTRUM.txt", "summary.json", "results.json",
                        "csv/range.csv", "csv/vacancies.csv", "csv/exits.csv"}) {
    EXPECT_TRUE(std::filesystem::exists(dir / f)) << f;
  }
  std::ifstream range(dir / "RANGE.txt");
  std::stringstream text;
  text << range.rdbuf();
  EXPECT_NE(text.str().find("Ion Average Range"), std::string::npos);
  const auto json = s.summary_json();
  EXPECT_NE(json.find("\"vacancies\""), std::string::npos);
  EXPECT_NE(json.find("\"mean_A\""), std::string::npos);
}

TEST(TrimEngine, SnapshotTalliesFeedTheGui) {
  auto c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 20);
  c.trajectory_ions = 3;
  TrimSimulation s(c);
  s.run();
  const auto live = s.snapshot_tallies();
  EXPECT_GT(live.total_displacements, 0u);
  std::vector<CollisionPoint> path;
  EXPECT_TRUE(s.try_pop_trajectory(path));
  EXPECT_FALSE(path.empty());
  EXPECT_GT(in_units(s.depth_extent(), angstrom), 0.0);
}

// ---------------------------------------------------------------------------
// Statistical regression (built-in stopping; deterministic but compiler
// dependent at the ulp level, hence the tolerances)
// ---------------------------------------------------------------------------

TEST(TrimRegression, FeOnFe100keVQuickCalculation) {
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 400);
  c.threads = 0;
  TrimSimulation s(c);
  s.run();
  const auto &t = s.tallies();
  // TRIM (SRIM-2013 stopping): R = 346.5 A, 658.3 vacancies per ion. The
  // built-in stopping is within ~25 % of SRIM's here, so allow 15 %.
  EXPECT_LT(rel(t.range.mean(), 346.5), 0.15);
  EXPECT_LT(rel(t.total_vacancies / static_cast<double>(t.ions), 658.3), 0.15);
}

// With the real SRIM-2013 tables (AMCSET_SRIM_STOPPING_TABLE) and pair
// potentials (AMCSET_SRIM_DATA = the SRIM "Data" directory) the SRIM mode
// reproduces TRIM. Reference values: SRIM-2013.00 TRIM, 1000-5000 ions.
TEST(TrimRegression, AgainstTrimWithSrimTablesIfAvailable) {
  const char *table = std::getenv("AMCSET_SRIM_STOPPING_TABLE");
  if (!table) GTEST_SKIP() << "AMCSET_SRIM_STOPPING_TABLE not set";
  struct Case {
    double energy_keV;
    DamageCalculation damage;
    double range_A;
    double vacancies;
    double tolerance;
  };
  for (const Case &k : {Case{100.0, DamageCalculation::kinchin_pease, 346.5, 658.29, 0.03},
                        Case{10.0, DamageCalculation::kinchin_pease, 54.16, 75.15, 0.03},
                        Case{100.0, DamageCalculation::full_cascade, 350.7, 1313.1, 0.04}}) {
    auto c = fe_on_fe(k.energy_keV, k.damage, 600, k.energy_keV == 10.0 ? 1200.0 : 2500.0);
    c.stopping_model = physics::ElectronicStoppingModel::srim_tables;
    c.stopping_table = table;
    if (const char *data = std::getenv("AMCSET_SRIM_DATA")) c.srim_data_directory = data;
    c.threads = 0;
    TrimSimulation s(c);
    s.run();
    const auto &t = s.tallies();
    // At least three standard errors of the mean range of this run (TRIM's
    // own value carries a similar error).
    const double se = t.range.stddev() / std::sqrt(t.range.n) / k.range_A;
    EXPECT_LT(rel(t.range.mean(), k.range_A), std::max(k.tolerance, 3.0 * std::sqrt(2.0) * se))
        << k.energy_keV;
    EXPECT_LT(rel(t.total_vacancies / static_cast<double>(t.ions), k.vacancies), k.tolerance)
        << k.energy_keV;
  }
}

// ---------------------------------------------------------------------------
// TRIM.DAT particle sources, EXYZ, continuing and pausing
// ---------------------------------------------------------------------------

namespace {

ParticleList sample_particles(unsigned z, std::size_t n, double emax_eV, double depth_A) {
  ParticleList list;
  list.title = "test particles";
  for (std::size_t i = 0; i < n; ++i) {
    SourceParticle p;
    p.name = "P" + std::to_string(i);
    p.z = z;
    p.energy = (0.2 + 0.8 * static_cast<double>(i % 5) / 4.0) * emax_eV * electron_volt;
    p.x = (i % 3 == 0 ? 0.0 : depth_A * static_cast<double>(i % 7) / 7.0) * angstrom;
    const double c = i % 4 == 0 ? 1.0 : 0.8;
    p.cos_x = c;
    p.cos_y = std::sqrt(1.0 - c * c);
    p.cos_z = 0.0;
    list.particles.push_back(p);
  }
  return list;
}

} // namespace

TEST(TrimDat, ReadsWritesAndSkipsBadLines) {
  const auto list = sample_particles(26, 6, 1e5, 1000.0);
  std::stringstream ss;
  write_trim_dat(ss, list);
  std::string text = ss.str();
  EXPECT_NE(text.find("\r\n"), std::string::npos);
  EXPECT_EQ(text.back(), '\x1a');
  text.insert(text.size() - 1, "BAD   26  abc 0 0 0 1 0 0\r\n");
  std::istringstream in(text);
  const auto back = parse_trim_dat(in);
  EXPECT_EQ(back.title, "test particles");
  ASSERT_EQ(back.particles.size(), 6u);
  EXPECT_EQ(back.warnings.size(), 1u);
  EXPECT_EQ(back.particles[2].name, "P2");
  EXPECT_NEAR(in_units(back.particles[3].energy, electron_volt),
              in_units(list.particles[3].energy, electron_volt), 1e-3);
  EXPECT_NEAR(back.particles[1].cos_x, 0.8, 1e-7);
  // TRIM's own number styles.
  std::istringstream odd("1\n2\n3\n4\n5\n6\n7\ntitle line\n9\n10\n"
                         " A-1  18 12345  0 0 0  1.00000 .000000 -.000000\n"
                         " AA#1, 18, 1.31E4, 123, 0, -154, 0.34234, -.336437, -.017437\n\x1a");
  const auto o = parse_trim_dat(odd);
  ASSERT_EQ(o.particles.size(), 2u);
  EXPECT_NEAR(in_units(o.particles[1].energy, electron_volt), 13100.0, 1e-9);
  EXPECT_NEAR(o.particles[1].cos_y, -0.336437, 1e-12);
}

TEST(TrimDat, IonSourceRunsEveryValidLine) {
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 0, 3000.0);
  c.source = TrimConfig::Source::particle_ions;
  auto list = sample_particles(26, 40, 1e5, 2000.0);
  list.particles[5].energy = 2e5 * electron_volt; // above the ion energy: skipped
  list.particles[6].x = 5000.0 * angstrom;        // outside the target: skipped
  c.particles = list.particles;
  c.trim_code();
  EXPECT_EQ(c.trim_code(), 4);
  TrimSimulation s(c);
  EXPECT_EQ(s.total_count(), 38u);
  EXPECT_EQ(s.warnings().size(), 2u);
  s.run();
  EXPECT_EQ(s.tallies().ions, 38u);
  EXPECT_GT(s.tallies().range.n, 20.0);
  // Ions starting deep inside stop deeper than the shallow ones on average.
  EXPECT_GT(s.tallies().range.mean(), 300.0);
}

TEST(TrimDat, RecoilSourceMakesCascadesWithoutAnIon) {
  for (auto damage : {DamageCalculation::kinchin_pease, DamageCalculation::full_cascade,
                      DamageCalculation::monolayer}) {
    const bool as_ions = damage == DamageCalculation::monolayer; // TRIM's code 7
    auto c = fe_on_fe(50.0, damage, 0, 5000.0);
    c.source = TrimConfig::Source::particle_recoils;
    auto list = sample_particles(26, 60, 2e4, 4000.0);
    for (auto &p : list.particles) p.x += 500.0 * angstrom;
    list.particles[3].z = 29; // not a target element: skipped
    c.particles = list.particles;
    EXPECT_EQ(c.trim_code(), as_ions ? 7 : 6);
    TrimSimulation s(c);
    EXPECT_EQ(s.total_count(), 59u);
    s.run();
    const auto &t = s.tallies();
    EXPECT_EQ(t.ions, 59u);
    // As in TRIM, code 6 books each event as an ion stopping where the atom
    // starts; code 7 moves the atom as an ion.
    EXPECT_EQ(t.range.n, 59.0);
    if (as_ions) {
      EXPECT_GT(t.total_ionization_ion, 0.0);
    } else {
      EXPECT_EQ(t.total_ionization_ion, 0.0); // no ion energy loss
      double sum = 0.0;
      for (std::size_t i = 0; i < c.particles.size(); ++i) {
        if (i != 3) sum += in_units(c.particles[i].x, angstrom);
      }
      EXPECT_NEAR(t.range.mean(), sum / 59.0, 1e-6);
    }
    // About NRT: 0.8 * 0.6 T / (2 Ed) for 4-20 keV Fe recoils.
    EXPECT_GT(t.total_vacancies / 59.0, 20.0);
    EXPECT_LT(t.total_vacancies / 59.0, 200.0);
  }
}

// source.repeat: every particle runs that many times in a row, exactly as
// a list with each line repeated; ions limits the runs.
TEST(TrimDat, RepeatRunsEachParticleSeveralTimes) {
  auto c = fe_on_fe(50.0, DamageCalculation::full_cascade, 0, 5000.0);
  c.source = TrimConfig::Source::particle_recoils;
  const auto list = sample_particles(26, 2, 2e4, 2000.0);
  auto repeated = c;
  repeated.particles = list.particles;
  repeated.particle_repeat = 3;
  EXPECT_EQ(repeated.particle_count(), 6u);
  auto listed = c;
  for (const auto &p : list.particles) {
    for (int k = 0; k < 3; ++k) listed.particles.push_back(p);
  }
  TrimSimulation a(repeated);
  TrimSimulation b(listed);
  EXPECT_EQ(a.total_count(), 6u);
  a.run();
  b.run();
  EXPECT_EQ(a.tallies().ions, 6u);
  EXPECT_EQ(a.tallies().total_vacancies, b.tallies().total_vacancies);
  EXPECT_EQ(a.tallies().total_phonons_recoils, b.tallies().total_phonons_recoils);
  repeated.ion_count = 4;
  EXPECT_EQ(TrimSimulation(repeated).total_count(), 4u);
  repeated.particle_repeat = 0;
  EXPECT_THROW(repeated.validate(), std::invalid_argument);
}

TEST(TrimEngine, ExyzRecordsTheIonEveryInterval) {
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 5);
  c.exyz_interval = 10.0 * kilo_electron_volt;
  TrimSimulation s(c);
  s.run();
  const auto &r = s.tallies().exyz;
  ASSERT_GT(r.size(), 5u * 10u);
  EXPECT_EQ(r.front().ion, 1u);
  EXPECT_DOUBLE_EQ(r.front().energy_keV, 100.0);
  EXPECT_DOUBLE_EQ(r.front().x, 0.0);
  // Energies fall within each ion; one record per 10 keV step at most.
  std::size_t per_ion = 0;
  for (std::size_t i = 1; i < r.size(); ++i) {
    if (r[i].ion == r[i - 1].ion) {
      EXPECT_LE(r[i].energy_keV, r[i - 1].energy_keV);
      ++per_ion;
    }
  }
  EXPECT_LE(per_ion, 5u * 12u);
  // As in TRIM, the record at the start sets the next one two steps down:
  // the second record of an ion is below 90 keV.
  for (std::size_t i = 1; i < r.size(); ++i) {
    if (r[i].ion != r[i - 1].ion || r[i - 1].x != 0.0 || r[i].energy_keV == 0.0) continue;
    EXPECT_LT(r[i].energy_keV, 90.0);
  }
}

TEST(TrimEngine, ContinuingGivesTheSameResultsAsOneRun) {
  auto c = fe_on_fe(50.0, DamageCalculation::full_cascade, 40);
  TrimSimulation whole(c);
  whole.run();
  auto first = c;
  first.ion_count = 24;
  TrimSimulation a(first);
  a.run();
  auto second = c;
  second.ion_count = 16;
  second.first_ion = 24;
  TrimSimulation b(second);
  b.continue_from(a.tallies());
  b.run();
  const auto &x = whole.tallies();
  const auto &y = b.tallies();
  EXPECT_EQ(x.ions, y.ions);
  EXPECT_DOUBLE_EQ(x.total_vacancies, y.total_vacancies);
  EXPECT_DOUBLE_EQ(x.range.s1, y.range.s1);
  EXPECT_EQ(x.exits.size(), y.exits.size());
}

TEST(TrimEngine, PauseHoldsTheWorkers) {
  auto c = fe_on_fe(200.0, DamageCalculation::kinchin_pease, 4000);
  TrimSimulation s(c);
  s.set_paused(true);
  std::thread worker([&] { s.run(); });
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  EXPECT_EQ(s.completed_count(), 0u);
  EXPECT_TRUE(s.is_paused());
  s.set_paused(false);
  while (s.completed_count() == 0) std::this_thread::sleep_for(std::chrono::milliseconds(5));
  s.set_paused(true);
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  const auto held = s.completed_count();
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  EXPECT_LE(s.completed_count(), held + 32u); // blocks in flight may finish
  s.request_cancel();                          // also releases the pause
  worker.join();
  EXPECT_LT(s.completed_count(), 4000u);
  EXPECT_EQ(s.snapshot().ions, s.completed_count());
}

TEST(TrimTables, DatasetMatchesTheTallies) {
  auto c = fe_on_fe(50.0, DamageCalculation::full_cascade, 30);
  TrimSimulation s(c);
  s.run();
  const auto d = make_dataset(s);
  for (const char *k : {"range", "lateral", "vacancies", "replacements", "damage", "ionization",
                        "phonons", "energy_to_recoils", "sputtering", "surface_energy",
                        "recoil_spectrum", "exits", "ion_positions"}) {
    EXPECT_NE(d.find(k), nullptr) << k;
  }
  // Integral of the ion distribution = fraction of ions stopped.
  const auto *range = d.find("range");
  double sum = 0.0;
  for (double v : range->find("ions")->values) sum += v;
  const double width_cm = in_units(s.bin_width() * angstrom, centimeter);
  EXPECT_NEAR(sum * width_cm, s.tallies().range.n / 30.0, 1e-9);
  // Vacancies per ion: table integral and scalar agree (TRIM's total also
  // counts the vacancies of sub-threshold atoms that were sputtered).
  double vac = 0.0;
  for (double v : d.find("vacancies")->find("total")->values) vac += v * s.bin_width();
  EXPECT_NEAR(vac, s.tallies().total_vacancies / 30.0, 0.02 * vac);
  // Surface arrivals include every sputtered atom.
  const auto *se = d.find("surface_energy");
  ASSERT_NE(se->find("yield_above_Fe"), nullptr);
  EXPECT_GE(se->find("yield_above_Fe")->values.front() + 1e-12, s.tallies().sputtered[0] / 30.0);
}

// ---------------------------------------------------------------------------
// Fast approximations, partial reruns, engine statistics
// ---------------------------------------------------------------------------

// physics.fast_math changes the arithmetic, not the physics: with the same
// seed the histories mostly stay the same, so the results agree far more
// closely than the statistics would require (a systematic error of the
// approximations would show here).
TEST(TrimEngine, FastMathAgreesWithTheExactTransport) {
  for (const auto damage : {DamageCalculation::full_cascade, DamageCalculation::kinchin_pease}) {
    auto c = fe_on_fe(20.0, damage, 200);
    TrimSimulation exact(c);
    exact.run();
    c.fast_math = true;
    TrimSimulation fast(c);
    fast.run();
    const auto &a = exact.tallies();
    const auto &b = fast.tallies();
    EXPECT_EQ(b.ions, 200u);
    EXPECT_LT(rel(b.total_vacancies, a.total_vacancies), 0.02) << to_string(damage);
    EXPECT_LT(rel(b.range.mean(), a.range.mean()), 0.02) << to_string(damage);
    EXPECT_LT(rel(b.total_ionization_ion + b.total_ionization_recoils,
                  a.total_ionization_ion + a.total_ionization_recoils),
              0.01)
        << to_string(damage);
  }
}

TEST(TrimEngine, FastMathIsRecordedInTheOutputs) {
  auto c = fe_on_fe(5.0, DamageCalculation::kinchin_pease, 16);
  c.fast_math = true;
  const auto dir = temp_dir("fast_math_outputs");
  c.output_directory = dir.string();
  c.outputs.formats = {output::Format::trim};
  TrimSimulation s(c);
  s.run();
  std::ifstream tdata(dir / "TDATA.txt");
  std::stringstream text;
  text << tdata.rdbuf();
  EXPECT_NE(text.str().find("Fast approximations    = on"), std::string::npos);
  EXPECT_NE(s.summary_json().find("\"fast_math\": true"), std::string::npos);
  // ... and only then: the default outputs are unchanged.
  c.fast_math = false;
  TrimSimulation e(c);
  EXPECT_EQ(e.summary_json().find("fast_math"), std::string::npos);
}

// A crash report names the first ion of the failing ion's block of 16: a run
// that starts there (first_ion) transports the following ions exactly as the
// full run did, because the recoil slots carried from ion to ion within a
// block start the same way.
TEST(TrimEngine, FirstIonRerunsTheSameIons) {
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 48);
  TrimSimulation full(c);
  full.run();
  c.first_ion = 16;
  c.ion_count = 16;
  TrimSimulation part(c);
  part.run();
  std::map<std::size_t, IonPosition> want;
  for (const auto &p : full.tallies().positions) want[p.ion] = p;
  ASSERT_FALSE(part.tallies().positions.empty());
  for (const auto &p : part.tallies().positions) {
    ASSERT_GE(p.ion, 17u);
    ASSERT_LE(p.ion, 32u);
    ASSERT_EQ(want.count(p.ion), 1u) << p.ion;
    EXPECT_EQ(p.x, want[p.ion].x) << p.ion;
    EXPECT_EQ(p.y, want[p.ion].y) << p.ion;
    EXPECT_EQ(p.z, want[p.ion].z) << p.ion;
  }
}

TEST(TrimEngine, SnapshotTalliesSummarizeTheTrimTallies) {
  auto c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 24);
  TrimSimulation s(c);
  s.run();
  const auto live = s.snapshot_tallies();
  const auto &t = s.tallies();
  EXPECT_EQ(live.incidents, 24u);
  EXPECT_EQ(live.total_displacements, static_cast<std::size_t>(std::llround(t.total_displacements)));
  double stopped = 0.0;
  for (double v : live.range_histogram.bins) stopped += v;
  double tally = 0.0;
  for (double v : t.ion_stopped.values) tally += v;
  EXPECT_DOUBLE_EQ(stopped, tally);
  EXPECT_NEAR(static_cast<double>(live.range_stats.mean), t.range.mean() * 1e-10, 1e-18);
  EXPECT_EQ(live.range_histogram.bins.size(), c.depth_bins);
}

// A recoil that reaches the front surface is held back by its own surface
// binding energy. The element entry of a recoil belongs to the layer it was
// struck in; it was once looked up in the layer it moved through last, so
// atoms from below a thin top layer were judged by a top-layer element's Es
// (or read beyond the top layer's element list when it has fewer elements).
TEST(TrimEngine, SputteredAtomsUseTheirOwnSurfaceBindingEnergy) {
  TrimConfig c;
  c.ion = IonBeam::defaults(18, 3.0 * kilo_electron_volt);
  TargetLayer top;
  top.name = "Cr";
  top.width = 4.0 * angstrom; // about two atomic layers
  top.density = 7.19 * gram_per_cubic_centimeter;
  top.elements.push_back(TargetElement::defaults(24));
  c.damage = DamageCalculation::full_cascade;
  c.ion_count = 300;
  c.threads = 2;
  // Fe and Ni below the Cr film: their sputtering yield must follow their
  // own Es. (With the old lookup they were judged by Cr's Es, or by memory
  // beyond the one-element Cr list, and the yields would not change.)
  const auto yield_below = [&](double es_eV) {
    TargetLayer below;
    below.name = "FeNi";
    below.width = 300.0 * angstrom;
    below.density = 8.0 * gram_per_cubic_centimeter;
    for (unsigned z : {26u, 28u}) {
      auto e = TargetElement::defaults(z, 0.5);
      e.surface_binding_energy = es_eV * electron_volt;
      below.elements.push_back(e);
    }
    c.layers = {top, below};
    TrimSimulation s(c);
    s.run();
    const auto &names = s.element_names();
    EXPECT_EQ(names.size(), 3u);
    double y = 0.0;
    for (std::size_t i = 0; i < names.size(); ++i) {
      if (names[i] != "Cr") y += s.tallies().sputtered[i];
    }
    return y;
  };
  const double easy = yield_below(1.0);
  const double hard = yield_below(30.0);
  EXPECT_GT(easy, 0.0);
  EXPECT_GT(easy, hard * 1.2) << "easy " << easy << " hard " << hard;
}

// Inputs outside what SRIM mode can simulate, or that make no physical
// sense, are refused with a message that names the field and the limit.
TEST(TrimConfig, RejectsInputsBeyondThePhysicalLimits) {
  const auto c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 10);
  const auto message = [](const TrimConfig &t) {
    try {
      t.validate();
    } catch (const std::invalid_argument &e) {
      return std::string(e.what());
    }
    return std::string();
  };
  auto bad = c;
  bad.ion.energy = 500.0 * giga_electron_volt;
  EXPECT_NE(message(bad).find("use at most 100 GeV"), std::string::npos) << message(bad);
  bad = c;
  bad.ion.energy = 0.5 * electron_volt; // below the 1 eV at which ions stop
  EXPECT_NE(message(bad).find("final energy"), std::string::npos) << message(bad);
  bad = c;
  bad.ion.energy = std::numeric_limits<double>::infinity() * electron_volt;
  EXPECT_NE(message(bad).find("finite"), std::string::npos);
  bad = c;
  bad.ion.mass = 2.0 * atomic_mass_unit; // not an iron atom
  EXPECT_NE(message(bad).find("is not a mass of Fe"), std::string::npos) << message(bad);
  bad = c;
  bad.layers[0].width = 1.0 * angstrom; // Fe atoms are 2.27 A apart
  EXPECT_NE(message(bad).find("thinner than one atomic layer"), std::string::npos) << message(bad);
  EXPECT_NE(message(bad).find("use at least 2.28 A"), std::string::npos) << message(bad);
  bad = c;
  bad.layers[0].width = -5.0 * angstrom;
  EXPECT_NE(message(bad).find("thickness must be positive"), std::string::npos);
  bad = c;
  bad.layers[0].density = 80.0 * gram_per_cubic_centimeter;
  EXPECT_NE(message(bad).find("denser than any material"), std::string::npos);
  bad = c;
  bad.layers[0].elements[0].displacement_energy = 50000.0 * electron_volt;
  EXPECT_NE(message(bad).find("is not physical"), std::string::npos);
  bad = c;
  bad.layers[0].elements[0].surface_binding_energy = -1.0 * electron_volt;
  EXPECT_NE(message(bad).find(">= 0"), std::string::npos);
  bad = c;
  bad.layers[0].compound_correction = 0.0;
  EXPECT_FALSE(message(bad).empty());
  bad = c;
  bad.depth_bins = 1000000;
  EXPECT_NE(message(bad).find("depth bins"), std::string::npos);
  bad = c;
  bad.window_start = 5000.0 * angstrom; // the target is 700 A thick
  EXPECT_NE(message(bad).find("beyond the target"), std::string::npos);
  bad = c;
  bad.final_energy = 5000.0 * electron_volt;
  EXPECT_NE(message(bad).find("too high"), std::string::npos);
  bad = c;
  bad.ion.angle = 90.0 * common::degree;
  EXPECT_FALSE(message(bad).empty());
  bad = c;
  bad.stopping_model = physics::ElectronicStoppingModel::srim_tables;
  bad.stopping_table = "no/such/stopping_table.amcset";
  EXPECT_NE(message(bad).find("does not exist"), std::string::npos) << message(bad);
  EXPECT_NO_THROW(bad.validate(false)); // (while a path is being edited)
  bad = c;
  bad.srim_data_directory = "no/such/SRIM/Data";
  EXPECT_NE(message(bad).find("is not a directory"), std::string::npos) << message(bad);
  // At the limits themselves: accepted.
  auto ok = c;
  ok.ion.energy = 100.0 * giga_electron_volt;
  ok.layers[0].width = 2.3 * angstrom;
  EXPECT_EQ(message(ok), "");
}

// Malformed TRIM.IN files fail with a message, never by indexing past the
// numbers they give (a negative element count once wrapped around).
TEST(TrimIn, RejectsMalformedCountsAndCodes) {
  const auto with = [](const std::string &from, const std::string &to) {
    std::string t = trim_in_text;
    const auto p = t.find(from);
    EXPECT_NE(p, std::string::npos) << from;
    t.replace(p, from.size(), to);
    return t;
  };
  const auto fails = [](const std::string &text) {
    std::istringstream in(text);
    try {
      parse_trim_in(in);
    } catch (const std::runtime_error &e) {
      return std::string(e.what());
    }
    return std::string();
  };
  EXPECT_NE(fails(with("\"FeCr on Ni\"       3", "\"FeCr on Ni\"       -1")).find("number of target elements"),
            std::string::npos);
  EXPECT_NE(fails(with("\"FeCr on Ni\"       3", "\"FeCr on Ni\"       1e30")).find("whole number"),
            std::string::npos);
  EXPECT_NE(fails(with("     2   4.003   2000", "     nan   4.003   2000")), "");
  EXPECT_NE(fails(with("     2   4.003   2000", "     200   4.003   2000")).find("Z"), std::string::npos);
  EXPECT_NE(fails(with("                      2   ", "                      9   ")).find("cascade type"),
            std::string::npos);
  EXPECT_NE(fails(with("123       0", "-5       0")).find("seed"), std::string::npos);
  EXPECT_NE(fails(with("Atom 3 = Ni =      28", "Atom 3 = Ni =      -28")).find("Z"), std::string::npos);
  // Truncated anywhere: an error, not a crash.
  const std::string full = trim_in_text;
  for (std::size_t cut = 0; cut < full.size(); cut += 37) {
    std::istringstream in(full.substr(0, cut));
    try {
      parse_trim_in(in);
    } catch (const std::exception &) {
    }
  }
}

TEST(TrimDat, SkipsLinesThatCannotBeTransported) {
  std::istringstream in(
      "TRIM.DAT\n1\n2\n3\n4\n5\n6\nlist\nEvent Atom Energy Depth Y Z cosx cosy cosz\n"
      "Name Numb (eV) (A) (A) (A)\n"
      "A0001 26 1000 10 0 0 1 0 0\n"
      "A0002 nan 1000 10 0 0 1 0 0\n"
      "A0003 26 1000 10 0 0 nan 0 0\n"
      "A0004 26 1000 10 0 0 0 0 0\n"
      "A0005 26 inf 10 0 0 1 0 0\n"
      "A0006 26.5 1000 10 0 0 1 0 0\n"
      "A0007 1e30 1000 10 0 0 1 0 0\n");
  const auto list = parse_trim_dat(in);
  ASSERT_EQ(list.particles.size(), 1u);
  EXPECT_EQ(list.particles[0].name, "A0001");
  EXPECT_EQ(list.warnings.size(), 6u);
}

// ---------------------------------------------------------------------------
// Automatic target depth and plausibility warnings
// ---------------------------------------------------------------------------

namespace {

std::string joined(const std::vector<std::string> &v) {
  std::string s;
  for (const auto &x : v) s += x + "\n";
  return s;
}

//! A thickness the automatic depth may choose: at most one decimal in the
//! unit of its size (A, nm, um, mm, m), whole numbers from 100 of it.
bool clean_thickness(double a_A) {
  double unit = 1.0;
  if (a_A >= 1e10) unit = 1e10;
  else if (a_A >= 1e7) unit = 1e7;
  else if (a_A >= 1e4) unit = 1e4;
  else if (a_A >= 10.0) unit = 10.0;
  const double v = a_A / unit * (a_A / unit < 100.0 ? 10.0 : 1.0);
  return std::abs(v - std::round(v)) < 1e-6 * std::max(1.0, v);
}

} // namespace

// The path length exceeds the mean projected range; the first guess of the
// automatic depth reaches the deepest ions without being many times deeper.
TEST(TrimAutoDepth, EstimateCoversTheDeepestIons) {
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 300, 1e5);
  const double path = in_units(estimate_penetration(c), angstrom);
  const double estimate = in_units(estimate_deepest(c), angstrom);
  TrimSimulation sim(c);
  sim.run();
  const auto &p = sim.tallies().positions;
  ASSERT_EQ(p.size(), 300u);
  double deepest = 0.0;
  double mean = 0.0;
  for (const auto &x : p) {
    deepest = std::max(deepest, x.x);
    mean += x.x / static_cast<double>(p.size());
  }
  EXPECT_GT(path, mean) << "path " << path << " A, mean " << mean << " A";
  EXPECT_GT(estimate, deepest) << "estimate " << estimate << " A, deepest " << deepest << " A";
  EXPECT_LT(estimate, 1.5 * deepest) << "estimate " << estimate << " A, deepest " << deepest << " A";
}

// From a far too thin and a far too thick last layer, the sized target stops
// every test ion, the deepest beyond auto_depth.shrink_below of the target
// (a too thick one is made thinner); only the thickness changes, it is a
// clean number, and it does not depend on the starting thickness.
TEST(TrimAutoDepth, DeepestTestIonStopsInTheLastMargin) {
  std::vector<double> sized;
  for (const double start_A : {50.0, 1e6}) {
    auto c = fe_on_fe(100.0, DamageCalculation::full_cascade, 10, start_A);
    c.auto_depth.enabled = true;
    c.auto_depth.test_ions = 300;
    const auto r = size_target(c);
    EXPECT_TRUE(r.changed);
    EXPECT_EQ(r.transmitted, 0u);
    EXPECT_GE(r.runs, 1u);
    EXPECT_LE(r.runs, 12u);
    const double total = in_units(c.target_thickness(), angstrom);
    const double deepest = in_units(r.deepest, angstrom);
    EXPECT_LE(deepest, total);
    EXPECT_GE(deepest, c.auto_depth.shrink_below * total * (1.0 - 1e-9)) << joined(r.notes);
    EXPECT_TRUE(clean_thickness(total)) << total << " A";
    EXPECT_EQ(in_units(r.last_layer_width, angstrom), total);
    EXPECT_EQ(c.damage, DamageCalculation::full_cascade);
    EXPECT_EQ(c.ion_count, 10u);
    EXPECT_FALSE(c.fast_math);
    ASSERT_FALSE(r.notes.empty());
    EXPECT_NE(r.notes.front().find("last layer ('Fe')"), std::string::npos) << r.notes.front();
    EXPECT_NO_THROW(c.validate());
    sized.push_back(total);
  }
  EXPECT_EQ(sized[0], sized[1]);
}

// Ions that stop above the last layer: it is made as thin as it can be (the
// domain shrinks to the ions), the layers above are kept.
TEST(TrimAutoDepth, IonsStoppingAboveTheLastLayerMakeItThin) {
  auto c = fe_on_fe(10.0, DamageCalculation::kinchin_pease, 10, 1e5); // 10 um; 10 keV Fe stops in nm
  auto back = c.layers.front();
  back.name = "Fe back";
  back.width = 1000.0 * angstrom;
  c.layers.push_back(back);
  const auto r = size_target(c);
  EXPECT_TRUE(r.changed);
  EXPECT_GE(r.runs, 1u);
  EXPECT_EQ(in_units(c.layers.front().width, angstrom), 1e5);
  const double last = in_units(c.layers.back().width, angstrom);
  EXPECT_LT(last, 5.0); // about a monolayer
  EXPECT_TRUE(clean_thickness(last)) << last;
  EXPECT_NE(joined(r.notes).find("as thin as it can be"), std::string::npos) << joined(r.notes);
  EXPECT_NO_THROW(c.validate());
}

// The threshold: a deepest ion between shrink_below and the back keeps the
// target; a higher threshold sizes it tighter.
TEST(TrimAutoDepth, ShrinkThresholdDecides) {
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 10, 1e6);
  c.auto_depth.enabled = true;
  c.auto_depth.test_ions = 200;
  c.auto_depth.shrink_below = 0.85;
  auto r = size_target(c);
  const double total = in_units(c.target_thickness(), angstrom);
  EXPECT_GE(in_units(r.deepest, angstrom), 0.85 * total * (1.0 - 1e-9)) << joined(r.notes);
  // A sized target is kept: sizing it again changes nothing.
  r = size_target(c);
  EXPECT_EQ(in_units(c.target_thickness(), angstrom), total) << joined(r.notes);
  auto bad = c;
  bad.auto_depth.shrink_below = 1.0;
  EXPECT_THROW(size_target(bad), std::invalid_argument);
}

TEST(TrimAutoDepth, RecoilSourcesAndCancellationKeepTheTarget) {
  auto c = fe_on_fe(100.0, DamageCalculation::full_cascade, 10, 1e6);
  c.source = TrimConfig::Source::particle_recoils;
  auto r = size_target(c);
  EXPECT_FALSE(r.changed);
  EXPECT_NE(joined(r.notes).find("recoil source"), std::string::npos);
  c.source = TrimConfig::Source::beam;
  const std::atomic<bool> cancel{true};
  r = size_target(c, &cancel);
  EXPECT_FALSE(r.changed);
  EXPECT_EQ(in_units(c.layers.back().width, angstrom), 1e6);
  EXPECT_NE(joined(r.notes).find("cancelled"), std::string::npos);
}

// A run on a thin target reports its transmitted ions while it runs (the
// front ends warn about them).
TEST(TrimAutoDepth, TransmittedIonsAreCountedDuringTheRun) {
  auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 64, 30.0);
  TrimSimulation sim(c);
  EXPECT_EQ(sim.transmitted_so_far(), 0u);
  sim.run();
  EXPECT_GT(sim.transmitted_so_far(), 0u);
  EXPECT_EQ(sim.transmitted_so_far(), sim.tallies().transmitted);
}

TEST(TrimConfig, PlausibilityWarnings) {
  const auto c = fe_on_fe(100.0, DamageCalculation::kinchin_pease, 10);
  EXPECT_EQ(joined(c.plausibility_warnings()), "");
  // 1 keV ions into 1 m of iron: they stop within nanometres.
  auto deep = fe_on_fe(1.0, DamageCalculation::kinchin_pease, 10, 1e10);
  const auto w = joined(deep.plausibility_warnings());
  EXPECT_NE(w.find("the target is 1 m thick"), std::string::npos) << w;
  EXPECT_NE(w.find("automatic target depth"), std::string::npos) << w;
  // Densities far from iron's, and gas flags that do not fit the density.
  auto light = c;
  light.layers[0].density = 0.5 * gram_per_cubic_centimeter;
  EXPECT_NE(joined(light.plausibility_warnings()).find("is far from that of Fe"), std::string::npos);
  auto gas = c;
  gas.layers[0].gas = true;
  EXPECT_NE(joined(gas.plausibility_warnings()).find("marked as a gas"), std::string::npos);
  auto thin_solid = c;
  thin_solid.layers[0].density = 0.001 * gram_per_cubic_centimeter;
  EXPECT_NE(joined(thin_solid.plausibility_warnings()).find("density of a gas"), std::string::npos);
  // A target much thinner than the reach: a note (transmission studies are fine).
  const auto thin = fe_on_fe(1000.0, DamageCalculation::kinchin_pease, 10, 100.0);
  EXPECT_NE(joined(thin.plausibility_warnings()).find("most will be transmitted"), std::string::npos);
  // Relativistic protons.
  TrimConfig p = c;
  p.ion = IonBeam::defaults(1, 10.0 * giga_electron_volt);
  EXPECT_NE(joined(p.plausibility_warnings()).find("relativistic"), std::string::npos);
}

// Found by tools/perf/fuzz_inputs.py: iron of gas density marked as a solid
// (spacing ~1000 A, impact parameters of hundreds of Angstrom) once ended
// the run with "the ion's energy or depth is not a number". It is warned
// about, and it runs.
TEST(TrimEngine, SolidOfGasDensityRuns) {
  for (const auto damage : {DamageCalculation::kinchin_pease, DamageCalculation::full_cascade}) {
    auto c = fe_on_fe(1000.0, damage, 8, 12000.0);
    c.layers[0].density = 1e-7 * gram_per_cubic_centimeter;
    EXPECT_NE(joined(c.plausibility_warnings()).find("density of a gas"), std::string::npos);
    TrimSimulation sim(c);
    EXPECT_NO_THROW(sim.run());
    EXPECT_EQ(sim.tallies().ions, 8u);
    EXPECT_EQ(sim.tallies().transmitted, 8u);
  }
}

// ---------------------------------------------------------------------------
// Trajectory bookkeeping, the trajectory file and hybrid BCA-MD events
// ---------------------------------------------------------------------------

namespace {

std::string read_file(const std::filesystem::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

//! Fe -> Fe with every TRIM bug fixed (the energy then adds up) and Eb = 0.
TrimConfig md_case(double energy_keV, std::size_t ions) {
  auto c = fe_on_fe(energy_keV, DamageCalculation::full_cascade, ions, 2.0 * energy_keV + 3000.0);
  c.bugs.set_all(false);
  c.threads = 1;
  return c;
}

//! The totals the TRIM tables are made of.
std::vector<double> totals(const TrimTallies &t) {
  std::vector<double> v = {t.total_displacements, t.total_vacancies, t.total_replacements, t.total_interstitials,
                           t.total_ionization_ion, t.total_ionization_recoils, t.total_phonons_ion,
                           t.total_phonons_recoils, t.total_energy_to_recoils, t.total_binding_energy_loss,
                           static_cast<double>(t.backscattered), t.backscattered_energy, t.range.mean()};
  for (double x : t.ion_stopped.values) v.push_back(x);
  for (const auto &h : t.vacancies) v.insert(v.end(), h.values.begin(), h.values.end());
  for (double x : t.sputtered) v.push_back(x);
  return v;
}

} // namespace

// Recording trajectories, the trajectory file and MD events (record mode)
// changes no tally, for any number of threads, and the files do not depend
// on the thread count either.
TEST(TrimMdEvents, RecordModeChangesNoTallyAndFilesAreThreadIndependent) {
  const auto dir = temp_dir("md_record");
  auto plain = md_case(20.0, 40);
  plain.bugs = TrimBugs{}; // TRIM's behaviour, bugs included
  TrimSimulation a(plain);
  a.run();
  std::string first_events;
  for (const unsigned threads : {1u, 4u}) {
    auto c = plain;
    c.threads = threads;
    c.output_directory = (dir / ("t" + std::to_string(threads))).string();
    c.outputs.formats.clear(); // only the export files
    c.record_ion_positions = false;
    c.trajectory_ions = TrimConfig::all_ions;
    c.trajectory_export.file = "traj.csv";
    c.md_events.file = "events.csv";
    c.md_events.transition_energies = {500.0 * electron_volt, 2000.0 * electron_volt};
    TrimSimulation b(c);
    b.set_trajectory_polling(false);
    b.run();
    EXPECT_EQ(totals(a.tallies()), totals(b.tallies())) << threads << " threads";
    ASSERT_EQ(b.tallies().md.size(), 2u);
    EXPECT_GT(b.tallies().md[0].events, 0.0);
    // A lower E_c hands more, smaller pieces to MD.
    EXPECT_GT(b.tallies().md[0].events, b.tallies().md[1].events);
    const std::string ev = read_file(std::filesystem::path(c.output_directory) / "events.csv");
    ASSERT_GT(ev.size(), 100u);
    EXPECT_EQ(ev.rfind("ion_id,event_id,E_c_eV,kind,t_ps,x_A,y_A,z_A,E_eV,ux,uy,uz,gen,Z,mass_amu,track,"
                       "parent_track,parent_event\n", 0), 0u);
    if (first_events.empty()) first_events = ev;
    else EXPECT_EQ(ev, first_events) << "the event file depends on the thread count";
    EXPECT_GT(std::filesystem::file_size(std::filesystem::path(c.output_directory) / "traj.csv"), 1000u);
  }
}

// Handoff: what MD takes over plus what the BCA deposits and what leaves
// the target is the ion energy (with TRIM's energy-creating bugs fixed, as
// the no-handoff run shows), and fewer atoms are followed.
TEST(TrimMdEvents, HandoffConservesTheEnergy) {
  const double e0 = 20000.0;
  auto c = md_case(20.0, 60);
  TrimSimulation full(c);
  full.run();
  const double reference = energy_balance(full, e0);
  c.md_events.file = "events.csv";
  c.output_directory = temp_dir("md_handoff").string();
  c.md_events.transition_energies = {1000.0 * electron_volt};
  c.md_events.handoff = true;
  TrimSimulation h(c);
  h.run();
  const auto &t = h.tallies();
  ASSERT_EQ(t.md.size(), 1u);
  const double md_share = t.md[0].energy_eV / (e0 * static_cast<double>(t.ions));
  EXPECT_GT(md_share, 0.05);
  // trim2013 books a stopping recoil's last electronic loss twice (TRIM's
  // bookkeeping, kept with the bugs fixed): the reference is ~1.0006. Handing
  // the low-energy part to MD avoids most of that.
  const double total = energy_balance(h, e0) + md_share;
  EXPECT_NEAR(total, 1.0, 1e-3) << "BCA " << energy_balance(h, e0) << " + MD " << md_share;
  EXPECT_LE(std::abs(total - 1.0), std::abs(reference - 1.0) + 1e-4) << "reference " << reference;
  EXPECT_LT(t.recoils_followed, full.tallies().recoils_followed);
  EXPECT_LT(t.total_vacancies, full.tallies().total_vacancies);
}

// The recorded points: one birth per track, parents that exist and have
// the child's position (or the surface the child was reflected at), times
// that never decrease, and the collision energies.
TEST(TrimTrajectories, PointsFormATreeWithTimes) {
  auto c = md_case(10.0, 6);
  c.trajectory_ions = TrimConfig::all_ions;
  TrimSimulation s(c);
  s.run();
  std::vector<common::CollisionPoint> pts;
  std::size_t ions = 0;
  while (s.try_pop_trajectory(pts)) {
    ++ions;
    ASSERT_FALSE(pts.empty());
    EXPECT_EQ(pts.front().kind, common::CollisionPoint::Kind::birth);
    EXPECT_EQ(pts.front().track, 0u);
    EXPECT_EQ(pts.front().parent_track, common::CollisionPoint::no_parent);
    std::map<std::uint32_t, double> last_time;
    std::set<std::uint32_t> born;
    for (std::size_t i = 0; i < pts.size(); ++i) {
      const auto &p = pts[i];
      EXPECT_EQ(p.point_index, i);
      EXPECT_EQ(p.ion, pts.front().ion);
      if (p.kind == common::CollisionPoint::Kind::birth) {
        EXPECT_TRUE(born.insert(p.track).second) << "two births of track " << p.track;
        if (p.track != 0) {
          ASSERT_LT(p.parent_point, i);
          const auto &q = pts[p.parent_point];
          EXPECT_EQ(q.track, p.parent_track);
          EXPECT_LE(q.time, p.time + 1e-21);
          EXPECT_NEAR(in_units(q.coordinate.x_, angstrom), in_units(p.coordinate.x_, angstrom), 1e-9);
          EXPECT_NEAR(in_units(q.coordinate.z_, angstrom), in_units(p.coordinate.z_, angstrom), 1e-3);
          EXPECT_GT(q.energy_transferred_eV, 0.0);
        }
      } else {
        EXPECT_TRUE(born.count(p.track)) << "point before the birth of track " << p.track;
      }
      if (last_time.count(p.track)) {
        EXPECT_GE(p.time, last_time[p.track]);
      }
      last_time[p.track] = p.time;
      if (p.kind == common::CollisionPoint::Kind::collision && p.track != 0 && p.energy_before_eV > 0.0) {
        EXPECT_NEAR(p.energy_before_eV - p.energy_transferred_eV - p.electronic_loss_eV, p.energy_eV,
                    1e-9 * p.energy_before_eV + 1e-12);
      }
    }
    EXPECT_GT(pts.back().time, 0.0);
  }
  EXPECT_EQ(ions, 6u);
}

// A target atom that never moves (TRIM's stale energy 0) leaves no points.
TEST(TrimTrajectories, AtomsThatNeverMoveAreNotRecorded) {
  auto c = fe_on_fe(5000.0, DamageCalculation::full_cascade, 2, 3e4);
  c.trajectory_ions = 2;
  TrimSimulation s(c);
  s.run();
  std::vector<common::CollisionPoint> pts;
  while (s.try_pop_trajectory(pts)) {
    std::map<std::uint32_t, bool> moving;
    for (const auto &p : pts) moving[p.track] = moving[p.track] || p.energy_eV > 0.0;
    for (const auto &[track, any] : moving) EXPECT_TRUE(any) << "track " << track << " never moves";
  }
}

TEST(TrimMdEvents, SettingsAreValidated) {
  const auto message = [](const TrimConfig &t) {
    try {
      t.validate();
    } catch (const std::invalid_argument &e) {
      return std::string(e.what());
    }
    return std::string();
  };
  auto c = md_case(10.0, 5);
  c.md_events.file = "e.csv";
  EXPECT_NE(message(c).find("1 to 32 transition energies"), std::string::npos);
  c.md_events.transition_energies = {1000.0 * electron_volt};
  EXPECT_EQ(message(c), "");
  auto bad = c;
  bad.md_events.file = "e.txt";
  EXPECT_NE(message(bad).find(".csv or .parquet"), std::string::npos);
  bad = c;
  bad.damage = DamageCalculation::kinchin_pease;
  EXPECT_NE(message(bad).find("full-cascade"), std::string::npos);
  bad = c;
  bad.transport = TrimConfig::Transport::refined; // (every transport writes the files)
  EXPECT_EQ(message(bad), "");
  bad = c;
  bad.md_events.transition_energies = {2000.0 * electron_volt, 1000.0 * electron_volt};
  EXPECT_NE(message(bad).find("increasing order"), std::string::npos);
  bad.md_events.transition_energies = {1000.0 * electron_volt, 2000.0 * electron_volt};
  bad.md_events.handoff = true;
  EXPECT_NE(message(bad).find("exactly one E_c"), std::string::npos);
  bad = c;
  bad.md_events.min_energy = 5000.0 * electron_volt;
  EXPECT_NE(message(bad).find("below the smallest E_c"), std::string::npos);
  bad = c;
  bad.trajectory_export.file = "t.parquet";
  EXPECT_NE(message(bad).find("output.trajectories"), std::string::npos);
  bad.trajectory_ions = TrimConfig::all_ions;
  EXPECT_EQ(message(bad), "");
}

// TRIM.DAT recoils (cascade code 6): a PKA below E_c is one event; above it
// the cascade's pieces are.
TEST(TrimMdEvents, RecoilSourceEvents) {
  auto c = md_case(30.0, 3); // (ion.energy limits the TRIM.DAT energies)
  c.source = TrimConfig::Source::particle_recoils;
  for (int i = 0; i < 3; ++i) {
    SourceParticle sp;
    sp.z = 26;
    sp.energy = (i == 0 ? 500.0 : 20000.0) * electron_volt;
    sp.x = 1500.0 * angstrom;
    sp.cos_x = 1.0;
    c.particles.push_back(sp);
  }
  c.ion_count = 3;
  c.md_events.file = "e.csv";
  c.output_directory = temp_dir("md_source").string();
  c.outputs.formats.clear();
  c.record_ion_positions = false;
  c.md_events.transition_energies = {1000.0 * electron_volt};
  TrimSimulation s(c);
  s.run();
  const std::string ev = read_file(std::filesystem::path(c.output_directory) / "e.csv");
  // The 500 eV PKA: exactly one event (the PKA itself, generation 1, no parent).
  std::istringstream in(ev);
  std::string line;
  std::getline(in, line);
  std::size_t ion1 = 0;
  bool root = false;
  while (std::getline(in, line)) {
    if (line.rfind("1,", 0) == 0) {
      ++ion1;
      root = line.find(",1,26,") != std::string::npos && line.size() > 3 &&
             line.compare(line.size() - 6, 6, ",-1,-1") == 0;
    }
  }
  EXPECT_EQ(ion1, 1u);
  EXPECT_TRUE(root) << ev.substr(0, 400);
  EXPECT_GT(s.tallies().md[0].events / 3.0, 5.0);
}

TEST(TrimMdEvents, SummaryTableAndJson) {
  auto c = md_case(10.0, 8);
  c.md_events.file = "e.parquet";
  c.output_directory = temp_dir("md_summary").string();
  c.outputs.formats.clear();
  c.record_ion_positions = false;
  c.md_events.transition_energies = {500.0 * electron_volt, 1000.0 * electron_volt};
  TrimSimulation s(c);
  s.run();
  const auto d = make_dataset(s);
  const auto *t = d.find("md_events_summary");
  ASSERT_NE(t, nullptr);
  EXPECT_EQ(t->rows(), 2u);
  EXPECT_EQ(t->find("E_c")->values[1], 1000.0);
  EXPECT_GT(t->find("events")->values[0], t->find("events")->values[1]);
  const std::string j = s.summary_json();
  EXPECT_NE(j.find("\"md_events\": {\"mode\": \"record\""), std::string::npos);
  EXPECT_TRUE(std::filesystem::exists(std::filesystem::path(c.output_directory) / "e.parquet"));
}

// ---------------------------------------------------------------------------
// Damage models (dpa table)
// ---------------------------------------------------------------------------

TEST(TrimDamageModels, PublishedConstantsAreNordlund2018TableOne) {
  const auto fe = published_damage_efficiency(26);
  ASSERT_TRUE(fe.has_value());
  EXPECT_EQ(fe->arc_b, -0.568);
  EXPECT_EQ(fe->arc_c, 0.286);
  EXPECT_NEAR(in_units(fe->rpa_b, electron_volt), 1018.0, 1e-9);
  EXPECT_EQ(fe->rpa_c, 0.95);
  EXPECT_NEAR(in_units(fe->fitted_ed, electron_volt), 40.0, 1e-9);
  const auto w = published_damage_efficiency(74);
  ASSERT_TRUE(w.has_value());
  EXPECT_NEAR(in_units(w->fitted_ed, electron_volt), 70.0, 1e-9);
  EXPECT_EQ(w->arc_c, 0.12);
  for (unsigned z : {29u, 28u, 46u, 78u}) EXPECT_TRUE(published_damage_efficiency(z).has_value()) << z;
  EXPECT_FALSE(published_damage_efficiency(14).has_value());
  EXPECT_FALSE(published_damage_efficiency(79).has_value());
}

TEST(TrimDamageModels, OnePrimaryFollowsTheFormulas) {
  // One 10 keV Fe recoil in Fe (TRIM.DAT recoil source, quick calculation):
  // exactly one primary, whose damage energy is Robinson's partition.
  auto c = fe_on_fe(50.0, DamageCalculation::kinchin_pease, 0, 5000.0);
  c.source = TrimConfig::Source::particle_recoils;
  auto list = sample_particles(26, 1, 1e4, 0.0);
  list.particles[0].energy = 10.0 * kilo_electron_volt;
  list.particles[0].x = 1000.0 * angstrom;
  c.particles = list.particles;
  TrimSimulation s(c);
  s.run();
  const auto &t = s.tallies();
  const double td = t.total_pka_damage_energy;
  EXPECT_LT(rel(td, in_units(robinson_damage_energy(26, 55.847 * atomic_mass_unit, 10.0 * kilo_electron_volt),
                             electron_volt)),
            2e-3);
  // Hand calculation with Nordlund et al. 2018 Eqs. 1, 7 and 9, Fe constants.
  const double ed = 40.0, knee = 2.0 * ed / 0.8;
  const double nrt = 0.8 * td / (2.0 * ed);
  const double xi_arc = (1.0 - 0.286) / std::pow(knee, -0.568) * std::pow(td, -0.568) + 0.286;
  const double bc = std::pow(1018.0, 0.95);
  const double xi_rpa = (bc / std::pow(knee, 0.95) + 1.0) * std::pow(td, 0.95) / (bc + std::pow(td, 0.95));
  EXPECT_NEAR(t.total_nrt, nrt, 1e-9 * nrt);
  EXPECT_NEAR(t.total_kp, td / (2.0 * ed), 1e-9 * nrt);
  EXPECT_NEAR(t.total_arc, nrt * xi_arc, 1e-9 * nrt);
  EXPECT_NEAR(t.total_rpa, nrt * xi_rpa, 1e-9 * nrt);
  EXPECT_EQ(t.pka_without_efficiency, 0.0);
  // The quick calculation's vacancies are NRT with TRIM's damage energy.
  EXPECT_NEAR(t.total_vacancies, nrt, 1e-6 * nrt);
  // Booked in the bin where the primary starts.
  const auto &h = t.nrt_displacements.values;
  const auto bin = static_cast<std::size_t>(1000.0 / (5000.0 / static_cast<double>(h.size())));
  EXPECT_NEAR(h[bin], nrt, 1e-9 * nrt);
}

TEST(TrimDamageModels, DpaTableScalesWithFluenceAndDensity) {
  for (auto damage : {DamageCalculation::kinchin_pease, DamageCalculation::full_cascade}) {
    auto c = fe_on_fe(50.0, damage, 40);
    c.dpa_fluence_per_cm2 = 2e14;
    TrimSimulation s(c);
    s.run();
    const auto &t = s.tallies();
    if (damage == DamageCalculation::kinchin_pease) {
      // Quick calculation: NRT over all primaries is TRIM's vacancy total.
      EXPECT_LT(rel(t.total_nrt, t.total_vacancies), 1e-9);
    }
    // Ordering of the models: arc <= NRT <= rpa, NRT < KP.
    EXPECT_LT(t.total_arc, t.total_nrt);
    EXPECT_GT(t.total_rpa, t.total_nrt);
    EXPECT_LT(t.total_nrt, t.total_kp);
    const auto d = make_dataset(s);
    const auto *tb = d.find("dpa");
    ASSERT_NE(tb, nullptr);
    for (const char *k : {"depth", "trim_vacancies", "nrt", "kinchin_pease", "arc_dpa", "rpa", "damage_energy"}) {
      EXPECT_NE(tb->find(k), nullptr) << k;
    }
    // dpa x N x bin width / fluence summed over depth = displacements per ion.
    const double n_A3 = in_units(c.layers[0].density / (55.847 * atomic_mass_unit), per_cubic_angstrom);
    const auto &dep = tb->find("depth")->values;
    const double w = dep[1] - dep[0];
    double sum_nrt = 0.0, sum_vac = 0.0, sum_td = 0.0;
    for (std::size_t i = 0; i < tb->rows(); ++i) {
      sum_nrt += tb->find("nrt")->values[i];
      sum_vac += tb->find("trim_vacancies")->values[i];
      sum_td += tb->find("damage_energy")->values[i] * w;
    }
    const double per_ion = n_A3 * w / (2e14 * 1e-16);
    const double ions = static_cast<double>(t.ions);
    EXPECT_LT(rel(sum_nrt * per_ion, t.total_nrt / ions), 1e-3) << to_string(damage);
    // trim_vacancies is TRIM's vacancy table (VACANCY.txt); in full cascades
    // TRIM's total also counts vacancies that table leaves out.
    double table_vac = 0.0;
    for (const auto &h : t.vacancies) {
      for (double v : h.values) table_vac += v;
    }
    EXPECT_LT(rel(sum_vac * per_ion, table_vac / ions), 1e-3) << to_string(damage);
    EXPECT_LT(rel(sum_td, t.total_pka_damage_energy / ions), 1e-3) << to_string(damage);
    EXPECT_NE(s.summary_json().find("\"damage_models\": {\"damage_energy_eV_per_ion\": "), std::string::npos);
  }
}

TEST(TrimDamageModels, ElementsWithoutConstantsLeaveArcDpaOut) {
  auto c = fe_on_fe(20.0, DamageCalculation::kinchin_pease, 10);
  c.layers[0].name = "Si";
  c.layers[0].density = 2.329 * gram_per_cubic_centimeter;
  c.layers[0].elements[0] = TargetElement::defaults(14);
  {
    TrimSimulation s(c);
    s.run();
    EXPECT_GT(s.tallies().pka_without_efficiency, 0.0);
    const auto d = make_dataset(s);
    const auto *tb = d.find("dpa");
    ASSERT_NE(tb, nullptr);
    EXPECT_NE(tb->find("nrt"), nullptr);
    EXPECT_EQ(tb->find("arc_dpa"), nullptr);
    EXPECT_EQ(tb->find("rpa"), nullptr);
    EXPECT_EQ(s.summary_json().find("arc_dpa_per_ion"), std::string::npos);
  }
  // Constants given by the user.
  DamageEfficiency de;
  de.arc_b = -0.6;
  de.arc_c = 0.3;
  de.rpa_b = 2000.0 * electron_volt;
  de.rpa_c = 0.9;
  c.layers[0].elements[0].damage_efficiency = de;
  TrimSimulation s(c);
  s.run();
  EXPECT_EQ(s.tallies().pka_without_efficiency, 0.0);
  EXPECT_NE(make_dataset(s).find("dpa")->find("arc_dpa"), nullptr);
  // Out-of-range constants and fluences are rejected.
  auto bad = c;
  bad.layers[0].elements[0].damage_efficiency->arc_c = 1.5;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.dpa_fluence_per_cm2 = 0.0;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Crystal transport (physics.transport: crystal; not TRIM)
// ---------------------------------------------------------------------------

namespace {

//! A layer with a built-in crystal structure of one element per sublattice.
TargetLayer crystal_layer(const std::string &structure, const std::vector<unsigned> &z, double a_A, double width_A,
                          double debye_K, std::vector<double> surface = {0, 0, 1},
                          std::vector<double> reference = {1, 0, 0}) {
  crystal::StructureSpec spec;
  spec.name = structure;
  spec.a_A = a_A;
  for (unsigned zi : z) spec.sublattices.push_back({crystal::Occupant{zi, 1.0}});
  crystal::Crystal c;
  c.cell = crystal::build_structure(spec);
  c.orientation.surface = std::move(surface);
  c.orientation.reference = std::move(reference);
  if (debye_K > 0.0) {
    for (unsigned zi : z) c.vibrations.push_back({zi, debye_K * kelvin, std::nullopt});
  }
  TargetLayer layer;
  layer.name = structure;
  layer.width = width_A * angstrom;
  layer.crystal = std::move(c);
  apply_crystal(layer);
  return layer;
}

TrimConfig crystal_config(unsigned ion_z, double energy_keV, TargetLayer layer, std::size_t ions,
                          DamageCalculation damage = DamageCalculation::kinchin_pease) {
  TrimConfig c;
  c.ion = IonBeam::defaults(ion_z, energy_keV * kilo_electron_volt);
  c.layers.push_back(std::move(layer));
  c.damage = damage;
  c.ion_count = ions;
  c.threads = 0; // every core: the results do not depend on it
  c.transport = TrimConfig::Transport::crystal;
  return c;
}

//! Silicon (diamond, a = 5.431 A, Debye temperature 645 K) with a (0 0 1) surface.
TargetLayer silicon(double width_A = 20000.0) {
  return crystal_layer("diamond", {14}, 5.431, width_A, 645.0, {0, 0, 1}, {1, 1, 0});
}

TrimTallies run_tallies(const TrimConfig &c) {
  TrimSimulation s(c);
  s.run();
  return s.tallies();
}

double per_ion(double total, const TrimTallies &t) { return total / static_cast<double>(t.ions); }

} // namespace

TEST(TrimCrystal, NamesAndValidation) {
  EXPECT_EQ(parse_transport("crystal"), TrimConfig::Transport::crystal);
  EXPECT_EQ(parse_transport("TRIM"), TrimConfig::Transport::trim2013);
  for (auto t : {TrimConfig::Transport::trim2013, TrimConfig::Transport::refined, TrimConfig::Transport::crystal}) {
    EXPECT_EQ(parse_transport(to_string(t)), t);
  }
  EXPECT_THROW(parse_transport("marlowe"), std::invalid_argument);
  for (auto s : {BeamSpot::Shape::point, BeamSpot::Shape::rectangle, BeamSpot::Shape::circle,
                 BeamSpot::Shape::gaussian}) {
    EXPECT_EQ(parse_beam_spot_shape(to_string(s)), s);
  }
  EXPECT_EQ(parse_beam_spot_shape("uniform"), BeamSpot::Shape::rectangle);
  EXPECT_THROW(parse_beam_spot_shape("ring"), std::invalid_argument);

  const auto c = crystal_config(5, 30.0, silicon(), 10);
  EXPECT_NO_THROW(c.validate());
  auto bad = c;
  bad.crystal.local_stopping_fraction = 1.5;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.crystal.local_stopping_decay = 0.0;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.crystal.max_impact_parameter = -1.0 * angstrom;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.crystal.simultaneous_distance = 9.0 * angstrom;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad = c;
  bad.ion.spot.width = 0.0 * angstrom;
  EXPECT_THROW(bad.validate(), std::invalid_argument);
  bad.ion.spot.shape = BeamSpot::Shape::point; // a point has no size
  EXPECT_NO_THROW(bad.validate());
  bad = c;
  bad.auto_depth.enabled = true; // sized in a deep box (TrimCrystal.AutomaticDepthHoldsTheChanneledIons)
  EXPECT_NO_THROW(bad.validate());
  bad = c;
  bad.source = TrimConfig::Source::particle_recoils;
  bad.particles = sample_particles(14, 3, 1e3, 100.0).particles;
  EXPECT_NO_THROW(bad.validate()); // particle sources run (TrimCrystal.RecoilStartsFromItsSite)
  bad.transport = TrimConfig::Transport::refined;
  EXPECT_THROW(bad.validate(), std::invalid_argument); // the refined transport runs a beam only

  // Warnings: the transport is not calibrated; atoms without thermal data stand still.
  const auto has = [](const TrimConfig &cfg, const char *text) {
    for (const auto &w : cfg.plausibility_warnings()) {
      if (w.find(text) != std::string::npos) return true;
    }
    return false;
  };
  EXPECT_TRUE(has(c, "not calibrated"));
  EXPECT_FALSE(has(c, "thermal vibrations"));
  auto still = c;
  still.layers[0].crystal->vibrations.clear();
  // A pure element of Gemmell's table takes its tabulated Debye temperature (and says so) ...
  EXPECT_TRUE(has(still, "no thermal vibrations are given for Si"));
  EXPECT_TRUE(has(still, "Gemmell"));
  EXPECT_TRUE(has(still, "543 K"));
  EXPECT_FALSE(has(still, "stand exactly on their sites"));
  // ... anything else stands still.
  auto gaas = crystal_config(5, 30.0, crystal_layer("zincblende", {31, 33}, 5.654, 20000.0, 0.0), 10);
  EXPECT_TRUE(has(gaas, "stand exactly on their sites"));
  still.crystal.thermal_vibrations = false; // asked for: nothing to warn about
  EXPECT_FALSE(has(still, "no thermal vibrations"));
  auto other = c;
  other.transport = TrimConfig::Transport::refined;
  EXPECT_TRUE(has(other, "is not used by this calculation"));
  EXPECT_TRUE(has(other, "refined"));
  auto none = c;
  none.layers[0].crystal.reset();
  EXPECT_TRUE(has(none, "no layer has a crystal"));
}

TEST(TrimCrystal, ResultsDoNotDependOnThreadCount) {
  auto c = crystal_config(26, 5.0, crystal_layer("bcc", {26}, 2.8665, 800.0, 470.0), 24,
                          DamageCalculation::full_cascade);
  c.ion.angle = 9.0 * amcset::common::degree;
  c.threads = 1;
  const auto a = run_tallies(c);
  c.threads = 4;
  const auto b = run_tallies(c);
  EXPECT_EQ(a.ions, 24u);
  EXPECT_GT(a.total_vacancies, 0.0);
  EXPECT_DOUBLE_EQ(a.range.s1, b.range.s1);
  EXPECT_DOUBLE_EQ(a.range.s2, b.range.s2);
  EXPECT_DOUBLE_EQ(a.total_vacancies, b.total_vacancies);
  EXPECT_DOUBLE_EQ(a.total_replacements, b.total_replacements);
  EXPECT_DOUBLE_EQ(a.total_ionization_recoils, b.total_ionization_recoils);
  EXPECT_DOUBLE_EQ(a.total_phonons_recoils, b.total_phonons_recoils);
  EXPECT_EQ(a.backscattered, b.backscattered);
  c.seed += 1;
  EXPECT_NE(run_tallies(c).range.s1, a.range.s1);
}

TEST(TrimCrystal, EnergyIsConservedAndDefectsAddUp) {
  for (auto damage : {DamageCalculation::kinchin_pease, DamageCalculation::full_cascade}) {
    auto c = crystal_config(26, 20.0, crystal_layer("bcc", {26}, 2.8665, 3000.0, 470.0), 40, damage);
    c.ion.angle = 7.0 * amcset::common::degree;
    TrimSimulation s(c);
    s.run();
    EXPECT_NEAR(energy_balance(s, 20000.0), 1.0, 0.01) << to_string(damage);
    const auto &t = s.tallies();
    EXPECT_GT(t.total_vacancies, 100.0 * 40) << to_string(damage); // some hundred per 20 keV ion
    if (damage == DamageCalculation::full_cascade) {
      // Every displaced atom leaves a vacancy or is replaced; none is lost.
      EXPECT_DOUBLE_EQ(t.total_displacements, t.total_vacancies + t.total_replacements);
      EXPECT_GT(t.total_replacements, 0.0);
      EXPECT_GT(t.recoils_followed, 0.0);
      double sputtered = 0.0, transmitted = 0.0;
      for (double v : t.sputtered) sputtered += v;
      for (double v : t.transmitted_recoils) transmitted += v;
      // Displaced atoms end as interstitials, on another atom's site, or outside.
      EXPECT_NEAR(t.total_interstitials + t.total_replacements + sputtered + transmitted, t.total_displacements,
                  0.02 * t.total_displacements);
    }
  }
}

TEST(TrimCrystal, RandomDirectionIsTheAmorphousTarget) {
  // 30 keV B into Si, 7 degrees off <001> with a 22 amcset::common::degree twist (the usual
  // way to avoid channeling): ranges and damage are the amorphous target's.
  auto c = crystal_config(5, 30.0, silicon(), 500);
  c.ion.angle = 7.0 * amcset::common::degree;
  c.layers[0].crystal->orientation.twist = 22.0 * amcset::common::degree;
  const auto tilted = run_tallies(c);
  auto amorphous = c;
  amorphous.transport = TrimConfig::Transport::trim2013; // ignores the crystal
  const auto ref = run_tallies(amorphous);
  EXPECT_LT(rel(tilted.range.mean(), ref.range.mean()), 0.10);
  EXPECT_LT(rel(tilted.total_vacancies, ref.total_vacancies), 0.08);
  // How the electronic loss is divided between the path and the collisions
  // does not matter in a random direction (the division is normalised).
  c.crystal.local_stopping_fraction = 0.0;
  const auto path_only = run_tallies(c);
  c.crystal.local_stopping_fraction = 1.0;
  const auto local_only = run_tallies(c);
  // (For a fast ion on a straight path the two agree to 0.2 %; for this slow
  // one, whose path bends in the close collisions, to a few per cent.)
  EXPECT_LT(rel(path_only.range.mean(), local_only.range.mean()), 0.10);
  EXPECT_LT(rel(path_only.total_ionization_ion, local_only.total_ionization_ion), 0.08);
  // Nor does a larger maximum impact parameter (the far collisions are weak).
  c.crystal.local_stopping_fraction = 0.5;
  c.crystal.max_impact_parameter = 3.0 * angstrom;
  EXPECT_LT(rel(run_tallies(c).range.mean(), tilted.range.mean()), 0.10);
}

TEST(TrimCrystal, ChannelsAreDeeperAndCleaner) {
  // Along <001> the ions are steered between the rows: deeper, straighter,
  // with fewer displacements than 7 degrees off; <110>, the widest channel
  // of the diamond structure, is deeper still.
  auto c = crystal_config(5, 30.0, silicon(40000.0), 300);
  const auto axial = run_tallies(c);
  auto off = c;
  off.ion.angle = 7.0 * amcset::common::degree;
  off.layers[0].crystal->orientation.twist = 22.0 * amcset::common::degree;
  const auto tilted = run_tallies(off);
  EXPECT_GT(axial.range.mean(), 1.8 * tilted.range.mean());
  EXPECT_LT(axial.lateral_abs_y.mean(), 0.6 * tilted.lateral_abs_y.mean());
  EXPECT_LT(axial.total_vacancies, 0.5 * tilted.total_vacancies);
  auto c110 = crystal_config(5, 30.0, crystal_layer("diamond", {14}, 5.431, 40000.0, 645.0, {1, 1, 0}, {0, 0, 1}), 300);
  EXPECT_GT(run_tallies(c110).range.mean(), 1.2 * axial.range.mean());
  // A thin amorphous layer on top scatters the beam out of the channel.
  auto oxide = c;
  TargetLayer sio2;
  sio2.name = "SiO2";
  sio2.width = 100.0 * angstrom;
  sio2.density = 2.2 * gram_per_cubic_centimeter;
  sio2.elements = {TargetElement::defaults(14, 1.0), TargetElement::defaults(8, 2.0)};
  oxide.layers.insert(oxide.layers.begin(), sio2);
  const auto screened = run_tallies(oxide);
  EXPECT_LT(screened.range.mean(), 0.8 * axial.range.mean());
  EXPECT_GT(screened.range.mean(), 1.1 * tilted.range.mean());
}

// The automatic target depth with the crystal transport: the box holds the
// channeled ions (far beyond the amorphous estimate), the deepest of them in
// its last margin, from a far too thin and a far too thick start alike.
TEST(TrimCrystal, AutomaticDepthHoldsTheChanneledIons) {
  std::vector<double> sized;
  for (const double start_A : {50.0, 1e7}) {
    auto c = crystal_config(5, 30.0, silicon(start_A), 10, DamageCalculation::full_cascade);
    c.auto_depth.enabled = true;
    c.auto_depth.test_ions = 200;
    const auto r = size_target(c);
    EXPECT_TRUE(r.changed);
    EXPECT_EQ(r.transmitted, 0u);
    EXPECT_LE(r.runs, 3u) << joined(r.notes); // the deep box, the trimmed one
    const double total = in_units(c.target_thickness(), angstrom);
    const double deepest = in_units(r.deepest, angstrom);
    EXPECT_LE(deepest, total);
    EXPECT_GE(deepest, c.auto_depth.shrink_below * total * (1.0 - 1e-9)) << joined(r.notes);
    EXPECT_TRUE(clean_thickness(total)) << total << " A";
    // Channeled along <001>: deeper than the estimate that knows no lattice.
    EXPECT_GT(deepest, 1.5 * in_units(estimate_deepest(c), angstrom)) << joined(r.notes);
    EXPECT_EQ(c.damage, DamageCalculation::full_cascade);
    EXPECT_EQ(c.ion_count, 10u);
    EXPECT_NE(joined(r.notes).find("crystal transport"), std::string::npos) << joined(r.notes);
    EXPECT_NO_THROW(c.validate());
    // The calculation itself, in the sized box: the same ions, none transmitted.
    auto run = c;
    run.auto_depth.enabled = false;
    run.ion_count = 200;
    run.damage = DamageCalculation::kinchin_pease;
    EXPECT_EQ(run_tallies(run).transmitted, 0u);
    sized.push_back(total);
  }
  EXPECT_EQ(sized[0], sized[1]);
}

// A recoil atom (PKA) that starts inside a crystal: it starts from the
// nearest lattice site of its element (not the given point), sets off a
// cascade, every run is booked as an event, and the results do not depend
// on the thread count. TRIM.DAT ions start at their points as well.
TEST(TrimCrystal, RecoilStartsFromItsSite) {
  auto c = crystal_config(26, 20.0, crystal_layer("bcc", {26}, 2.8665, 3000.0, 420.0), 0,
                          DamageCalculation::full_cascade);
  c.layers[0].elements[0].displacement_energy = 40.0 * electron_volt;
  SourceParticle p;
  p.name = "PKA";
  p.z = 26;
  p.energy = 5000.0 * electron_volt;
  p.x = 1000.0 * angstrom;
  p.y = 0.3 * angstrom;
  p.z_position = 0.2 * angstrom;
  p.cos_x = 0.6;
  p.cos_y = 0.8;
  p.cos_z = 0.0;
  c.source = TrimConfig::Source::particle_recoils;
  c.particles = {p};
  c.particle_repeat = 12;
  c.trajectory_ions = TrimConfig::all_ions;
  c.threads = 1;
  TrimSimulation one(c);
  EXPECT_EQ(one.total_count(), 12u);
  one.run();
  const auto &t = one.tallies();
  EXPECT_EQ(t.ions, 12u);
  EXPECT_EQ(t.range.n, 12.0);                 // booked where it starts (as TRIM's code 6)
  EXPECT_EQ(t.total_ionization_ion, 0.0);     // no ion
  EXPECT_GT(t.total_displacements, 12.0 * 5.0); // the PKA and its cascade
  std::vector<common::CollisionPoint> pts;
  std::size_t runs = 0;
  std::set<double> starts;
  while (one.try_pop_trajectory(pts)) {
    ++runs;
    ASSERT_FALSE(pts.empty());
    const auto &b = pts.front();
    EXPECT_EQ(b.kind, common::CollisionPoint::Kind::birth);
    EXPECT_EQ(b.parent_track, common::CollisionPoint::no_parent);
    const double dx = in_units(b.coordinate.z_, angstrom) - 1000.0; // (output axes: z = depth)
    const double dy = in_units(b.coordinate.x_, angstrom) - 0.3;
    const double dz = in_units(b.coordinate.y_, angstrom) - 0.2;
    const double d = std::sqrt(dx * dx + dy * dy + dz * dz);
    EXPECT_GT(d, 1e-6) << "the PKA did not move to a lattice site";
    EXPECT_LT(d, 2.0) << "the nearest site of bcc Fe is within 2 A of any point";
    starts.insert(in_units(b.coordinate.z_, angstrom));
    EXPECT_GT(pts.size(), 3u);
  }
  EXPECT_EQ(runs, 12u);
  EXPECT_GT(starts.size(), 1u); // thermal displacements differ between histories
  auto four = c;
  four.threads = 4;
  four.trajectory_ions = 0;
  const auto t4 = run_tallies(four);
  EXPECT_EQ(t4.total_displacements, t.total_displacements);
  EXPECT_EQ(t4.total_ionization_recoils, t.total_ionization_recoils);

  // TRIM.DAT ions in the crystal transport.
  auto ions = c;
  ions.source = TrimConfig::Source::particle_ions;
  ions.particles = sample_particles(26, 6, 2e4, 1000.0).particles;
  ions.particle_repeat = 1;
  ions.trajectory_ions = 0;
  const auto ti = run_tallies(ions);
  EXPECT_EQ(ti.ions, 6u);
  EXPECT_GT(ti.total_ionization_ion, 0.0);
}

// physics.crystal.auto_crystal: an amorphous layer of a database crystal's
// composition runs as that crystal (exactly as if it had been given it),
// says so, and nothing changes with the switch off or another transport.
TEST(TrimCrystal, AutoCrystalFromTheDatabase) {
  auto amorphous = md_case(20.0, 20); // amorphous Fe
  amorphous.transport = TrimConfig::Transport::crystal;
  amorphous.damage = DamageCalculation::kinchin_pease;
  auto automatic = amorphous;
  automatic.crystal.auto_crystal = true;
  TrimSimulation a(automatic);
  bool said = false;
  for (const auto &w : a.warnings()) said = said || w.find("followed as the crystal Fe (bcc") != std::string::npos;
  EXPECT_TRUE(said);
  ASSERT_TRUE(a.config().layers[0].crystal.has_value());
  auto given = amorphous;
  given.layers[0].crystal = crystal::make_database_crystal(crystal::find_database_crystal("Fe"));
  apply_crystal(given.layers[0]);
  a.run();
  const auto t_given = run_tallies(given);
  EXPECT_EQ(a.tallies().range.s1, t_given.range.s1);
  EXPECT_EQ(a.tallies().total_vacancies, t_given.total_vacancies);
  // Off, or with another transport: the layer stays amorphous.
  EXPECT_FALSE(TrimSimulation(amorphous).config().layers[0].crystal.has_value());
  auto trim = automatic;
  trim.transport = TrimConfig::Transport::trim2013;
  EXPECT_FALSE(TrimSimulation(trim).config().layers[0].crystal.has_value());
  // An alloy matches nothing.
  auto alloy = automatic;
  alloy.layers[0].elements.push_back(TargetElement::defaults(24, 0.2));
  EXPECT_FALSE(TrimSimulation(alloy).config().layers[0].crystal.has_value());
}

TEST(TrimCrystal, ThermalVibrationsDechannel) {
  // The hotter the crystal, the larger the vibrations and the sooner the
  // ions leave the <001> channel; a static lattice channels best.
  auto c = crystal_config(5, 30.0, silicon(40000.0), 300);
  const auto room = run_tallies(c);
  c.layers[0].crystal->temperature = 1200.0 * kelvin;
  const auto hot = run_tallies(c);
  c.crystal.thermal_vibrations = false;
  const auto frozen = run_tallies(c);
  EXPECT_GT(room.range.mean(), 1.15 * hot.range.mean());
  EXPECT_GT(frozen.range.mean(), 1.15 * room.range.mean());
  EXPECT_LT(frozen.total_vacancies, room.total_vacancies);
  EXPECT_LT(room.total_vacancies, hot.total_vacancies);
}

TEST(TrimCrystal, AxialDipOfHeliumInTungsten) {
  // 1 MeV He through 600 A of W along <001>. Atoms displaced by the ion
  // need impact parameters below about 0.03 A (T > Ed), less than the
  // thermal amplitude, so their number follows the close-encounter
  // probability of a backscattering experiment. Lindhard's critical angle:
  // psi_1 = sqrt(2 Z1 Z2 e^2 / (E d)) = sqrt(2 * 2 * 74 * 14.4 eV A /
  // (1e6 eV * 3.165 A)) = 0.0367 rad = 2.10 degrees.
  const auto yield = [](double tilt_deg, bool crystal_transport) {
    auto layer = crystal_layer("bcc", {74}, 3.165, 600.0, 310.0);
    layer.crystal->orientation.twist = 17.0 * amcset::common::degree; // tilt away from the low-index planes
    auto c = crystal_config(2, 1000.0, layer, 1200);
    c.ion.angle = tilt_deg * amcset::common::degree;
    if (!crystal_transport) c.transport = TrimConfig::Transport::trim2013;
    const auto t = run_tallies(c);
    return per_ion(t.total_displacements, t);
  };
  const double random = yield(0.0, false); // the amorphous target: about 0.28 knock-ons per ion
  EXPECT_GT(random, 0.2);
  EXPECT_LT(yield(0.0, true), 0.05 * random);  // on the axis: a few per cent
  EXPECT_LT(yield(1.0, true), 0.30 * random);  // inside the dip
  EXPECT_GT(yield(2.1, true), 0.60 * random);  // at psi_1: out of it
  EXPECT_GT(yield(3.0, true), 0.90 * random);  // the shoulder
}

TEST(TrimCrystal, BeamSpotDecidesWhereTheLatticeIsMet) {
  // A static lattice, no straggling, every ion at the same point: one
  // trajectory, repeated. (75As is heavier than Si: never turned back in
  // one collision.) Spread over the surface, the ions differ.
  auto c = crystal_config(33, 50.0, silicon(), 12);
  c.crystal.thermal_vibrations = false;
  c.electronic_straggling = false;
  c.ion.spot.shape = BeamSpot::Shape::point;
  const auto point = run_tallies(c);
  ASSERT_EQ(point.range.n, 12.0);
  EXPECT_NEAR(point.range.stddev(), 0.0, 1e-3); // (the moments' rounding)
  for (auto shape : {BeamSpot::Shape::rectangle, BeamSpot::Shape::circle, BeamSpot::Shape::gaussian}) {
    c.ion.spot.shape = shape;
    c.ion.spot.width = 200.0 * angstrom;
    c.ion.spot.height = 100.0 * angstrom;
    const auto spread = run_tallies(c);
    EXPECT_GT(spread.range.stddev(), 10.0) << to_string(shape);
    // Positions stay relative to each ion's entry point: no 100 A offsets.
    EXPECT_LT(std::abs(spread.lateral_y.mean()), 60.0) << to_string(shape);
  }
}

TEST(TrimCrystal, AmorphousLayersTakeRandomPartners) {
  // Without a crystal the transport is a monolayer-step amorphous transport:
  // TRIM's ranges and damage within the statistics.
  auto c = crystal_config(5, 30.0, silicon(), 600);
  c.layers[0].crystal.reset();
  const auto t = run_tallies(c);
  auto trim = c;
  trim.transport = TrimConfig::Transport::trim2013;
  const auto ref = run_tallies(trim);
  EXPECT_LT(rel(t.range.mean(), ref.range.mean()), 0.06);
  EXPECT_LT(rel(t.total_vacancies, ref.total_vacancies), 0.06);
  EXPECT_LT(rel(t.total_ionization_ion, ref.total_ionization_ion), 0.03);
}

TEST(TrimCrystal, HeavySlowIonsKeepTheAmorphousRange) {
  // 30 keV W into amorphous W: large-angle collisions, whose deflection
  // point lies x1 = p tan(theta/2) before the closest approach. One partner
  // per monolayer of path and the electronic loss charged once: the range
  // is TRIM's. (Charging the last x1 again and starting the next monolayer
  // at the deflection point made it 4.9 % short: 50.3 against 52.9 A with
  // 20000 ions.) 6000 ions: the ratio's standard error is about 1 %.
  auto layer = crystal_layer("bcc", {74}, 3.165, 2000.0, 310.0);
  layer.crystal.reset();
  auto c = crystal_config(74, 30.0, layer, 6000);
  const auto t = run_tallies(c);
  c.transport = TrimConfig::Transport::trim2013;
  const auto ref = run_tallies(c);
  EXPECT_NEAR(ref.range.mean(), 52.9, 1.5);
  EXPECT_LT(rel(t.range.mean(), ref.range.mean()), 0.025);
}

TEST(TrimCrystal, AlloysInterstitialAtomsAndSputtering) {
  // Fe-Cr on the bcc sites with carbon on a fifth of the octahedral sites:
  // atoms of all three elements are struck and displaced.
  crystal::StructureSpec spec;
  spec.name = "bcc";
  spec.a_A = 2.87;
  spec.sublattices = {{crystal::Occupant{26, 0.5}, crystal::Occupant{24, 0.5}}};
  spec.interstitials = {{"octahedral", {crystal::Occupant{6, 0.2}}}};
  TargetLayer layer;
  layer.width = 1500.0 * angstrom;
  layer.crystal.emplace();
  layer.crystal->cell = crystal::build_structure(spec);
  apply_crystal(layer);
  auto c = crystal_config(26, 10.0, layer, 40, DamageCalculation::full_cascade);
  c.ion.angle = 12.0 * amcset::common::degree;
  TrimSimulation s(c);
  s.run();
  const auto &t = s.tallies();
  ASSERT_EQ(t.vacancies.size(), 3u);
  std::vector<double> vacancies;
  for (const auto &h : t.vacancies) {
    double sum = 0.0;
    for (double v : h.values) sum += v;
    vacancies.push_back(sum);
    EXPECT_GT(sum, 0.0);
  }
  // Fe and Cr share the sites equally and have nearly the same mass.
  const auto names = s.element_names();
  ASSERT_EQ(names.size(), 3u);
  std::map<std::string, double> by_name;
  for (std::size_t i = 0; i < names.size(); ++i) by_name[names[i]] = vacancies[i];
  EXPECT_LT(rel(by_name.at("Fe"), by_name.at("Cr")), 0.25);
  EXPECT_LT(by_name.at("C"), by_name.at("Fe")); // 1.2 C per cell against 1 Fe, but small and light
  EXPECT_NEAR(energy_balance(s, 10000.0), 1.0, 0.01);

  // 1 keV Ar onto Cu (0 0 1): atoms leave through the surface barrier.
  auto cu = crystal_config(18, 1.0, crystal_layer("fcc", {29}, 3.615, 300.0, 343.0), 150, DamageCalculation::full_cascade);
  cu.ion.angle = 15.0 * amcset::common::degree;
  TrimSimulation sp(cu);
  sp.run();
  double sputtered = 0.0;
  for (double v : sp.tallies().sputtered) sputtered += v;
  EXPECT_GT(sputtered / 150.0, 0.3);
  EXPECT_LT(sputtered / 150.0, 10.0);
  EXPECT_NEAR(energy_balance(sp, 1000.0), 1.0, 0.01);
}

TEST(TrimCrystal, CutOffAtEdGivesTheNrtCount) {
  // The displacement model of Robinson and Torrens (Sec. V): a projectile
  // that keeps less than min(E_c, Ed) after displacing an atom of its own
  // element takes its site. With the cut-off energy E_c at Ed this is
  // Kinchin and Pease's model, and the vacancy count is the NRT standard's
  // (which was fitted to such calculations); with TRIM's 1 eV cut-off far
  // fewer collisions are replacements and the count is TRIM's full-cascade
  // one, well above NRT.
  auto layer = silicon();
  layer.elements[0].displacement_energy = 15.0 * electron_volt;
  layer.elements[0].lattice_binding_energy = 0.0 * electron_volt;
  auto c = crystal_config(5, 30.0, layer, 300, DamageCalculation::full_cascade);
  c.ion.angle = 7.0 * amcset::common::degree;
  c.layers[0].crystal->orientation.twist = 22.0 * amcset::common::degree;
  const auto trim_like = run_tallies(c);
  c.ion_stop_energy = 15.0 * electron_volt;
  c.recoil_stop_energy = 15.0 * electron_volt;
  const auto kp_like = run_tallies(c);
  auto quick = c;
  quick.transport = TrimConfig::Transport::trim2013;
  quick.damage = DamageCalculation::kinchin_pease; // TRIM's quick calculation: NRT
  quick.ion_stop_energy = 0.0 * electron_volt;
  quick.recoil_stop_energy = 0.0 * electron_volt;
  const auto nrt = run_tallies(quick);
  // (11 % above NRT with 2000 ions, tools/crystal/channeling_checks.py.)
  EXPECT_LT(rel(kp_like.total_vacancies, nrt.total_vacancies), 0.20);
  EXPECT_GT(trim_like.total_vacancies, 1.25 * kp_like.total_vacancies);
  EXPECT_GT(kp_like.total_replacements, 2.0 * trim_like.total_replacements);
  // Either way every displaced atom is accounted for.
  EXPECT_DOUBLE_EQ(kp_like.total_displacements, kp_like.total_vacancies + kp_like.total_replacements);
  EXPECT_NEAR(kp_like.total_displacements, trim_like.total_displacements, 0.1 * trim_like.total_displacements);
}

TEST(TrimCrystal, TakesNoneOfTrimsBugs) {
  // The crystal transport is not TRIM: whatever the bug switches say, a
  // crystal run takes none of TRIM.exe's bugs. The scenario keeps its
  // switches for the other transports.
  auto c = crystal_config(5, 30.0, silicon(), 400, DamageCalculation::full_cascade);
  c.ion.angle = 7.0 * amcset::common::degree;
  ASSERT_TRUE(c.bugs.all_reproduced()); // the default: as TRIM
  ASSERT_TRUE(c.single_precision_tables);
  for (const auto &m : TrimBugs::members()) EXPECT_FALSE(c.effective_bugs().*(m.second)) << m.first;
  EXPECT_FALSE(c.effective_single_precision_tables());
  for (auto t : {TrimConfig::Transport::trim2013, TrimConfig::Transport::refined}) {
    auto other = c;
    other.transport = t;
    EXPECT_TRUE(other.effective_bugs().all_reproduced()) << to_string(t);
    EXPECT_TRUE(other.effective_single_precision_tables()) << to_string(t);
    other.bugs.struck_atom_direction = false;
    other.single_precision_tables = false;
    EXPECT_FALSE(other.effective_bugs().struck_atom_direction) << to_string(t);
    EXPECT_TRUE(other.effective_bugs().stale_primary_energy) << to_string(t);
    EXPECT_FALSE(other.effective_single_precision_tables()) << to_string(t);
  }
  // The engine runs with the effective switches ...
  const TrimSimulation as_given(c);
  EXPECT_FALSE(as_given.config().single_precision_tables);
  for (const auto &m : TrimBugs::members()) EXPECT_FALSE(as_given.config().bugs.*(m.second)) << m.first;
  // ... and the results are those of a scenario with every bug fixed, number
  // for number.
  auto fixed = c;
  fixed.bugs.set_all(false);
  fixed.single_precision_tables = false;
  const auto a = run_tallies(c);
  const auto b = run_tallies(fixed);
  EXPECT_EQ(a.ions, 400u);
  EXPECT_EQ(a.range.mean(), b.range.mean());
  EXPECT_EQ(a.total_vacancies, b.total_vacancies);
  EXPECT_EQ(a.total_replacements, b.total_replacements);
  EXPECT_EQ(a.total_interstitials, b.total_interstitials);
  EXPECT_EQ(a.total_ionization_ion, b.total_ionization_ion);
  EXPECT_EQ(a.total_phonons_recoils, b.total_phonons_recoils);
}

TEST(TrimCrystal, PolycrystalGrains) {
  // 30 keV B along the surface normal of (0 0 1) Si, a channel of the
  // single crystal, into grains of that silicon.
  const std::size_t ions = 1500;
  const auto single = run_tallies(crystal_config(5, 30.0, silicon(), ions));
  auto amorphous = crystal_config(5, 30.0, silicon(), ions);
  amorphous.layers[0].crystal.reset();
  const double channeled = single.range.mean();
  const double random = run_tallies(amorphous).range.mean();
  ASSERT_GT(channeled, 2.0 * random); // (2838 against 1120 A)
  const auto grains = [&](double size_A, crystal::Polycrystal::Texture texture, double width_A) {
    auto c = crystal_config(5, 30.0, silicon(), ions);
    crystal::Polycrystal p;
    p.grain_size = size_A * angstrom;
    p.texture = texture;
    p.boundary_width = width_A * angstrom;
    c.layers[0].crystal->polycrystal = p;
    // The beam over one period of the grains.
    c.ion.spot.width = c.ion.spot.height = 8.0 * size_A * angstrom;
    return c;
  };
  // Large grains of the crystal's own orientation (a mosaic crystal without
  // spread): the ions still channel. The grains' lattices are not in
  // registry, so each boundary dechannels some of them (about 15 % of the
  // range here).
  const auto mosaic = run_tallies(grains(4000.0, crystal::Polycrystal::Texture::mosaic, 0.0));
  EXPECT_GT(mosaic.range.mean(), 0.7 * channeled);
  EXPECT_LT(mosaic.range.mean(), 1.05 * channeled);
  // Small grains of random orientations: most are not channels. The ions go
  // little deeper than in the amorphous target.
  const auto small = run_tallies(grains(100.0, crystal::Polycrystal::Texture::random, 0.0));
  EXPECT_LT(small.range.mean(), 0.65 * channeled);
  EXPECT_GT(small.range.mean(), 0.95 * random);
  // A fiber texture keeps the channel along the normal in every grain.
  const auto fiber = run_tallies(grains(4000.0, crystal::Polycrystal::Texture::fiber, 0.0));
  EXPECT_GT(fiber.range.mean(), 0.8 * channeled);
  // Wide disordered boundaries around small grains leave little crystal:
  // grains of 60 A with 50 A of boundary are nearly the amorphous target.
  const auto ordered = run_tallies(grains(60.0, crystal::Polycrystal::Texture::mosaic, 0.0));
  const auto disordered = run_tallies(grains(60.0, crystal::Polycrystal::Texture::mosaic, 50.0));
  EXPECT_LT(disordered.range.mean(), ordered.range.mean());
  EXPECT_LT(disordered.range.mean(), 1.3 * random);

  // The grains do not depend on the number of threads, nor does the run;
  // every displaced atom is accounted for and the energy is conserved.
  auto c = grains(100.0, crystal::Polycrystal::Texture::random, 6.0);
  c.ion_count = 40;
  c.damage = DamageCalculation::full_cascade;
  c.threads = 1;
  TrimSimulation one(c);
  one.run();
  c.threads = 4;
  const auto four = run_tallies(c);
  const auto &t = one.tallies();
  EXPECT_DOUBLE_EQ(t.range.s1, four.range.s1);
  EXPECT_DOUBLE_EQ(t.total_vacancies, four.total_vacancies);
  EXPECT_DOUBLE_EQ(t.total_replacements, four.total_replacements);
  EXPECT_DOUBLE_EQ(t.total_phonons_recoils, four.total_phonons_recoils);
  EXPECT_GT(t.total_vacancies, 100.0 * 40);
  EXPECT_DOUBLE_EQ(t.total_displacements, t.total_vacancies + t.total_replacements);
  EXPECT_NEAR(energy_balance(one, 30000.0), 1.0, 0.01);
  // Another set of grains, another result.
  c.layers[0].crystal->polycrystal->seed = 2;
  EXPECT_NE(run_tallies(c).range.s1, four.range.s1);

  // A polycrystal whose file has not been read cannot run.
  auto unread = grains(100.0, crystal::Polycrystal::Texture::random, 0.0);
  unread.layers[0].crystal->polycrystal->source = crystal::Polycrystal::Source::ebsd;
  unread.layers[0].crystal->polycrystal->file = "map.ang";
  EXPECT_THROW(TrimSimulation{unread}, std::invalid_argument);
  auto bad = grains(1.0, crystal::Polycrystal::Texture::random, 0.0);
  EXPECT_THROW(bad.validate(), std::invalid_argument); // grains of 1 A
}

// The trajectory and MD event files in the refined and crystal transports:
// the recorded points form a tree with times, recording changes no tally,
// and in the handoff mode what MD takes plus what the BCA books is the ion
// energy.
TEST(TrimTrajectories, EveryTransportRecordsTreesAndEvents) {
  struct Case {
    const char *name;
    TrimConfig config;
  };
  std::vector<Case> cases;
  {
    auto c = md_case(20.0, 30);
    c.transport = TrimConfig::Transport::refined;
    cases.push_back({"refined", c});
  }
  {
    auto c = crystal_config(26, 20.0, crystal_layer("bcc", {26}, 2.8665, 3000.0, 420.0), 30,
                            DamageCalculation::full_cascade);
    c.ion.angle = 7.0 * amcset::common::degree;
    c.layers[0].elements[0].displacement_energy = 40.0 * electron_volt;
    c.threads = 1;
    cases.push_back({"crystal", c});
  }
  {
    // The crystal transport through an amorphous layer (random partners).
    auto c = md_case(20.0, 30);
    c.transport = TrimConfig::Transport::crystal;
    cases.push_back({"crystal, amorphous layer", c});
  }
  for (auto &[name, plain] : cases) {
    SCOPED_TRACE(name);
    TrimSimulation reference(plain);
    reference.run();

    // Points: one birth per track, the parent's collision point recorded
    // before the child's birth and near it (the struck atom sits within
    // an impact parameter of the path), times that never decrease along a
    // track, and the energies of each collision.
    auto c = plain;
    c.trajectory_ions = TrimConfig::all_ions;
    c.md_events.file = "events.csv";
    c.md_events.transition_energies = {300.0 * electron_volt, 2000.0 * electron_volt};
    c.trajectory_export.file = "traj.csv";
    c.output_directory = temp_dir(std::string("traj_") + (name[0] == 'r' ? "r" : name[8] == ',' ? "ca" : "c")).string();
    c.outputs.formats.clear();
    TrimSimulation s(c);
    s.run();
    EXPECT_EQ(totals(reference.tallies()), totals(s.tallies())) << "recording changed a tally";
    ASSERT_EQ(s.tallies().md.size(), 2u);
    EXPECT_GT(s.tallies().md[0].events, s.tallies().md[1].events);
    EXPECT_GT(s.tallies().md[1].events, 0.0);
    std::vector<common::CollisionPoint> pts;
    std::size_t ions = 0;
    std::size_t recoils = 0;
    while (s.try_pop_trajectory(pts)) {
      ++ions;
      ASSERT_FALSE(pts.empty());
      EXPECT_EQ(pts.front().kind, common::CollisionPoint::Kind::birth);
      EXPECT_EQ(pts.front().track, 0u);
      std::map<std::uint32_t, double> last_time;
      std::set<std::uint32_t> born;
      std::set<std::uint32_t> ended;
      for (std::size_t i = 0; i < pts.size(); ++i) {
        const auto &p = pts[i];
        EXPECT_EQ(p.point_index, i);
        if (p.kind == common::CollisionPoint::Kind::birth) {
          EXPECT_TRUE(born.insert(p.track).second) << "two births of track " << p.track;
          if (p.track != 0) {
            ++recoils;
            ASSERT_LT(p.parent_point, i);
            const auto &q = pts[p.parent_point];
            EXPECT_EQ(q.track, p.parent_track);
            EXPECT_LE(q.time, p.time + 1e-21);
            const double dx = in_units(q.coordinate.x_ - p.coordinate.x_, angstrom);
            const double dy = in_units(q.coordinate.y_ - p.coordinate.y_, angstrom);
            const double dz = in_units(q.coordinate.z_ - p.coordinate.z_, angstrom);
            EXPECT_LT(std::sqrt(dx * dx + dy * dy + dz * dz), 15.0) << "track " << p.track;
            EXPECT_GT(q.energy_transferred_eV, 0.0);
          }
        } else {
          EXPECT_TRUE(born.count(p.track)) << "point before the birth of track " << p.track;
          EXPECT_FALSE(ended.count(p.track)) << "point after the end of track " << p.track;
        }
        if (p.kind == common::CollisionPoint::Kind::stop || p.kind == common::CollisionPoint::Kind::exit) {
          ended.insert(p.track);
        }
        if (last_time.count(p.track)) {
          EXPECT_GE(p.time, last_time[p.track]);
        }
        last_time[p.track] = p.time;
        if (p.kind == common::CollisionPoint::Kind::collision && p.energy_before_eV > 0.0) {
          EXPECT_NEAR(p.energy_before_eV - p.energy_transferred_eV - p.electronic_loss_eV, p.energy_eV,
                      1e-9 * p.energy_before_eV + 1e-9);
        }
      }
      EXPECT_GT(pts.back().time, 0.0);
    }
    EXPECT_EQ(ions, 30u);
    EXPECT_GT(recoils, 30u);
    const auto events = read_file(std::filesystem::path(c.output_directory) / "events.csv");
    EXPECT_GT(events.size(), 200u);
    EXPECT_GT(std::filesystem::file_size(std::filesystem::path(c.output_directory) / "traj.csv"), 1000u);

    // Handoff: what MD takes over plus what the BCA deposits and what
    // leaves the target is the ion energy; fewer atoms are followed.
    auto h = plain;
    h.md_events.file = "events.csv";
    h.md_events.transition_energies = {1000.0 * electron_volt};
    h.md_events.handoff = true;
    h.output_directory = c.output_directory;
    h.outputs.formats.clear();
    TrimSimulation handed(h);
    handed.run();
    const auto &t = handed.tallies();
    ASSERT_EQ(t.md.size(), 1u);
    const double e0 = 20000.0;
    const double md_share = t.md[0].energy_eV / (e0 * static_cast<double>(t.ions));
    EXPECT_GT(md_share, 0.05);
    EXPECT_NEAR(energy_balance(handed, e0) + md_share, energy_balance(reference, e0), 0.01);
    EXPECT_LT(t.recoils_followed, reference.tallies().recoils_followed);
  }
}

namespace {

//! The rows of a CSV file without comment lines, by column name.
std::vector<std::map<std::string, std::string>> csv_rows(const std::filesystem::path &p) {
  std::istringstream in(read_file(p));
  std::vector<std::string> header;
  std::vector<std::map<std::string, std::string>> out;
  for (std::string line; std::getline(in, line);) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> cells;
    std::istringstream ls(line);
    for (std::string cell; std::getline(ls, cell, ',');) cells.push_back(cell);
    if (header.empty()) {
      header = cells;
      continue;
    }
    std::map<std::string, std::string> row;
    for (std::size_t i = 0; i < header.size() && i < cells.size(); ++i) row[header[i]] = cells[i];
    out.push_back(std::move(row));
  }
  return out;
}

} // namespace

// MD events with several E_c in a channeling crystal (a user's report,
// 5 MeV Fe -> bcc Fe (1 1 1)): the electronic loss of the steps with no atom
// within reach (in a channel) took the ion below some E_c without a crossing
// event, so those levels kept it as the BCA's. Every level must get one
// crossing per atom that falls below it, and the levels' events must nest:
// each event of a finer E_c lies in the subtree of exactly one event of the
// next coarser E_c, and no coarser event lies inside a finer event's subtree.
TEST(TrimMdEvents, CrystalCrossesEveryLevelOnceAndLevelsNest) {
  auto c = crystal_config(26, 100.0, crystal_layer("bcc", {26}, 2.859, 10000.0, 420.0, {1, 1, 1}, {1, -1, 0}), 20,
                          DamageCalculation::full_cascade);
  c.layers[0].elements[0].displacement_energy = 40.0 * electron_volt;
  c.seed = 20261126;
  c.threads = 1;
  c.output_directory = temp_dir("md_crystal_levels").string();
  c.outputs.formats.clear();
  c.trajectory_ions = TrimConfig::all_ions;
  c.trajectory_export.file = "births.csv";
  c.trajectory_export.births_only = true;
  c.md_events.file = "events.csv";
  c.md_events.min_energy = 40.0 * electron_volt;
  std::vector<double> levels = {500.0, 1000.0, 2000.0, 5000.0, 10000.0, 20000.0, 50000.0};
  for (double e : levels) c.md_events.transition_energies.push_back(e * electron_volt);
  TrimSimulation s(c);
  s.set_trajectory_polling(false);
  s.run();

  // Lineage: (ion, track) -> (parent track, birth time in ps).
  std::map<std::pair<long, long>, std::pair<long, double>> parent;
  for (const auto &r : csv_rows(std::filesystem::path(c.output_directory) / "births.csv")) {
    parent[{std::stol(r.at("ion")), std::stol(r.at("track"))}] = {std::stol(r.at("parent_track")),
                                                                   std::stod(r.at("t_ps"))};
  }
  struct Event {
    long ion, track;
    double t_ps;
    bool crossing;
  };
  std::map<double, std::vector<Event>> events;
  for (const auto &r : csv_rows(std::filesystem::path(c.output_directory) / "events.csv")) {
    events[std::stod(r.at("E_c_eV"))].push_back(
        {std::stol(r.at("ion_id")), std::stol(r.at("track")), std::stod(r.at("t_ps")), r.at("kind") == "crossing"});
  }
  ASSERT_EQ(events.size(), levels.size());

  // One crossing per track and level; the ion crosses every level once
  // when it stops in the target (it does not leave this thick layer).
  std::size_t ion_crossings = 0;
  for (double L : levels) {
    SCOPED_TRACE(L);
    std::map<std::pair<long, long>, int> per_track;
    for (const auto &e : events[L]) {
      if (e.crossing) ++per_track[{e.ion, e.track}];
    }
    for (const auto &[key, n] : per_track) EXPECT_EQ(n, 1) << "ion " << key.first << " track " << key.second;
    for (long ion = 1; ion <= 20; ++ion) {
      EXPECT_EQ(per_track.count({ion, 0}), 1u) << "no crossing of ion " << ion;
      ion_crossings += per_track.count({ion, 0});
    }
  }
  EXPECT_EQ(ion_crossings, 20u * levels.size());

  // Nesting. An event of level L holds a point (track, t) in its subtree if
  // it is on that track at or before t, or on an ancestor track before the
  // branch towards it was born. At equal times (a collision right after the
  // event that took no flight time, or the event's own collision) the files
  // do not give the order: the branch is held if its track has no level-L
  // event of its own (it inherited the event's mark).
  const auto has_event = [&](long ion, long track, double L) {
    return std::any_of(events[L].begin(), events[L].end(), [&](const Event &e) { return e.ion == ion && e.track == track; });
  };
  const auto holders = [&](long ion, long track, double t, double L) {
    int n = 0;
    long tr = track;
    double cut = t;
    long child = -1;
    for (;;) {
      for (const auto &e : events[L]) {
        if (e.ion != ion || e.track != tr) continue;
        if (child < 0 ? e.t_ps <= cut : (e.t_ps < cut || (e.t_ps == cut && !has_event(ion, child, L)))) ++n;
      }
      const auto it = parent.find({ion, tr});
      if (tr == 0 || it == parent.end() || it->second.first < 0) return n;
      child = tr;
      tr = it->second.first;
      cut = it->second.second;
    }
  };
  for (std::size_t i = 0; i + 1 < levels.size(); ++i) {
    const double fine = levels[i];
    const double coarse = levels[i + 1];
    SCOPED_TRACE(std::to_string(fine) + " in " + std::to_string(coarse));
    for (const auto &e : events[fine]) {
      EXPECT_EQ(holders(e.ion, e.track, e.t_ps, coarse), 1)
          << "ion " << e.ion << " track " << e.track << " t " << e.t_ps << " ps";
    }
    for (const auto &e : events[coarse]) {
      // (One collision that crosses both levels gives both events at one point.)
      int same = 0;
      for (const auto &f : events[fine]) same += f.ion == e.ion && f.track == e.track && f.t_ps == e.t_ps;
      EXPECT_EQ(holders(e.ion, e.track, e.t_ps, fine), same)
          << "coarser event inside a finer one: ion " << e.ion << " track " << e.track << " t " << e.t_ps << " ps";
    }
  }
}

TEST(TrimCrystal, OrientationIsReported) {
  // Si (0 0 1) with [1 1 0] along TRIM's y: the beam at 7 degrees tilt,
  // twist 0, travels 7 degrees from [0 0 1].
  auto c = crystal_config(5, 30.0, silicon(), 30, DamageCalculation::full_cascade);
  c.ion.angle = 7.0 * amcset::common::degree;
  auto reports = crystal_reports(c);
  ASSERT_EQ(reports.size(), 1u);
  auto r = reports.front();
  EXPECT_TRUE(r.followed);
  EXPECT_EQ(r.depth_uvw, (std::array<int, 3>{0, 0, 1}));
  EXPECT_NEAR(r.depth_off_deg, 0.0, 1e-9);
  EXPECT_EQ(r.y_uvw, (std::array<int, 3>{1, 1, 0}));
  // Bunge angles that give back the rotation.
  crystal::Orientation o;
  o.kind = crystal::OrientationKind::euler;
  o.euler_deg = r.euler_deg;
  const auto again = crystal::rotation(c.layers[0].crystal->cell.lattice, o);
  for (int i = 0; i < 3; ++i) {
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(again[i][k], r.to_trim[i][k], 1e-12);
  }
  // The output axes are TRIM's (y, z, x).
  for (int k = 0; k < 3; ++k) {
    EXPECT_EQ(r.to_output[0][k], r.to_trim[1][k]);
    EXPECT_EQ(r.to_output[2][k], r.to_trim[0][k]);
  }
  // 7 degrees from [0 0 1]; Lindhard's angle along <001> of Si (spacing
  // a = 5.431 A) for 30 keV B: 0.307 sqrt(5 * 14 / (5.431 * 0.030)) degrees.
  ASSERT_FALSE(r.axes.empty());
  EXPECT_EQ(r.axes.front().uvw, (std::array<int, 3>{0, 0, 1}));
  EXPECT_NEAR(r.axes.front().angle_deg, 7.0, 1e-9);
  EXPECT_NEAR(r.axes.front().spacing_A, 5.431, 1e-9);
  EXPECT_NEAR(r.axes.front().psi1_deg, 0.307 * std::sqrt(5.0 * 14.0 / (5.431 * 0.030)), 0.01);
  // Below E' psi1 exceeds a/d (a = 0.8853 a0 (5^2/3 + 14^2/3)^-1/2 =
  // 0.1585 A: 1.67 degrees), and the critical angle is psi2 =
  // (sqrt 3 a psi1 / (d sqrt 2))^1/2 = 3.6 degrees.
  {
    const double psi1 = r.axes.front().psi1_deg * pi_v / 180.0;
    const double a = 0.8853 * 0.529177 / std::sqrt(std::pow(5.0, 2.0 / 3.0) + std::pow(14.0, 2.0 / 3.0));
    const double psi2 = std::sqrt(std::sqrt(3.0) * a * psi1 / (5.431 * std::sqrt(2.0))) * 180.0 / pi_v;
    EXPECT_NEAR(r.axes.front().psi2_deg, psi2, 1e-9);
    EXPECT_NEAR(r.axes.front().critical_deg, psi2, 1e-9);
    EXPECT_NEAR(psi2, 3.6, 0.1);
  }
  // The mean spacing of the atoms along a row: a / sqrt 2 along <1 1 0> of
  // diamond; along <1 1 1> two atoms per a sqrt 3 (a sqrt(3) / 4 and
  // 3 a sqrt(3) / 4 apart alternately), a sqrt(3) / 2 on average. At 54.74
  // degrees the beam runs along [1 1 1].
  c.ion.angle = std::atan(std::sqrt(2.0)) * amcset::common::radian;
  r = crystal_reports(c).front();
  EXPECT_EQ(r.beam_uvw, (std::array<int, 3>{1, 1, 1}));
  EXPECT_NEAR(r.beam_off_deg, 0.0, 1e-6);
  ASSERT_FALSE(r.axes.empty());
  EXPECT_EQ(r.axes.front().uvw, (std::array<int, 3>{1, 1, 1}));
  EXPECT_NEAR(r.axes.front().spacing_A, 5.431 * std::sqrt(3.0) / 2.0, 1e-9);
  {
    auto c110 = crystal_config(5, 30.0, crystal_layer("diamond", {14}, 5.431, 20000.0, 645.0, {1, 1, 0}, {0, 0, 1}), 1);
    const auto r110 = crystal_reports(c110).front();
    EXPECT_EQ(r110.beam_uvw, (std::array<int, 3>{1, 1, 0}));
    ASSERT_FALSE(r110.axes.empty());
    EXPECT_EQ(r110.axes.front().uvw, (std::array<int, 3>{1, 1, 0}));
    EXPECT_NEAR(r110.axes.front().spacing_A, 5.431 / std::sqrt(2.0), 1e-9);
  }
  // In summary.json, the result tables and the export files.
  c.ion.angle = 7.0 * amcset::common::degree;
  c.trajectory_ions = 3;
  c.trajectory_export.file = "traj.csv";
  c.md_events.file = "events.csv";
  c.md_events.transition_energies = {1000.0 * electron_volt};
  c.output_directory = temp_dir("crystal_report").string();
  c.outputs.formats = {output::Format::json};
  c.ion.spot.width = c.ion.spot.height = 50.0 * angstrom;
  TrimSimulation s(c);
  s.run();
  const std::string summary = s.summary_json();
  EXPECT_NE(summary.find("\"crystals\": ["), std::string::npos);
  EXPECT_NE(summary.find("\"rotation_crystal_to_output_axes\""), std::string::npos);
  EXPECT_NE(summary.find("\"axes_near_beam\": [{\"uvw\": \"[0 0 1]\""), std::string::npos);
  const auto d = make_dataset(s);
  const auto *table = d.find("crystals");
  ASSERT_NE(table, nullptr);
  ASSERT_EQ(table->rows(), 1u);
  const auto traj = read_file(std::filesystem::path(c.output_directory) / "traj.csv");
  ASSERT_GT(traj.size(), 100u);
  EXPECT_NE(traj.substr(0, traj.find('\n')).find(",lattice_dx_A,lattice_dy_A"), std::string::npos);
  const auto ev = read_file(std::filesystem::path(c.output_directory) / "events.csv");
  EXPECT_NE(ev.substr(0, ev.find('\n')).find(",lattice_dx_A,lattice_dy_A"), std::string::npos);
  // Offsets within the spot, and the same for every row of an ion.
  std::istringstream lines(traj);
  std::string line;
  std::getline(lines, line);
  std::map<std::string, std::string> per_ion;
  while (std::getline(lines, line)) {
    const auto ion = line.substr(0, line.find(','));
    const auto tail = line.substr(line.rfind(',', line.rfind(',') - 1) + 1);
    const double dx = std::stod(tail.substr(0, tail.find(',')));
    EXPECT_LE(std::abs(dx), 25.0);
    if (per_ion.count(ion)) {
      EXPECT_EQ(per_ion[ion], tail);
    }
    per_ion[ion] = tail;
  }
  EXPECT_EQ(per_ion.size(), 3u);
  // Amorphous runs keep the columns they had.
  auto amorphous = md_case(10.0, 2);
  amorphous.trajectory_ions = 2;
  amorphous.trajectory_export.file = "traj.csv";
  amorphous.output_directory = temp_dir("crystal_report_amorphous").string();
  amorphous.outputs.formats.clear();
  TrimSimulation a(amorphous);
  a.run();
  const auto plain = read_file(std::filesystem::path(amorphous.output_directory) / "traj.csv");
  EXPECT_EQ(plain.substr(0, plain.find('\n')).find("lattice_dx_A"), std::string::npos);
  EXPECT_TRUE(crystal_reports(amorphous).empty());
}
