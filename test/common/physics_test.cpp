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
 * \file physics_test.cpp
 *
 * \brief Tests of the units layer, the element database, electronic stopping
 * (Lindhard-Scharff, Bethe-Bloch, built-in, SRIM tables) and binary
 * collision kinematics (screening functions, magic formula, scattering
 * integral, nuclear stopping).
 *
 * Reference values are quoted next to each check with their source.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "amcset_units.h"
#include "element_data.h"
#include "physics/scattering.h"
#include "physics/stopping.h"

using namespace amcset;
using namespace amcset::common;
using namespace amcset::physics;

namespace {

//! Relative difference |a - b| / |b|.
double rel(double a, double b) { return std::abs(a - b) / std::abs(b); }

//! Stopping power in MeV cm^2/g for a cross section per molecule of mass m.
double mev_cm2_per_g(stopping_cross_section_quantity s, mass_quantity per_target) {
  const auto mass_stopping = s / per_target; // energy * area / mass
  return in_units(mass_stopping,
                  mega_electron_volt * square_centimeter / gram);
}

} // namespace

// ---------------------------------------------------------------------------
// Units and constants
// ---------------------------------------------------------------------------

TEST(Units, LengthAndEnergyConversions) {
  EXPECT_DOUBLE_EQ(in_units(1.0 * angstrom, meter), 1e-10);
  EXPECT_DOUBLE_EQ(in_units(1.0 * nanometer, angstrom), 10.0);
  EXPECT_NEAR(in_units(1.0 * micrometer, angstrom), 1e4, 1e-9);
  EXPECT_NEAR(in_units(1.0 * kilo_electron_volt, electron_volt), 1000.0, 1e-9);
  // Boost.Units uses CODATA 2006 (1.602176487e-19 J).
  EXPECT_NEAR(in_units(1.0 * electron_volt, joule) / 1.602176634e-19, 1.0, 1e-6);
  EXPECT_NEAR(in_units(1.0 * square_centimeter, angstrom * angstrom), 1e16, 1.0);
  EXPECT_NEAR(in_units(1.0 * degree, radian), pi_v / 180.0, 1e-15);
}

TEST(Units, PhysicalConstants) {
  // CODATA values (Boost.Units uses CODATA 2006; differences ~1e-7).
  EXPECT_NEAR(in_units(coulomb_energy_length, electron_volt * angstrom), 14.399645, 1e-5);
  EXPECT_NEAR(in_units(bohr_radius, angstrom), 0.529177, 1e-6);
  EXPECT_NEAR(in_units(electron_rest_energy, kilo_electron_volt), 510.99895, 1e-4);
  EXPECT_NEAR(in_units(atomic_mass_unit_energy, mega_electron_volt), 931.494102, 1e-3);
  EXPECT_NEAR(1.0 / fine_structure_constant, 137.035999, 1e-5);
  EXPECT_NEAR(in_units(rydberg_energy, electron_volt), 13.605693, 1e-5);
  EXPECT_NEAR(in_units(bohr_velocity, meters_per_second), 2.18769126e6, 10.0);
}

TEST(Units, SrimStoppingUnit) {
  // 1 eV/(1e15 atoms/cm^2) = 1e-15 eV cm^2 = 10 eV Angstrom^2.
  EXPECT_NEAR(in_units(1.0 * srim_stopping_unit, electron_volt * angstrom * angstrom), 10.0, 1e-9);
}

// ---------------------------------------------------------------------------
// Element data
// ---------------------------------------------------------------------------

TEST(Elements, LookupBySymbolAndNumber) {
  const auto &fe = Element::by_symbol("Fe");
  EXPECT_EQ(fe.z(), 26u);
  EXPECT_EQ(Element::by_z(26).symbol(), "Fe");
  EXPECT_EQ(Element::by_symbol("fe").z(), 26u);   // case-insensitive
  EXPECT_EQ(Element::by_symbol("Iron").z(), 26u); // English name
  EXPECT_EQ(Element::find_z("Cr"), std::optional<unsigned>(24));
  EXPECT_FALSE(Element::find_z("Xx").has_value());
  EXPECT_THROW(Element::by_symbol("Xx"), std::invalid_argument);
  EXPECT_THROW(Element::by_z(0), std::out_of_range);
  EXPECT_EQ(Element::all().size(), 92u);
  for (const auto &e : Element::all()) {
    EXPECT_GT(e.natural_mass(), 0.0 * atomic_mass_unit) << e.symbol();
    EXPECT_GT(e.mean_excitation_energy(), 0.0 * electron_volt) << e.symbol();
  }
}

TEST(Elements, TrimAtomicData) {
  // SRIM-2013 ATOMDATA / TRIM defaults for Fe.
  const auto &fe = Element::by_z(26);
  EXPECT_NEAR(in_units(fe.natural_mass(), atomic_mass_unit), 55.847, 0.01);
  EXPECT_NEAR(in_units(fe.most_abundant_isotope_mass(), atomic_mass_unit), 55.935, 0.01);
  EXPECT_NEAR(in_units(fe.solid_density(), gram_per_cubic_centimeter), 7.8658, 1e-3);
  EXPECT_NEAR(in_units(fe.sublimation_energy(), electron_volt), 4.34, 1e-6);
  // Atomic density 8.48e22 atoms/cm^3 (TRIM prints 8.481E22).
  EXPECT_NEAR(in_units(fe.solid_atomic_density(), per_cubic_centimeter), 8.481e22, 0.003e22);
  // NIST ICRU-37 I-value of Fe: 286 eV.
  EXPECT_NEAR(in_units(fe.mean_excitation_energy(), electron_volt), 286.0, 1.0);
  EXPECT_NEAR(in_units(Element::by_z(1).mean_excitation_energy(), electron_volt), 19.2, 0.1);
}

// ---------------------------------------------------------------------------
// Kinematics
// ---------------------------------------------------------------------------

TEST(Stopping, RelativisticKinematics) {
  const mass_quantity mp = 1.007276 * atomic_mass_unit;
  // A proton with T = m c^2 has gamma = 2 and beta = sqrt(3)/2.
  const energy_quantity t = 938.272 * mega_electron_volt;
  EXPECT_NEAR(lorentz_gamma(t, mp), 2.0, 1e-4);
  EXPECT_NEAR(relativistic_beta(t, mp), std::sqrt(3.0) / 2.0, 1e-4);
  // Non-relativistic limit: v = sqrt(2E/m).
  const energy_quantity e = 1.0 * kilo_electron_volt;
  const double v = in_units(relativistic_speed(e, mp), meters_per_second);
  const double v_classical = std::sqrt(2.0 * in_units(e, joule) / in_units(mp, kilogram));
  EXPECT_LT(rel(v, v_classical), 1e-5);
  EXPECT_LT(relativistic_beta(100.0 * giga_electron_volt, mp), 1.0);
}

// ---------------------------------------------------------------------------
// Electronic stopping models
// ---------------------------------------------------------------------------

TEST(Stopping, LindhardScharffIsVelocityProportional) {
  const mass_quantity m = 55.935 * atomic_mass_unit;
  const auto s1 = lindhard_scharff_stopping(26, m, 26, 10.0 * kilo_electron_volt);
  const auto s4 = lindhard_scharff_stopping(26, m, 26, 40.0 * kilo_electron_volt);
  EXPECT_NEAR(s4 / s1, 2.0, 1e-5); // (relativistic speed)
  // Hand evaluation for Fe in Fe at 10 keV: xi_e = 26^(1/6) = 1.7213,
  // 8 pi e^2 a0 = 191.5 eV A^2, Z1 Z2 / (Z1^(2/3)+Z2^(2/3))^(3/2) = 9.193,
  // v/v0 = 0.0849 -> 257 eV A^2 = 25.7 eV/(1e15 atoms/cm^2).
  EXPECT_NEAR(in_units(s1, srim_stopping_unit), 25.7, 0.1);
}

TEST(Stopping, BetheBlochMatchesPstarForWater) {
  // NIST PSTAR, liquid water (I = 75 eV), electronic stopping of 100 MeV
  // protons: 7.289 MeV cm^2/g. Evaluated per H2O molecule (Z = 10).
  const mass_quantity mp = 1.007276 * atomic_mass_unit;
  const mass_quantity water = 18.015 * atomic_mass_unit;
  const number_density_quantity n_molecules =
      (1.0 * gram_per_cubic_centimeter) / water;
  BetheBlochOptions options;
  options.effective_charge = false;
  const auto s = bethe_bloch_stopping(1, mp, 10, 75.0 * electron_volt,
                                      100.0 * mega_electron_volt,
                                      10.0 * n_molecules, options);
  EXPECT_NEAR(mev_cm2_per_g(s, water), 7.289, 0.05);
  // 10 MeV: PSTAR 45.67 MeV cm^2/g (no shell correction here, so ~1 %).
  const auto s10 = bethe_bloch_stopping(1, mp, 10, 75.0 * electron_volt,
                                        10.0 * mega_electron_volt,
                                        10.0 * n_molecules, options);
  EXPECT_LT(rel(mev_cm2_per_g(s10, water), 45.67), 0.02);
}

TEST(Stopping, BetheBlochMinimumIonizationAndDensityEffect) {
  // Minimum ionization near beta*gamma ~ 3.5; the density effect lowers the
  // relativistic rise.
  const mass_quantity mp = 1.007276 * atomic_mass_unit;
  const auto &fe = Element::by_z(26);
  const auto n = fe.solid_atomic_density();
  BetheBlochOptions with;
  BetheBlochOptions without;
  without.density_effect = false;
  const auto i = fe.mean_excitation_energy();
  const auto s_min = bethe_bloch_stopping(1, mp, 26, i, 2.5 * giga_electron_volt, 26.0 * n, with);
  const auto s_low = bethe_bloch_stopping(1, mp, 26, i, 0.3 * giga_electron_volt, 26.0 * n, with);
  const auto s_high = bethe_bloch_stopping(1, mp, 26, i, 100.0 * giga_electron_volt, 26.0 * n, with);
  const auto s_high_nd = bethe_bloch_stopping(1, mp, 26, i, 100.0 * giga_electron_volt, 26.0 * n, without);
  EXPECT_LT(s_min, s_low);
  EXPECT_LT(s_min, s_high);
  EXPECT_LT(s_high, s_high_nd);
  // PDG: minimum ionization of iron 1.451 MeV cm^2/g.
  EXPECT_NEAR(mev_cm2_per_g(s_min, fe.natural_mass()), 1.451, 0.05);
}

TEST(Stopping, BuiltinIsCloseToSrimAboveThePeak) {
  // SRIM-2013 electronic stopping of H in Fe (eV/Angstrom, density 7.8658):
  // 1 MeV: 10.38, 10 MeV: 2.260, 100 MeV: 0.3975.
  const mass_quantity mh = 1.008 * atomic_mass_unit;
  const auto &fe = Element::by_z(26);
  const auto n = fe.solid_atomic_density();
  const auto sp = [&](double mev) {
    return in_units(builtin_stopping(1, mh, 26, mev * mega_electron_volt, 26.0 * n) * n,
                    electron_volt_per_angstrom);
  };
  EXPECT_LT(rel(sp(1.0), 10.38), 0.12);
  EXPECT_LT(rel(sp(10.0), 2.260), 0.06);
  EXPECT_LT(rel(sp(100.0), 0.3975), 0.03);
  // At low energy it is Lindhard-Scharff-like.
  EXPECT_LT(rel(sp(0.01), 9.20), 0.30);
}

TEST(Stopping, StragglingBohrLimit) {
  // Bohr: Omega^2/L = 4 pi Z1^2 Z2 e^4 N; 3 MeV H in Fe is close to it.
  const mass_quantity mh = 1.008 * atomic_mass_unit;
  const auto n = Element::by_z(26).solid_atomic_density();
  const auto w = trim_electronic_straggling(1, mh, 26, 3.0 * mega_electron_volt, n);
  const double bohr = 4.0 * pi_v * 26.0 *
                      std::pow(in_units(coulomb_energy_length, electron_volt * angstrom), 2) *
                      in_units(n, per_cubic_angstrom);
  EXPECT_LT(rel(in_units(w, electron_volt * electron_volt / angstrom), bohr), 0.01);
  // The low-energy branch takes over at low energy (smaller value).
  const auto w_low = trim_electronic_straggling(1, mh, 26, 1.0 * kilo_electron_volt, n);
  EXPECT_LT(w_low, w);
}

TEST(Stopping, ModelNamesRoundTrip) {
  for (auto m : {ElectronicStoppingModel::srim_tables, ElectronicStoppingModel::builtin,
                 ElectronicStoppingModel::bethe_bloch, ElectronicStoppingModel::lindhard_scharff}) {
    EXPECT_EQ(parse_electronic_stopping_model(to_string(m)), m);
  }
  EXPECT_THROW(parse_electronic_stopping_model("nonsense"), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// Stopping tables
// ---------------------------------------------------------------------------

namespace {

//! A synthetic table: S_e = 10 * sqrt(E/keV), S_n = 1/E, from 1 keV to 1 MeV.
StoppingTable synthetic_table(unsigned z1, double m1_amu, unsigned z2) {
  StoppingTable t;
  t.z1 = z1;
  t.m1 = m1_amu * atomic_mass_unit;
  t.z2 = z2;
  for (int i = 0; i <= 30; ++i) {
    const double e = std::pow(10.0, i / 10.0);
    t.log_energy.push_back(std::log(e));
    t.log_electronic.push_back(std::log(10.0 * std::sqrt(e)));
    t.log_nuclear.push_back(std::log(1.0 / e));
  }
  return t;
}

} // namespace

TEST(StoppingTable, LogLogInterpolationAndExtrapolation) {
  const auto t = synthetic_table(1, 1.0, 26);
  const mass_quantity m = 1.0 * atomic_mass_unit;
  // Power laws are reproduced exactly by log-log interpolation.
  EXPECT_NEAR(in_units(t.electronic(50.0 * kilo_electron_volt, m), srim_stopping_unit),
              10.0 * std::sqrt(50.0), 1e-9);
  EXPECT_NEAR(in_units(t.nuclear(50.0 * kilo_electron_volt), srim_stopping_unit), 1.0 / 50.0, 1e-12);
  // Velocity scaling: a projectile of twice the mass at twice the energy has
  // the same velocity and so the same S_e.
  EXPECT_NEAR(in_units(t.electronic(100.0 * kilo_electron_volt, 2.0 * m), srim_stopping_unit),
              10.0 * std::sqrt(50.0), 1e-9);
  // Below the table: proportional to the velocity.
  EXPECT_NEAR(in_units(t.electronic(0.25 * kilo_electron_volt, m), srim_stopping_unit),
              10.0 * 0.5, 1e-9);
  EXPECT_NEAR(in_units(t.min_energy(), kilo_electron_volt), 1.0, 1e-9);
  EXPECT_NEAR(in_units(t.max_energy(), kilo_electron_volt), 1000.0, 1e-6);
}

TEST(StoppingTable, LoadsAmcsetFormat) {
  const auto dir = std::filesystem::temp_directory_path() / "amcset_stopping_test";
  std::filesystem::create_directories(dir);
  const auto file = dir / "tiny.amcset";
  {
    std::ofstream out(file);
    out << "amcset-stopping-table 1\n"
           "# comment\n"
           "source unit test\n"
           "energy_unit keV\n"
           "stopping_unit eV/(1e15 atoms/cm2)\n"
           "energies 3 1 10 100\n"
           "table 1 1.008 26 s\n"
           "Se 10 20 30\n"
           "Sn 1 0.5 0.1\n";
  }
  std::shared_ptr<const SrimStoppingDatabase> db;
  ASSERT_NO_THROW(db = SrimStoppingDatabase::load(file));
  ASSERT_EQ(db->size(), 1u);
  const auto *t = db->find(1, 26, false);
  ASSERT_NE(t, nullptr);
  EXPECT_NEAR(in_units(t->electronic(10.0 * kilo_electron_volt, t->m1), srim_stopping_unit), 20.0, 1e-9);
  // A gas lookup falls back to the solid table.
  EXPECT_EQ(db->find(1, 26, true), t);
  EXPECT_EQ(db->find(2, 26, false), nullptr);
  EXPECT_THROW(SrimStoppingDatabase::load(dir / "missing.amcset"), std::runtime_error);
}

// Packed files are indexed at load time and each table is parsed when it is
// first asked for: loading must not trip over a table nobody uses, a later
// table of the same pair replaces an earlier one, and the values equal those
// of eager parsing (the same text through the same parser).
TEST(StoppingTable, PackedTablesAreParsedOnDemand) {
  const auto dir = std::filesystem::temp_directory_path() / "amcset_stopping_lazy_test";
  std::filesystem::create_directories(dir);
  const auto file = dir / "lazy.amcset";
  {
    std::ofstream out(file, std::ios::binary);
    out << "amcset-stopping-table 1\r\n" // CRLF lines are accepted
           "energy_unit keV\n"
           "stopping_unit eV/(1e15 atoms/cm2)\n"
           "energies 3 1 10 100\n"
           "table 1 1.008 26 s\n"
           "Se 1 2 3\n"
           "Sn 1 1 1\n"
           "table 2 4.0026 26 s\n"
           "Se 10 20\n" // malformed: two values for three energies
           "Sn 1 1 1\n"
           "table 1 1.008 26 s\n" // replaces the first table
           "Se 10 20 30\n"
           "Sn 1 0.5 0.1\n";
  }
  std::shared_ptr<const SrimStoppingDatabase> db;
  ASSERT_NO_THROW(db = SrimStoppingDatabase::load(file));
  EXPECT_EQ(db->size(), 2u);
  const auto *t = db->find(1, 26, false);
  ASSERT_NE(t, nullptr);
  EXPECT_NEAR(in_units(t->electronic(10.0 * kilo_electron_volt, t->m1), srim_stopping_unit), 20.0, 1e-9);
  EXPECT_EQ(db->find(1, 26, false), t); // parsed once, then kept
  EXPECT_THROW(db->find(2, 26, false), std::runtime_error);
  EXPECT_EQ(db->size(), 2u);
  // A truncated file is still an error when loading.
  const auto cut = dir / "cut.amcset";
  {
    std::ofstream out(cut, std::ios::binary);
    out << "amcset-stopping-table 1\nenergies 2 1 10\ntable 1 1.008 26 s\nSe 1 2\n";
  }
  EXPECT_THROW(SrimStoppingDatabase::load(cut), std::runtime_error);
}

TEST(StoppingTable, LayerStoppingUsesBraggsRuleAndCorrection) {
  SrimStoppingDatabase db;
  db.add(synthetic_table(1, 1.0, 26));
  auto cr = synthetic_table(1, 1.0, 24);
  for (auto &v : cr.log_electronic) v += std::log(2.0); // Cr: twice the stopping
  db.add(cr);
  const mass_quantity m = 1.0 * atomic_mass_unit;
  const auto n = 8.0e22 * per_cubic_centimeter;
  const std::vector<StoppingComponent> alloy{{26, 0.75}, {24, 0.25}};
  const auto e = 10.0 * kilo_electron_volt;
  const double s_fe = 10.0 * std::sqrt(10.0);
  const auto s = layer_stopping(ElectronicStoppingModel::srim_tables, &db, 1, m, alloy,
                                false, 0.9, e, n);
  EXPECT_NEAR(in_units(s, srim_stopping_unit), 0.9 * (0.75 * s_fe + 0.25 * 2.0 * s_fe), 1e-9);
  // A missing pair is an error in table mode.
  const std::vector<StoppingComponent> ni{{28, 1.0}};
  EXPECT_THROW(layer_stopping(ElectronicStoppingModel::srim_tables, &db, 1, m, ni, false, 1.0, e, n),
               std::exception);
}

TEST(StoppingTable, CurveContinuesWithBetheBlochAboveTheTable) {
  SrimStoppingDatabase db;
  db.add(synthetic_table(1, 1.0, 26));
  const mass_quantity m = 1.0 * atomic_mass_unit;
  const auto n = Element::by_z(26).solid_atomic_density();
  const std::vector<StoppingComponent> fe{{26, 1.0}};
  StoppingCurve curve(ElectronicStoppingModel::srim_tables, &db, 1, m, fe, false, 1.0, n);
  ASSERT_FALSE(curve.empty());
  // Inside the table the curve reproduces it.
  EXPECT_LT(rel(in_units(curve.cross_section(100.0 * kilo_electron_volt), srim_stopping_unit),
                10.0 * std::sqrt(100.0)), 2e-3);
  // Continuous at the upper end of the table.
  const double below = in_units(curve.cross_section(0.99 * mega_electron_volt), srim_stopping_unit);
  const double above = in_units(curve.cross_section(1.01 * mega_electron_volt), srim_stopping_unit);
  EXPECT_LT(rel(above, below), 0.02);
  // Far above, it falls like Bethe-Bloch (roughly 1/E).
  const double s1 = in_units(curve.cross_section(10.0 * mega_electron_volt), srim_stopping_unit);
  const double s2 = in_units(curve.cross_section(100.0 * mega_electron_volt), srim_stopping_unit);
  EXPECT_GT(s1 / s2, 4.0);
  EXPECT_LT(s1 / s2, 10.0);
  // stopping_power = N * S.
  EXPECT_LT(rel(in_units(curve.stopping_power(1.0 * kilo_electron_volt), electron_volt_per_angstrom),
                in_units(curve.cross_section(1.0 * kilo_electron_volt) * n, electron_volt_per_angstrom)),
            1e-12);
}

// Optional: checks against the real SRIM-2013 tables when they are available
// (set AMCSET_SRIM_STOPPING_TABLE to the file made by tools/srim).
TEST(StoppingTable, RealSrimTablesIfAvailable) {
  const char *path = std::getenv("AMCSET_SRIM_STOPPING_TABLE");
  if (!path) GTEST_SKIP() << "AMCSET_SRIM_STOPPING_TABLE not set";
  const auto db = SrimStoppingDatabase::load(path);
  const auto *t = db->find(1, 26, false);
  ASSERT_NE(t, nullptr);
  const auto n = Element::by_z(26).solid_atomic_density();
  // SR Module, H in Fe at 100 keV: 24.21 eV/Angstrom.
  EXPECT_LT(rel(in_units(t->electronic(100.0 * kilo_electron_volt, 1.008 * atomic_mass_unit) * n,
                         electron_volt_per_angstrom), 24.21), 0.01);
}

// ---------------------------------------------------------------------------
// Screening functions and scattering
// ---------------------------------------------------------------------------

TEST(Scattering, UniversalScreeningFunction) {
  const auto f = ScreeningFunction::universal();
  // phi(0) = 1: the coefficients sum to one (the 0.02817 coefficient; the
  // native engine used 0.2817 before physics fix P4).
  EXPECT_NEAR(f.phi(0.0), 1.0, 1e-4);
  EXPECT_NEAR(f.c[3], 0.02817, 1e-5);
  // dphi is the derivative of phi.
  const double x = 0.7, h = 1e-6;
  EXPECT_NEAR(f.dphi(x), (f.phi(x + h) - f.phi(x - h)) / (2 * h), 1e-6);
  // ZBL screening length, Fe-Fe: 0.8854 a0 / (2 * 26^0.23) = 0.1107 A.
  EXPECT_NEAR(in_units(universal_screening_length(26, 26), angstrom), 0.11066, 2e-4);
}

TEST(Scattering, MagicFormulaAgreesWithTheScatteringIntegral) {
  const auto f = ScreeningFunction::universal();
  for (double eps : {0.01, 0.1, 1.0, 10.0}) {
    for (double b : {0.05, 0.3, 1.0, 3.0}) {
      const double magic = magic_sin2_half_theta(f, eps, b);
      const double exact = integral_sin2_half_theta(f, eps, b, 200);
      if (exact < 1e-3) continue; // tiny transfers (below TRIM's cut-off)
      EXPECT_LT(rel(magic, exact), 0.05) << "eps " << eps << " b " << b;
    }
  }
}

TEST(Scattering, RutherfordLimitAtHighEnergy) {
  // Unscreened Rutherford: sin^2(theta/2) = 1/(1 + (2 eps b)^2).
  const auto f = ScreeningFunction::universal();
  const double eps = 1000.0, b = 0.002;
  const double ruth = 1.0 / (1.0 + std::pow(2 * eps * b, 2));
  EXPECT_LT(rel(integral_sin2_half_theta(f, eps, b, 200), ruth), 0.02);
  EXPECT_LT(rel(rutherford_sin2_half_theta(eps, b), ruth), 0.01);
  // TRIM switches to the Rutherford form above eps = 10.
  EXPECT_DOUBLE_EQ(trim_sin2_half_theta(f, 20.0, 0.1), rutherford_sin2_half_theta(20.0, 0.1));
  EXPECT_DOUBLE_EQ(trim_sin2_half_theta(f, 5.0, 0.1), magic_sin2_half_theta(f, 5.0, 0.1));
}

// Extreme impact parameters and energies (a solid of gas density gives
// b of thousands, where the potential underflows to zero): finite angles in
// [0, 1], no deflection far outside the screening radius. A NaN here once
// stopped a run ("the ion's energy or depth is not a number").
TEST(Scattering, ExtremeImpactParametersGiveFiniteAngles) {
  const auto f = ScreeningFunction::universal();
  const CollisionPair pair(26, 55.935 * atomic_mass_unit, 26, 55.847 * atomic_mass_unit);
  for (double eps : {1e-6, 1e-3, 0.1, 1.0, 5.7, 9.99, 10.5, 1e3, 1e6}) {
    for (double b : {1e-12, 1e-6, 1e-3, 0.5, 10.0, 255.0, 300.0, 3000.0, 3495.35, 3600.0, 3700.0, 5000.0, 1e5, 1e9}) {
      for (const auto method : {ScatteringMethod::magic, ScatteringMethod::integral}) {
        const double s2 = trim_sin2_half_theta(f, eps, b, method);
        EXPECT_TRUE(std::isfinite(s2) && s2 >= 0.0 && s2 <= 1.0) << "eps " << eps << " b " << b << ": " << s2;
        double sh = -1.0, ch = -1.0;
        fast_half_angle(pair, eps, b, method, sh, ch);
        EXPECT_TRUE(std::isfinite(sh) && std::isfinite(ch) && sh >= 0.0 && sh <= 1.0 && ch >= 0.0 && ch <= 1.0)
            << "fast: eps " << eps << " b " << b << ": " << sh << " " << ch;
      }
    }
  }
  EXPECT_EQ(magic_sin2_half_theta(f, 5.7, 5000.0), 0.0);
  EXPECT_EQ(magic_sin2_half_theta(f, 5.7, 3495.35), 0.0); // denormal potential
  EXPECT_EQ(magic_sin2_half_theta(f, 1e-3, 1e9), 0.0);
}

TEST(Scattering, HeadOnCollisionTransfersTheMaximum) {
  const auto f = ScreeningFunction::universal();
  EXPECT_NEAR(magic_sin2_half_theta(f, 1.0, 0.0), 1.0, 1e-9);
  EXPECT_NEAR(integral_sin2_half_theta(f, 1.0, 1e-8, 100), 1.0, 1e-4);
  // Closest approach: at b = 0 it is the root of 1 - phi(x)/(x eps) = 0.
  const double r0 = closest_approach(f, 1.0, 0.0);
  EXPECT_NEAR(1.0 - f.phi(r0) / (r0 * 1.0), 0.0, 1e-8);
}

TEST(Scattering, ReducedNuclearStoppingMatchesTheScatteringIntegral) {
  // S_n(eps) = (1/eps) int_0^inf sin^2(theta/2) d(b^2) (reduced units),
  // which the ZBL universal fit reproduces to a few percent.
  const auto f = ScreeningFunction::universal();
  for (double eps : {0.01, 0.1, 1.0, 10.0}) {
    double sum = 0.0;
    const int n = 4000;
    const double bmax = 40.0;
    for (int i = 0; i < n; ++i) {
      // integrate in u = ln b (b from 1e-5 to bmax)
      const double u0 = std::log(1e-5), u1 = std::log(bmax);
      const double u = u0 + (u1 - u0) * (i + 0.5) / n;
      const double b = std::exp(u);
      sum += integral_sin2_half_theta(f, eps, b, 64) * 2.0 * b * b * (u1 - u0) / n;
    }
    const double sn = eps * sum; // ZBL normalization: S_n = eps int sin^2 d(b^2)
    EXPECT_LT(rel(sn, reduced_nuclear_stopping(eps)), 0.06) << "eps " << eps;
  }
}

TEST(Scattering, CollisionKinematicsConserveEnergyAndMomentum) {
  CollisionPair pair(1, 1.008 * atomic_mass_unit, 26, 55.847 * atomic_mass_unit);
  const energy_quantity e = 100.0 * kilo_electron_volt;
  for (double s2 : {1e-4, 0.1, 0.5, 0.9, 1.0}) {
    const auto k = collision_kinematics(pair, e, s2);
    const double t = in_units(k.transferred, kilo_electron_volt);
    EXPECT_NEAR(t, pair.gamma * s2 * 100.0, 1e-9);
    // Momentum along and across the incident direction (p ~ sqrt(2 m E)).
    const double m1 = 1.008, m2 = 55.847;
    const double p0 = std::sqrt(2 * m1 * 100.0);
    const double p1 = std::sqrt(2 * m1 * (100.0 - t));
    const double p2 = std::sqrt(2 * m2 * t);
    EXPECT_NEAR(p1 * k.cos_projectile + p2 * k.cos_recoil, p0, 1e-6 * p0) << s2;
    EXPECT_NEAR(p1 * k.sin_projectile - p2 * k.sin_recoil, 0.0, 1e-6 * p0) << s2;
  }
  // Equal masses: gamma = 1 and the two directions are perpendicular.
  CollisionPair same(26, 55.847 * atomic_mass_unit, 26, 55.847 * atomic_mass_unit);
  EXPECT_NEAR(same.gamma, 1.0, 1e-12);
  const auto k = collision_kinematics(same, e, 0.3);
  EXPECT_NEAR(k.cos_projectile * k.cos_recoil - k.sin_projectile * k.sin_recoil, 0.0, 1e-9);
}

TEST(Scattering, FreeFlightImpactParameterDecreasesWithEnergy) {
  const double gamma = 1.0, eps_min = 1e-4;
  double prev = 1e9;
  for (double eps : {0.01, 0.1, 1.0, 10.0, 100.0}) {
    const double b = trim_max_reduced_impact_parameter(eps, eps_min, gamma);
    EXPECT_GT(b, 0.0);
    EXPECT_LT(b, prev);
    prev = b;
  }
}
