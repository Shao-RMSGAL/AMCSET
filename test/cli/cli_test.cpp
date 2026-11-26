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
 * \file cli_test.cpp
 *
 * \brief Tests of the JSON run configuration (run_config.h) and of the
 * command-line layer (cli.h): unit parsing, schema checks, overrides,
 * round trips, engine construction and argument handling.
 */

#include <gtest/gtest.h>

#include <atomic>
#include <csignal>

#include <filesystem>
#include <fstream>
#include <sstream>

#include "cli.h"
#include "interrupts.h"
#include "element_data.h"
#include "run_config.h"
#include "trim/trim.h"

using namespace amcset;
using namespace amcset::common;
using namespace amcset::config;

// ---------------------------------------------------------------------------
// Quantities with units
// ---------------------------------------------------------------------------

TEST(RunConfigUnits, ParsesQuantitiesWithUnits) {
  EXPECT_NEAR(in_units(parse_energy("2 MeV"), kilo_electron_volt), 2000.0, 1e-9);
  EXPECT_NEAR(in_units(parse_energy("150keV"), electron_volt), 150000.0, 1e-6);
  EXPECT_NEAR(in_units(parse_energy("40 eV"), electron_volt), 40.0, 1e-12);
  EXPECT_NEAR(in_units(parse_energy("1.5 GeV"), mega_electron_volt), 1500.0, 1e-9);
  EXPECT_NEAR(in_units(parse_length("1 um"), angstrom), 1e4, 1e-6);
  EXPECT_NEAR(in_units(parse_length("250 nm"), angstrom), 2500.0, 1e-9);
  EXPECT_NEAR(in_units(parse_length("3.5 A"), angstrom), 3.5, 1e-12);
  EXPECT_NEAR(in_units(parse_length("2 mm"), micrometer), 2000.0, 1e-6);
  EXPECT_NEAR(in_units(parse_mass("55.847 amu"), atomic_mass_unit), 55.847, 1e-12);
  EXPECT_NEAR(in_units(parse_density("7.8658 g/cm3"), gram_per_cubic_centimeter), 7.8658, 1e-12);
  EXPECT_NEAR(in_units(parse_angle("7 deg"), degree), 7.0, 1e-12);
  EXPECT_NEAR(in_units(parse_angle("0.1 rad"), radian), 0.1, 1e-12);
}

TEST(RunConfigUnits, RejectsMissingOrWrongUnits) {
  EXPECT_THROW(parse_energy("100"), std::invalid_argument);       // no unit
  EXPECT_THROW(parse_energy("100 um"), std::invalid_argument);    // wrong dimension
  EXPECT_THROW(parse_length("abc nm"), std::invalid_argument);
  EXPECT_THROW(parse_density("7.8"), std::invalid_argument);
}

TEST(RunConfigUnits, FormatRoundTrips) {
  for (double kev : {0.001, 1.0, 150.0, 2000.0, 5e6}) {
    const auto e = kev * kilo_electron_volt;
    EXPECT_NEAR(in_units(parse_energy(format_energy(e)), kilo_electron_volt), kev, 1e-9 * kev);
  }
  const auto l = 12345.6 * angstrom;
  EXPECT_NEAR(in_units(parse_length(format_length(l)), angstrom), 12345.6, 1e-6);
  EXPECT_NEAR(in_units(parse_density(format_density(7.8658 * gram_per_cubic_centimeter)),
                       gram_per_cubic_centimeter), 7.8658, 1e-9);
}

// ---------------------------------------------------------------------------
// JSON configuration
// ---------------------------------------------------------------------------

namespace {

const char *alloy_json = R"({
  "mode": "srim",
  "title": "Fe on FeCr",
  "ion": {"element": "Fe", "energy": "1 MeV", "angle": "7 deg"},
  "layers": [
    {"name": "FeCr", "thickness": "1.5 um", "density": "7.8 g/cm3",
     "elements": [
       {"element": "Fe", "fraction": 0.91, "Ed": "40 eV", "Eb": "0 eV", "Es": "4.34 eV"},
       {"element": "Cr", "fraction": 0.09, "Ed": "28 eV"}
     ]},
    {"name": "Substrate", "thickness": "2 um", "density": "8.9 g/cm3", "compound_correction": 0.97,
     "elements": [{"element": "Ni"}]}
  ],
  "ions": 250,
  "seed": 42,
  "threads": 2,
  "damage": "full",
  "stopping": {"model": "builtin"},
  "output": {"directory": "out", "bins": 50}
})";

} // namespace

TEST(RunConfig, ParsesAFullSrimModeScenario) {
  const auto c = parse_config(alloy_json, "/base");
  EXPECT_EQ(c.mode, Mode::srim);
  const auto &t = c.trim;
  EXPECT_EQ(t.title, "Fe on FeCr");
  EXPECT_EQ(t.ion.z, 26u);
  EXPECT_NEAR(in_units(t.ion.mass, atomic_mass_unit), 55.935, 1e-3); // default isotope
  EXPECT_NEAR(in_units(t.ion.energy, kilo_electron_volt), 1000.0, 1e-9);
  EXPECT_NEAR(in_units(t.ion.angle, degree), 7.0, 1e-12);
  ASSERT_EQ(t.layers.size(), 2u);
  EXPECT_NEAR(in_units(t.layers[0].width, angstrom), 15000.0, 1e-6);
  ASSERT_EQ(t.layers[0].elements.size(), 2u);
  EXPECT_NEAR(in_units(t.layers[0].elements[1].displacement_energy, electron_volt), 28.0, 1e-12);
  // Unspecified per-element energies take the element defaults.
  EXPECT_NEAR(in_units(t.layers[0].elements[1].surface_binding_energy, electron_volt),
              in_units(trim::TargetElement::defaults(24).surface_binding_energy, electron_volt),
              1e-12);
  EXPECT_NEAR(t.layers[1].compound_correction, 0.97, 1e-12);
  EXPECT_EQ(t.ion_count, 250u);
  EXPECT_EQ(t.seed, 42u);
  EXPECT_EQ(t.damage, trim::DamageCalculation::full_cascade);
  EXPECT_EQ(t.depth_bins, 50u);
  // Relative paths are resolved against the configuration's directory.
  EXPECT_EQ(std::filesystem::path(t.output_directory), std::filesystem::path("/base/out"));
}

TEST(RunConfig, RejectsUnknownKeysAndBadValues) {
  EXPECT_THROW(parse_config(R"({"ion": {"element": "Fe", "energy": "1 MeV"}, "layerz": []})"),
               std::exception);
  EXPECT_THROW(parse_config(R"({"ion": {"element": "Qq", "energy": "1 MeV"}, "layers": []})"),
               std::exception);
  EXPECT_THROW(parse_config("{not json"), std::exception);
  EXPECT_THROW(parse_config(R"({"mode": "warp"})"), std::exception);
}

TEST(RunConfig, JsonRoundTrip) {
  const auto a = parse_config(alloy_json, "/base");
  const auto text = to_json(a);
  const auto b = parse_config(text, "/base");
  EXPECT_EQ(b.trim.ion.z, a.trim.ion.z);
  EXPECT_NEAR(in_units(b.trim.ion.energy, electron_volt), in_units(a.trim.ion.energy, electron_volt), 1e-6);
  ASSERT_EQ(b.trim.layers.size(), a.trim.layers.size());
  EXPECT_NEAR(b.trim.layers[0].elements[0].stoichiometry, a.trim.layers[0].elements[0].stoichiometry, 1e-12);
  EXPECT_NEAR(in_units(b.trim.layers[0].elements[1].displacement_energy, electron_volt), 28.0, 1e-9);
  EXPECT_EQ(b.trim.damage, a.trim.damage);
  EXPECT_EQ(b.trim.seed, a.trim.seed);
  EXPECT_EQ(to_json(b), text);
}

TEST(RunConfig, OverridesUseDottedPaths) {
  auto c = parse_config(alloy_json);
  apply_override(c, "ion.energy=2 MeV");
  apply_override(c, "layers.0.elements.1.Ed=33 eV");
  apply_override(c, "ions=10");
  apply_override(c, "damage=kp");
  EXPECT_NEAR(in_units(c.trim.ion.energy, kilo_electron_volt), 2000.0, 1e-9);
  EXPECT_NEAR(in_units(c.trim.layers[0].elements[1].displacement_energy, electron_volt), 33.0, 1e-12);
  EXPECT_EQ(c.trim.ion_count, 10u);
  EXPECT_EQ(c.trim.damage, trim::DamageCalculation::kinchin_pease);
  EXPECT_THROW(apply_override(c, "no_equals_sign"), std::exception);
  EXPECT_THROW(apply_override(c, "layers.9.name=x"), std::exception);
}

TEST(RunConfig, TransportOption) {
  auto c = parse_config(alloy_json);
  EXPECT_EQ(c.trim.transport, trim::TrimConfig::Transport::trim2013); // default
  apply_override(c, "physics.transport=refined");
  EXPECT_EQ(c.trim.transport, trim::TrimConfig::Transport::refined);
  // Written out and read back.
  const auto again = parse_config(to_json(c));
  EXPECT_EQ(again.trim.transport, trim::TrimConfig::Transport::refined);
  apply_override(c, "physics.transport=trim");
  EXPECT_EQ(c.trim.transport, trim::TrimConfig::Transport::trim2013);
  EXPECT_THROW(apply_override(c, "physics.transport=fast"), std::exception);
}

TEST(RunConfig, ExamplesAreValidAndBuildEngines) {
  for (auto mode : {Mode::srim, Mode::native}) {
    auto c = example_config(mode);
    EXPECT_EQ(c.mode, mode);
    const auto again = parse_config(to_json(c));
    EXPECT_EQ(again.mode, mode);
    c.trim.ion_count = 3;
    c.trim.threads = 1;
    const auto dir = std::filesystem::temp_directory_path() / ("amcset_cli_example_" + to_string(mode));
    std::filesystem::remove_all(dir);
    c.trim.output_directory = dir.string();
    auto engine = make_engine(c);
    ASSERT_NE(engine, nullptr);
    EXPECT_EQ(engine->total_count(), 3u);
    EXPECT_FALSE(engine->describe().empty());
  }
}

TEST(RunConfig, SaveAndLoad) {
  const auto dir = std::filesystem::temp_directory_path() / "amcset_cli_save";
  std::filesystem::create_directories(dir);
  const auto file = dir / "scenario.json";
  auto c = parse_config(alloy_json, dir);
  save_config(c, file);
  const auto d = load_config(file);
  EXPECT_EQ(to_json(d), to_json(c));
}

// ---------------------------------------------------------------------------
// Command line
// ---------------------------------------------------------------------------

TEST(Cli, InterruptsSetTheFlagAndGiveTheSignalsStatus) {
  // 128 + the signal number, as POSIX shells report a program ended by it.
  EXPECT_EQ(cli::interrupted_exit_status(SIGINT), 130);
  EXPECT_EQ(cli::interrupted_exit_status(SIGTERM), 143);
  std::atomic<bool> flag{true};
  {
    const cli::InterruptGuard guard(flag);
    EXPECT_FALSE(flag.load()); // reset by the guard
    EXPECT_EQ(cli::interrupt_signal(), 0);
    EXPECT_EQ(cli::exit_after_interrupt(0), 0); // not interrupted: the status as it is
#ifndef _WIN32
    // The handler runs in this thread: it only sets the flag and the signal.
    std::raise(SIGTERM);
    EXPECT_TRUE(flag.load());
    EXPECT_EQ(cli::interrupt_signal(), SIGTERM);
    // Another status than the signal's (e.g. an error after the interrupt) is kept.
    EXPECT_EQ(cli::exit_after_interrupt(1), 1);
#endif
  }
}

TEST(Cli, ParsesCommandsAndOptions) {
  const auto o = cli::parse_arguments(std::vector<std::string>{
      "run", "scenario.json", "-n", "500", "-j", "4", "--seed", "7", "--damage", "full",
      "--stopping", "srim", "--stopping-table", "t.amcset", "--set", "ion.energy=3 MeV",
      "--set", "ions=9", "-o", "out", "--json", "-q"});
  EXPECT_EQ(o.command, "run");
  ASSERT_EQ(o.positional.size(), 1u);
  EXPECT_EQ(o.positional[0], "scenario.json");
  EXPECT_EQ(o.ions, std::optional<std::size_t>(500));
  EXPECT_EQ(o.threads, std::optional<unsigned>(4));
  EXPECT_EQ(o.seed, std::optional<std::uint32_t>(7));
  EXPECT_EQ(o.damage, std::optional<std::string>("full"));
  EXPECT_EQ(o.stopping_table, std::optional<std::string>("t.amcset"));
  EXPECT_EQ(o.overrides.size(), 2u);
  EXPECT_TRUE(o.json);
  EXPECT_TRUE(o.quiet);
  EXPECT_THROW(cli::parse_arguments(std::vector<std::string>{"run", "--bogus"}), cli::UsageError);
  EXPECT_THROW(cli::parse_arguments(std::vector<std::string>{"run", "-n"}), cli::UsageError);
  EXPECT_THROW(cli::parse_arguments(std::vector<std::string>{"run", "-n", "-3"}), cli::UsageError);
}

TEST(Cli, OptionsOverrideTheConfiguration) {
  auto c = parse_config(alloy_json);
  const auto o = cli::parse_arguments(std::vector<std::string>{
      "run", "x.json", "-n", "12", "--seed", "99", "--damage", "kp", "--stopping", "bethe-bloch",
      "--set", "ion.energy=5 MeV", "-o", "elsewhere"});
  cli::apply_options(c, o);
  EXPECT_EQ(c.trim.ion_count, 12u);
  EXPECT_EQ(c.trim.seed, 99u);
  EXPECT_EQ(c.trim.damage, trim::DamageCalculation::kinchin_pease);
  EXPECT_EQ(c.trim.stopping_model, physics::ElectronicStoppingModel::bethe_bloch);
  EXPECT_NEAR(in_units(c.trim.ion.energy, mega_electron_volt), 5.0, 1e-12);
  EXPECT_EQ(c.trim.output_directory, "elsewhere");
}

TEST(Cli, CompositionSpecifications) {
  const auto a = cli::parse_composition("Fe:0.91,Cr:0.09");
  ASSERT_EQ(a.size(), 2u);
  EXPECT_EQ(a[0].first, 26u);
  EXPECT_NEAR(a[0].second, 0.91, 1e-12);
  const auto b = cli::parse_composition("Fe:3,Ni:1"); // normalized
  EXPECT_NEAR(b[0].second, 0.75, 1e-12);
  const auto c = cli::parse_composition("W");
  ASSERT_EQ(c.size(), 1u);
  EXPECT_EQ(c[0].first, 74u);
  EXPECT_THROW(cli::parse_composition("Fe:-1"), std::invalid_argument);
  EXPECT_THROW(cli::parse_composition("Zz"), std::invalid_argument);
}

TEST(Cli, RecognizesTrimInFiles) {
  const auto dir = std::filesystem::temp_directory_path() / "amcset_cli_trimin";
  std::filesystem::create_directories(dir);
  {
    std::ofstream(dir / "TRIM.IN") << "==> SRIM-2013.00 This file controls TRIM Calculations.\r\n";
    std::ofstream(dir / "a.json") << "{}";
  }
  EXPECT_TRUE(cli::looks_like_trim_in((dir / "TRIM.IN").string()));
  EXPECT_FALSE(cli::looks_like_trim_in((dir / "a.json").string()));
}

TEST(Cli, StoppingTableOutput) {
  const auto o = cli::parse_arguments(std::vector<std::string>{
      "stopping", "--ion", "H", "--target", "Fe", "--from", "100 keV", "--to", "10 MeV",
      "--per-decade", "1", "--csv"});
  std::ostringstream out;
  EXPECT_EQ(cli::print_stopping_table(o, out), 0);
  const auto text = out.str();
  // Header plus three energies (100 keV, 1 MeV, 10 MeV).
  std::size_t rows = 0;
  std::istringstream lines(text);
  for (std::string line; std::getline(lines, line);) {
    if (!line.empty() && line[0] != '#' && std::isdigit(static_cast<unsigned char>(line[0]))) ++rows;
  }
  EXPECT_EQ(rows, 3u);
  // Units option.
  const auto p = cli::parse_arguments(std::vector<std::string>{
      "stopping", "--ion", "He", "--target", "Fe:0.9,Cr:0.1", "--density", "7.8 g/cm3",
      "--units", "MeV/(mg/cm2)", "--from", "1 MeV", "--to", "1 MeV"});
  std::ostringstream out2;
  EXPECT_EQ(cli::print_stopping_table(p, out2), 0);
  EXPECT_NE(out2.str().find("MeV/(mg/cm2)"), std::string::npos);
}

TEST(Cli, FastMathFirstIonAndLogOptions) {
  auto c = parse_config(alloy_json);
  EXPECT_FALSE(c.trim.fast_math);
  const auto o = cli::parse_arguments(std::vector<std::string>{
      "run", "x.json", "--fast", "--log", "run.log", "--set", "first_ion=32"});
  EXPECT_TRUE(o.fast_math);
  EXPECT_EQ(o.log_file, "run.log");
  cli::apply_options(c, o);
  EXPECT_TRUE(c.trim.fast_math);
  EXPECT_EQ(c.trim.first_ion, 32u);
  // Both keys survive the JSON round trip; first_ion is written only when set.
  const auto back = parse_config(to_json(c));
  EXPECT_TRUE(back.trim.fast_math);
  EXPECT_EQ(back.trim.first_ion, 32u);
  c.trim.first_ion = 0;
  EXPECT_EQ(to_json(c).find("first_ion"), std::string::npos);
  EXPECT_NE(to_json(c).find("\"fast_math\": true"), std::string::npos);
  EXPECT_THROW(parse_config(R"({"physics": {"fast_math": "fast"}})"), std::exception);
}

TEST(Cli, RandomSeedKeyAndOption) {
  auto c = parse_config(alloy_json);
  EXPECT_FALSE(c.trim.random_seed);
  apply_override(c, "seed=random");
  EXPECT_TRUE(c.trim.random_seed);
  EXPECT_NE(to_json(c).find("\"seed\": \"random\""), std::string::npos);
  EXPECT_TRUE(parse_config(to_json(c)).trim.random_seed);
  // A number sets a fixed seed again; --seed takes both.
  cli::apply_options(c, cli::parse_arguments(std::vector<std::string>{"run", "x.json", "--seed", "42"}));
  EXPECT_FALSE(c.trim.random_seed);
  EXPECT_EQ(c.trim.seed, 42u);
  cli::apply_options(c, cli::parse_arguments(std::vector<std::string>{"run", "x.json", "-s", "random"}));
  EXPECT_TRUE(c.trim.random_seed);
  EXPECT_THROW(cli::parse_arguments(std::vector<std::string>{"run", "x.json", "-s", "sometimes"}), cli::UsageError);
  EXPECT_THROW(parse_config(R"({"seed": "often"})"), std::invalid_argument);
}

TEST(Cli, RecoilOptionAndRepeatKey) {
  auto c = parse_config(alloy_json);
  const auto o = cli::parse_arguments(
      std::vector<std::string>{"run", "x.json", "-n", "25", "--recoil", "Fe,10keV,50nm,0,0.6,0.8"});
  cli::apply_options(c, o);
  EXPECT_EQ(c.trim.source, trim::TrimConfig::Source::particle_recoils);
  ASSERT_EQ(c.trim.particles.size(), 1u);
  const auto &p = c.trim.particles.front();
  EXPECT_EQ(p.z, 26u);
  EXPECT_DOUBLE_EQ(in_units(p.energy, electron_volt), 10000.0);
  EXPECT_DOUBLE_EQ(in_units(p.x, angstrom), 500.0);
  EXPECT_DOUBLE_EQ(p.cos_y, 0.6);
  EXPECT_EQ(c.trim.particle_repeat, 25u);
  EXPECT_EQ(c.trim.ion_count, 0u);
  EXPECT_EQ(c.trim.particle_count(), 25u);
  EXPECT_GE(in_units(c.trim.ion.energy, electron_volt), 10000.0);
  // The repeat survives the JSON round trip (written only when not 1).
  const auto back = parse_config(to_json(c));
  EXPECT_EQ(back.trim.particle_repeat, 25u);
  ASSERT_EQ(back.trim.particles.size(), 1u);
  EXPECT_EQ(back.trim.particles.front().name, "PKA");
  EXPECT_EQ(to_json(parse_config(alloy_json)).find("repeat"), std::string::npos);
  for (const char *bad : {"Fe,10keV", "Xx,10keV,5nm", "Fe,10keV,5nm,1,0", "Fe,10keV,5nm,a,b,c"}) {
    auto d = parse_config(alloy_json);
    EXPECT_ANY_THROW(cli::apply_options(d, cli::parse_arguments(std::vector<std::string>{"run", "x.json", "--recoil", bad})))
        << bad;
  }
  const std::string base = R"({"ion": {"element": "Fe", "energy": "20 keV"},
    "layers": [{"thickness": "1 um", "elements": [{"element": "Fe"}]}], )";
  EXPECT_THROW(parse_config(base + R"("source": {"type": "recoils", "repeat": 0, "particles": []}})"),
               std::invalid_argument);
}

TEST(Cli, AutoDepthOptionAndKeys) {
  auto c = parse_config(alloy_json);
  EXPECT_FALSE(c.trim.auto_depth.enabled);
  // Defaults are not written.
  EXPECT_EQ(to_json(c).find("auto_depth"), std::string::npos);
  const auto o = cli::parse_arguments(std::vector<std::string>{"run", "x.json", "--auto-depth"});
  EXPECT_TRUE(o.auto_depth);
  cli::apply_options(c, o);
  EXPECT_TRUE(c.trim.auto_depth.enabled);
  apply_override(c, "auto_depth.test_ions=500");
  apply_override(c, "auto_depth.margin=0.2");
  apply_override(c, "auto_depth.shrink_below=0.6");
  apply_override(c, "auto_depth.warn_transmitted=false");
  const auto back = parse_config(to_json(c));
  EXPECT_TRUE(back.trim.auto_depth.enabled);
  EXPECT_EQ(back.trim.auto_depth.test_ions, 500u);
  EXPECT_EQ(back.trim.auto_depth.margin, 0.2);
  EXPECT_EQ(back.trim.auto_depth.shrink_below, 0.6);
  EXPECT_FALSE(back.trim.auto_depth.warn_transmitted);
  const std::string base = R"({"ion": {"element": "Fe", "energy": "1 MeV"},
    "layers": [{"thickness": "1 um", "elements": [{"element": "Fe"}]}], )";
  EXPECT_NO_THROW(parse_config(base + R"("auto_depth": {"enabled": true}})"));
  EXPECT_THROW(parse_config(base + R"("auto_depth": {"margin": 0.9}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("auto_depth": {"test_ions": 0}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("auto_depth": {"shrink_below": 1.2}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("auto_depth": {"test_ions": 1e9}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("auto_depth": {"depth": "1 um"}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("auto_depth": {"enabled": "maybe"}})"), std::exception);
}

TEST(Cli, NickelStoppingCorrection) {
  auto c = parse_config(alloy_json);
  EXPECT_FALSE(c.trim.stopping_corrections.nickel_ions);
  EXPECT_NE(to_json(c).find("\"nickel_ions\": false"), std::string::npos) << to_json(c);
  const auto o = cli::parse_arguments(std::vector<std::string>{"run", "x.json", "--correct-nickel"});
  EXPECT_TRUE(o.correct_nickel);
  cli::apply_options(c, o);
  EXPECT_TRUE(c.trim.stopping_corrections.nickel_ions);
  EXPECT_TRUE(parse_config(to_json(c)).trim.stopping_corrections.nickel_ions);
  auto d = parse_config(alloy_json);
  apply_override(d, "physics.stopping_corrections.nickel_ions=true");
  EXPECT_TRUE(d.trim.stopping_corrections.nickel_ions);
  const std::string base = R"({"ion": {"element": "Ni", "energy": "5 MeV"},
    "layers": [{"thickness": "3 um", "elements": [{"element": "Fe"}]}], )";
  EXPECT_THROW(parse_config(base + R"("physics": {"stopping_corrections": {"cobalt_ions": true}}})"),
               std::invalid_argument);
  // The stopping command: Ni in Fe at 3.5 MeV (0.06 MeV/u) stops much more with it.
  const auto table = [](bool fix) {
    std::vector<std::string> args{"stopping", "--ion", "Ni", "--target", "Fe", "--from", "3.5 MeV", "--to",
                                  "3.5 MeV", "--csv"};
    if (fix) args.push_back("--correct-nickel");
    std::ostringstream out;
    cli::print_stopping_table(cli::parse_arguments(args), out);
    std::istringstream in(out.str());
    std::string line;
    while (std::getline(in, line)) {
      if (!line.empty() && line[0] != '#') return std::stod(line.substr(line.find(',') + 1));
    }
    return 0.0;
  };
  EXPECT_GT(table(true), 1.4 * table(false));
}

// Malformed numbers in a scenario are errors, never silently truncated or
// wrapped values.
TEST(RunConfig, RejectsMalformedNumbers) {
  const std::string base = R"({"ion": {"element": "Fe", "energy": "1 MeV"},
    "layers": [{"thickness": "1 um", "elements": [{"element": "Fe"}]}])";
  const auto with = [&](const std::string &extra) { return base + "," + extra + "}"; };
  EXPECT_NO_THROW(parse_config(base + "}"));
  // 2^32 + 26 would wrap to Z = 26 in a 32-bit integer.
  EXPECT_THROW(parse_config(R"({"ion": {"z": 4294967322, "energy": "1 MeV"},
    "layers": [{"thickness": "1 um", "elements": [{"element": "Fe"}]}]})"), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("ions": "12abc")")), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("ions": "nan")")), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("ions": 1e30)")), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("ions": -3)")), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("seed": 4294967296)")), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("threads": 100000)")), std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("output": {"bins": 1e12})")), std::invalid_argument);
  EXPECT_THROW(parse_config(R"({"ion": {"element": "Fe", "energy": "inf MeV"},
    "layers": [{"thickness": "1 um", "elements": [{"element": "Fe"}]}]})"), std::invalid_argument);
  // A table's file name cannot leave the output folder.
  EXPECT_THROW(parse_config(with(R"("output": {"files": {"tables": {"range": {"file_name": "../../evil"}}}})")),
               std::invalid_argument);
  EXPECT_THROW(parse_config(with(R"("output": {"files": {"tables": {"range": {"file_name": "C:\x"}}}})")),
               std::invalid_argument);
  EXPECT_NO_THROW(parse_config(with(R"("output": {"files": {"tables": {"range": {"file_name": "my_range"}}}})")));
  EXPECT_THROW(parse_config(with(R"("output": {"files": {"number_format": {"digits": 1e9}}})")), std::invalid_argument);
}

TEST(Cli, TrajectoryAndMdEventOptionsAndKeys) {
  auto c = parse_config(alloy_json);
  const auto o = cli::parse_arguments(std::vector<std::string>{
      "run", "x.json", "--trajectory-file", "traj.parquet", "--trajectory-points", "births", "--md-events",
      "events.csv", "--md-ec", "0.5,1,2 keV", "--md-emin", "30 eV", "--md-handoff"});
  cli::apply_options(c, o);
  EXPECT_EQ(c.trim.trajectory_export.file, "traj.parquet");
  EXPECT_TRUE(c.trim.trajectory_export.births_only);
  EXPECT_EQ(c.trim.trajectory_ions, trim::TrimConfig::all_ions); // the file implies recorded paths
  EXPECT_EQ(c.trim.md_events.file, "events.csv");
  ASSERT_EQ(c.trim.md_events.transition_energies.size(), 3u);
  EXPECT_NEAR(in_units(c.trim.md_events.transition_energies[0], electron_volt), 500.0, 1e-9);
  EXPECT_NEAR(in_units(c.trim.md_events.transition_energies[2], electron_volt), 2000.0, 1e-9);
  EXPECT_NEAR(in_units(c.trim.md_events.min_energy, electron_volt), 30.0, 1e-12);
  EXPECT_TRUE(c.trim.md_events.handoff);
  // Every key survives the JSON round trip ("all", "births", "Ed", the list).
  c.trim.md_events.min_energy = 0.0 * electron_volt;
  const std::string j = to_json(c);
  EXPECT_NE(j.find("\"trajectories\": \"all\""), std::string::npos) << j.substr(0, 2000);
  EXPECT_NE(j.find("\"E_min\": \"Ed\""), std::string::npos);
  const auto back = parse_config(j);
  EXPECT_EQ(back.trim.trajectory_ions, trim::TrimConfig::all_ions);
  EXPECT_TRUE(back.trim.trajectory_export.births_only);
  EXPECT_EQ(back.trim.md_events.transition_energies.size(), 3u);
  EXPECT_TRUE(back.trim.md_events.handoff);
  EXPECT_EQ(in_units(back.trim.md_events.min_energy, electron_volt), 0.0);
  // Defaults are not written.
  EXPECT_EQ(to_json(parse_config(alloy_json)).find("md_events"), std::string::npos);
  // Malformed values.
  const std::string base = R"({"ion": {"element": "Fe", "energy": "1 MeV"},
    "layers": [{"thickness": "1 um", "elements": [{"element": "Fe"}]}], "damage": "full", )";
  EXPECT_NO_THROW(parse_config(base + R"("output": {"md_events": {"file": "e.csv", "E_c": "1 keV"}}})"));
  EXPECT_THROW(parse_config(base + R"("output": {"md_events": {"file": "e.csv", "E_c": [1000]}}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("output": {"md_events": {"file": "e.csv", "mode": "both"}}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("output": {"md_events": {"file": "e.csv", "Ec": "1 keV"}}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("output": {"trajectories": "some"}})"), std::invalid_argument);
  EXPECT_THROW(parse_config(base + R"("output": {"trajectory_file": {"file": "t.csv", "points": "few"}}})"),
               std::invalid_argument);
  EXPECT_THROW(cli::apply_options(c, cli::parse_arguments(std::vector<std::string>{"run", "x", "--trajectory-points", "x"})),
               cli::UsageError);
}

TEST(RunConfig, EnergyLists) {
  const auto v = parse_energy_list("0.5, 1, 2 keV", "E_c");
  ASSERT_EQ(v.size(), 3u);
  EXPECT_NEAR(in_units(v[0], electron_volt), 500.0, 1e-9);
  EXPECT_NEAR(in_units(v[1], electron_volt), 1000.0, 1e-9);
  const auto w = parse_energy_list("500 eV,1e3 eV, 2 keV", "E_c");
  ASSERT_EQ(w.size(), 3u);
  EXPECT_NEAR(in_units(w[1], electron_volt), 1000.0, 1e-9);
  EXPECT_NEAR(in_units(w[2], electron_volt), 2000.0, 1e-9);
  EXPECT_THROW(parse_energy_list("1, 2", "E_c"), std::invalid_argument);   // no unit at all
  EXPECT_THROW(parse_energy_list("1 keV,,2 keV", "E_c"), std::invalid_argument);
  EXPECT_THROW(parse_energy_list("1 nm", "E_c"), std::invalid_argument);
}

// The table for LAMMPS fix electron/stopping: energy (eV) and one N*S_e
// column (eV/A) per atom type, energies increasing.
TEST(Cli, LammpsStoppingTable) {
  const auto o = cli::parse_arguments(std::vector<std::string>{
      "stopping", "--lammps", "--ion", "Fe,Cr", "--target", "Fe:0.9,Cr:0.1", "--density", "7.8 g/cm3",
      "--from", "100 eV", "--to", "1 MeV", "--per-decade", "5", "--model", "builtin"});
  std::ostringstream out;
  EXPECT_EQ(cli::print_stopping_table(o, out), 0);
  std::istringstream in(out.str());
  std::string line;
  int rows = 0;
  double last = 0.0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream ls(line);
    double e = 0, a = 0, b = 0, extra = 0;
    ls >> e >> a >> b;
    EXPECT_FALSE(ls.fail()) << line;
    EXPECT_FALSE(ls >> extra) << "three columns: " << line;
    EXPECT_GT(e, last);
    EXPECT_GT(a, 0.0);
    EXPECT_GT(b, 0.0);
    last = e;
    ++rows;
  }
  EXPECT_EQ(rows, 21); // 4 decades x 5 + 1
  EXPECT_NEAR(last, 1e6, 1e-3);
  EXPECT_NE(out.str().find("metal units"), std::string::npos);
  // --ntypes repeats a single ion's column.
  const auto r = cli::parse_arguments(std::vector<std::string>{
      "stopping", "--lammps", "--ion", "Fe", "--ntypes", "3", "--target", "Fe", "--from", "1 keV", "--to", "1 keV",
      "--model", "builtin"});
  std::ostringstream o3;
  cli::print_stopping_table(r, o3);
  std::string data;
  std::istringstream i3(o3.str());
  while (std::getline(i3, data)) {
    if (!data.empty() && data[0] != '#') break;
  }
  std::istringstream f(data);
  double e = 0, s1 = 0, s2 = 0, s3 = 0;
  f >> e >> s1 >> s2 >> s3;
  EXPECT_EQ(s1, s2);
  EXPECT_EQ(s2, s3);
  EXPECT_THROW(cli::print_stopping_table(cli::parse_arguments(std::vector<std::string>{
                   "stopping", "--ion", "Fe,Cr", "--target", "Fe"}), out), cli::UsageError);
}

// --fix-trim-bugs, then --set re-enables one bug (Eb = 0 for MD coupling).
TEST(Cli, SetAppliesAfterFixTrimBugs) {
  auto c = parse_config(alloy_json);
  cli::apply_options(c, cli::parse_arguments(std::vector<std::string>{
                            "run", "x.json", "--fix-trim-bugs", "--set",
                            "physics.trim_bugs.zero_energy_defaults_skip_first_atom=true"}));
  EXPECT_TRUE(c.trim.bugs.zero_energy_defaults_skip_first_atom);
  EXPECT_FALSE(c.trim.bugs.stale_primary_energy);
  EXPECT_FALSE(c.trim.bugs.recoil_energy_created);
}

TEST(RunConfig, DamageEfficiencyAndDpaFluence) {
  auto c = parse_config(R"({
    "ion": {"element": "Fe", "energy": "100 keV"},
    "layers": [{"name": "Si", "thickness": "1 um", "density": "2.33 g/cm3",
                "elements": [{"element": "Si",
                              "damage_efficiency": {"arc_b": -0.6, "arc_c": 0.3, "rpa_b": "2 keV", "rpa_c": 0.9}}]}],
    "output": {"dpa_fluence": "1e19 /m2"}
  })");
  const auto &de = c.trim.layers[0].elements[0].damage_efficiency;
  ASSERT_TRUE(de.has_value());
  EXPECT_EQ(de->arc_b, -0.6);
  EXPECT_NEAR(in_units(de->rpa_b, electron_volt), 2000.0, 1e-9);
  EXPECT_NEAR(c.trim.dpa_fluence_per_cm2, 1e15, 1.0);
  const auto b = parse_config(to_json(c));
  ASSERT_TRUE(b.trim.layers[0].elements[0].damage_efficiency.has_value());
  EXPECT_EQ(b.trim.layers[0].elements[0].damage_efficiency->rpa_c, 0.9);
  EXPECT_NEAR(b.trim.dpa_fluence_per_cm2, 1e15, 1.0);
  EXPECT_EQ(to_json(b), to_json(c));
  // All four constants are needed; the unit of the fluence is checked.
  EXPECT_THROW(parse_config(R"({"ion": {"element": "Fe", "energy": "100 keV"},
    "layers": [{"name": "Si", "thickness": "1 um", "density": "2.33 g/cm3",
                "elements": [{"element": "Si", "damage_efficiency": {"arc_b": -0.6}}]}]})"),
               std::invalid_argument);
  EXPECT_THROW(parse_config(R"({"ion": {"element": "Fe", "energy": "100 keV"},
    "layers": [{"name": "Si", "thickness": "1 um", "density": "2.33 g/cm3", "elements": [{"element": "Si"}]}],
    "output": {"dpa_fluence": "1e15 /mm2"}})"),
               std::invalid_argument);
  apply_override(c, "output.dpa_fluence=3e14");
  EXPECT_NEAR(c.trim.dpa_fluence_per_cm2, 3e14, 1.0);
}

// ---------------------------------------------------------------------------
// Crystals: the layers' "crystal" object and the cell command
// ---------------------------------------------------------------------------

namespace {

const char *si_crystal_json = R"({
  "ion": {"element": "B", "energy": "30 keV", "angle": "7 deg"},
  "layers": [{"thickness": "1 um",
              "crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
                          "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0], "twist": "22 deg"},
                          "temperature": "27 C",
                          "vibrations": [{"element": "Si", "debye_temperature": "645 K"}]}}],
  "ions": 20, "threads": 2, "stopping": {"model": "builtin"}
})";

//! A scenario with one layer whose members are \p layer (JSON text).
std::string one_layer(const std::string &layer) {
  return R"({"ion": {"element": "He", "energy": "100 keV"}, "layers": [{"thickness": "1 um", )" + layer + "}]}";
}

} // namespace

TEST(RunConfig, CrystalFromABuiltInStructure) {
  const auto c = parse_config(si_crystal_json);
  const auto &l = c.trim.layers[0];
  ASSERT_TRUE(l.crystal.has_value());
  EXPECT_EQ(l.name, "Si");
  // The layer's elements, stoichiometry and density are the unit cell's.
  ASSERT_EQ(l.elements.size(), 1u);
  EXPECT_EQ(l.elements[0].z, 14u);
  EXPECT_DOUBLE_EQ(l.elements[0].stoichiometry, 8.0);
  EXPECT_NEAR(in_units(l.density, gram_per_cubic_centimeter), 2.329, 0.005); // silicon
  EXPECT_EQ(l.crystal->cell.sites.size(), 8u);
  EXPECT_NEAR(in_units(l.crystal->cell.lattice.a(), angstrom), 5.431, 1e-12);
  EXPECT_EQ(l.crystal->orientation.reference, (std::vector<double>{1, 1, 0}));
  EXPECT_NEAR(in_units(l.crystal->orientation.twist, degree), 22.0, 1e-12);
  EXPECT_NEAR(in_units(l.crystal->temperature, kelvin), 300.15, 1e-9);
  const auto u = l.crystal->rms_displacement(14, l.elements[0].mass);
  ASSERT_TRUE(u.has_value());
  EXPECT_EQ(*u, crystal::debye_rms_displacement(l.elements[0].mass, 645.0 * kelvin, 300.15 * kelvin));
  EXPECT_NO_THROW(c.trim.validate());
  // The canonical form (lattice vectors and sites written out) reads back to itself.
  const std::string text = to_json(c);
  EXPECT_NE(text.find("\"vectors_A\""), std::string::npos);
  EXPECT_NE(text.find("\"surface\": [0,0,1]"), std::string::npos);
  EXPECT_NE(text.find("\"element\": \"Si\""), std::string::npos); // one line per fully occupied site
  const auto b = parse_config(text);
  EXPECT_EQ(to_json(b), text);
  ASSERT_TRUE(b.trim.layers[0].crystal.has_value());
  EXPECT_TRUE(*b.trim.layers[0].crystal == *l.crystal);
  EXPECT_EQ(b.trim.layers[0].density, l.density);

  // No transport follows the lattice yet: the run is the amorphous layer's,
  // number for number, and a warning says so.
  auto amorphous = c;
  amorphous.trim.layers[0].crystal.reset();
  bool warned = false;
  for (const auto &w : c.trim.plausibility_warnings()) warned = warned || w.find("crystal") != std::string::npos;
  EXPECT_TRUE(warned);
  for (const auto &w : amorphous.trim.plausibility_warnings()) EXPECT_EQ(w.find("crystal"), std::string::npos) << w;
  auto e1 = make_engine(c);
  auto e2 = make_engine(amorphous);
  e1->run();
  e2->run();
  const auto &t1 = dynamic_cast<const trim::TrimSimulation &>(*e1).tallies();
  const auto &t2 = dynamic_cast<const trim::TrimSimulation &>(*e2).tallies();
  EXPECT_EQ(t1.ions, 20u);
  EXPECT_EQ(t1.range.mean(), t2.range.mean());
  EXPECT_EQ(t1.total_vacancies, t2.total_vacancies);
  EXPECT_EQ(t1.total_ionization_ion, t2.total_ionization_ion);
}

TEST(RunConfig, CrystalLatticeConstantFromTheDensity) {
  // A substitutional alloy with interstitial carbon; no lattice constant.
  const auto c = parse_config(one_layer(R"("density": "7.8 g/cm3",
      "elements": [{"element": "Fe", "Ed": "40 eV"}, {"element": "Cr", "Ed": "28 eV"}, {"element": "C"}],
      "crystal": {"structure": "bcc", "sublattices": [{"Fe": 0.9, "Cr": 0.1}],
                  "interstitials": [{"sites": "octahedral", "occupants": {"C": 0.01}}]})"));
  const auto &l = c.trim.layers[0];
  EXPECT_NEAR(in_units(l.density, gram_per_cubic_centimeter), 7.8, 1e-9);
  ASSERT_EQ(l.elements.size(), 3u);
  EXPECT_NEAR(l.elements[0].stoichiometry, 1.8, 1e-12);  // 2 sites x 0.9
  EXPECT_NEAR(l.elements[1].stoichiometry, 0.2, 1e-12);
  EXPECT_NEAR(l.elements[2].stoichiometry, 0.06, 1e-12); // 6 octahedral sites x 0.01
  EXPECT_NEAR(in_units(l.elements[0].displacement_energy, electron_volt), 40.0, 1e-12); // kept
  EXPECT_NEAR(in_units(l.elements[1].displacement_energy, electron_volt), 28.0, 1e-12);
  // a^3 = (mass of the cell) / density, with 1 amu = 1.66053907e-24 g.
  const double mass = 1.8 * in_units(Element::by_z(26).natural_mass(), atomic_mass_unit) +
                      0.2 * in_units(Element::by_z(24).natural_mass(), atomic_mass_unit) +
                      0.06 * in_units(Element::by_z(6).natural_mass(), atomic_mass_unit);
  EXPECT_NEAR(in_units(l.crystal->cell.lattice.a(), angstrom), std::cbrt(mass * 1.66053907 / 7.8), 1e-6);
  EXPECT_EQ(l.crystal->cell.sites.size(), 8u); // 2 lattice + 6 interstitial
  EXPECT_NO_THROW(c.trim.validate());
  EXPECT_EQ(to_json(parse_config(to_json(c))), to_json(c));
  // A pure element without a density takes the element's: bcc W, a = 3.16 A.
  const auto w = parse_config(one_layer(R"("crystal": {"structure": "bcc", "sublattices": ["W"]})"));
  EXPECT_NEAR(in_units(w.trim.layers[0].crystal->cell.lattice.a(), angstrom), 3.16, 0.01);
  EXPECT_NEAR(in_units(w.trim.layers[0].density, gram_per_cubic_centimeter),
              in_units(Element::by_z(74).solid_density(), gram_per_cubic_centimeter), 1e-9);
  // A given lattice constant decides the density; a density within 2 % is accepted and replaced.
  const auto fe = parse_config(one_layer(R"("density": "7.8658 g/cm3",
      "crystal": {"structure": "bcc", "a": "2.8665 A", "sublattices": ["Fe"]})"));
  const double rho = 2.0 * in_units(Element::by_z(26).natural_mass(), atomic_mass_unit) * 1.66053907 /
                     std::pow(2.8665, 3);
  EXPECT_NEAR(in_units(fe.trim.layers[0].density, gram_per_cubic_centimeter), rho, 1e-6);
  EXPECT_NEAR(rho, 7.874, 0.005);
}

TEST(RunConfig, CrystalFromSitesAndFiles) {
  // A lattice with sites; the elements are left out and come from the cell.
  const auto c = parse_config(one_layer(R"("crystal": {
      "name": "B2 NiAl", "lattice": {"a": "2.88 A"},
      "sites": [{"position": [0, 0, 0], "element": "Ni"},
                {"label": "Al", "position": [1.5, 0.5, -0.5], "occupants": [{"element": "Al", "fraction": 0.95}]},
                {"position": [0.5, 0.5, 0], "kind": "interstitial", "occupants": []}],
      "orientation": {"euler_deg": [0, 45, 0]}})"));
  const auto &l = c.trim.layers[0];
  ASSERT_EQ(l.elements.size(), 2u);
  EXPECT_EQ(l.name, "AlNi"); // in the order of Z
  EXPECT_EQ(l.crystal->cell.name, "B2 NiAl");
  EXPECT_DOUBLE_EQ(l.crystal->cell.atoms_per_cell(), 1.95);
  EXPECT_EQ(l.crystal->cell.sites[0].label, "Ni1");
  EXPECT_EQ(l.crystal->cell.sites[1].position, (crystal::Vec3{0.5, 0.5, 0.5})); // wrapped into the cell
  EXPECT_EQ(l.crystal->cell.sites[2].kind, crystal::SiteKind::interstitial);
  EXPECT_TRUE(l.crystal->cell.sites[2].occupants.empty());
  EXPECT_EQ(l.crystal->orientation.kind, crystal::OrientationKind::euler);
  EXPECT_NO_THROW(c.trim.validate());
  EXPECT_EQ(to_json(parse_config(to_json(c))), to_json(c));
  // Overrides reach the crystal.
  auto o = c;
  apply_override(o, "layers.0.crystal.temperature=77 K");
  apply_override(o, "layers.0.crystal.orientation.twist=15 deg");
  EXPECT_NEAR(in_units(o.trim.layers[0].crystal->temperature, kelvin), 77.0, 1e-12);
  EXPECT_NEAR(in_units(o.trim.layers[0].crystal->orientation.twist, degree), 15.0, 1e-12);

  // A cell file, relative to the scenario's directory.
  const auto dir = std::filesystem::temp_directory_path() / "amcset_test_crystal_scenario";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  {
    std::ofstream f(dir / "POSCAR", std::ios::binary);
    f << "GaAs\n5.6533\n 1 0 0\n 0 1 0\n 0 0 1\n Ga As\n 4 4\nDirect\n0 0 0\n0.5 0.5 0\n0.5 0 0.5\n0 0.5 0.5\n"
         "0.25 0.25 0.25\n0.75 0.75 0.25\n0.75 0.25 0.75\n0.25 0.75 0.75\n";
  }
  const auto g = parse_config(one_layer(R"("crystal": {"file": "POSCAR", "format": "poscar"})"), dir);
  ASSERT_EQ(g.trim.layers[0].elements.size(), 2u);
  EXPECT_EQ(g.trim.layers[0].crystal->cell.sites.size(), 8u);
  EXPECT_NEAR(in_units(g.trim.layers[0].density, gram_per_cubic_centimeter), 5.32, 0.01); // GaAs
  EXPECT_THROW(parse_config(one_layer(R"("crystal": {"file": "missing.cif"})"), dir), std::runtime_error);
  EXPECT_THROW(parse_config(one_layer(R"("crystal": {"file": "POSCAR", "format": "pdb"})"), dir),
               std::invalid_argument);
  std::filesystem::remove_all(dir);

  // The crystal object alone (as the cell command prints it).
  const auto d = parse_crystal(R"({"structure": "fcc", "sublattices": ["Cu"]})");
  EXPECT_TRUE(d.lattice_from_density);
  EXPECT_EQ(d.crystal.cell.sites.size(), 4u);
  const auto again = parse_crystal(crystal_to_json(d.crystal));
  EXPECT_FALSE(again.lattice_from_density);
  EXPECT_TRUE(again.crystal == d.crystal);
  EXPECT_THROW(parse_crystal("[1, 2]"), std::invalid_argument);
  EXPECT_THROW(parse_crystal("{oops"), std::invalid_argument);
}

TEST(RunConfig, CrystalErrorsNameTheKey) {
  const auto message = [](const std::string &layer) -> std::string {
    try {
      const auto c = parse_config(one_layer(layer));
      c.trim.validate();
    } catch (const std::invalid_argument &e) {
      return e.what();
    }
    return "";
  };
  const auto expect = [&](const std::string &layer, const std::string &part) {
    const std::string m = message(layer);
    EXPECT_NE(m.find(part), std::string::npos) << layer << "\n -> " << (m.empty() ? "(accepted)" : m);
  };
  // Which description is meant.
  expect(R"("crystal": {})", "exactly one of");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "file": "x.cif"})", "exactly one of");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "format": "cif"})", "layers.0.crystal.format");
  expect(R"("crystal": {"file": "x.cif", "a": "3 A"})", "layers.0.crystal.a");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "colour": 1})", "layers.0.crystal.colour");
  // Built-in structures.
  expect(R"("crystal": {"structure": "bbc", "sublattices": ["Fe"]})", "layers.0.crystal: unknown crystal structure");
  expect(R"("crystal": {"structure": "bcc"})", "layers.0.crystal.sublattices");
  expect(R"("crystal": {"structure": "rocksalt", "sublattices": ["Mg"]})", "2 sublattice");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Qq"]})", "layers.0.crystal.sublattices.0");
  expect(R"("crystal": {"structure": "bcc", "sublattices": [{"Fe": 0.7, "Cr": 0.5}]})", "more than 1");
  expect(R"("crystal": {"structure": "bcc", "sublattices": [["Fe", "Fe"]]})", "listed twice");
  expect(R"("crystal": {"structure": "bcc", "sublattices": [{"Fe": 0}]})", "above 0");
  expect(R"("crystal": {"structure": "bcc", "a": "-1 A", "sublattices": ["Fe"]})", "layers.0.crystal.a");
  expect(R"("crystal": {"structure": "hcp", "c_over_a": 0, "sublattices": ["Ti"]})", "layers.0.crystal.c_over_a");
  expect(R"("crystal": {"structure": "sc", "sublattices": ["Po"], "interstitials": [{"sites": "octahedral"}]})",
         "interstitial");
  // Density and composition must agree with the cell.
  expect(R"("crystal": {"structure": "rocksalt", "sublattices": ["Mg", "O"]})", "lattice constant");
  expect(R"("density": "5 g/cm3", "crystal": {"structure": "bcc", "a": "2.8665 A", "sublattices": ["Fe"]})",
         "layers.0.density");
  expect(R"("elements": [{"element": "Fe"}, {"element": "Ni"}],
            "crystal": {"structure": "bcc", "sublattices": ["Fe"]})", "layers.0.elements.1");
  expect(R"("density": "7.8 g/cm3", "elements": [{"element": "Fe", "fraction": 0.5}, {"element": "Cr", "fraction": 0.5}],
            "crystal": {"structure": "bcc", "sublattices": [{"Fe": 0.9, "Cr": 0.1}]})",
         "layers.0.elements.0.stoichiometry");
  expect(R"("gas": true, "crystal": {"structure": "fcc", "sublattices": ["Ar"]})", "gas");
  // Lattice and sites.
  expect(R"("crystal": {"sites": [{"position": [0, 0, 0], "element": "Fe"}]})", "layers.0.crystal.lattice");
  expect(R"("crystal": {"lattice": {"a": "3 A"}})", "layers.0.crystal.sites");
  expect(R"("crystal": {"lattice": {"a": "3 A", "gamma": "200 deg"}, "sites": [{"position": [0, 0, 0], "element": "Fe"}]})",
         "layers.0.crystal.lattice");
  expect(R"("crystal": {"lattice": {"vectors_A": [[3, 0, 0], [0, 3, 0]]}, "sites": [{"position": [0, 0, 0], "element": "Fe"}]})",
         "three vectors");
  expect(R"("crystal": {"lattice": {"vectors_A": [[3, 0, 0], [0, 3, 0], [3, 3, 0]]}, "sites": [{"position": [0, 0, 0], "element": "Fe"}]})",
         "degenerate");
  expect(R"("crystal": {"lattice": {"a": "3 A"}, "sites": [{"position": [0, 0], "element": "Fe"}]})",
         "layers.0.crystal.sites.0.position");
  expect(R"("crystal": {"lattice": {"a": "3 A"}, "sites": [{"position": [0, 0, 0], "element": "Fe", "kind": "vacancy"}]})",
         "layers.0.crystal.sites.0.kind");
  expect(R"("crystal": {"lattice": {"a": "3 A"}, "sites": [{"position": [0, 0, 0], "element": "Fe"},
                                                           {"position": [0.01, 0, 0], "element": "Cr"}]})",
         "apart");
  expect(R"("crystal": {"lattice": {"a": "3 A"}, "sites": [{"position": [0, 0, 0], "occupants": []}]})", "no atoms");
  // Orientation, temperature, vibrations (checked by validate).
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "orientation": {"surface": [0, 0, 1], "reference": [0, 0, 1]}})",
         "reference");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "orientation": {"surface": [0, 0, 0]}})", "surface");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "orientation": {"surface": [0, 0, 0, 1]}})", "hexagonal");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "orientation": {"euler_deg": [0, 0, 0], "surface": [1, 0, 0]}})",
         "not both");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "temperature": "300"})", "needs a unit");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "temperature": "-5 K"})", "temperature");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "vibrations": [{"element": "Cr", "rms": "0.05 A"}]})",
         "not in the unit cell");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "vibrations": [{"element": "Fe"}]})", "debye_temperature");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "vibrations": [{"element": "Fe", "debye_temperature": "0 K"}]})",
         "debye_temperature");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"], "vibrations": [{"element": "Fe", "rms": "2 A"}]})", "rms");
  expect(R"("crystal": {"structure": "bcc", "sublattices": ["Fe"],
                        "vibrations": [{"element": "Fe", "rms": "0.05 A"}, {"element": "Fe", "rms": "0.06 A"}]})", "twice");
  // The checks hold for a crystal changed in memory (as the GUI will do).
  auto c = parse_config(si_crystal_json);
  c.trim.layers[0].crystal->cell.sites[0].occupants = {crystal::Occupant{32, 1.0}}; // Si -> Ge
  EXPECT_THROW(c.trim.validate(), std::invalid_argument); // the layer still says pure Si
  trim::apply_crystal(c.trim.layers[0]);
  EXPECT_NO_THROW(c.trim.validate());
  ASSERT_EQ(c.trim.layers[0].elements.size(), 2u);
  EXPECT_DOUBLE_EQ(c.trim.layers[0].elements[0].stoichiometry, 7.0);
  EXPECT_EQ(c.trim.layers[0].elements[1].z, 32u);
  c.trim.layers[0].density = 2.0 * gram_per_cubic_centimeter;
  EXPECT_THROW(c.trim.validate(), std::invalid_argument);
}

TEST(Cli, CellCommand) {
  const auto run = [](std::vector<std::string> args) {
    args.insert(args.begin(), "cell");
    std::ostringstream out;
    EXPECT_EQ(cli::print_cell(cli::parse_arguments(args), out), 0);
    return out.str();
  };
  // A built-in structure: bcc Fe, 8 neighbours at a sqrt(3)/2 = 2.4825 A, 6 at a.
  std::string text = run({"--structure", "bcc", "--sublattice", "Fe", "--a", "2.8665 A", "--beam", "1 1 1"});
  EXPECT_NE(text.find("2 atoms per cell"), std::string::npos) << text;
  EXPECT_NE(text.find("8 at 2.4825 A, 6 at 2.8665 A"), std::string::npos) << text;
  EXPECT_NE(text.find("depth axis (lab x) = [0 0 1]"), std::string::npos) << text;
  // [1 1 1] is acos(1/sqrt 3) = 54.735610 degrees from the [0 0 1] surface normal.
  EXPECT_NE(text.find("Beam along [1 1 1]: ion.angle = 54.73561 deg"), std::string::npos) << text;
  // The lattice constant from the element's density; a beam direction from a tilt.
  text = run({"--structure", "diamond", "--sublattice", "Si", "--surface", "0 0 1", "--reference", "0 1 0", "--tilt",
              "45 deg"});
  EXPECT_NE(text.find("lattice constant from the density"), std::string::npos) << text;
  EXPECT_NE(text.find("travels along [0 1 1]"), std::string::npos) << text; // 45 deg from [0 0 1] towards [0 1 0]
  // An alloy with interstitials, as the scenario's "crystal" object.
  text = run({"--structure", "bcc", "--sublattice", "Fe:0.9,Cr:0.1", "--interstitial", "octahedral=C:0.01", "--density",
              "7.8 g/cm3", "--json"});
  const auto d = parse_crystal(text);
  EXPECT_EQ(d.crystal.cell.sites.size(), 8u);
  EXPECT_NEAR(d.crystal.cell.atoms_per_cell(), 2.06, 1e-12);
  EXPECT_NE(run({"--list"}).find("wurtzite"), std::string::npos);
  // Conversion: POSCAR -> CIF -> summary.
  const auto dir = std::filesystem::temp_directory_path() / "amcset_test_cell_command";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  run({"--structure", "rocksalt", "--sublattice", "Mg", "--sublattice", "O", "--a", "4.212 A", "-o",
       (dir / "MgO.vasp").string()});
  text = run({(dir / "MgO.vasp").string(), "-o", (dir / "MgO.cif").string()});
  EXPECT_NE(text.find("8 atoms per cell"), std::string::npos) << text;
  EXPECT_NE(text.find("6 at 2.1060 A"), std::string::npos) << text; // Mg-O = a/2
  text = run({(dir / "MgO.cif").string(), "--from-format", "cif"});
  EXPECT_NE(text.find("Mg 4 (50.000 at.%)"), std::string::npos) << text;
  std::filesystem::remove_all(dir);
  // Errors.
  const auto fails = [](std::vector<std::string> args) {
    args.insert(args.begin(), "cell");
    std::ostringstream out;
    cli::print_cell(cli::parse_arguments(args), out);
  };
  EXPECT_THROW(fails({}), cli::UsageError);
  EXPECT_THROW(fails({"x.cif", "--structure", "bcc"}), cli::UsageError);
  EXPECT_THROW(fails({"--structure", "bcc", "--sublattice", "Fe", "--interstitial", "octahedral"}), cli::UsageError);
  EXPECT_THROW(fails({"--structure", "bbc", "--sublattice", "Fe"}), std::invalid_argument);
  EXPECT_THROW(fails({"--structure", "rocksalt", "--sublattice", "Mg", "--sublattice", "O"}), std::invalid_argument); // no size
  EXPECT_THROW(fails({"--structure", "bcc", "--sublattice", "Fe", "--beam", "0 0 -1"}), std::invalid_argument);
  EXPECT_THROW(fails({"--structure", "bcc", "--sublattice", "Fe", "--surface", "1 1"}), std::invalid_argument);
  EXPECT_THROW(fails({"missing_file.cif"}), std::runtime_error);
}

TEST(RunConfig, CrystalTransportKeys) {
  // The transport, its settings and the beam spot.
  const auto c = parse_config(R"({
    "ion": {"element": "B", "energy": "30 keV", "angle": "7 deg",
            "spot": {"shape": "circle", "diameter": "50 nm"}},
    "layers": [{"thickness": "1 um", "crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"]}}],
    "physics": {"transport": "crystal",
                "crystal": {"max_impact_parameter": "2.5 A", "simultaneous_distance": "0.3 A",
                            "local_stopping_fraction": 0.4, "local_stopping_decay": 0.35,
                            "thermal_vibrations": false}}
  })");
  const auto &t = c.trim;
  EXPECT_EQ(t.transport, trim::TrimConfig::Transport::crystal);
  EXPECT_EQ(t.ion.spot.shape, trim::BeamSpot::Shape::circle);
  EXPECT_NEAR(in_units(t.ion.spot.width, angstrom), 500.0, 1e-9);
  EXPECT_NEAR(in_units(t.ion.spot.height, angstrom), 500.0, 1e-9);
  EXPECT_NEAR(in_units(t.crystal.max_impact_parameter, angstrom), 2.5, 1e-12);
  EXPECT_NEAR(in_units(t.crystal.simultaneous_distance, angstrom), 0.3, 1e-12);
  EXPECT_DOUBLE_EQ(t.crystal.local_stopping_fraction, 0.4);
  EXPECT_DOUBLE_EQ(t.crystal.local_stopping_decay, 0.35);
  EXPECT_FALSE(t.crystal.thermal_vibrations);
  EXPECT_NO_THROW(t.validate());
  const std::string text = to_json(c);
  const auto b = parse_config(text);
  EXPECT_EQ(to_json(b), text);
  EXPECT_TRUE(b.trim.crystal == t.crystal);
  EXPECT_TRUE(b.trim.ion.spot == t.ion.spot);
  // "auto" is the automatic maximum impact parameter; width and height separately.
  auto o = c;
  apply_override(o, "physics.crystal.max_impact_parameter=auto");
  apply_override(o, "ion.spot.shape=rectangle");
  apply_override(o, "ion.spot.height=20 nm");
  EXPECT_EQ(in_units(o.trim.crystal.max_impact_parameter, angstrom), 0.0);
  EXPECT_EQ(o.trim.ion.spot.shape, trim::BeamSpot::Shape::rectangle);
  EXPECT_NEAR(in_units(o.trim.ion.spot.height, angstrom), 200.0, 1e-9);
  EXPECT_NEAR(in_units(o.trim.ion.spot.width, angstrom), 500.0, 1e-9);
  // Lists can be set where a list belongs: the surface plane, a new reference direction.
  apply_override(o, "layers.0.crystal.orientation.surface=[1, 1, 0]");
  apply_override(o, "layers.0.crystal.orientation.reference=[0,0,1]");
  EXPECT_EQ(o.trim.layers[0].crystal->orientation.surface, (std::vector<double>{1, 1, 0}));
  EXPECT_EQ(o.trim.layers[0].crystal->orientation.reference, (std::vector<double>{0, 0, 1}));
  // Scenarios that do not use them do not mention them (their files are unchanged).
  const std::string plain = to_json(example_config(Mode::srim));
  EXPECT_EQ(plain.find("\"spot\""), std::string::npos);
  EXPECT_EQ(plain.find("\"crystal\""), std::string::npos);
  EXPECT_NE(plain.find("\"transport\": \"trim2013\""), std::string::npos);
  auto refined = example_config(Mode::srim);
  apply_override(refined, "physics.transport=refined");
  EXPECT_EQ(refined.trim.transport, trim::TrimConfig::Transport::refined);
  // Errors name the key.
  const auto message = [](const std::string &json) -> std::string {
    try {
      parse_config(json).trim.validate();
    } catch (const std::invalid_argument &e) {
      return e.what();
    }
    return "";
  };
  const std::string head = R"({"ion": {"element": "B", "energy": "30 keV")";
  const std::string layers = R"(, "layers": [{"thickness": "1 um", "elements": [{"element": "Si"}]}])";
  EXPECT_NE(message(head + "}" + layers + R"(, "physics": {"transport": "marlowe"}})").find("physics.transport"),
            std::string::npos);
  EXPECT_NE(message(head + "}" + layers + R"(, "physics": {"crystal": {"local_fraction": 1}}})")
                .find("physics.crystal.local_fraction"),
            std::string::npos);
  EXPECT_NE(message(head + "}" + layers + R"(, "physics": {"crystal": {"local_stopping_fraction": 2}}})")
                .find("physics.crystal.local_stopping_fraction"),
            std::string::npos);
  EXPECT_NE(message(head + R"(, "spot": {"shape": "ring"}})" + layers + "}").find("ion.spot.shape"), std::string::npos);
  EXPECT_NE(message(head + R"(, "spot": {"width": "0 A"}})" + layers + "}").find("ion.spot.width"), std::string::npos);
  EXPECT_NE(message(head + R"(, "spot": {"radius": "1 um"}})" + layers + "}").find("ion.spot.radius"), std::string::npos);
}

TEST(RunConfig, CrystalRunsRecordEveryTrimBugAsFixed) {
  // The scenario keeps its bug switches (they apply to the other
  // transports); what a crystal run records in results.json is what it
  // does: no TRIM bug, exact sums in the TRIM-format tables.
  auto c = parse_config(one_layer(R"("crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"]})"));
  c.trim.ion_count = 4;
  ASSERT_TRUE(c.trim.bugs.all_reproduced());
  const auto recorded = [](const RunConfig &config) {
    auto engine = make_engine(config);
    return parse_config(dynamic_cast<const trim::TrimSimulation &>(*engine).config().configuration_json);
  };
  const auto as_trim = recorded(c);
  EXPECT_TRUE(as_trim.trim.bugs.all_reproduced());
  EXPECT_TRUE(as_trim.trim.single_precision_tables);
  c.trim.transport = trim::TrimConfig::Transport::crystal;
  EXPECT_NE(to_json(c).find("\"trim_bugs\": \"trim\""), std::string::npos); // the scenario is unchanged
  const auto crystal = recorded(c);
  EXPECT_EQ(crystal.trim.transport, trim::TrimConfig::Transport::crystal);
  for (const auto &m : trim::TrimBugs::members()) EXPECT_FALSE(crystal.trim.bugs.*(m.second)) << m.first;
  EXPECT_FALSE(crystal.trim.single_precision_tables);
}

TEST(RunConfig, PolycrystalKeys) {
  const auto c = parse_config(one_layer(R"("crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
      "polycrystal": {"grain_size": "20 nm", "texture": "fiber", "spread": "2 deg", "columnar": true, "seed": 7,
                      "period": ["100 nm", "80 nm"], "boundary_width": "5 A"}})"));
  ASSERT_TRUE(c.trim.layers[0].crystal.has_value());
  ASSERT_TRUE(c.trim.layers[0].crystal->polycrystal.has_value());
  const auto &p = *c.trim.layers[0].crystal->polycrystal;
  EXPECT_EQ(p.source, crystal::Polycrystal::Source::voronoi);
  EXPECT_DOUBLE_EQ(in_units(p.grain_size, angstrom), 200.0);
  EXPECT_EQ(p.texture, crystal::Polycrystal::Texture::fiber);
  EXPECT_DOUBLE_EQ(in_units(p.spread, amcset::common::degree), 2.0);
  EXPECT_TRUE(p.columnar);
  EXPECT_EQ(p.seed, 7u);
  EXPECT_DOUBLE_EQ(in_units(p.period_y, angstrom), 1000.0);
  EXPECT_DOUBLE_EQ(in_units(p.period_z, angstrom), 800.0);
  EXPECT_DOUBLE_EQ(in_units(p.boundary_width, angstrom), 5.0);
  EXPECT_NO_THROW(c.trim.validate());
  // The scenario written out reads back to itself.
  const std::string text = to_json(c);
  EXPECT_NE(text.find("\"polycrystal\""), std::string::npos);
  const auto b = parse_config(text);
  EXPECT_EQ(to_json(b), text);
  EXPECT_TRUE(*b.trim.layers[0].crystal == *c.trim.layers[0].crystal);

  // Grains from a file: the source follows from its name, and it is read.
  const auto file = std::filesystem::temp_directory_path() / "amcset_polycrystal_keys_test.txt";
  {
    std::ofstream out(file);
    out << "unit nm\nperiod 40 20\ngrain 5 10 10 0 0 0\ngrain 5 30 10 30 0 0\n";
  }
  const auto f = parse_config(one_layer(R"("crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
      "polycrystal": {"file": ")" + file.generic_string() + R"("}})"));
  const auto &fp = *f.trim.layers[0].crystal->polycrystal;
  EXPECT_EQ(fp.source, crystal::Polycrystal::Source::grain_file);
  ASSERT_TRUE(fp.map != nullptr);
  EXPECT_EQ(fp.map->grain_count(), 2u);
  EXPECT_EQ(parse_config(to_json(f)).trim.layers[0].crystal->polycrystal->map->grain_count(), 2u);
  std::filesystem::remove(file);

  // Mistakes name the key.
  const auto message = [](const std::string &polycrystal) {
    try {
      parse_config(one_layer(R"("crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
          "polycrystal": )" + polycrystal + "}"));
    } catch (const std::exception &e) {
      return std::string(e.what());
    }
    return std::string();
  };
  EXPECT_NE(message(R"({"texture": "cube"})").find("polycrystal.texture"), std::string::npos);
  EXPECT_NE(message(R"({"grain_size": "1 A"})").find("polycrystal.grain_size"), std::string::npos);
  EXPECT_NE(message(R"({"source": "ebsd"})").find("polycrystal.file"), std::string::npos);
  EXPECT_NE(message(R"({"period": ["1 um"]})").find("polycrystal.period"), std::string::npos);
  EXPECT_NE(message(R"({"grains": 5})").find("grains"), std::string::npos);
  EXPECT_NE(message(R"({"file": "no_such_grain_file.txt"})").find("no_such_grain_file.txt"), std::string::npos);
}
