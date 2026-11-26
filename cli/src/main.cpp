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

/*! \file main.cpp
 *  \brief The AMCSET command-line program (no Qt required).
 *
 * Subcommands: run, trim, convert, stopping, elements, template, validate,
 * version, help. See docs/CLI.md or `amcset help`.
 */

#include "amcset_common.h"
#include "cli.h"
#include "config.h"
#include "diagnostics.h"
#include "element_data.h"
#include "interrupts.h"
#include "physics/scattering.h"
#include "physics/stopping.h"
#include "run_config.h"
#include "trim/trim.h"
#include "trim/trim_in.h"

#include <absl/log/globals.h>
#include <absl/log/initialize.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <csignal>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>

using namespace amcset;
using namespace amcset::common;

namespace {

// Set by a signal or console event during a run (interrupts.h); the progress
// loop forwards it to the engine (request_cancel() locks, a handler must not).
std::atomic<bool> g_interrupted{false};

//! The exit status of an interrupted run: 128 + the signal (130 for Ctrl+C).
int interrupted_status() {
  const int sig = cli::interrupt_signal();
  return cli::interrupted_exit_status(sig != 0 ? sig : SIGINT);
}

const char *usage = R"(AMCSET - Monte Carlo transport of ions and electrons in matter

Usage:
  amcset run <scenario.json> [options]      Run a scenario (SRIM mode or native)
  amcset trim <TRIM.IN> [options]           Run a TRIM batch input in SRIM mode
  amcset convert <in> <out>                 Convert TRIM.IN <-> scenario JSON
  amcset stopping [options]                 Stopping powers and ranges (like SR Module)
  amcset elements [symbol|Z ...]            Element data (masses, densities, Es, I)
  amcset cell <file> | --structure NAME     A crystal's unit cell: summary, conversion,
                                            scenario "crystal" object, beam directions
  amcset template [srim|native] [-o FILE]   Print (or save) an example scenario
  amcset validate <scenario.json|TRIM.IN>   Check a scenario without running it
  amcset version | help

Options for run/trim/validate:
  -o, --output DIR          Output directory (TRIM-format files, summary.json)
  -n, --ions N              Number of incident ions
  -j, --threads N           Worker threads (default: all cores)
  -s, --seed S              Random number seed
      --mode srim|native    Engine (scenario files only)
      --damage kp|full|monolayer
      --stopping MODEL      srim | builtin | bethe-bloch | lss
      --stopping-table PATH SRIM stopping tables (file or SR Module output dir;
                            default: the SRIM-2013 tables in data/srim2013)
      --srim-data DIR       SRIM Data dir: SNUC03 nuclear stopping, and pair
                            potentials if physics.pair_potentials asks for them
                            (TRIM.exe itself uses the universal potential)
      --fix-trim-bugs       Correct every bug of TRIM.exe that SRIM mode reproduces
      --fix-trim-bug NAME   Correct one of them (repeatable); see docs/SRIM_MODE.md
                            "Bugs of TRIM.exe". Results then depart from TRIM's.
      --fast                Fast approximations (physics.fast_math): tabulated
                            scattering angle, about 2x faster; statistically the
                            same results, but not number for number (not TRIM
                            step by step). Leave off for TRIM comparisons.
      --auto-depth          Size the last layer automatically: quick test runs until
                            the deepest of auto_depth.test_ions ions (1000) stops in
                            the last auto_depth.margin (10 %) of the target
      --log FILE            Write a diagnostic log (run steps, warnings, errors;
                            crash reports go next to it)
      --set KEY=VALUE       Override any scenario key, e.g. --set "ion.energy=2 MeV"
                            or --set "layers.0.elements.0.Ed=25 eV" (repeatable)
      --save-config FILE    Write the final scenario (after overrides) as JSON
      --json                Print the summary as JSON
      --particles FILE      Particle source: a TRIM.DAT (ions, or recoils with
                            --source recoils), or BACKSCAT/TRANSMIT/SPUTTER.txt
      --source beam|ions|recoils
      --exyz ENERGY         Write EXYZ.txt (ion position every ENERGY of loss)

Trajectory and MD event files (SRIM mode, every transport; not TRIM outputs;
docs/OUTPUT.md):
      --trajectory-file F   Every recorded point (ion, track, parent, kind, position,
                            direction, energies, time) to F (.csv or .parquet; a
                            relative name goes into the output directory). Records
                            all ions unless output.trajectories gives a number.
      --trajectory-points P all (default) or births (first point of every track)
      --md-events F         Hybrid BCA-MD recoil events to F (.csv or .parquet);
                            needs --damage full (or monolayer) and --md-ec
      --md-ec LIST          Transition energies E_c, e.g. "1 keV" or "0.5,1,2,5 keV"
      --md-emin E           Recoils below E are not events (default "Ed": the struck
                            element's displacement energy)
      --md-handoff          Stop following events (MD takes them; one E_c; not
                            TRIM; tallies then cover the BCA part only)

Output files (SRIM mode):
      --formats LIST        Any of trim,csv,json,parquet (default trim,csv,json)
      --output-units LIST   Units per dimension, e.g. length=nm,energy=keV,
                            energy_per_length=keV/um (or just: nm,keV)
      --number-format F     shortest | general | fixed | scientific | engineering,
                            optionally with digits, e.g. scientific:4
      --output-spec FILE    JSON file choosing tables, columns, rows, units and
                            formats (see docs/OUTPUT.md)
  -q, --quiet               No progress bar
  -v, --verbose             Log engine messages

Options for stopping:
      --ion X               Projectile (symbol or Z), e.g. --ion He
      --ion-mass M          Projectile mass (default: most abundant isotope)
      --target SPEC         Element or composition, e.g. Fe or "Fe:0.91,Cr:0.09"
      --density D           Target density (default: element density)
      --gas                 Gas-phase target
      --from E --to E       Energy range (default 10 keV .. 10 MeV)
      --per-decade N        Points per decade (default 10)
      --correct-nickel      Ni projectiles: Co/Cu-interpolated stopping (also for run/trim;
                            physics.stopping_corrections.nickel_ions, docs/VALIDATION.md)
      --model MODEL         srim | builtin | bethe-bloch | lss (default srim with
                            the bundled SRIM-2013 tables; builtin without them)
      --stopping-table PATH, --srim-data DIR
      --units U             eV/A (default), keV/um, MeV/(mg/cm2), eV/(1e15 atoms/cm2)
      --csv                 Comma-separated output
      --lammps              Electronic stopping for LAMMPS fix electron/stopping
                            (metal units: E in eV, N*S_e in eV/A, one column per
                            atom type: --ion Fe,Cr lists the types in order)
      --ntypes N            With --lammps and one --ion: repeat its column N times

Options for cell (unit cells of crystalline layers; docs/CLI.md "Crystals"):
      <file>                A cell file: CIF, VASP POSCAR/CONTCAR, LAMMPS data,
                            extended XYZ (format from the name or --from-format)
      --list                The built-in structures and their sublattices
      --structure NAME      A built-in structure, e.g. bcc, fcc, diamond, rocksalt
      --sublattice SPEC     The elements of one sublattice, e.g. Fe or
                            "Fe:0.9,Cr:0.1" (repeat in the structure's order)
      --interstitial S=SPEC Atoms on built-in interstitial sites, e.g.
                            octahedral=C:0.01 (repeatable)
      --a LENGTH            Lattice constant; without it --density (or the
                            element's density) sizes the cell
      --c-over-a X          c/a of hexagonal and tetragonal structures
      --elements LIST       Elements of a POSCAR without a species line, e.g. Ga,As
      --surface "h k l"     Surface plane (default "0 0 1"); --reference "u v w":
                            the direction in the surface along the lab y axis;
                            --twist ANGLE: rotation about the surface normal
      --beam "u v w"        Print the ion.angle (tilt) and orientation.twist that
                            send the beam along this crystal direction
      --tilt ANGLE          Print the crystal direction of a beam with this tilt
  -o, --output FILE         Write the cell (format from the name or --to-format)
      --json                Print the "crystal" object for a scenario layer

Physical values take units: "1 MeV", "1.2 um", "7.87 g/cm3", "40 eV", "7 deg".
)";

void print_summary(const config::RunConfig &cfg, const SimulationEngine &engine, bool as_json) {
  if (const auto *t = dynamic_cast<const trim::TrimSimulation *>(&engine)) {
    if (as_json) {
      std::cout << t->summary_json() << "\n";
      return;
    }
    const auto &r = t->tallies();
    const double n = std::max<double>(1.0, static_cast<double>(r.ions));
    std::cout << std::fixed << std::setprecision(1);
    std::cout << "\nSRIM mode results (" << r.ions << " ions, " << std::setprecision(1)
              << t->elapsed_seconds() << " s)\n";
    std::cout << "  Projected range      " << r.range.mean() << " A   straggling "
              << r.range.stddev() << " A\n";
    std::cout << "  Lateral range        " << r.lateral_abs_y.mean() << " A   straggling "
              << r.lateral_y.rms() << " A\n";
    std::cout << "  Radial range         " << r.radial.mean() << " A   straggling "
              << r.radial.stddev() << " A\n";
    std::cout << std::setprecision(3) << "  Skewness / kurtosis  " << r.range.skewness() << " / "
              << r.range.kurtosis() << "\n";
    std::cout << std::setprecision(2) << "  Vacancies per ion    " << r.total_vacancies / n << "\n";
    if (cfg.trim.damage != trim::DamageCalculation::kinchin_pease) {
      std::cout << "  Displacements / ion  " << r.total_displacements / n << "\n";
      std::cout << "  Replacements / ion   " << r.total_replacements / n << "\n";
    }
    std::cout << std::setprecision(4) << "  Backscattered        " << r.backscattered / n * 100.0
              << " %   transmitted " << r.transmitted / n * 100.0 << " %\n";
    for (std::size_t i = 0; i < t->element_names().size(); ++i) {
      if (r.sputtered[i] > 0) {
        std::cout << "  Sputtering yield " << t->element_names()[i] << "  " << r.sputtered[i] / n
                  << " atoms/ion\n";
      }
    }
    for (const auto &m : r.md) {
      std::cout << std::setprecision(4) << "  MD events, E_c " << m.transition_energy_eV << " eV: "
                << m.events / n << " per ion (" << m.birth_events / n << " born, " << m.crossing_events / n
                << " crossing), " << m.energy_eV / n << " eV per ion to MD"
                << (cfg.trim.md_events.handoff ? " (handed off)" : " (recorded)") << "\n";
    }
    if (!cfg.trim.output_directory.empty()) {
      std::cout << "  Output written to " << cfg.trim.output_directory << "\n";
    }
    return;
  }
  const auto tallies = engine.snapshot_tallies();
  if (as_json) {
    std::cout << tallies.to_json() << "\n";
    return;
  }
  std::cout << std::fixed << std::setprecision(1) << "\nNative AMCSET results ("
            << tallies.incidents << " histories)\n";
  std::cout << "  Mean range           " << in_units(static_cast<double>(tallies.range_stats.mean) * meter, angstrom)
            << " A   straggling "
            << in_units(std::sqrt(static_cast<double>(tallies.range_stats.variance())) * meter, angstrom)
            << " A\n";
  std::cout << "  Displacements        " << tallies.total_displacements << "\n";
  std::cout << "  Left through surface " << tallies.sputtered << "\n";
}

int run_engine(config::RunConfig &cfg, const cli::Options &opt) {
  // Ctrl+C, SIGTERM, SIGHUP (Windows: the console events) set g_interrupted
  // from here on (the sizing test runs and the run stop on it); the default
  // handlers are back on every way out.
  const cli::InterruptGuard interrupt_guard(g_interrupted);
  if (cfg.mode == config::Mode::srim && cfg.trim.auto_depth.enabled) {
    // Automatic target depth: quick test runs size the last layer first.
    cfg.trim.validate();
    if (!opt.quiet) std::cerr << "Sizing the target (automatic depth)...\n";
    const auto r = trim::size_target(cfg.trim, &g_interrupted);
    for (const auto &n : r.notes) {
      std::cerr << "Note: " << n << "\n";
      AMCSET_LOG(info) << n;
    }
    if (g_interrupted.load()) {
      std::cerr << "Interrupted while sizing the target; nothing was run.\n";
      return interrupted_status();
    }
  }
  // The scenario as it runs (with the sized target).
  if (!opt.save_config.empty()) config::save_config(cfg, opt.save_config);
  auto engine = config::make_engine(cfg);
  if (opt.verbose) std::cerr << engine->describe() << "\n";
  if (auto *t = dynamic_cast<trim::TrimSimulation *>(engine.get())) {
    for (const auto &w : t->warnings()) std::cerr << "Warning: " << w << "\n";
    // Nothing polls trajectories here: recorded ones only go to the
    // trajectory file (output.trajectory_file), not into memory as well.
    t->set_trajectory_polling(false);
  }
  AMCSET_LOG(info) << "engine ready: " << engine->total_count() << " particles";
  std::exception_ptr error;
  std::mutex done_mutex;
  std::condition_variable done_cv;
  bool done = false; // guarded by done_mutex
  std::thread worker([&] {
    try {
      engine->run();
    } catch (...) {
      error = std::current_exception();
    }
    {
      const std::lock_guard<std::mutex> lock(done_mutex);
      done = true;
    }
    done_cv.notify_all();
  });
  // The worker is joined on every way out of this function (destroying a
  // joinable std::thread would terminate the program).
  struct Joiner {
    std::thread &thread;
    common::SimulationEngine &engine;
    ~Joiner() {
      if (thread.joinable()) {
        engine.request_cancel();
        thread.join();
      }
    }
  } joiner{worker, *engine};
  // Wakes as soon as the run ends (not up to a poll interval later), and
  // every 200 ms for the progress bar.
  const auto wait_done = [&] {
    std::unique_lock<std::mutex> lock(done_mutex);
    return done_cv.wait_for(lock, std::chrono::milliseconds(200), [&] { return done; });
  };
  cli::ProgressBar bar(engine->total_count(), !opt.quiet);
  bool cancel_sent = false;
  const auto *trim_engine = dynamic_cast<const trim::TrimSimulation *>(engine.get());
  const bool warn_escape = trim_engine != nullptr && cfg.trim.auto_depth.warn_transmitted &&
                           cfg.trim.source == trim::TrimConfig::Source::beam;
  bool escape_warned = false;
  while (true) {
    if (g_interrupted.load() && !cancel_sent) {
      AMCSET_LOG(info) << "interrupted (signal " << cli::interrupt_signal()
                       << "): stopping after the running particles";
      engine->request_cancel();
      cancel_sent = true;
      bar.clear_line();
      std::cerr << "Stopping after the ions in progress; the results so far will be written "
                   "(interrupt again to quit at once).\n";
    }
    if (warn_escape && !escape_warned && trim_engine->transmitted_so_far() > 0) {
      // The run goes on; the user decides afterwards.
      escape_warned = true;
      bar.clear_line();
      std::cerr << "Warning: ions are leaving through the back of the target (transmitted); the run "
                   "continues. For a target deep enough to stop every ion, run again with --auto-depth "
                   "(or ignore this for a transmission study; --set auto_depth.warn_transmitted=false "
                   "turns it off).\n";
      AMCSET_LOG(warning) << "ions are transmitted through the target";
    }
    bar.update(engine->completed_count());
    if (wait_done()) break;
  }
  worker.join();
  bar.finish(engine->completed_count());
  if (error) std::rethrow_exception(error);
  if (engine->is_cancelled()) std::cerr << "Interrupted: results cover the finished ions only.\n";
  if (warn_escape && trim_engine->transmitted_so_far() > 0) {
    std::cerr << "Warning: " << trim_engine->transmitted_so_far() << " of " << engine->completed_count()
              << " ions escaped through the back of the target.\n";
  }
  print_summary(cfg, *engine, opt.json);
  return engine->is_cancelled() ? interrupted_status() : 0;
}

config::RunConfig load_any(const std::string &path) {
  if (cli::looks_like_trim_in(path)) {
    config::RunConfig c;
    c.mode = config::Mode::srim;
    c.trim = trim::read_trim_in(path);
    return c;
  }
  return config::load_config(path);
}

int cmd_run(const cli::Options &opt, bool trim_in) {
  if (opt.positional.size() != 1) throw cli::UsageError("expected one input file");
  config::RunConfig cfg = trim_in ? config::RunConfig{} : load_any(opt.positional[0]);
  if (trim_in) {
    cfg.mode = config::Mode::srim;
    cfg.trim = trim::read_trim_in(opt.positional[0]);
    cfg.trim.output_directory = "output";
  }
  cli::apply_options(cfg, opt);
  return run_engine(cfg, opt);
}

int cmd_validate(const cli::Options &opt) {
  if (opt.positional.size() != 1) throw cli::UsageError("expected one input file");
  auto cfg = load_any(opt.positional[0]);
  cli::apply_options(cfg, opt);
  if (cfg.mode == config::Mode::srim) {
    cfg.trim.validate();
    std::cout << cfg.trim.describe();
    for (const auto &w : cfg.trim.plausibility_warnings()) std::cout << "Warning: " << w << "\n";
  } else {
    auto engine = config::make_engine(cfg);
    std::cout << engine->describe();
  }
  std::cout << "OK\n";
  return 0;
}

int cmd_convert(const cli::Options &opt) {
  if (opt.positional.size() != 2) throw cli::UsageError("expected <input> <output>");
  const std::string &in = opt.positional[0];
  const std::string &out = opt.positional[1];
  auto cfg = load_any(in);
  cli::apply_options(cfg, opt);
  if (cli::looks_like_trim_in(out) || out.ends_with(".IN") || out.ends_with(".in")) {
    std::ofstream o(out, std::ios::binary);
    if (!o) throw std::runtime_error("Cannot write " + out);
    trim::write_trim_in(o, cfg.trim, trim::extras_for(cfg.trim));
  } else {
    config::save_config(cfg, out);
  }
  std::cout << "Wrote " << out << "\n";
  return 0;
}

int cmd_template(const cli::Options &opt) {
  const auto mode = opt.positional.empty() ? config::Mode::srim
                                           : config::parse_mode(opt.positional[0]);
  const auto example = config::example_config(mode);
  if (opt.output) {
    const std::filesystem::path file(*opt.output);
    if (file.has_parent_path()) std::filesystem::create_directories(file.parent_path());
    config::save_config(example, file);
    if (!opt.quiet) std::cout << "Wrote " << file.string() << "\n";
  } else {
    std::cout << config::to_json(example);
  }
  return 0;
}

int cmd_elements(const cli::Options &opt) {
  std::vector<unsigned> zs;
  for (const auto &p : opt.positional) {
    try {
      zs.push_back(static_cast<unsigned>(std::stoul(p)));
    } catch (const std::exception &) {
      zs.push_back(Element::by_symbol(p).z());
    }
  }
  if (zs.empty()) {
    for (unsigned z = 1; z <= Element::max_z; ++z) zs.push_back(z);
  }
  std::cout << "  Z  Sym  Name            A_MAI   M_MAI(amu) M_nat(amu)  rho(g/cm3)  N(1/cm3)   Es(eV)  I(eV)\n";
  for (unsigned z : zs) {
    const auto &e = Element::by_z(z);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%3u  %-3s  %-14s %5u  %10.3f %10.3f  %10.4f  %10.3E  %6.2f  %6.1f\n",
                  z, std::string(e.symbol()).c_str(), std::string(e.name()).c_str(),
                  e.most_abundant_mass_number(), (e.most_abundant_isotope_mass() / atomic_mass_unit).value(),
                  (e.natural_mass() / atomic_mass_unit).value(),
                  in_units(e.solid_density(), gram_per_cubic_centimeter),
                  in_units(e.solid_atomic_density(), per_cubic_centimeter),
                  in_units(e.sublimation_energy(), electron_volt),
                  in_units(e.mean_excitation_energy(), electron_volt));
    std::cout << buf;
  }
  return 0;
}

int cmd_stopping(const cli::Options &opt) {
  return cli::print_stopping_table(opt, std::cout);
}

int cmd_cell(const cli::Options &opt) {
  return cli::print_cell(opt, std::cout);
}

// --- diag-crash: check that crash reports work (not in the help text) -------

volatile int g_depth_guard = 1;

int overflow_stack(int n) {
  volatile char frame[4096];
  frame[0] = static_cast<char>(n);
  if (g_depth_guard != 0) return overflow_stack(n + 1) + frame[0];
  return frame[0];
}

/*!
 * `amcset diag-crash KIND` crashes on purpose (KIND: segv, abort, terminate,
 * thread-exception, stack-overflow) so that the crash report (stderr and the
 * crash directory, see diagnostics.h) can be checked on a machine.
 */
int cmd_diag_crash(const cli::Options &opt) {
  const std::string kind = opt.positional.empty() ? "segv" : opt.positional[0];
  diag::set_context("diag-crash " + kind + " (a deliberate crash; no calculation)");
  diag::Activity activity("diag-crash");
  activity.set_particle(41, 32, 716381);
  AMCSET_LOG(info) << "diag-crash: crashing on purpose (" << kind << ")";
  if (kind == "segv") {
    volatile int *p = nullptr;
    *p = 1;
  } else if (kind == "abort") {
    std::abort();
  } else if (kind == "terminate") {
    std::terminate();
  } else if (kind == "thread-exception") {
    std::thread([] { throw std::runtime_error("test exception escaping a worker thread"); }).join();
  } else if (kind == "stack-overflow") {
    return overflow_stack(0);
  }
  throw cli::UsageError("diag-crash: expected segv, abort, terminate, thread-exception or stack-overflow");
}

} // namespace

int main(int argc, char **argv) {
  absl::InitializeLog();
  // Crash reports from the start (diagnostics.h); the log is set up below.
  diag::install_crash_handlers(argc > 0 ? argv[0] : "amcset");
  try {
    cli::Options opt = cli::parse_arguments(argc, argv);
    absl::SetMinLogLevel(opt.verbose ? absl::LogSeverityAtLeast::kInfo
                                     : absl::LogSeverityAtLeast::kError);
    absl::SetStderrThreshold(opt.verbose ? absl::LogSeverityAtLeast::kInfo
                                         : absl::LogSeverityAtLeast::kError);
    {
      // Run-level events are always kept in memory (for crash reports);
      // --log writes them to a file, -v also to stderr.
      diag::Options d;
      d.program = "amcset";
      d.log_file = opt.log_file;
      d.level = opt.verbose ? diag::Level::debug : diag::Level::info;
      d.echo_to_stderr = opt.verbose;
      diag::initialize(argc > 0 ? argv[0] : "amcset", d);
      std::string line = "amcset " AMCSET_VERSION;
      for (int i = 1; i < argc; ++i) line += std::string(" ") + argv[i];
      AMCSET_LOG(info) << line;
    }
    if (opt.command == "help" || opt.command.empty() || opt.help) {
      std::cout << usage;
      return 0;
    }
    if (opt.command == "version") {
      std::cout << "AMCSET " << AMCSET_VERSION << "\n";
      return 0;
    }
    // An interrupted run ends by its signal (POSIX), after its results are written.
    if (opt.command == "run") return cli::exit_after_interrupt(cmd_run(opt, false));
    if (opt.command == "trim") return cli::exit_after_interrupt(cmd_run(opt, true));
    if (opt.command == "validate") return cmd_validate(opt);
    if (opt.command == "convert") return cmd_convert(opt);
    if (opt.command == "template") return cmd_template(opt);
    if (opt.command == "elements") return cmd_elements(opt);
    if (opt.command == "stopping") return cmd_stopping(opt);
    if (opt.command == "cell") return cmd_cell(opt);
    if (opt.command == "diag-crash") return cmd_diag_crash(opt);
    throw cli::UsageError("unknown command '" + opt.command + "'");
  } catch (const cli::UsageError &e) {
    std::cerr << "amcset: " << e.what() << "\n(run 'amcset help' for usage)\n";
    return 2;
  } catch (const std::exception &e) {
    AMCSET_LOG(error) << "amcset failed: " << e.what();
    std::cerr << "amcset: error: ";
    print_exception(e);
    return 1;
  }
}
