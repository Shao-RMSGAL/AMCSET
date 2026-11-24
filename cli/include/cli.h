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

/*! \file cli.h
 *  \brief Argument parsing and helpers of the command-line program (kept in
 *  a library so that the test suite can exercise them).
 */

#pragma once

#include "run_config.h"

#include <cstddef>
#include <iosfwd>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace amcset {
namespace cli {

//! Bad command line (exit status 2).
class UsageError : public std::runtime_error {
public:
  using std::runtime_error::runtime_error;
};

//! Parsed command line.
struct Options {
  std::string command;
  std::vector<std::string> positional;
  bool help = false;
  bool quiet = false;
  bool verbose = false;
  bool json = false;
  bool csv = false;
  bool gas = false;
  std::optional<std::string> output;
  std::optional<std::size_t> ions;
  std::optional<unsigned> threads;
  std::optional<std::uint32_t> seed;
  std::optional<std::string> mode;
  std::optional<std::string> damage;
  std::optional<std::string> stopping;
  std::optional<std::string> stopping_table;
  std::optional<std::string> srim_data;
  std::vector<std::string> overrides;
  //! --fix-trim-bug NAME (repeatable; "all" = every bug of TRIM.exe).
  std::vector<std::string> fixed_bugs;
  bool fast_math = false;       //!< --fast (physics.fast_math)
  bool auto_depth = false;      //!< --auto-depth (auto_depth.enabled)
  std::string log_file;         //!< --log FILE (diagnostic log)
  std::string save_config;
  // output files
  std::optional<std::string> formats;       //!< --formats trim,csv,json,parquet
  std::optional<std::string> output_spec;   //!< --output-spec FILE (JSON)
  std::optional<std::string> output_units;  //!< --output-units length=nm,energy=keV
  std::optional<std::string> number_format; //!< --number-format scientific:4
  // particle source
  std::optional<std::string> particles;     //!< --particles TRIM.DAT
  std::optional<std::string> source;        //!< --source beam|ions|recoils
  std::optional<std::string> exyz;          //!< --exyz "10 keV"
  // trajectory and MD event files (output.trajectory_file, output.md_events)
  std::optional<std::string> trajectory_file;   //!< --trajectory-file FILE
  std::optional<std::string> trajectory_points; //!< --trajectory-points all|births
  std::optional<std::string> md_events;         //!< --md-events FILE
  std::vector<std::string> md_ec;               //!< --md-ec LIST (repeatable)
  std::optional<std::string> md_emin;           //!< --md-emin ENERGY|Ed
  bool md_handoff = false;                      //!< --md-handoff
  // stopping command
  std::string ion;
  std::string ion_mass;
  std::string target;
  std::string density;
  std::string from = "10 keV";
  std::string to = "10 MeV";
  unsigned per_decade = 10;
  std::string units = "eV/A";
  bool lammps = false;  //!< --lammps: a table for LAMMPS fix electron/stopping
  unsigned ntypes = 0;  //!< --ntypes N: repeat a single ion's column N times (LAMMPS)
  // cell command (unit cells: built-in structures and cell files)
  bool list = false;                      //!< --list: the built-in structures
  std::string structure;                  //!< --structure NAME
  std::vector<std::string> sublattices;   //!< --sublattice SPEC (repeatable), e.g. "Fe:0.9,Cr:0.1"
  std::vector<std::string> interstitials; //!< --interstitial SITES=SPEC (repeatable), e.g. octahedral=C:0.01
  std::string lattice_a;                  //!< --a LENGTH
  std::string c_over_a;                   //!< --c-over-a X
  std::string cell_from_format;           //!< --from-format F
  std::string cell_to_format;             //!< --to-format F
  std::string cell_elements;              //!< --elements LIST (POSCAR without a species line)
  std::string surface;                    //!< --surface "h k l"
  std::string reference;                  //!< --reference "u v w"
  std::string twist;                      //!< --twist ANGLE
  std::string beam;                       //!< --beam "u v w"
  std::string tilt;                       //!< --tilt ANGLE
};


//! Parse argv. \throws UsageError
Options parse_arguments(int argc, const char *const *argv);
Options parse_arguments(const std::vector<std::string> &args);

//! Apply the command-line options to a scenario (options win over the file).
void apply_options(config::RunConfig &config, const Options &options);

//! Whether \p path names a TRIM.IN-style file (by name or content).
bool looks_like_trim_in(const std::string &path);

/*!
 * \brief Parse a target composition "Fe:0.91,Cr:0.09" (or a single
 * element "Fe") into (Z, atom fraction) pairs normalized to 1.
 * \throws std::invalid_argument
 */
std::vector<std::pair<unsigned, double>> parse_composition(const std::string &spec);

//! The `stopping` command: SR-Module-like table to \p out.
int print_stopping_table(const Options &options, std::ostream &out);

/*!
 * \brief The `cell` command: a unit cell from a file or a built-in
 * structure; prints its summary (or with --json the scenario's "crystal"
 * object) to \p out, converts it with -o, and gives the tilt and twist of a
 * beam direction.
 * \throws UsageError, std::invalid_argument, std::runtime_error
 */
int print_cell(const Options &options, std::ostream &out);

//! A text progress bar on stderr (only when enabled and stderr is a terminal).
class ProgressBar {
public:
  ProgressBar(std::size_t total, bool enabled);
  void update(std::size_t done);
  void finish(std::size_t done);
  //! Erase the bar (before printing a message); the next update redraws it.
  void clear_line();

private:
  std::size_t total_;
  bool enabled_;
  double start_;
  int last_permille_ = -1;
};

} // namespace cli
} // namespace amcset
