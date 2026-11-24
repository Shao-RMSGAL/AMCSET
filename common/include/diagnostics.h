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

/*! \file diagnostics.h
 *  \brief Diagnostic log, crash reports and per-thread breadcrumbs.
 *
 * - A diagnostic log (off unless a front end enables it) with levels debug,
 *   info, warning and error, written to a file and kept in memory (the last
 *   lines go into crash reports). Use AMCSET_LOG(info) << "..."; a disabled
 *   level costs one relaxed atomic load. The log is for run-level events
 *   (configuration, start, threads, files, errors), never for the transport
 *   loops.
 * - Crash handlers (installed by install_crash_handlers()): on a fatal
 *   signal (SIGSEGV, SIGFPE, SIGILL, SIGBUS, SIGABRT), an unhandled Windows
 *   exception, std::terminate (e.g. an exception escaping a thread) or a
 *   pure-virtual call, a crash report is written to the crash directory and
 *   to stderr: what happened, the stack trace (symbolised from the PDB on
 *   Windows), the run context, what every worker thread was doing (which
 *   ion of which seed, and how to rerun exactly that ion) and the last log
 *   lines. On Windows a minidump (.dmp, for Visual Studio or WinDbg) is
 *   written next to it.
 * - Breadcrumbs: a worker thread registers an Activity and updates the
 *   particle it is transporting (one relaxed store per particle).
 *
 * Library code only records; front ends (CLI, GUI) call initialize().
 */

#pragma once

#include <cstdint>
#include <filesystem>
#include <sstream>
#include <string>
#include <string_view>

namespace amcset {
namespace diag {

enum class Level { debug = 0, info = 1, warning = 2, error = 3, off = 4 };

//! Name of a level ("debug", ...).
const char *to_string(Level level) noexcept;
//! Parse "debug", "info", "warning", "error", "off". \throws std::invalid_argument
Level parse_level(std::string_view text);

//! Front-end settings.
struct Options {
  std::string program = "amcset";        //!< Name used in reports and file names
  std::filesystem::path log_file;        //!< Log file (empty: none; memory only)
  Level level = Level::warning;          //!< Lowest level recorded
  bool echo_to_stderr = false;           //!< Also print log lines to stderr
  //! Where crash reports go (empty: the log file's directory, else the
  //! system temporary directory / amcset-crashes).
  std::filesystem::path crash_directory;
};

/*!
 * \brief Set up logging (and the abseil log sink that forwards abseil
 * messages into the log) and install the crash handlers. Call once from
 * main() before starting threads; calling it again replaces the options.
 * \p argv0 helps symbolisation. Failures to open the log file are reported
 * on stderr and leave logging to memory.
 */
void initialize(const char *argv0, const Options &options);

//! Install only the crash handlers (initialize() does this too).
void install_crash_handlers(const char *argv0);

//! Change the level at run time (e.g. the GUI's "detailed log" switch).
void set_level(Level level) noexcept;
Level level() noexcept;
//! Whether messages of \p level are recorded.
bool enabled(Level level) noexcept;

//! Record one message (thread-safe; a trailing newline is added).
void write(Level level, std::string_view message);

//! The log file in use (empty if none).
std::filesystem::path log_file();
//! Directory of the crash reports.
std::filesystem::path crash_directory();

/*!
 * \brief A short description of the running calculation, printed in crash
 * reports (e.g. engine, seed, ions, first ion, title). Truncated to 2 kB.
 */
void set_context(std::string_view text);

/*!
 * \brief Breadcrumb of the calling thread: registers it under \p what for
 * crash reports until destruction. Not copyable; at most 256 threads are
 * tracked at a time (further ones are not listed).
 */
class Activity {
public:
  explicit Activity(const char *what) noexcept;
  ~Activity();
  Activity(const Activity &) = delete;
  Activity &operator=(const Activity &) = delete;

  /*!
   * \brief The particle this thread works on: \p index (0-based) of a run
   * with \p seed; \p block_first is the first particle of the group the
   * engine processes together (SRIM mode: blocks of 16), from which the
   * particle can be rerun exactly.
   */
  void set_particle(std::uint64_t index, std::uint64_t block_first, std::uint32_t seed) noexcept;
  //! Done with particles (idle, merging, writing files).
  void clear_particle() noexcept;

private:
  int slot_ = -1;
};

//! Builds one log line; written when it goes out of scope.
class Line {
public:
  explicit Line(Level level) : level_(level) {}
  ~Line() { write(level_, stream_.str()); }
  Line(const Line &) = delete;
  Line &operator=(const Line &) = delete;
  template <class T> Line &operator<<(const T &v) {
    stream_ << v;
    return *this;
  }

private:
  Level level_;
  std::ostringstream stream_;
};

} // namespace diag
} // namespace amcset

//! AMCSET_LOG(info) << "text " << value; (the stream is not evaluated when
//! the level is disabled).
#define AMCSET_LOG(level)                                                    \
  if (!::amcset::diag::enabled(::amcset::diag::Level::level)) {             \
  } else                                                                     \
    ::amcset::diag::Line(::amcset::diag::Level::level)
