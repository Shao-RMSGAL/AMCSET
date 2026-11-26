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

/*! \file interrupts.h
 *  \brief Stopping a run on signals and console events.
 *
 * POSIX: SIGINT (Ctrl+C), SIGTERM (kill, systemd, a batch scheduler's time
 * limit) and SIGHUP (the terminal went away). Windows: Ctrl+C, Ctrl+Break and
 * the console's close, log-off and shutdown events. The first one asks the
 * run to stop after the ions in progress, so that the partial results are
 * written; a second one ends the program at once. After a stopped run the
 * program exits as POSIX programs do: by the same signal with its default
 * action (the shell then reports 128 + the signal number, e.g. 130 for SIGINT
 * and 143 for SIGTERM), on Windows with exit code 130. Other signals keep
 * their default action (SIGPIPE ends the program quietly, SIGQUIT dumps core).
 */

#pragma once

#include <atomic>

namespace amcset {
namespace cli {

/*!
 * \brief While it exists, interrupts set \p flag instead of ending the
 * program (one guard at a time). Only flags are set in the handler (a
 * handler must not lock mutexes); the caller polls \p flag.
 */
class InterruptGuard {
public:
  explicit InterruptGuard(std::atomic<bool> &flag);
  ~InterruptGuard();
  InterruptGuard(const InterruptGuard &) = delete;
  InterruptGuard &operator=(const InterruptGuard &) = delete;
};

//! The signal that interrupted the run (SIGINT for every Windows console event), or 0.
int interrupt_signal() noexcept;

//! The conventional exit status after signal \p signal_number: 128 + the number.
constexpr int interrupted_exit_status(int signal_number) noexcept { return 128 + signal_number; }

/*!
 * \brief The results are written. On Windows a console close, log-off or
 * shutdown handler waits for this (up to 4.5 s) before it lets the system end
 * the process.
 */
void interrupt_handled() noexcept;

/*!
 * \brief End the program after an interrupted run: POSIX re-raises the
 * signal with its default action, so that a calling shell or script sees
 * that the program was stopped by it (a script then stops too); returns
 * \p status where that is not possible (Windows, or no signal).
 */
int exit_after_interrupt(int status);

} // namespace cli
} // namespace amcset
