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

/*! \file interrupts.cpp
 *  \brief Signal and console-event handlers of the command-line program.
 */

#include "interrupts.h"

#include <csignal>
#include <cstdio>
#include <iostream>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <signal.h>
#include <unistd.h>
#endif

namespace amcset {
namespace cli {
namespace {

// The handlers touch only these lock-free atomics (and, POSIX, _exit;
// Windows, TerminateProcess and Sleep): nothing that may lock.
static_assert(std::atomic<bool>::is_always_lock_free && std::atomic<int>::is_always_lock_free &&
                  std::atomic<std::atomic<bool> *>::is_always_lock_free,
              "the signal handlers need lock-free atomics");
std::atomic<std::atomic<bool> *> g_flag{nullptr};
std::atomic<int> g_signal{0};
std::atomic<bool> g_handled{false};

#ifdef _WIN32

BOOL WINAPI on_console_event(DWORD event) {
  std::atomic<bool> *flag = g_flag.load();
  if (flag == nullptr) return FALSE; // the default: end the program
  switch (event) {
  case CTRL_C_EVENT:
  case CTRL_BREAK_EVENT:
    if (flag->exchange(true)) TerminateProcess(GetCurrentProcess(), interrupted_exit_status(SIGINT));
    g_signal.store(SIGINT);
    return TRUE;
  case CTRL_CLOSE_EVENT:
  case CTRL_LOGOFF_EVENT:
  case CTRL_SHUTDOWN_EVENT:
    // The system ends the process when this returns (or after 5 s for a
    // closed console window): until then the run stops and writes its
    // partial results.
    g_signal.store(SIGINT);
    flag->store(true);
    for (int i = 0; i < 45 && !g_handled.load(); ++i) Sleep(100);
    return TRUE;
  default: return FALSE;
  }
}

#else

struct sigaction g_previous[3];
constexpr int g_signals[3] = {SIGINT, SIGTERM, SIGHUP};

extern "C" void on_signal(int sig) {
  std::atomic<bool> *flag = g_flag.load();
  if (flag == nullptr || flag->exchange(true)) {
    // The second signal (or one between runs): quit at once. Only
    // async-signal-safe calls here.
    _exit(interrupted_exit_status(sig));
  }
  g_signal.store(sig);
}

#endif

} // namespace

InterruptGuard::InterruptGuard(std::atomic<bool> &flag) {
  flag.store(false);
  g_signal.store(0);
  g_handled.store(false);
  g_flag.store(&flag);
#ifdef _WIN32
  SetConsoleCtrlHandler(on_console_event, TRUE);
#else
  for (int i = 0; i < 3; ++i) {
    struct sigaction sa {};
    sa.sa_handler = on_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART; // (interrupted reads and writes resume)
    sigaction(g_signals[i], &sa, &g_previous[i]);
  }
#endif
}

InterruptGuard::~InterruptGuard() {
#ifdef _WIN32
  SetConsoleCtrlHandler(on_console_event, FALSE);
#else
  for (int i = 0; i < 3; ++i) sigaction(g_signals[i], &g_previous[i], nullptr);
#endif
  g_flag.store(nullptr);
}

int interrupt_signal() noexcept { return g_signal.load(); }

void interrupt_handled() noexcept { g_handled.store(true); }

int exit_after_interrupt(int status) {
  const int sig = g_signal.load();
  interrupt_handled();
  if (sig == 0 || status != interrupted_exit_status(sig)) return status;
#ifndef _WIN32
  std::cout.flush();
  std::cerr.flush();
  std::fflush(nullptr);
  struct sigaction sa {};
  sa.sa_handler = SIG_DFL;
  sigemptyset(&sa.sa_mask);
  sigaction(sig, &sa, nullptr);
  sigset_t set;
  sigemptyset(&set);
  sigaddset(&set, sig);
  sigprocmask(SIG_UNBLOCK, &set, nullptr);
  raise(sig);
#endif
  return status;
}

} // namespace cli
} // namespace amcset
