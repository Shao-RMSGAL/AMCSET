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

#include "app/unix_signals.h"

#ifndef _WIN32
#include <QSocketNotifier>

#include <atomic>
#include <csignal>
#include <fcntl.h>
#include <signal.h>
#include <sys/socket.h>
#include <unistd.h>
#endif

namespace amcset {
namespace gui {

#ifdef _WIN32

void install_unix_signal_handlers(QObject *, std::function<void(int)>) {}

#else

namespace {

int g_fds[2] = {-1, -1};
std::atomic<int> g_received{0};

extern "C" void on_unix_signal(int sig) {
  // Async-signal-safe only: an atomic, write and _exit.
  if (g_received.fetch_add(1) > 0) _exit(128 + sig);
  const char c = static_cast<char>(sig);
  [[maybe_unused]] const auto n = ::write(g_fds[0], &c, 1);
}

} // namespace

void install_unix_signal_handlers(QObject *parent, std::function<void(int)> on_signal) {
  if (g_fds[0] >= 0) return; // once
  if (::socketpair(AF_UNIX, SOCK_STREAM, 0, g_fds) != 0) return;
  for (int fd : g_fds) ::fcntl(fd, F_SETFD, FD_CLOEXEC);
  auto *notifier = new QSocketNotifier(g_fds[1], QSocketNotifier::Read, parent);
  QObject::connect(notifier, &QSocketNotifier::activated, parent, [notifier, on_signal = std::move(on_signal)] {
    notifier->setEnabled(false);
    char c = 0;
    if (::read(g_fds[1], &c, 1) == 1) on_signal(static_cast<int>(c));
  });
  for (int sig : {SIGINT, SIGTERM, SIGHUP}) {
    struct sigaction sa {};
    sa.sa_handler = on_unix_signal;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    ::sigaction(sig, &sa, nullptr);
  }
}

#endif

} // namespace gui
} // namespace amcset
