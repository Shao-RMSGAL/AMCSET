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

/*! \file unix_signals.h
 *  \brief SIGINT, SIGTERM and SIGHUP in the GUI (POSIX).
 */

#pragma once

#include <QObject>

#include <functional>

namespace amcset {
namespace gui {

/*!
 * \brief Calls \p on_signal(signal number) on the GUI thread when the
 * program receives SIGINT, SIGTERM or SIGHUP (a desktop session ending,
 * `kill`, Ctrl+C in the terminal that started it). The handler only writes
 * the number to a socket pair that a QSocketNotifier reads (the pattern of
 * Qt's "Calling Qt Functions From Unix Signal Handlers"); a second signal
 * ends the program at once. Nothing happens on Windows, whose GUI programs
 * get no signals (the session's end arrives as QGuiApplication::commitDataRequest).
 */
void install_unix_signal_handlers(QObject *parent, std::function<void(int)> on_signal);

} // namespace gui
} // namespace amcset
