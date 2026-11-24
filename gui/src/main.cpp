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
 *  \brief Entry point of the AMCSET GUI.
 *
 * Usage: amcset-gui [scenario.json | TRIM.IN | results.json]
 * The environment variable AMCSET_RHI (opengl, vulkan, d3d11, d3d12, metal)
 * selects the graphics API of the path views.
 */

#include "app/guard.h"
#include "app/main_window.h"
#include "config.h"
#include "diagnostics.h"

#include <QApplication>
#include <QCommandLineParser>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QMessageBox>
#include <QSettings>
#include <QStandardPaths>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <thread>

namespace {

//! Qt's own warnings and errors go into the diagnostic log (and, as before,
//! to stderr); a qFatal ends in a crash report.
void qt_message(QtMsgType type, const QMessageLogContext &, const QString &msg) {
  using amcset::diag::Level;
  Level level = Level::info;
  switch (type) {
  case QtDebugMsg: level = Level::debug; break;
  case QtInfoMsg: level = Level::info; break;
  case QtWarningMsg: level = Level::warning; break;
  case QtCriticalMsg:
  case QtFatalMsg: level = Level::error; break;
  }
  const QByteArray text = msg.toLocal8Bit();
  if (amcset::diag::enabled(level)) amcset::diag::write(level, "[Qt] " + text.toStdString());
  if (type != QtDebugMsg) std::fprintf(stderr, "%s\n", text.constData());
  if (type == QtFatalMsg) std::abort();
}

//! The GUI's log: <app data>/logs/amcset-gui.log, the previous one kept as
//! amcset-gui.old.log once it passes 5 MB.
std::filesystem::path gui_log_file() {
  const QString dir = QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) +
                      QStringLiteral("/logs");
  QDir().mkpath(dir);
  const QString file = dir + QStringLiteral("/amcset-gui.log");
  const QFileInfo info(file);
  if (info.exists() && info.size() > 5 * 1024 * 1024) {
    const QString old = dir + QStringLiteral("/amcset-gui.old.log");
    QFile::remove(old);
    QFile::rename(file, old);
  }
  return std::filesystem::path(file.toStdWString());
}

} // namespace

int main(int argc, char *argv[]) {
  // Crash reports from the very start (the log file follows below).
  amcset::diag::install_crash_handlers(argc > 0 ? argv[0] : "amcset-gui");
  // QApplication with an exception-catching notify() (app/guard.h).
  amcset::gui::Application app(argc, argv);
  QCoreApplication::setOrganizationName(QStringLiteral("AMCSET"));
  QCoreApplication::setApplicationName(QStringLiteral("AMCSET"));
  QCoreApplication::setApplicationVersion(QStringLiteral(AMCSET_VERSION));
  {
    amcset::diag::Options d;
    d.program = "amcset-gui";
    d.log_file = gui_log_file();
    d.level = QSettings().value(QStringLiteral("diagnostics/detailed"), false).toBool()
                  ? amcset::diag::Level::debug
                  : amcset::diag::Level::info;
    amcset::diag::initialize(argc > 0 ? argv[0] : "amcset-gui", d);
    qInstallMessageHandler(qt_message);
    QStringList args;
    for (int i = 0; i < argc; ++i) args << QString::fromLocal8Bit(argv[i]);
    AMCSET_LOG(info) << "amcset-gui " AMCSET_VERSION " started (Qt " << qVersion() << "): "
                     << args.join(QLatin1Char(' ')).toStdString();
  }

  QCommandLineParser parser;
  parser.setApplicationDescription(QStringLiteral("AMCSET: Monte Carlo transport of ions and electrons in matter"));
  parser.addHelpOption();
  parser.addVersionOption();
  parser.addPositionalArgument(QStringLiteral("file"),
                               QStringLiteral("Scenario (.json), TRIM.IN or results.json to open"));
  const QCommandLineOption run_opt(QStringLiteral("run"), QStringLiteral("Start the calculation at once."));
  const QCommandLineOption shots_opt(QStringLiteral("screenshots"),
                                     QStringLiteral("After the calculation, save images of the views to <dir>."),
                                     QStringLiteral("dir"));
  const QCommandLineOption quit_opt(QStringLiteral("quit"), QStringLiteral("Quit when the calculation ends."));
  const QCommandLineOption timeout_opt(QStringLiteral("timeout"),
                                       QStringLiteral("With --run: exit with code 3 after <seconds>."),
                                       QStringLiteral("seconds"), QStringLiteral("0"));
  parser.addOption(run_opt);
  parser.addOption(timeout_opt);
  parser.addOption(shots_opt);
  parser.addOption(quit_opt);
  parser.process(app);

  const bool unattended = parser.isSet(run_opt) || parser.isSet(shots_opt);
  const QString shots = parser.value(shots_opt);
  int rc = 0;
  // Exceptions in event handlers are caught by Application::notify and the
  // slot guards; this catches what escapes building or closing the window.
  try {
    amcset::gui::MainWindow window;
    window.show();
    const QStringList files = parser.positionalArguments();
    if (!files.isEmpty()) window.open_path(files.front());
    if (unattended) window.run_unattended(shots, parser.isSet(quit_opt), parser.value(timeout_opt).toInt());
    rc = app.exec();
    if (unattended) {
      qInfo().noquote() << QStringLiteral("Event loop ended (exit code %1); closing.").arg(rc);
      // An unattended run must end: if closing the window and the graphics
      // device takes longer than 60 s, say so and leave with code 10 + rc.
      std::thread([rc] {
        std::this_thread::sleep_for(std::chrono::seconds(60));
        std::fprintf(stderr, "Shutdown did not finish within 60 s.\n");
        std::fflush(stderr);
        std::_Exit(10 + rc);
      }).detach();
    }
  } catch (...) {
    const QString text = amcset::gui::exception_text(std::current_exception());
    AMCSET_LOG(error) << "amcset-gui failed: " << text.toStdString();
    std::fprintf(stderr, "amcset-gui: %s\n", text.toLocal8Bit().constData());
    if (!unattended) QMessageBox::critical(nullptr, QStringLiteral("AMCSET"), text);
    rc = 2;
  }
  if (unattended) {
    qInfo().noquote() << QStringLiteral("Window closed.");
    // The exit code, for scripts that cannot read it from the process.
    if (!shots.isEmpty()) {
      QFile f(QDir(shots).filePath(QStringLiteral("exit_code.txt")));
      if (f.open(QIODevice::WriteOnly | QIODevice::Text)) f.write(QByteArray::number(rc) + '\n');
    }
  }
  AMCSET_LOG(info) << "amcset-gui exits with code " << rc;
  return rc;
}
