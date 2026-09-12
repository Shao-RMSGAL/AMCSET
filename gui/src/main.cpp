// Copyright 2024, Texas A&M University
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

#include <iostream>
#include <qboxlayout.h>
#include <qlabel.h>

#include <thread>

#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/math/constants/constants.hpp>

#include <QApplication>
#include <QDebug>
#include <QFile>
#include <QFrame>
#include <QLabel>
#include <QMainWindow>
#include <QTemporaryDir>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "amcset_common.h"
#include "amcset_utilities.h"
#include "application.h"
#include "config.h"
#include "mainwindow.h"

DEFINE_bool(gui, true, "Enable or disable the graphical user interface (GUI)");

using namespace amcset::common;
// Wow
int main(int argc, char *argv[]) {

  // TODO: See how this behaves under Windows
  QTemporaryDir tempDir;
  if (tempDir.isValid()) {
    // Extract fonts.conf from resources to temporary directory
    QFile resFile(":resources/fonts.conf");
    if (resFile.open(QIODevice::ReadOnly)) {
      QString configPath = tempDir.filePath("fonts.conf");
      QFile outFile(configPath);
      if (outFile.open(QIODevice::WriteOnly)) {
        outFile.write(resFile.readAll());
        outFile.close();

        // Set Fontconfig environment variables
        qputenv("FONTCONFIG_FILE", configPath.toUtf8());
        qputenv("FONTCONFIG_PATH", tempDir.path().toUtf8());
      } else {
        qWarning() << "Failed to write fonts.conf to temporary directory";
      }
    } else {
      qWarning() << "Failed to open embedded fonts.conf resource";
    }
  } else {
    qWarning() << "Failed to create temporary directory for Fontconfig";
  }

  try { // {{{
    setlocale(LC_ALL, "C.UTF-8");
    qputenv("LC_ALL", "C.UTF-8");
    qputenv("LANG", "C.UTF-8");

    // Initialize logger
    google::InitGoogleLogging(argv[0]);
    FLAGS_alsologtostderr = 1;

    std::string usage_message = "AMCSET. Call " + std::string(argv[0]);
    gflags::SetUsageMessage(usage_message);

    std::ostringstream version;
    version << amcset_VERSION_MAJOR << '.' << amcset_VERSION_MINOR << '.'
            << amcset_VERSION_PATCH;
    gflags::SetVersionString(version.str());
    // Process flags
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    // PERF: Add an option to use GPU (Use Kokkos for cross-platform
    // compatability. Note that it is not available in conan. Either make a
    // conan recipe for it or use another installation method.

    if (FLAGS_gui) {
      LOG(INFO) << "Starting AMCSET GUI...";
      QApplication app(argc, argv);
      amcset::gui::Application application;
      VLOG(1) << "Initializing window";
      amcset::gui::MainWindow mainWindow(&application);
      VLOG(1) << "Executing application";
      return app.exec();
    } else {
      LOG(INFO) << "Starting AMCSET in CLI mode...";
      LOG(ERROR) << "CLI implemementation not complete";
      LOG(INFO) << "...Quitting AMCSET CLI.";
      return 0;
    }
  } catch (const std::exception &e) {
    LOG(FATAL) << (e.what());
  } // }}}
}
