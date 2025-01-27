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

// Author note: I am not a UI dev. I just did what I needed to in order to get
// this to work. Forgive any jankyness, and feel free to reimplement this.

/*! \file amcset_server.h
 * \brief Header file that declares the UI elements for the user interface.
 *
 *  This file contains Qt definitions for UI elements that are presented
 *  to the user when they open AMCSET.
 */

#pragma once

// Qt libraries
#include <QApplication>
#include <QBoxLayout>
#include <QCheckBox>
#include <QCursor>
#include <QFrame>
#include <QLabel>
#include <QListWidget>
#include <QMainWindow>
#include <QMessageBox>
#include <QPushButton>
#include <QRadioButton>
#include <QScreen>

// Logging libraries
#include <glog/logging.h>

// Local libraries
#include "amcset_common.h"

namespace amcset {
namespace server {

//! The application controls the state of the entire application.
class Application : public QApplication {

public:
  /*!
   * \brief Construct an application with commandline arguments;
   *
   * This effectively serves as a main function for a Qt application. It starts
   * an application, which can then be used to open windows.
   *
   * \param argc Number of command line arguments.
   * \param argv Strings of the command line arguments.
   */
  Application(int &argc, char **argv, amcset::common::Simulation *simulation)
      : QApplication(argc, argv), simulation_(simulation) {
    VLOG(1) << "Constructing application";
    dark_mode = false;
    toggle_dark_mode();
  }

  /*!
   * \brief Toggle dark mode.
   *
   * This function sets background colors to a dark color and text colors to
   * white when enabling dark mode. Vice versa when enabling light mode.
   */
  void toggle_dark_mode();

  /*!
   * \brief return whether dark mode is enabled or not.
   */
  bool get_dark_mode() { return dark_mode; };

  /*!
   * \brief Start a simulation.
   *
   * Once a simulation is initialized, this will run the simulation. After
   * a period of time, the simulation will complete and the results will
   * become available inside the 'simulation' class member.
   */
  void runSimulation() {
    LOG(INFO) << "Starting simulation";
    simulation_->run_simulation();
  };

private:
  bool dark_mode;
  amcset::common::Simulation *simulation_; // HACK Remove raw pointer
};

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  Application *app;
  MainWindow(Application *app_in);

private:
  QFrame *base_frame_;
  // QBoxLayout *base_layout_;
  // QFrame testZone{&frame};
  // QPushButton generateHandledExceptionButton{&frame};
  // QPushButton generateExceptionButton{&frame};
  // QPushButton generateUnknownExceptionButton{&frame};
  // int button1Clicked = 0;
  // int button2Clicked = 0;
  // QLabel label1{&frame};
  // QLabel label{&frame};
  QPushButton *dark_mode_button_;
  // QListWidget listWidget{&frame};
  // QRadioButton radioButton1{&frame};
  // QRadioButton radioButton2{&frame};
  // QCheckBox checkBox1{&frame};
  // QCheckBox checkBox2{&frame};
  QPushButton *start_simulation_button_;
};
// common::Simulation start_simulation();

} // namespace server
} // namespace amcset
