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

// Author note: I am not a UI dev. I just did what I needed to in order to get
// this to work. Forgive any jankyness, and feel free to reimplement this.

#include <QHBoxLayout>
#include <QString>
#include <QVBoxLayout>

#include "amcset_server.h"
#include "config.h"

namespace amcset {
namespace server {

// Dark mode toggle{{{
void Application::toggle_dark_mode() {
#ifndef Q_OS_MACOS
  qApp->setStyle("Fusion");
  auto palette = QPalette();
  auto disabledColor = QColor(127, 127, 127);

  // Slightly inefficient, but easier than manually declaring each variable
  QColor color;
  QColor color_1;
  Qt::GlobalColor primaryColor;
  Qt::GlobalColor secondaryColor;

  if (dark_mode) {
    LOG(INFO) << "Enabling Light Mode";
    color_1 = QColor(200, 200, 200);
    color = QColor(255, 255, 255);
    primaryColor = Qt::black;
    secondaryColor = Qt::white;
    dark_mode = false;
  } else {
    LOG(INFO) << "Enabling Dark Mode";
    color_1 = QColor(18, 18, 18);
    color = QColor(38, 38, 38);
    primaryColor = Qt::white;
    secondaryColor = Qt::black;
    dark_mode = true;
  }
  palette.setColor(QPalette::Window, color);
  palette.setColor(QPalette::WindowText, primaryColor);
  palette.setColor(QPalette::Base, color);
  palette.setColor(QPalette::AlternateBase, color);
  palette.setColor(QPalette::ToolTipBase, secondaryColor);
  palette.setColor(QPalette::ToolTipText, secondaryColor);
  palette.setColor(QPalette::Text, primaryColor);
  palette.setColor(QPalette::Disabled, QPalette::Text, disabledColor);
  palette.setColor(QPalette::Button, color);
  palette.setColor(QPalette::ButtonText, secondaryColor);
  palette.setColor(QPalette::Disabled, QPalette::ButtonText, disabledColor);
  palette.setColor(QPalette::BrightText, Qt::red);
  palette.setColor(QPalette::Link, QColor(42, 130, 218));
  palette.setColor(QPalette::Highlight, QColor(42, 130, 218));
  palette.setColor(QPalette::HighlightedText, primaryColor);
  palette.setColor(QPalette::Disabled, QPalette::HighlightedText,
                   disabledColor);
  // palette.setColor(QPalette::, QPalette::HighlightedText, color_1);

  palette.setColor(QPalette::ButtonText, primaryColor);

  qApp->setPalette(palette);
  qApp->setStyleSheet("QToolTip { color: #000000; background-color: #ffffff; "
                      "border: 1px solid black; }");
#else
  LOG(INFO) << "Dark mode cannot be set. You are on MacOS."
#endif
} // }}}

// MainWindow {{{
MainWindow::MainWindow(Application *app_in) {
  app = app_in;

  // Layout code
  base_frame_ = new QFrame;
  auto base_vbox = new QVBoxLayout;
  base_frame_->setLayout(base_vbox);

  dark_mode_button_ = new QPushButton;
  start_simulation_button_ = new QPushButton;
  base_vbox->addWidget(dark_mode_button_);
  base_vbox->addWidget(start_simulation_button_);

  if (app->get_dark_mode()) {
    dark_mode_button_->setText("Switch to Light Mode");
  } else {
    dark_mode_button_->setText("Switch to Dark Mode");
  }
  connect(dark_mode_button_, &QPushButton::clicked, [&] {
    try {
      app->toggle_dark_mode();
      if (app->get_dark_mode()) {
        dark_mode_button_->setText("Switch to Light Mode");
      } else {
        dark_mode_button_->setText("Switch to Dark Mode");
      }
    } catch (const std::exception &e) {
      QMessageBox::information(this, "Exception handled", e.what());
    }
  });

  start_simulation_button_->setText("Run Simulation");
  connect(start_simulation_button_, &QPushButton::clicked,
          [&] { app->runSimulation(); });

  setCentralWidget(base_frame_);
  setWindowTitle(QString("AMCSET version %1.%2.%3")
                     .arg(amcset_VERSION_MAJOR)
                     .arg(amcset_VERSION_MINOR)
                     .arg(amcset_VERSION_PATCH));
  // resize(label1.sizeHint());
} // }}}

} // namespace server
} // namespace amcset
