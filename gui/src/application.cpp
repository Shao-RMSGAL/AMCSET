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

// Qt includes
#include <QPalette>

// Logging libraries
#include <glog/logging.h>

// Local library includes
#include "application.h"

namespace amcset {
namespace gui {

Application::Application(amcset::common::Simulation *simulation)
    : QObject(nullptr), simulation_(simulation) {
  // Initialize settings
  darkMode_ = true;
  settings_.beginGroup("MainWindow");
}

Application::~Application() {
  // Save settings on destruction
  settings_.sync();
}

void Application::toggleDarkMode() {
#ifndef Q_OS_MACOS
  qApp->setStyle("Fusion");
  auto palette = QPalette();
  auto disabledColor = QColor(127, 127, 127);

  // Slightly inefficient, but easier than manually declaring each variable
  QColor color;
  QColor color_1;
  Qt::GlobalColor primaryColor;
  Qt::GlobalColor secondaryColor;

  if (darkMode_) {
    LOG(INFO) << "Enabling Light Mode";
    color_1 = QColor(200, 200, 200);
    color = QColor(255, 255, 255);
    primaryColor = Qt::black;
    secondaryColor = Qt::white;
  } else {
    LOG(INFO) << "Enabling Dark Mode";
    color_1 = QColor(18, 18, 18);
    color = QColor(38, 38, 38);
    primaryColor = Qt::white;
    secondaryColor = Qt::black;
  }

  darkMode_ = !darkMode_;

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
  qApp->setStyleSheet("QToolTip { color: #000000; background-color:#ffffff;"
                      "border: 1px solid black; }");
#else
  LOG(INFO) << "Dark mode cannot be set. You are on MacOS."
#endif
}
// }}}

bool Application::getDarkMode() const noexcept { return darkMode_; };

void Application::runSimulation() { simulation_->run_simulation(); }

void Application::saveWindowSize(const QSize &size) {
  settings_.setValue("size", size);
}

void Application::saveWindowGeometry(const QByteArray &geometry) {
  settings_.setValue("geometry", geometry);
}

void Application::saveFullscreenState(bool isFullscreen) {
  settings_.setValue("fullscreen", isFullscreen);
}

QSize Application::loadWindowSize() const {
  return settings_.value("size", QSize(1920, 1080)).toSize();
}

QByteArray Application::loadWindowGeometry() const {
  return settings_.value("geometry", QByteArray()).toByteArray();
}

bool Application::loadFullscreenState() const {
  return settings_.value("fullscreen", false).toBool();
}
} // namespace gui
} // namespace amcset
