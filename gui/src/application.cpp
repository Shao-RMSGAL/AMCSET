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

// TODO: Modify to do proper window management (need to figure out what to do
// with this)

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

// void Application::runSimulation() { simulation_->run_simulation(); }

void Application::saveWindowSize(const QSize &size) {
  settings_.setValue("size", size);
}

void Application::saveWindowGeometry(const QByteArray &geometry) {
  settings_.setValue("geometry", geometry);
}

void Application::saveFullscreenState(bool isFullscreen) {
  settings_.setValue("fullscreen", isFullscreen);
}

void Application::saveMainWindowState(const QVariantMap &state) {
  settings_.beginGroup("mainWindowStateData");
  settings_.remove("");

  for (auto it = state.cbegin(); it != state.cend(); ++it) {
    if (it.key() != "layers") {
      settings_.setValue("fields/" + it.key(), it.value());
    }
  }

  const auto layers = state.value("layers").toList();
  settings_.beginWriteArray("layers", layers.size());
  for (qsizetype layerIndex = 0; layerIndex < layers.size(); ++layerIndex) {
    settings_.setArrayIndex(static_cast<int>(layerIndex));
    const auto layer = layers.at(layerIndex).toMap();
    for (auto it = layer.cbegin(); it != layer.cend(); ++it) {
      if (it.key() != "composition") {
        settings_.setValue(it.key(), it.value());
      }
    }

    const auto composition = layer.value("composition").toList();
    settings_.beginWriteArray("composition", composition.size());
        for (qsizetype elementIndex = 0; elementIndex < composition.size();
          ++elementIndex) {
      settings_.setArrayIndex(static_cast<int>(elementIndex));
      const auto element = composition.at(elementIndex).toMap();
      for (auto it = element.cbegin(); it != element.cend(); ++it) {
        settings_.setValue(it.key(), it.value());
      }
    }
    settings_.endArray();
  }
  settings_.endArray();
  settings_.endGroup();

  settings_.remove("mainWindowState");
  settings_.sync();
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

QVariantMap Application::loadMainWindowState() const {
  auto &settings = settings_;
  settings.beginGroup("mainWindowStateData");
  const auto hasStructuredState = settings.childGroups().contains("layers") ||
                                  settings.childGroups().contains("fields");
  if (!hasStructuredState) {
    settings.endGroup();
    return settings_.value("mainWindowState", QVariantMap()).toMap();
  }

  QVariantMap state;
  settings.beginGroup("fields");
  for (const auto &key : settings.childKeys()) {
    state.insert(key, settings.value(key));
  }
  settings.endGroup();

  QVariantList layers;
  const auto layerCount = settings.beginReadArray("layers");
  for (int layerIndex = 0; layerIndex < layerCount; ++layerIndex) {
    settings.setArrayIndex(layerIndex);
    QVariantMap layer;
    for (const auto &key : settings.childKeys()) {
      layer.insert(key, settings.value(key));
    }

    QVariantList composition;
    const auto elementCount = settings.beginReadArray("composition");
    for (int elementIndex = 0; elementIndex < elementCount; ++elementIndex) {
      settings.setArrayIndex(elementIndex);
      QVariantMap element;
      for (const auto &key : settings.childKeys()) {
        element.insert(key, settings.value(key));
      }
      composition.append(element);
    }
    settings.endArray();
    layer.insert("composition", composition);
    layers.append(layer);
  }
  settings.endArray();
  settings.endGroup();

  state.insert("layers", layers);
  return state;
}
} // namespace gui
} // namespace amcset
