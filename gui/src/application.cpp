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
#include <QStyle>
#include <QStyleFactory>

// Logging libraries
#include <glog/logging.h>

// Local library includes
#include "application.h"

namespace amcset {
namespace gui {

Application::Application(amcset::common::Simulation *simulation)
    : QObject(nullptr), simulation_(simulation) {
  systemStyleName_ = qApp->style()->objectName();
  systemPalette_ = qApp->palette();
  settings_.beginGroup("MainWindow");
  const auto storedMode = settings_.value("themeMode", 0).toInt();
  themeMode_ = storedMode >= static_cast<int>(ThemeMode::System) &&
                       storedMode <= static_cast<int>(ThemeMode::Dark)
                   ? static_cast<ThemeMode>(storedMode)
                   : ThemeMode::System;
  setThemeMode(themeMode_);
}

Application::~Application() {
  // Save settings on destruction
  settings_.sync();
}

void Application::toggleDarkMode() {
  setThemeMode(getDarkMode() ? ThemeMode::Light : ThemeMode::Dark);
}

void Application::setThemeMode(ThemeMode mode) {
  themeMode_ = mode;
  settings_.setValue("themeMode", static_cast<int>(themeMode_));

  if (themeMode_ == ThemeMode::System) {
    qApp->setStyle(QStyleFactory::create(systemStyleName_));
    qApp->setPalette(systemPalette_);
    qApp->setStyleSheet(QString());
    return;
  }

  qApp->setStyle(QStyleFactory::create(QStringLiteral("Fusion")));
  auto palette = QPalette();
  const auto disabledColor = QColor(127, 127, 127);
  const bool dark = themeMode_ == ThemeMode::Dark;
  const auto windowColor = dark ? QColor(38, 38, 38) : QColor(242, 242, 242);
  const auto baseColor = dark ? QColor(25, 25, 25) : QColor(255, 255, 255);
  const auto alternateColor = dark ? QColor(45, 45, 45) : QColor(248, 248, 248);
  const auto textColor = dark ? QColor(238, 238, 238) : QColor(32, 32, 32);
  const auto buttonColor = dark ? QColor(52, 52, 52) : QColor(232, 232, 232);
  const auto tooltipColor = dark ? QColor(48, 48, 48) : QColor(255, 255, 220);

  palette.setColor(QPalette::Window, windowColor);
  palette.setColor(QPalette::WindowText, textColor);
  palette.setColor(QPalette::Base, baseColor);
  palette.setColor(QPalette::AlternateBase, alternateColor);
  palette.setColor(QPalette::ToolTipBase, tooltipColor);
  palette.setColor(QPalette::ToolTipText, textColor);
  palette.setColor(QPalette::Text, textColor);
  palette.setColor(QPalette::Disabled, QPalette::Text, disabledColor);
  palette.setColor(QPalette::Button, buttonColor);
  palette.setColor(QPalette::ButtonText, textColor);
  palette.setColor(QPalette::Disabled, QPalette::ButtonText, disabledColor);
  palette.setColor(QPalette::BrightText, Qt::red);
  palette.setColor(QPalette::Link, QColor(42, 130, 218));
  palette.setColor(QPalette::Highlight, QColor(42, 130, 218));
  palette.setColor(QPalette::HighlightedText, Qt::white);
  palette.setColor(QPalette::Disabled, QPalette::HighlightedText,
                   disabledColor);

  qApp->setPalette(palette);
  qApp->setStyleSheet(QString());
}

bool Application::getDarkMode() const noexcept {
  if (themeMode_ == ThemeMode::Dark) {
    return true;
  }
  if (themeMode_ == ThemeMode::Light) {
    return false;
  }
  return qApp->palette().color(QPalette::Window).lightness() < 128;
}

Application::ThemeMode Application::getThemeMode() const noexcept {
  return themeMode_;
}

bool Application::showLargeCsvWarning() const {
  return settings_.value("showLargeCsvWarning", true).toBool();
}

void Application::setShowLargeCsvWarning(bool show) {
  settings_.setValue("showLargeCsvWarning", show);
}

bool Application::showOutputOverwriteWarning() const {
  return settings_.value("showOutputOverwriteWarning", true).toBool();
}

void Application::setShowOutputOverwriteWarning(bool show) {
  settings_.setValue("showOutputOverwriteWarning", show);
}

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
