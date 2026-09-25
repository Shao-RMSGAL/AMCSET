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

#pragma once

// Qt library includes
#include <QApplication>
#include <QObject>
#include <QSettings>
#include <QSize>
#include <QVariant>

// Local includes
#include "amcset_common.h"

namespace amcset {
namespace gui {

/*!
 * \brief A custom object for storing Application information.
 */
class Application : public QObject {
  Q_OBJECT

public:
  explicit Application(amcset::common::Simulation *simulation = nullptr);
  ~Application();

  void toggleDarkMode();
  bool getDarkMode() const noexcept;
  // void runSimulation();

  void saveWindowSize(const QSize &size);
  void saveWindowGeometry(const QByteArray &geometry);
  void saveFullscreenState(bool isFullscreen);
  void saveMainWindowState(const QVariantMap &state);

  QSize loadWindowSize() const;
  QByteArray loadWindowGeometry() const;
  bool loadFullscreenState() const;
  QVariantMap loadMainWindowState() const;

private:
  bool darkMode_;
  amcset::common::Simulation *simulation_;
  mutable QSettings settings_;
};
} // namespace gui
} // namespace amcset
