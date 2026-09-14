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
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

/*! \file canvas.cpp
 *  \brief A class for drawing particle trajectories in the user interface
 */

#include <QPainter>
#include <QPainterPath>
#include <QTimer>

// Local includes
#include "canvas.h"

Canvas::Canvas(QWidget *parent) : QWidget(parent) {
  QTimer *timer = new QTimer(this);
  connect(timer, &QTimer::timeout, this, [this]() {
    update();
  });
  timer->start(16); // ~60 FPS
}

void Canvas::addLine(const QLine &line) {
  std::lock_guard<std::mutex> lock(mutex);
  lines.emplace_back(line);
}

void Canvas::clear() {
  std::lock_guard<std::mutex> lock(mutex);
  lines.clear();
  update();
}

void Canvas::paintEvent(QPaintEvent *event) {
  Q_UNUSED(event);
  std::lock_guard<std::mutex> lock(mutex);
  QPainter painter(this);
  for (auto &line : lines) {
    painter.drawLine(line);
  }
}
