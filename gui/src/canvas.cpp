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

// Qt includes
#include <QPainter>
#include <QPainterPath>
#include <QTimer>

// Local includes
#include "canvas.h"

Canvas::Canvas(QWidget *parent) : QWidget(parent) {
  renderCache = QPixmap(size());
  renderCache.fill(Qt::transparent);
}

void Canvas::flushLinesToBitmap() {
  if (lines.isEmpty()) {
    return;
  }

  if (renderCache.size() != size()) {
    renderCache = QPixmap(size());
    renderCache.fill(Qt::transparent);
  }

  QPainter painter(&renderCache);
  for (const auto &line : lines) {
    painter.drawLine(line);
  }
  lines.clear();
}

void Canvas::addLine(const QLine &line) {
  {
    QMutexLocker locker(&mutex);
    lines.emplace_back(line);
    if (lines.size() >= kBitmapFlushThreshold) {
      flushLinesToBitmap();
    }
  }

  update();
}

void Canvas::clear() {
  QMutexLocker locker(&mutex);
  lines.clear();
  renderCache = QPixmap(size());
  renderCache.fill(Qt::transparent);
  update();
}

void Canvas::paintEvent(QPaintEvent *event) {
  Q_UNUSED(event);
  QMutexLocker locker(&mutex);

  QPainter painter(this);
  if (!renderCache.isNull()) {
    painter.drawPixmap(0, 0, renderCache);
  }

  for (const auto &line : lines) {
    painter.drawLine(line);
  }
}
