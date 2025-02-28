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
#include <QMutex>
#include <QPainter>
#include <QPainterPath>
#include <QTimer>

// Local includes
#include "canvas.h"

Canvas::Canvas(QWidget *parent) : QWidget(parent) {
  QTimer *timer = new QTimer(this);
  // connect(timer, &QTimer::timeout, this, &update);
  timer->start(16); // ~60 FPS
}

void Canvas::updateCanvas() { update(); } // Probably don't need

void Canvas::addLine(const QLine &line) {
  // mutex.lock();
  lines.emplace_back(line);
  update();
  // mutex.unlock();
}

void Canvas::clear() {
  lines.clear();
  update();
}

void Canvas::paintEvent(QPaintEvent *event) {
  QPainter painter(this);
  for (auto &line : lines) {
    painter.drawLine(line);
  }
}

// void Canvas::paintEvent(QPaintEvent *) {
//   // Draw a house

//   // Main house structure (rectangle)
//   painter->setPen(QPen(Qt::black, 2));
//   painter->setBrush(QBrush(QColor(255, 220, 180)));
//   painter->drawRect(300, 200, 200, 150);

//   // Roof (triangle)
//   QPolygon roof;
//   roof << QPoint(280, 200) << QPoint(400, 120) << QPoint(520, 200);
//   painter->setBrush(QBrush(QColor(139, 69, 19)));
//   painter->drawPolygon(roof);

//   // Door
//   painter->setBrush(QBrush(QColor(139, 69, 19)));
//   painter->drawRect(375, 270, 50, 80);

//   // Doorknob
//   painter->setBrush(QBrush(Qt::yellow));
//   painter->drawEllipse(410, 310, 8, 8);

//   // Windows
//   painter->setBrush(QBrush(Qt::lightGray));
//   painter->drawRect(320, 230, 40, 40); // Left window
//   painter->drawRect(440, 230, 40, 40); // Right window

//   // Window crossbars
//   painter->drawLine(320, 250, 360, 250); // Left window horizontal
//   painter->drawLine(340, 230, 340, 270); // Left window vertical
//   painter->drawLine(440, 250, 480, 250); // Right window horizontal
//   painter->drawLine(460, 230, 460, 270); // Right window vertical

//   // Sun
//   painter->setPen(Qt::NoPen);
//   painter->setBrush(QBrush(Qt::yellow));
//   painter->drawEllipse(50, 50, 60, 60);

//   // Ground
//   painter->setPen(QPen(Qt::black, 2));
//   painter->setBrush(QBrush(QColor(34, 139, 34)));
//   painter->drawRect(0, 350, width(), 250);
// }
