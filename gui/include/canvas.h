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

/*! \file canvas.h
 *  \brief A class for drawing particle trajectories in the user interface
 */
#pragma once

// Qt includes
#include <QMutex>
#include <QPen>
#include <QPixmap>
#include <QVector>
#include <QWidget>

class Canvas : public QWidget {
  Q_OBJECT

public:
  Canvas(QWidget *parent = nullptr);

  void addLine(const QLine &line);

  void clear();

protected:
  void updateCanvas();

  void paintEvent(QPaintEvent *event) override;
  QVector<QLine> lines;

  QMutex mutex;
};
