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
#include <QVector>
#include <QWidget>

class Canvas : public QWidget {
  Q_OBJECT

public:
  struct Point {
    double depthAngstrom;
    double lateralAngstrom;
    std::size_t generation;
  };

  Canvas(QWidget *parent = nullptr);

  void addPoint(double depthAngstrom, double lateralAngstrom,
                std::size_t generation);
  void addPoints(const QVector<Point> &points);
  void setRangeAngstrom(double rangeAngstrom);

  void clear();

protected:
  void paintEvent(QPaintEvent *event) override;

private:
  QVector<Point> points_;
  double rangeAngstrom_ = 1.0;
  QMutex mutex;
};
