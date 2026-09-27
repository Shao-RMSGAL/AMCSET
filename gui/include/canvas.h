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

/*! \file canvas.h
 *  \brief A class for drawing particle trajectories in the user interface
 */
#pragma once

// Qt includes
#include <QMutex>
#include <QSet>
#include <QVector>
#include <QWidget>

class Canvas : public QWidget {
  Q_OBJECT

public:
  enum class Projection { XY, YZ, XZ };

  struct Point {
    double xAngstrom;
    double yAngstrom;
    double zAngstrom;
    std::size_t generation;
  };

  Canvas(QWidget *parent = nullptr);

  void addPoint(double depthAngstrom, double lateralAngstrom,
                std::size_t generation);
  void addPoints(const QVector<Point> &points);
  void setProjection(Projection projection);
  void setGenerationFilter(const QSet<std::size_t> &generations,
                          bool showAllGenerations);
  void setRangeAngstrom(double rangeAngstrom);

  void clear();

protected:
  void paintEvent(QPaintEvent *event) override;

private:
  QVector<Point> points_;
  double rangeAngstrom_ = 1.0;
  Projection projection_ = Projection::XZ;
  QSet<std::size_t> visibleGenerations_;
  bool showAllGenerations_ = true;
  QMutex mutex;
};
