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
#include <QApplication>
#include <QFontDatabase>
#include <QMap>

// Local includes
#include "canvas.h"

Canvas::Canvas(QWidget *parent) : QWidget(parent) {}

void Canvas::addPoint(double depthAngstrom, double lateralAngstrom,
                      std::size_t generation) {
  addPoints({{depthAngstrom, lateralAngstrom, generation}});
}

void Canvas::addPoints(const QVector<Point> &points) {
  if (points.isEmpty()) {
    return;
  }
  {
    QMutexLocker locker(&mutex);
    points_ += points;
  }
  update();
}

void Canvas::setRangeAngstrom(double rangeAngstrom) {
  QMutexLocker locker(&mutex);
  rangeAngstrom_ = rangeAngstrom > 0.0 ? rangeAngstrom : 1.0;
  update();
}

void Canvas::clear() {
  QMutexLocker locker(&mutex);
  points_.clear();
  update();
}

void Canvas::paintEvent(QPaintEvent *event) {
  Q_UNUSED(event);
  QVector<Point> points;
  double rangeAngstrom = 1.0;
  {
    QMutexLocker locker(&mutex);
    points = points_;
    rangeAngstrom = rangeAngstrom_;
  }

  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing, true);
  const auto colors = palette();
  painter.fillRect(rect(), colors.color(QPalette::Base));

  constexpr double leftMargin = 72.0;
  constexpr double rightMargin = 112.0;
  constexpr double topMargin = 20.0;
  constexpr double bottomMargin = 58.0;
  const QRectF plotArea(leftMargin, topMargin,
                        std::max(1.0, width() - leftMargin - rightMargin),
                        std::max(1.0, height() - topMargin - bottomMargin));

  painter.setFont(QFontDatabase::systemFont(QFontDatabase::GeneralFont));
  const auto textColor = colors.color(QPalette::Text);
  const auto gridColor = colors.color(QPalette::AlternateBase);
  const double halfRange = rangeAngstrom / 2.0;
  for (int tick = 0; tick <= 4; ++tick) {
    const double fraction = tick / 4.0;
    const double x = plotArea.left() + fraction * plotArea.width();
    const double y = plotArea.bottom() - fraction * plotArea.height();
    painter.setPen(QPen(gridColor, 1.0));
    painter.drawLine(QPointF(x, plotArea.top()), QPointF(x, plotArea.bottom()));
    painter.drawLine(QPointF(plotArea.left(), y), QPointF(plotArea.right(), y));
    painter.setPen(textColor);
    painter.drawText(QRectF(x - 35.0, plotArea.bottom() + 4.0, 70.0, 18.0),
                     Qt::AlignHCenter | Qt::AlignTop,
                     QString::number(fraction * rangeAngstrom, 'g', 4));
    painter.drawText(QRectF(0.0, y - 9.0, leftMargin - 8.0, 18.0),
                     Qt::AlignRight | Qt::AlignVCenter,
                     QString::number(fraction * rangeAngstrom_ - halfRange, 'g', 4));
  }
  painter.setPen(QPen(colors.color(QPalette::Mid), 1.2));
  painter.drawRect(plotArea);
  painter.setPen(textColor);
  painter.drawText(QRectF(leftMargin, height() - 24.0,
                          width() - leftMargin - rightMargin, 20.0),
                   Qt::AlignCenter, tr("Depth (Å)"));
  painter.save();
  painter.rotate(-90.0);
  painter.drawText(QRectF(-height(), 0.0, height(), 20.0), Qt::AlignCenter,
                   tr("Lateral displacement (Å)"));
  painter.restore();

  QMap<std::size_t, QVector<QPointF>> generationPoints;
  const double inverseRange = 1.0 / rangeAngstrom;
  for (const auto &point : points) {
    const double xFraction = point.depthAngstrom * inverseRange;
    const double yFraction = (point.lateralAngstrom + halfRange) * inverseRange;
    if (xFraction < 0.0 || xFraction > 1.0 || yFraction < 0.0 || yFraction > 1.0) {
      continue;
    }
    generationPoints[point.generation].append(
        QPointF(plotArea.left() + xFraction * plotArea.width(),
                plotArea.bottom() - yFraction * plotArea.height()));
  }

  const bool dark = colors.color(QPalette::Window).lightness() < 128;
  const QVector<QColor> generationColors = dark
      ? QVector<QColor>{QColor("#56B4E9"), QColor("#F28E2B"), QColor("#4CC9A6"),
                        QColor("#FF77A8"), QColor("#FFE066"), QColor("#7AA2FF"),
                        QColor("#C792EA")}
      : QVector<QColor>{QColor("#0072B2"), QColor("#D55E00"), QColor("#009E73"),
                        QColor("#CC79A7"), QColor("#A66E00"), QColor("#007F9E"),
                        QColor("#6B5B95")};
  int legendRow = 0;
  painter.save();
  painter.setClipRect(plotArea);
  for (auto it = generationPoints.cbegin(); it != generationPoints.cend(); ++it) {
    const auto color = generationColors.at(static_cast<int>(it.key() % generationColors.size()));
    painter.setPen(QPen(color, 2.0));
    painter.setBrush(color);
    painter.drawPoints(it.value().constData(), it.value().size());
  }
  painter.restore();

  for (auto it = generationPoints.cbegin(); it != generationPoints.cend(); ++it) {
    const auto color = generationColors.at(static_cast<int>(it.key() % generationColors.size()));
    const double legendY = topMargin + legendRow * 20.0;
    painter.setPen(QPen(color, 2.0));
    painter.setBrush(color);
    painter.drawEllipse(QPointF(width() - rightMargin + 12.0, legendY + 8.0), 3.0, 3.0);
    painter.setPen(textColor);
    painter.drawText(QPointF(width() - rightMargin + 21.0, legendY + 12.0),
                     tr("Generation %1").arg(it.key()));
    ++legendRow;
  }
}
