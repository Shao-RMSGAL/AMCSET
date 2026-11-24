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

/*! \file plot_widget.h
 *  \brief A light 2D plot (lines, histogram steps, points) for result tables.
 *
 * Drag with the left button to zoom into a rectangle, use the wheel to zoom
 * the horizontal axis, double-click to reset; the context menu switches log
 * scales, copies the data or saves the image.
 */

#pragma once

#include <QColor>
#include <QString>
#include <QVector>
#include <QWidget>

#include <optional>

namespace amcset {
namespace gui {

struct PlotSeries {
  enum class Style { line, steps, points };
  QString name;
  QVector<double> x;
  QVector<double> y;
  QColor color;
  Style style = Style::steps;
  bool fill = false;
};

class PlotWidget : public QWidget {
  Q_OBJECT

public:
  explicit PlotWidget(QWidget *parent = nullptr);

  void set_title(const QString &title);
  void set_axis_labels(const QString &x, const QString &y);
  //! Replace the data; keeps the zoom if \p keep_zoom.
  void set_series(QVector<PlotSeries> series, bool keep_zoom = false);
  void set_log_x(bool on);
  void set_log_y(bool on);
  bool log_x() const noexcept { return log_x_; }
  bool log_y() const noexcept { return log_y_; }
  //! A vertical line with a label (e.g. the surface binding energy).
  void set_markers(QVector<QPair<double, QString>> markers);
  void reset_zoom();

  //! The data as tab-separated text (one block per series).
  QString data_text() const;
  //! Draw the plot into an image of \p size logical pixels.
  QImage render_image(const QSize &size, qreal device_pixel_ratio = 2.0) const;
  //! PNG/JPEG/BMP by extension. \return false on failure.
  bool save_image(const QString &path) const;
  const QVector<PlotSeries> &series() const noexcept { return series_; }

signals:
  void hover_text(const QString &text);
  void log_scale_changed(bool log_x, bool log_y);

protected:
  void paintEvent(QPaintEvent *) override;
  void mousePressEvent(QMouseEvent *e) override;
  void mouseMoveEvent(QMouseEvent *e) override;
  void mouseReleaseEvent(QMouseEvent *e) override;
  void mouseDoubleClickEvent(QMouseEvent *e) override;
  void wheelEvent(QWheelEvent *e) override;
  void contextMenuEvent(QContextMenuEvent *e) override;
  void leaveEvent(QEvent *e) override;

private:
  struct Range {
    double x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  };
  void paint(QPainter &p, const QRectF &area, const QPalette &pal) const;
  QRectF plot_rect(const QRectF &area) const;
  Range auto_range() const;
  Range current_range() const;
  double to_screen_x(double v, const Range &r, const QRectF &pr) const;
  double to_screen_y(double v, const Range &r, const QRectF &pr) const;
  double from_screen_x(double s, const Range &r, const QRectF &pr) const;
  double from_screen_y(double s, const Range &r, const QRectF &pr) const;

  QString title_;
  QString xlabel_;
  QString ylabel_;
  QVector<PlotSeries> series_;
  QVector<QPair<double, QString>> markers_;
  bool log_x_ = false;
  bool log_y_ = false;
  std::optional<Range> zoom_;
  // Interaction.
  bool dragging_ = false;
  QPointF drag_from_;
  QPointF drag_to_;
  QPointF hover_;
  bool hovering_ = false;
};

} // namespace gui
} // namespace amcset
