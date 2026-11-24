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

#include "app/guard.h"
#include "views/plot_widget.h"

#include "views/axis_ticks.h"

#include <QApplication>
#include <QClipboard>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QMenu>
#include <QMouseEvent>
#include <QPainter>
#include <QPainterPath>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>
#include <limits>

namespace amcset {
namespace gui {

PlotWidget::PlotWidget(QWidget *parent) : QWidget(parent) {
  setMouseTracking(true);
  setMinimumSize(240, 180);
  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
}

void PlotWidget::set_title(const QString &t) {
  title_ = t;
  update();
}

void PlotWidget::set_axis_labels(const QString &x, const QString &y) {
  xlabel_ = x;
  ylabel_ = y;
  update();
}

void PlotWidget::set_series(QVector<PlotSeries> s, bool keep_zoom) {
  series_ = std::move(s);
  if (!keep_zoom) zoom_.reset();
  update();
}

void PlotWidget::set_log_x(bool on) {
  if (on == log_x_) return;
  log_x_ = on;
  zoom_.reset();
  emit log_scale_changed(log_x_, log_y_);
  update();
}

void PlotWidget::set_log_y(bool on) {
  if (on == log_y_) return;
  log_y_ = on;
  zoom_.reset();
  emit log_scale_changed(log_x_, log_y_);
  update();
}

void PlotWidget::set_markers(QVector<QPair<double, QString>> m) {
  markers_ = std::move(m);
  update();
}

void PlotWidget::reset_zoom() {
  zoom_.reset();
  update();
}

PlotWidget::Range PlotWidget::auto_range() const {
  double x0 = std::numeric_limits<double>::max(), x1 = -x0;
  double y0 = x0, y1 = -x0;
  double ypos_min = x0;
  for (const auto &s : series_) {
    const int n = static_cast<int>(std::min(s.x.size(), s.y.size()));
    for (int i = 0; i < n; ++i) {
      const double x = s.x[i], y = s.y[i];
      if (!std::isfinite(x) || !std::isfinite(y)) continue;
      if (log_x_ && x <= 0) continue;
      x0 = std::min(x0, x);
      x1 = std::max(x1, x);
      y0 = std::min(y0, y);
      y1 = std::max(y1, y);
      if (y > 0) ypos_min = std::min(ypos_min, y);
    }
    // Histogram steps start one bin before the first label.
    if (s.style == PlotSeries::Style::steps && n >= 2 && !log_x_) x0 = std::min(x0, s.x[0] - (s.x[1] - s.x[0]));
  }
  Range r;
  if (x0 > x1) return r;
  if (x1 == x0) {
    x0 -= 0.5;
    x1 += 0.5;
  }
  r.x0 = x0;
  r.x1 = x1;
  if (log_y_) {
    if (ypos_min > y1 || y1 <= 0) {
      r.y0 = 0.1;
      r.y1 = 1;
    } else {
      r.y0 = std::max(ypos_min, y1 * 1e-8) / 1.5;
      r.y1 = y1 * 1.5;
    }
  } else {
    const double lo = std::min(0.0, y0);
    const double hi = std::max(0.0, y1);
    const double pad = (hi - lo) * 0.06;
    r.y0 = lo < 0 ? lo - pad : lo;
    r.y1 = hi > 0 ? hi + pad : hi + 1.0;
    if (r.y1 == r.y0) r.y1 = r.y0 + 1;
  }
  return r;
}

PlotWidget::Range PlotWidget::current_range() const { return zoom_ ? *zoom_ : auto_range(); }

QRectF PlotWidget::plot_rect(const QRectF &area) const {
  const QFontMetricsF fm(font());
  const double left = fm.horizontalAdvance(QStringLiteral("-8.88E-88")) + fm.height() + 10;
  const double top = fm.height() * (title_.isEmpty() ? 0.8 : 2.0);
  const double bottom = 2.4 * fm.height() + 6;
  return area.adjusted(left, top, -14, -bottom);
}

double PlotWidget::to_screen_x(double v, const Range &r, const QRectF &pr) const {
  const double t = log_x_ ? (std::log10(v) - std::log10(r.x0)) / (std::log10(r.x1) - std::log10(r.x0))
                          : (v - r.x0) / (r.x1 - r.x0);
  return pr.left() + t * pr.width();
}

double PlotWidget::to_screen_y(double v, const Range &r, const QRectF &pr) const {
  double t;
  if (log_y_) {
    if (v <= 0) return pr.bottom() + 10;
    t = (std::log10(v) - std::log10(r.y0)) / (std::log10(r.y1) - std::log10(r.y0));
  } else {
    t = (v - r.y0) / (r.y1 - r.y0);
  }
  return pr.bottom() - t * pr.height();
}

double PlotWidget::from_screen_x(double s, const Range &r, const QRectF &pr) const {
  const double t = (s - pr.left()) / pr.width();
  return log_x_ ? std::pow(10.0, std::log10(r.x0) + t * (std::log10(r.x1) - std::log10(r.x0)))
                : r.x0 + t * (r.x1 - r.x0);
}

double PlotWidget::from_screen_y(double s, const Range &r, const QRectF &pr) const {
  const double t = (pr.bottom() - s) / pr.height();
  return log_y_ ? std::pow(10.0, std::log10(r.y0) + t * (std::log10(r.y1) - std::log10(r.y0)))
                : r.y0 + t * (r.y1 - r.y0);
}

void PlotWidget::paint(QPainter &p, const QRectF &area, const QPalette &pal) const {
  const QColor fg = pal.color(QPalette::WindowText);
  QColor grid = fg;
  grid.setAlphaF(0.12f);
  QColor faint = fg;
  faint.setAlphaF(0.55f);
  p.fillRect(area, pal.color(QPalette::Base));
  const QFontMetricsF fm(p.font());
  const QRectF pr = plot_rect(area);
  const Range r = current_range();
  if (!title_.isEmpty()) {
    p.save();
    QFont bold = p.font();
    bold.setBold(true);
    p.setFont(bold);
    p.setPen(fg);
    p.drawText(QRectF(area.left(), area.top() + 2, area.width(), fm.height() * 1.6), Qt::AlignCenter, title_);
    p.restore();
  }
  // Grid and ticks.
  const auto xt = log_x_ ? log_ticks(r.x0, r.x1, 8) : nice_ticks(r.x0, r.x1, std::max(3, static_cast<int>(pr.width() / 80)));
  const auto yt = log_y_ ? log_ticks(r.y0, r.y1, 8) : nice_ticks(r.y0, r.y1, std::max(3, static_cast<int>(pr.height() / 45)));
  p.setPen(grid);
  for (double t : xt) {
    const double sx = to_screen_x(t, r, pr);
    p.drawLine(QPointF(sx, pr.top()), QPointF(sx, pr.bottom()));
  }
  for (double t : yt) {
    const double sy = to_screen_y(t, r, pr);
    p.drawLine(QPointF(pr.left(), sy), QPointF(pr.right(), sy));
  }
  p.setPen(fg);
  for (double t : xt) {
    const double sx = to_screen_x(t, r, pr);
    const QString s = format_tick(t);
    p.drawText(QPointF(sx - fm.horizontalAdvance(s) / 2, pr.bottom() + fm.ascent() + 4), s);
  }
  for (double t : yt) {
    const double sy = to_screen_y(t, r, pr);
    const QString s = format_tick(t);
    p.drawText(QPointF(pr.left() - fm.horizontalAdvance(s) - 5, sy + fm.ascent() / 2 - 1), s);
  }
  // Axis labels.
  p.drawText(QRectF(pr.left(), pr.bottom() + fm.height() + 6, pr.width(), fm.height() + 2), Qt::AlignCenter, xlabel_);
  p.save();
  p.translate(area.left() + 2, pr.center().y());
  p.rotate(-90);
  p.drawText(QRectF(-pr.height() / 2, 0, pr.height(), fm.height() + 2), Qt::AlignCenter, ylabel_);
  p.restore();

  // Data, clipped to the plot area.
  p.save();
  p.setClipRect(pr.adjusted(-1, -1, 1, 1));
  p.setRenderHint(QPainter::Antialiasing, true);
  const double base_y = log_y_ ? pr.bottom() : to_screen_y(std::max(r.y0, std::min(0.0, r.y1)), r, pr);
  for (const auto &s : series_) {
    const int n = static_cast<int>(std::min(s.x.size(), s.y.size()));
    if (n == 0) continue;
    QPolygonF poly;
    poly.reserve(s.style == PlotSeries::Style::steps ? 2 * n + 2 : n);
    if (s.style == PlotSeries::Style::steps) {
      for (int i = 0; i < n; ++i) {
        const double xa = i > 0 ? s.x[i - 1] : (n > 1 ? s.x[0] - (s.x[1] - s.x[0]) : s.x[0] - 1);
        if (log_x_ && xa <= 0) continue;
        const double sy = to_screen_y(s.y[i], r, pr);
        poly << QPointF(to_screen_x(xa, r, pr), sy) << QPointF(to_screen_x(s.x[i], r, pr), sy);
      }
    } else {
      for (int i = 0; i < n; ++i) {
        if (!std::isfinite(s.y[i]) || (log_x_ && s.x[i] <= 0)) continue;
        poly << QPointF(to_screen_x(s.x[i], r, pr), to_screen_y(s.y[i], r, pr));
      }
    }
    if (poly.isEmpty()) continue;
    if (s.fill && s.style != PlotSeries::Style::points) {
      QPolygonF area_poly = poly;
      area_poly << QPointF(poly.last().x(), base_y) << QPointF(poly.first().x(), base_y);
      QColor c = s.color;
      c.setAlphaF(0.18f);
      p.setPen(Qt::NoPen);
      p.setBrush(c);
      p.drawPolygon(area_poly);
      p.setBrush(Qt::NoBrush);
    }
    QPen pen(s.color, 1.6);
    p.setPen(pen);
    if (s.style == PlotSeries::Style::points) {
      pen.setWidthF(3.0);
      pen.setCapStyle(Qt::RoundCap);
      p.setPen(pen);
      p.drawPoints(poly);
    } else {
      p.drawPolyline(poly);
    }
  }
  // Markers.
  for (const auto &m : markers_) {
    if (log_x_ && m.first <= 0) continue;
    const double sx = to_screen_x(m.first, r, pr);
    QPen pen(faint, 1.0, Qt::DashLine);
    p.setPen(pen);
    p.drawLine(QPointF(sx, pr.top()), QPointF(sx, pr.bottom()));
    p.setPen(fg);
    p.drawText(QPointF(sx + 4, pr.top() + fm.ascent() + 2), m.second);
  }
  p.restore();
  p.setPen(faint);
  p.drawRect(pr);

  // Legend.
  if (series_.size() > 1 || (series_.size() == 1 && !series_.front().name.isEmpty())) {
    double w = 0;
    for (const auto &s : series_) w = std::max(w, fm.horizontalAdvance(s.name));
    const double h = fm.height() * series_.size() + 6;
    QRectF box(pr.right() - w - 34, pr.top() + 6, w + 28, h);
    QColor bg = pal.color(QPalette::Base);
    bg.setAlphaF(0.85f);
    p.fillRect(box, bg);
    double y = box.top() + 3;
    for (const auto &s : series_) {
      p.fillRect(QRectF(box.left() + 6, y + fm.height() / 2 - 2, 14, 4), s.color);
      p.setPen(fg);
      p.drawText(QPointF(box.left() + 24, y + fm.ascent()), s.name);
      y += fm.height();
    }
  }
}

void PlotWidget::paintEvent(QPaintEvent *) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    QPainter p(this);
    paint(p, QRectF(rect()), palette());
    const QRectF pr = plot_rect(QRectF(rect()));
    if (dragging_) {
      QColor c = palette().color(QPalette::Highlight);
      c.setAlphaF(0.2f);
      p.fillRect(QRectF(drag_from_, drag_to_).normalized().intersected(pr), c);
    } else if (hovering_ && pr.contains(hover_)) {
      QColor c = palette().color(QPalette::WindowText);
      c.setAlphaF(0.3f);
      p.setPen(QPen(c, 1, Qt::DotLine));
      p.drawLine(QPointF(hover_.x(), pr.top()), QPointF(hover_.x(), pr.bottom()));
    }
    if (series_.isEmpty()) {
      p.setPen(palette().color(QPalette::PlaceholderText));
      p.drawText(pr, Qt::AlignCenter, tr("No data yet"));
    }
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (paintEvent)"), std::current_exception());
  }
}

void PlotWidget::mousePressEvent(QMouseEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (e->button() != Qt::LeftButton) return;
    dragging_ = true;
    drag_from_ = drag_to_ = e->position();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (mousePressEvent)"), std::current_exception());
  }
}

void PlotWidget::mouseMoveEvent(QMouseEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    hover_ = e->position();
    hovering_ = true;
    if (dragging_) drag_to_ = e->position();
    const QRectF pr = plot_rect(QRectF(rect()));
    if (pr.contains(hover_) && !series_.isEmpty()) {
      const Range r = current_range();
      const double x = from_screen_x(hover_.x(), r, pr);
      // Value of each series at this x (nearest point).
      QStringList parts;
      parts << QStringLiteral("x = %1").arg(QString::number(x, 'g', 5));
      for (const auto &s : series_) {
        if (s.x.isEmpty()) continue;
        int best = 0;
        double bd = std::numeric_limits<double>::max();
        for (int i = 0; i < s.x.size(); ++i) {
          const double d = std::abs(s.x[i] - x);
          if (d < bd) {
            bd = d;
            best = i;
          }
        }
        if (s.style == PlotSeries::Style::steps) {
          // The bin containing x.
          for (int i = 0; i < s.x.size(); ++i) {
            if (s.x[i] >= x) {
              best = i;
              break;
            }
          }
        }
        parts << QStringLiteral("%1 = %2").arg(s.name.isEmpty() ? tr("y") : s.name,
                                                QString::number(s.y.value(best), 'g', 5));
      }
      emit hover_text(parts.join(QStringLiteral(",  ")));
    }
    update();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (mouseMoveEvent)"), std::current_exception());
  }
}

void PlotWidget::mouseReleaseEvent(QMouseEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (!dragging_ || e->button() != Qt::LeftButton) return;
    dragging_ = false;
    const QRectF pr = plot_rect(QRectF(rect()));
    const QRectF sel = QRectF(drag_from_, drag_to_).normalized().intersected(pr);
    if (sel.width() > 6 && sel.height() > 6) {
      const Range r = current_range();
      Range z;
      z.x0 = from_screen_x(sel.left(), r, pr);
      z.x1 = from_screen_x(sel.right(), r, pr);
      z.y0 = from_screen_y(sel.bottom(), r, pr);
      z.y1 = from_screen_y(sel.top(), r, pr);
      zoom_ = z;
    }
    update();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (mouseReleaseEvent)"), std::current_exception());
  }
}

void PlotWidget::mouseDoubleClickEvent(QMouseEvent *) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    reset_zoom();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (mouseDoubleClickEvent)"), std::current_exception());
  }
}

void PlotWidget::wheelEvent(QWheelEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    const QRectF pr = plot_rect(QRectF(rect()));
    if (!pr.contains(e->position())) return;
    Range r = current_range();
    const double f = std::pow(1.0015, -e->angleDelta().y());
    const double c = from_screen_x(e->position().x(), r, pr);
    if (log_x_) {
      const double lc = std::log10(c);
      r.x0 = std::pow(10.0, lc - (lc - std::log10(r.x0)) * f);
      r.x1 = std::pow(10.0, lc + (std::log10(r.x1) - lc) * f);
    } else {
      r.x0 = c - (c - r.x0) * f;
      r.x1 = c + (r.x1 - c) * f;
    }
    zoom_ = r;
    update();
    e->accept();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (wheelEvent)"), std::current_exception());
  }
}

void PlotWidget::leaveEvent(QEvent *) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    hovering_ = false;
    emit hover_text({});
    update();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (leaveEvent)"), std::current_exception());
  }
}

void PlotWidget::contextMenuEvent(QContextMenuEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    QMenu menu(this);
    auto *lx = menu.addAction(tr("Logarithmic x axis"));
    lx->setCheckable(true);
    lx->setChecked(log_x_);
    connect(lx, &QAction::toggled, this, guarded(this, tr("Plot (toggled)"), [this](auto &&...a) { set_log_x(std::forward<decltype(a)>(a)...); }));
    auto *ly = menu.addAction(tr("Logarithmic y axis"));
    ly->setCheckable(true);
    ly->setChecked(log_y_);
    connect(ly, &QAction::toggled, this, guarded(this, tr("Plot (toggled)"), [this](auto &&...a) { set_log_y(std::forward<decltype(a)>(a)...); }));
    menu.addAction(tr("Reset zoom"), this, &PlotWidget::reset_zoom);
    menu.addSeparator();
    menu.addAction(tr("Copy data"), this, guarded(this, tr("Plot menu"), [this] { QApplication::clipboard()->setText(data_text()); }));
    menu.addAction(tr("Copy image"), this, guarded(this, tr("Plot menu"), [this] { QApplication::clipboard()->setImage(render_image(size())); }));
    menu.addAction(tr("Save image..."), this, guarded(this, tr("Plot menu"), [this] {
      const QString f = QFileDialog::getSaveFileName(this, tr("Save plot"), title_ + QStringLiteral(".png"),
                                                     tr("Images (*.png *.jpg *.bmp)"));
      if (!f.isEmpty()) save_image(f);
    }));
    menu.exec(e->globalPos());
  } catch (...) {
    report_exception_quietly(QStringLiteral("Plot (contextMenuEvent)"), std::current_exception());
  }
}

QString PlotWidget::data_text() const {
  QString out;
  for (const auto &s : series_) {
    out += QStringLiteral("# %1\n%2\t%3\n").arg(s.name, xlabel_, ylabel_);
    const int n = static_cast<int>(std::min(s.x.size(), s.y.size()));
    for (int i = 0; i < n; ++i) {
      out += QString::number(s.x[i], 'g', 10) + QLatin1Char('\t') + QString::number(s.y[i], 'g', 10) +
             QLatin1Char('\n');
    }
    out += QLatin1Char('\n');
  }
  return out;
}

QImage PlotWidget::render_image(const QSize &size, qreal dpr) const {
  QImage img(size * dpr, QImage::Format_ARGB32_Premultiplied);
  img.setDevicePixelRatio(dpr);
  QPainter p(&img);
  p.setFont(font());
  paint(p, QRectF(QPointF(0, 0), QSizeF(size)), palette());
  return img;
}

bool PlotWidget::save_image(const QString &path) const {
  return render_image(QSize(std::max(width(), 900), std::max(height(), 600))).save(path);
}

} // namespace gui
} // namespace amcset
