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

#include "views/cell_view.h"

#include "app/guard.h"
#include "element_data.h"

#include <QMouseEvent>
#include <QPainter>
#include <QPolygonF>
#include <QRadialGradient>
#include <QToolTip>
#include <QWheelEvent>

#include <algorithm>
#include <cmath>

namespace amcset {
namespace gui {

using crystal::Mat3;
using crystal::Vec3;

namespace {

// More atoms than this are not drawn (a supercell read from a file).
constexpr std::size_t max_drawn_atoms = 6000;

QString symbol_of(unsigned z) {
  return QString::fromStdString(std::string(common::Element::by_z(z).symbol()));
}

QString number(double v, int decimals = 4) {
  QString s = QString::number(v, 'f', decimals);
  // 0.2500 -> 0.25, 1.0000 -> 1
  while (s.contains(QLatin1Char('.')) && (s.endsWith(QLatin1Char('0')) || s.endsWith(QLatin1Char('.')))) {
    const bool dot = s.endsWith(QLatin1Char('.'));
    s.chop(1);
    if (dot) break;
  }
  return s == QLatin1String("-0") ? QStringLiteral("0") : s;
}

} // namespace

QColor element_color(unsigned z) {
  // The usual colours of ball models for a few light elements; every other
  // element gets a hue of its own (golden-angle steps keep neighbours in the
  // periodic table apart).
  switch (z) {
  case 1: return QColor(235, 235, 235);
  case 6: return QColor(95, 95, 95);
  case 7: return QColor(48, 80, 248);
  case 8: return QColor(230, 40, 40);
  case 9: return QColor(120, 200, 80);
  case 15: return QColor(255, 128, 0);
  case 16: return QColor(235, 210, 50);
  case 17: return QColor(40, 200, 40);
  default: break;
  }
  const double hue = std::fmod(static_cast<double>(z) * 137.508, 360.0);
  return QColor::fromHsvF(static_cast<float>(hue / 360.0), 0.55f, 0.85f);
}

CellView::CellView(QWidget *parent) : QWidget(parent) {
  setMouseTracking(true);
  setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  // The standard oblique view (reset_view() without the signal).
  const Vec3 out = normalized(Vec3{0.80, -0.50, 0.33});
  const Vec3 up = normalized(Vec3{0, 0, 1} - out * dot(Vec3{0, 0, 1}, out));
  view_ = Mat3::from_rows(cross(up, out), up, out);
}

void CellView::set_crystal(const crystal::Crystal &crystal) {
  crystal_ = crystal;
  note_.clear();
  if (selected_ >= static_cast<int>(crystal.cell.sites.size())) selected_ = -1;
  hovered_ = -1;
  rebuild();
  update();
}

void CellView::clear(const QString &note) {
  crystal_.reset();
  note_ = note;
  selected_ = hovered_ = -1;
  rebuild();
  update();
}

void CellView::set_repeat(int n) {
  n = std::clamp(n, 1, 4);
  if (n == repeat_) return;
  repeat_ = n;
  rebuild();
  update();
}

void CellView::set_directions(const std::optional<Vec3> &beam, const std::optional<Vec3> &normal) {
  beam_ = beam;
  normal_ = normal;
  update();
}

void CellView::set_selected_site(int site) {
  if (site == selected_) return;
  selected_ = site;
  update();
}

Vec3 CellView::view_direction() const noexcept { return normalized(-view_[2]); }

Vec3 CellView::view_up() const noexcept { return normalized(view_[1]); }

void CellView::look_along(const Vec3 &direction) {
  const Vec3 f = normalized(direction);
  if (!(norm(f) > 0.0)) return;
  const Vec3 out = -f;
  // Keep the screen's "up" as far as the new direction allows.
  Vec3 up = view_[1] - out * dot(view_[1], out);
  if (norm(up) < 1e-6) up = view_[0] - out * dot(view_[0], out); // looking along the old "up"
  up = normalized(up);
  view_ = Mat3::from_rows(cross(up, out), up, out);
  update();
  emit view_changed();
}

void CellView::reset_view() {
  const Vec3 out = normalized(Vec3{0.80, -0.50, 0.33});
  const Vec3 up = normalized(Vec3{0, 0, 1} - out * dot(Vec3{0, 0, 1}, out));
  view_ = Mat3::from_rows(cross(up, out), up, out);
  zoom_ = 1.0;
  update();
  emit view_changed();
}

void CellView::orthonormalize() {
  const Vec3 out = normalized(view_[2]);
  const Vec3 up = normalized(view_[1] - out * dot(view_[1], out));
  view_ = Mat3::from_rows(cross(up, out), up, out);
}

void CellView::rebuild() {
  atoms_.clear();
  balls_.clear();
  truncated_ = false;
  drawn_repeat_ = 1;
  if (!crystal_) return;
  const auto &cell = crystal_->cell;
  const auto &lattice = cell.lattice;
  const std::size_t sites = std::min(cell.sites.size(), max_drawn_atoms);
  truncated_ = sites < cell.sites.size();
  int n = repeat_;
  // (The images on the far faces add up to (n + 1)^3 per site at most.)
  while (n > 1 && sites * static_cast<std::size_t>((n + 1) * (n + 1) * (n + 1)) > max_drawn_atoms) --n;
  drawn_repeat_ = n;
  for (std::size_t si = 0; si < sites; ++si) {
    const Vec3 f = cell.sites[si].position;
    // A site on a face, edge or corner (coordinate 0) is also drawn on the
    // opposite one, as in the usual picture of a cell.
    const int ni = f.x < 1e-6 ? n : n - 1;
    const int nj = f.y < 1e-6 ? n : n - 1;
    const int nk = f.z < 1e-6 ? n : n - 1;
    for (int i = 0; i <= ni; ++i) {
      for (int j = 0; j <= nj; ++j) {
        for (int k = 0; k <= nk; ++k) {
          atoms_.push_back({lattice.to_cartesian_A(Vec3{f.x + i, f.y + j, f.z + k}), static_cast<int>(si)});
        }
      }
    }
  }
  const double h = 0.5 * n;
  centre_ = lattice.to_cartesian_A(Vec3{h, h, h});
  extent_A_ = 0.0;
  for (int c = 0; c < 8; ++c) {
    const Vec3 corner = lattice.to_cartesian_A(Vec3{(c & 1) ? double(n) : 0.0, (c & 2) ? double(n) : 0.0, (c & 4) ? double(n) : 0.0});
    extent_A_ = std::max(extent_A_, norm(corner - centre_));
  }
  // Sphere size from the volume per lattice site (about a quarter of the
  // nearest-neighbour distance: the atoms do not hide each other).
  std::size_t lattice_sites = 0;
  for (const auto &s : cell.sites) lattice_sites += s.kind == crystal::SiteKind::lattice ? 1 : 0;
  if (lattice_sites == 0) lattice_sites = std::max<std::size_t>(cell.sites.size(), 1);
  ball_A_ = 0.21 * std::cbrt(std::abs(lattice.volume_A3()) / static_cast<double>(lattice_sites));
}

double CellView::scale() const {
  const double room = 0.5 * std::min(width(), height()) - 16.0;
  return std::max(room, 10.0) / std::max(extent_A_ + ball_A_, 1e-9) * zoom_;
}

QPointF CellView::project(const Vec3 &r, double *depth) const {
  const Vec3 d = r - centre_;
  if (depth != nullptr) *depth = dot(view_[2], d);
  const double s = scale();
  return {0.5 * width() + dot(view_[0], d) * s, 0.5 * height() - dot(view_[1], d) * s};
}

int CellView::site_at(const QPointF &pos) const {
  // balls_ is sorted far to near: the nearest atom under the pointer wins.
  for (auto it = balls_.rbegin(); it != balls_.rend(); ++it) {
    const double r = std::max(it->radius, 4.0);
    const QPointF d = pos - it->centre;
    if (d.x() * d.x() + d.y() * d.y() <= r * r) return it->site;
  }
  return -1;
}

void CellView::draw_atom(QPainter &p, const Ball &b) const {
  const auto &site = crystal_->cell.sites[static_cast<std::size_t>(b.site)];
  const double r = std::max(b.radius, 2.0);
  const QRectF box(b.centre.x() - r, b.centre.y() - r, 2.0 * r, 2.0 * r);
  // Far atoms are darker.
  const double nearness = std::clamp(0.5 + 0.5 * b.depth / std::max(extent_A_, 1e-9), 0.0, 1.0);
  const auto shaded = [nearness](const QColor &c) {
    return c.darker(static_cast<int>(100.0 + 45.0 * (1.0 - nearness)));
  };
  const QColor outline = palette().color(QPalette::Text);
  if (site.occupants.empty()) {
    // An empty site.
    p.setPen(QPen(outline, 1.2, Qt::DashLine));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(box);
  } else {
    p.setPen(Qt::NoPen);
    if (site.occupants.size() == 1 && site.occupancy() >= 1.0 - 1e-9) {
      p.setBrush(shaded(element_color(site.occupants.front().z)));
      p.drawEllipse(box);
    } else {
      // Shared or partly vacant: a pie chart of the occupants, clockwise
      // from the top; the vacant part stays empty (see-through when most
      // of the site is vacant).
      if (site.occupancy() >= 0.5) {
        p.setBrush(palette().base());
        p.drawEllipse(box);
      }
      double start = 90.0;
      for (const auto &o : site.occupants) {
        const double span = -360.0 * o.fraction;
        p.setBrush(shaded(element_color(o.z)));
        p.drawPie(box, static_cast<int>(std::lround(start * 16.0)), static_cast<int>(std::lround(span * 16.0)));
        start += span;
      }
    }
    // A highlight and a dark rim make it a sphere; a site that is mostly
    // vacant (an interstitial site with a few per cent of solute) stays a
    // flat ring with its sliver, so that it does not look like an atom.
    const bool mostly_vacant = site.occupancy() < 0.5;
    if (!mostly_vacant) {
      QRadialGradient g(b.centre + QPointF(-0.35 * r, -0.35 * r), 1.3 * r);
      g.setColorAt(0.0, QColor(255, 255, 255, 150));
      g.setColorAt(0.35, QColor(255, 255, 255, 25));
      g.setColorAt(1.0, QColor(0, 0, 0, 80));
      p.setBrush(g);
      p.drawEllipse(box);
    }
    QColor edge = outline;
    edge.setAlphaF(mostly_vacant ? 0.35f : 0.55f);
    p.setPen(QPen(edge, 1.0, mostly_vacant ? Qt::DashLine : Qt::SolidLine));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(box);
  }
  if (b.site == selected_ || b.site == hovered_) {
    p.setPen(QPen(palette().color(QPalette::Highlight), b.site == selected_ ? 2.6 : 1.4));
    p.setBrush(Qt::NoBrush);
    p.drawEllipse(box.adjusted(-2.5, -2.5, 2.5, 2.5));
  }
  // The element's symbol where there is room.
  if (r >= 9.0 && !site.occupants.empty() && site.occupancy() >= 0.5) {
    const auto major = std::max_element(site.occupants.begin(), site.occupants.end(),
                                       [](const auto &a, const auto &c) { return a.fraction < c.fraction; });
    QFont f = font();
    f.setPixelSize(std::clamp(static_cast<int>(r * 0.8), 8, 14));
    p.setFont(f);
    p.setPen(shaded(element_color(major->z)).lightnessF() > 0.6f ? QColor(20, 20, 20) : QColor(255, 255, 255));
    p.drawText(box, Qt::AlignCenter, symbol_of(major->z) + (site.occupants.size() > 1 ? QStringLiteral("+") : QString()));
  }
}

void CellView::draw_compass(QPainter &p) const {
  if (!beam_ && !normal_) return;
  const double radius = 24.0;
  const QColor text = palette().color(QPalette::Text);
  const QColor beam_color(230, 126, 34);
  QColor normal_color = text;
  normal_color.setAlphaF(0.6f);
  QFont f = font();
  f.setPixelSize(11);
  p.setFont(f);
  const QFontMetrics fm(f);
  const QPointF c(width() - radius - 12.0, height() - radius - 2.0 * fm.height() - 8.0);
  QColor ring = text;
  ring.setAlphaF(0.25f);
  p.setPen(QPen(ring, 1.0));
  p.setBrush(Qt::NoBrush);
  p.drawEllipse(c, radius, radius);
  const auto arrow = [&](const Vec3 &direction, const QColor &color) {
    const Vec3 d = normalized(direction);
    const double x = dot(view_[0], d);
    const double y = dot(view_[1], d);
    const double z = dot(view_[2], d);
    const double length = std::hypot(x, y);
    p.setPen(QPen(color, 2.0));
    if (length < 0.12) {
      // Along the view: into the screen (a cross) or out of it (a dot).
      p.setBrush(Qt::NoBrush);
      p.drawEllipse(c, 6.5, 6.5);
      if (z < 0.0) {
        p.drawLine(c + QPointF(-4.5, -4.5), c + QPointF(4.5, 4.5));
        p.drawLine(c + QPointF(-4.5, 4.5), c + QPointF(4.5, -4.5));
      } else {
        p.setBrush(color);
        p.drawEllipse(c, 2.0, 2.0);
      }
      return;
    }
    const QPointF u(x / length, -y / length);
    const QPointF v(-u.y(), u.x());
    const QPointF tip = c + u * (length * radius);
    p.drawLine(c, tip);
    p.setBrush(color);
    QPolygonF head;
    head << tip << tip - u * 7.0 + v * 3.5 << tip - u * 7.0 - v * 3.5;
    p.drawPolygon(head);
  };
  if (normal_) arrow(*normal_, normal_color);
  if (beam_) arrow(*beam_, beam_color);
  // The legend, against the right edge.
  const QRectF line1(width() - 168.0, c.y() + radius + 3.0, 160.0, fm.height());
  const QRectF line2 = line1.translated(0.0, fm.height());
  if (beam_) {
    p.setPen(beam_color);
    p.drawText(line1, Qt::AlignRight | Qt::AlignVCenter, tr("beam"));
  }
  if (normal_) {
    p.setPen(normal_color);
    p.drawText(line2, Qt::AlignRight | Qt::AlignVCenter, tr("depth (surface normal)"));
  }
}

void CellView::paint(QPainter &p) {
  p.setRenderHint(QPainter::Antialiasing, true);
  const QPalette &pal = palette();
  p.fillRect(rect(), pal.base());
  balls_.clear();
  const QColor text = pal.color(QPalette::Text);
  if (!crystal_) {
    p.setPen(pal.color(QPalette::Disabled, QPalette::Text));
    p.drawText(rect().adjusted(12, 12, -12, -12), Qt::AlignCenter | Qt::TextWordWrap,
               note_.isEmpty() ? tr("No crystal") : note_);
    return;
  }
  const auto &cell = crystal_->cell;
  const auto &lattice = cell.lattice;
  const int n = drawn_repeat_;
  const double s = scale();

  // The cells: the outline of the drawn block, the cells inside it fainter.
  QColor faint = text;
  faint.setAlphaF(0.16f);
  QColor line = text;
  line.setAlphaF(0.5f);
  const auto corner = [&](int i, int j, int k) {
    return project(lattice.to_cartesian_A(Vec3{double(i), double(j), double(k)}));
  };
  for (int axis = 0; axis < 3; ++axis) {
    for (int u = 0; u <= n; ++u) {
      for (int v = 0; v <= n; ++v) {
        const bool outer = (u == 0 || u == n) && (v == 0 || v == n);
        p.setPen(QPen(outer ? line : faint, outer ? 1.3 : 1.0));
        int a[3], b[3];
        a[axis] = 0;
        b[axis] = n;
        a[(axis + 1) % 3] = b[(axis + 1) % 3] = u;
        a[(axis + 2) % 3] = b[(axis + 2) % 3] = v;
        p.drawLine(corner(a[0], a[1], a[2]), corner(b[0], b[1], b[2]));
      }
    }
  }
  // The lattice vectors, from the origin of the cell.
  static const QColor axis_color[3] = {QColor(214, 69, 65), QColor(46, 160, 67), QColor(56, 120, 220)};
  static const char *const axis_name[3] = {"a", "b", "c"};
  QFont label_font = font();
  label_font.setBold(true);
  p.setFont(label_font);
  const QPointF origin = corner(0, 0, 0);
  for (int i = 0; i < 3; ++i) {
    const QPointF tip = project(lattice.vectors_A()[i]);
    const QPointF d = tip - origin;
    const double length = std::hypot(d.x(), d.y());
    p.setPen(QPen(axis_color[i], 2.0));
    p.drawLine(origin, tip);
    if (length > 10.0) {
      const QPointF u = d / length;
      const QPointF v(-u.y(), u.x());
      p.setBrush(axis_color[i]);
      QPolygonF head;
      head << tip << tip - u * 8.0 + v * 4.0 << tip - u * 8.0 - v * 4.0;
      p.drawPolygon(head);
      p.drawText(QRectF(tip.x() + u.x() * 12.0 - 8.0, tip.y() + u.y() * 12.0 - 8.0, 16.0, 16.0), Qt::AlignCenter,
                 QString::fromLatin1(axis_name[i]));
    }
  }

  // The atoms, far to near.
  balls_.reserve(atoms_.size());
  for (const auto &a : atoms_) {
    Ball b;
    b.centre = project(a.position, &b.depth);
    b.site = a.site;
    const bool interstitial = cell.sites[static_cast<std::size_t>(a.site)].kind == crystal::SiteKind::interstitial;
    b.radius = s * ball_A_ * (interstitial ? 0.6 : 1.0);
    balls_.push_back(b);
  }
  std::stable_sort(balls_.begin(), balls_.end(), [](const Ball &a, const Ball &b) { return a.depth < b.depth; });
  for (const auto &b : balls_) draw_atom(p, b);

  // Captions: the direction of view, and what is left out.
  p.setFont(font());
  p.setPen(text);
  const auto nearest = crystal::nearest_direction(lattice, view_direction(), 6);
  const double off = common::in_units(nearest.off, common::degree);
  const QString indices = QStringLiteral("[%1 %2 %3]").arg(nearest.uvw[0]).arg(nearest.uvw[1]).arg(nearest.uvw[2]);
  p.drawText(QRectF(8.0, 6.0, width() - 16.0, 20.0), Qt::AlignLeft | Qt::AlignVCenter,
             off < 0.05 ? tr("View along %1").arg(indices)
                        : tr("View %1 deg from %2").arg(QString::number(off, 'f', off < 10.0 ? 1 : 0), indices));
  QString below;
  if (truncated_) {
    below = tr("%1 of %2 sites shown").arg(max_drawn_atoms).arg(cell.sites.size());
  } else if (n != repeat_) {
    below = tr("One cell shown (too many atoms for more)");
  }
  if (!below.isEmpty()) {
    p.drawText(QRectF(8.0, height() - 24.0, width() - 16.0, 20.0), Qt::AlignLeft | Qt::AlignVCenter, below);
  }
  draw_compass(p);
}

void CellView::paintEvent(QPaintEvent *) {
  try {
    QPainter p(this);
    paint(p);
  } catch (...) {
    report_exception_quietly(QStringLiteral("Unit cell view (paintEvent)"), std::current_exception());
  }
}

void CellView::mousePressEvent(QMouseEvent *event) {
  if (event->button() != Qt::LeftButton) return;
  dragging_ = true;
  moved_ = false;
  press_ = last_ = event->pos();
}

void CellView::mouseMoveEvent(QMouseEvent *event) {
  try {
    if (dragging_) {
      if (!moved_ && (event->pos() - press_).manhattanLength() < 4) return;
      moved_ = true;
      const QPoint d = event->pos() - last_;
      last_ = event->pos();
      const double k = 0.01; // radians per pixel
      if (event->modifiers() & Qt::ShiftModifier) {
        // Roll about the direction of view.
        view_ = crystal::rotation_about(Vec3{0, 0, 1}, -d.x() * k) * view_;
      } else {
        // The front of the cell follows the pointer.
        view_ = crystal::rotation_about(Vec3{0, 1, 0}, d.x() * k) * crystal::rotation_about(Vec3{1, 0, 0}, d.y() * k) * view_;
      }
      orthonormalize();
      update();
      emit view_changed();
      return;
    }
    const int site = crystal_ ? site_at(event->position()) : -1;
    if (site != hovered_) {
      hovered_ = site;
      update();
    }
    if (site < 0) {
      QToolTip::hideText();
      return;
    }
    const auto &s = crystal_->cell.sites[static_cast<std::size_t>(site)];
    QString occupants;
    for (const auto &o : s.occupants) {
      occupants += (occupants.isEmpty() ? QString() : QStringLiteral(", ")) +
                   QStringLiteral("%1 %2 %").arg(symbol_of(o.z), number(100.0 * o.fraction, 2));
    }
    if (occupants.isEmpty()) occupants = tr("empty");
    QToolTip::showText(event->globalPosition().toPoint(),
                       tr("%1: %2 site at (%3, %4, %5)\n%6\nClick to edit")
                           .arg(QString::fromStdString(s.label),
                                s.kind == crystal::SiteKind::interstitial ? tr("interstitial") : tr("lattice"),
                                number(s.position.x), number(s.position.y), number(s.position.z), occupants),
                       this);
  } catch (...) {
    report_exception_quietly(QStringLiteral("Unit cell view (mouseMoveEvent)"), std::current_exception());
  }
}

void CellView::mouseReleaseEvent(QMouseEvent *event) {
  try {
    if (event->button() != Qt::LeftButton || !dragging_) return;
    dragging_ = false;
    if (moved_ || !crystal_) return;
    const int site = site_at(event->position());
    if (site >= 0) emit site_clicked(site, event->globalPosition().toPoint());
  } catch (...) {
    report_exception_quietly(QStringLiteral("Unit cell view (mouseReleaseEvent)"), std::current_exception());
  }
}

void CellView::wheelEvent(QWheelEvent *event) {
  zoom_ = std::clamp(zoom_ * std::pow(1.0015, event->angleDelta().y()), 0.3, 12.0);
  update();
  event->accept();
}

void CellView::leaveEvent(QEvent *) {
  if (hovered_ >= 0) {
    hovered_ = -1;
    update();
  }
}

} // namespace gui
} // namespace amcset
