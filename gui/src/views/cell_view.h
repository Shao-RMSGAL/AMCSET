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

/*! \file cell_view.h
 *  \brief A three-dimensional view of a crystal's unit cell.
 *
 * Atoms are drawn as shaded spheres (sites shared by several elements as
 * pie charts, interstitial sites smaller, empty sites as rings), the cell as
 * its twelve edges, with the lattice vectors, and the beam and the surface
 * normal in a corner. The projection is parallel, so rows of atoms line up
 * exactly when the view is along a crystal direction (the channels show).
 * Drag to turn, Shift+drag to roll, wheel to zoom, click an atom to edit its
 * site. The view direction is a crystal direction: the panel uses it to
 * point the beam the way the user looks. Painted with QPainter (a unit cell
 * has few atoms; no GPU needed).
 */

#pragma once

#include "crystal/crystal.h"

#include <QColor>
#include <QPoint>
#include <QPointF>
#include <QWidget>

#include <optional>
#include <vector>

namespace amcset {
namespace gui {

//! The colour an element is drawn with.
QColor element_color(unsigned z);

class CellView : public QWidget {
  Q_OBJECT

public:
  explicit CellView(QWidget *parent = nullptr);

  //! Show \p crystal (its cell); keeps the view direction.
  void set_crystal(const crystal::Crystal &crystal);
  //! Show nothing (an amorphous layer); \p note is written in the middle.
  void clear(const QString &note = {});
  bool has_crystal() const noexcept { return crystal_.has_value(); }
  //! Unit cells drawn along each lattice vector (1 to 4).
  void set_repeat(int n);
  int repeat() const noexcept { return repeat_; }
  /*!
   * \brief The beam's direction and the inward surface normal (the depth
   * axis), both in the crystal's Cartesian frame; drawn in a corner.
   */
  void set_directions(const std::optional<crystal::Vec3> &beam, const std::optional<crystal::Vec3> &normal);
  //! The site drawn highlighted (-1: none).
  void set_selected_site(int site);
  int selected_site() const noexcept { return selected_; }

  //! The direction the user looks along (into the screen), crystal Cartesian frame, unit length.
  crystal::Vec3 view_direction() const noexcept;
  //! The direction drawn towards the top of the view (crystal Cartesian frame, unit length).
  crystal::Vec3 view_up() const noexcept;
  //! Turn the view to look along \p direction (crystal Cartesian frame, any length).
  void look_along(const crystal::Vec3 &direction);
  //! The standard oblique view.
  void reset_view();

  QSize sizeHint() const override { return {340, 320}; }
  QSize minimumSizeHint() const override { return {220, 200}; }

signals:
  //! An atom (site index in UnitCell::sites) was clicked.
  void site_clicked(int site, const QPoint &global_position);
  //! The user turned the view.
  void view_changed();

protected:
  void paintEvent(QPaintEvent *event) override;
  void mousePressEvent(QMouseEvent *event) override;
  void mouseMoveEvent(QMouseEvent *event) override;
  void mouseReleaseEvent(QMouseEvent *event) override;
  void wheelEvent(QWheelEvent *event) override;
  void leaveEvent(QEvent *event) override;

private:
  //! One drawn atom (a site and one of its periodic images).
  struct Atom {
    crystal::Vec3 position; // Cartesian, Angstrom
    int site = -1;
  };
  //! An atom on the screen, as last painted (for picking).
  struct Ball {
    QPointF centre;
    double radius = 0.0;
    double depth = 0.0; // towards the viewer
    int site = -1;
  };
  void rebuild();
  void paint(QPainter &p);
  void orthonormalize();
  double scale() const;
  QPointF project(const crystal::Vec3 &r, double *depth = nullptr) const;
  int site_at(const QPointF &pos) const;
  void draw_atom(QPainter &p, const Ball &b) const;
  void draw_compass(QPainter &p) const;

  std::optional<crystal::Crystal> crystal_;
  QString note_;
  std::optional<crystal::Vec3> beam_;
  std::optional<crystal::Vec3> normal_;
  // Rows: screen right, screen up, towards the viewer (crystal Cartesian frame).
  crystal::Mat3 view_;
  double zoom_ = 1.0;
  int repeat_ = 1;
  int drawn_repeat_ = 1;   // repeat_, or 1 when that would be too many atoms
  int selected_ = -1;
  int hovered_ = -1;
  // Derived from the crystal.
  std::vector<Atom> atoms_;
  bool truncated_ = false;  // too many atoms: one cell only, or nothing
  crystal::Vec3 centre_;
  double extent_A_ = 1.0;   // radius of the drawn block
  double ball_A_ = 0.5;     // sphere radius of a lattice atom
  std::vector<Ball> balls_; // as last painted, far to near
  // Dragging.
  bool dragging_ = false;
  bool moved_ = false;
  QPoint press_;
  QPoint last_;
};

} // namespace gui
} // namespace amcset
