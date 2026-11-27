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

/*! \file track_view.h
 *  \brief GPU view of the particle paths: 2D projections and a 3D view.
 *
 * One QRhiWidget (Vulkan, Direct3D, Metal or OpenGL, chosen by Qt) draws the
 * paths of the TrackStore as line segments from a single vertex buffer that
 * only grows: new paths are uploaded as they arrive, so the cost of a frame
 * does not depend on how many paths were drawn before it, and millions of
 * segments stay interactive. The widget shows one view or four (depth-Y,
 * depth-Z, beam view and 3D) side by side, each with its own camera; axes and
 * labels are painted over the GPU image.
 */

#pragma once

#include <QColor>
#include <QMatrix4x4>
#include <QPointF>
#include <QRhiWidget>
#include <QVector3D>

#include <array>
#include <memory>
#include <vector>

class QRhi;
class QRhiBuffer;
class QRhiGraphicsPipeline;
class QRhiResourceUpdateBatch;
class QRhiShaderResourceBindings;

namespace amcset {
namespace gui {

class TrackStore;
class TrackOverlay;
struct SampleScene;

class TrackView : public QRhiWidget {
  Q_OBJECT

public:
  //! x = depth, y and z lateral (TRIM's axes).
  enum class Projection { depth_y, depth_z, beam, perspective };
  enum class Arrangement { single, quad };
  enum class ColorBy { element, generation };

  explicit TrackView(TrackStore *store, QWidget *parent = nullptr);
  ~TrackView() override;

  void set_arrangement(Arrangement a);
  Arrangement arrangement() const noexcept { return arrangement_; }
  //! The projection shown in the single arrangement.
  void set_single_projection(Projection p);
  Projection single_projection() const noexcept { return single_; }
  void set_color_by(ColorBy c);
  void set_show_incident(bool on);
  void set_show_recoils(bool on);
  //! Deepest recoil generation drawn (1 = primary recoils only).
  void set_max_generation(int g);
  //! Opacity of the paths (0.05-1; low values show the density of many paths).
  void set_opacity(double a);
  //! The target: layer boundaries from the surface (Angstrom, first = 0) and
  //! the atomic number of the incident particle (0 = electrons).
  void set_target(const std::vector<double> &boundaries, unsigned incident_z);
  //! Call after the store grew or was cleared.
  void data_changed();
  //! Fit every view to the target and the paths.
  void reset_views();
  /*!
   * Fit the views to the paths alone, not to the whole target depth (for
   * recoil sources: a cascade of a few nanometres deep inside the target
   * would be a dot in a view of the target). Off: the target and the paths.
   */
  void set_fit_to_paths(bool on);
  /*!
   * The whole beam spot: every path at its true place on the sample (its
   * ion's entry point added), the entry points marked on the surface, and
   * the sample (set_scene) drawn under them. Off: every path from the
   * origin, as the results give them.
   */
  void set_spot_mode(bool on);
  bool spot_mode() const noexcept { return spot_mode_; }
  //! The sample the spot mode draws (null: none yet).
  void set_scene(std::shared_ptr<const SampleScene> scene);
  const SampleScene *scene() const noexcept { return scene_.get(); }

  //! Colours of the paths (also for legends).
  static QColor incident_color();
  static QColor element_color(unsigned z);
  static QColor generation_color(int generation);

  // For the overlay.
  struct Viewport {
    Projection projection;
    QRectF rect; //!< Device-independent pixels within the widget
  };
  std::vector<Viewport> viewports() const;
  //! World rectangle (horizontal, vertical axis) shown in a 2D viewport.
  QRectF visible_world(const Viewport &v) const;
  ColorBy color_by() const noexcept { return color_by_; }
  unsigned incident_z() const noexcept { return incident_z_; }
  const TrackStore *store() const noexcept { return store_; }
  bool render_failed() const noexcept { return failed_; }

signals:
  //! Mouse position in target coordinates, or an empty text outside a 2D view.
  void hover_text(const QString &text);
  //! The graphics API and device the views draw with, once initialised.
  void graphics_ready(const QString &description);

protected:
  void initialize(QRhiCommandBuffer *cb) override;
  void render(QRhiCommandBuffer *cb) override;
  void releaseResources() override;
  void resizeEvent(QResizeEvent *e) override;
  void mousePressEvent(QMouseEvent *e) override;
  void mouseMoveEvent(QMouseEvent *e) override;
  void mouseReleaseEvent(QMouseEvent *e) override;
  void mouseDoubleClickEvent(QMouseEvent *e) override;
  void wheelEvent(QWheelEvent *e) override;

private:
  struct Camera {
    // 2D: world point at the viewport centre and pixels per Angstrom.
    QPointF center;
    double scale = 1.0;
    // 3D: orbit around target.
    QVector3D target;
    float yaw = -35.0f;
    float pitch = 25.0f;
    float distance = 1000.0f;
    bool user_moved = false;
  };
  struct GpuVertex {
    float x, y, z;
    quint32 rgba;
    float generation;
  };

  int viewport_at(const QPointF &p) const;
  Camera &camera(Projection p) { return cameras_[static_cast<std::size_t>(p)]; }
  const Camera &camera(Projection p) const { return cameras_[static_cast<std::size_t>(p)]; }
  void fit(Projection p, const QRectF &rect);
  QMatrix4x4 matrix(Projection p, const QRectF &rect, const QSize &pixel_size) const;
  void append_vertices(std::size_t first_segment, std::size_t last_segment, std::vector<GpuVertex> &out) const;
  quint32 color_of(unsigned z, int generation) const;
  void build_outline(std::vector<GpuVertex> &out) const;
  void upload_scene(QRhiResourceUpdateBatch *u);
  void upload_marks(QRhiResourceUpdateBatch *u);
  QRhiGraphicsPipeline *make_pipeline(bool triangles);
  // World (scene) position of a TRIM point for the 3D view: up = -depth.
  static QVector3D scene(float x, float y, float z) { return {y, -x, z}; }

  TrackStore *store_;
  TrackOverlay *overlay_;
  Arrangement arrangement_ = Arrangement::quad;
  Projection single_ = Projection::depth_y;
  ColorBy color_by_ = ColorBy::element;
  bool show_incident_ = true;
  bool show_recoils_ = true;
  int max_generation_ = 255;
  float opacity_ = 0.8f;
  std::vector<double> boundaries_{0.0, 1000.0};
  unsigned incident_z_ = 0;
  bool fit_to_paths_ = false;
  std::array<Camera, 4> cameras_;
  bool failed_ = false;

  // Interaction.
  int drag_view_ = -1;
  QPointF drag_last_;
  Qt::MouseButton drag_button_ = Qt::NoButton;

  // GPU resources.
  QRhi *rhi_ = nullptr;
  std::unique_ptr<QRhiBuffer> vbuf_;
  std::unique_ptr<QRhiBuffer> outline_buf_;
  // Per viewport: [0..3] recoil pass, [4..7] incident-particle pass (drawn on top).
  std::array<std::unique_ptr<QRhiBuffer>, 8> ubufs_;
  std::array<std::unique_ptr<QRhiShaderResourceBindings>, 8> srbs_;
  std::unique_ptr<QRhiGraphicsPipeline> pipeline_;
  std::size_t capacity_ = 0;       // segments the vertex buffer holds
  std::size_t uploaded_ = 0;       // segments uploaded
  quint64 uploaded_revision_ = 0;  // store revision of the upload
  bool reupload_ = true;           // colours or store changed
  bool outline_dirty_ = true;
  std::size_t outline_2d_ = 0;     // vertices of the depth outline (2D)
  std::size_t outline_3d_ = 0;     // vertices of the 3D outline

  // The whole beam spot (set_spot_mode).
  bool spot_mode_ = false;
  std::shared_ptr<const SampleScene> scene_;
  bool scene_dirty_ = true;
  std::unique_ptr<QRhiGraphicsPipeline> tri_pipeline_; // the sample's filled cells
  std::unique_ptr<QRhiBuffer> scene_tri_buf_;
  std::unique_ptr<QRhiBuffer> scene_line_buf_;
  std::array<quint32, 4> tri_first_{}, tri_count_{}, line_first_{}, line_count_{}; // per projection, in vertices
  quint32 spot_first_ = 0, spot_count_ = 0; // the spot's outline (in the line buffer, drawn over the paths)
  std::unique_ptr<QRhiBuffer> marks_buf_;  // the entry points (crosses on the surface)
  std::size_t marks_capacity_ = 0;         // entries the buffer holds
  std::size_t marks_uploaded_ = 0;
  quint64 marks_revision_ = 0;
};

} // namespace gui
} // namespace amcset
