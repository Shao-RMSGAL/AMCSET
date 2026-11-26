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
#include "views/track_view.h"

#include "model/quantity_text.h"
#include "run/track_store.h"
#include "views/axis_ticks.h"

#include <QFile>
#include <QMouseEvent>
#include <QPainter>
#include <QWheelEvent>
#include <rhi/qrhi.h>

#include <algorithm>
#include <cmath>

namespace amcset {
namespace gui {

// ---------------------------------------------------------------------------
// Overlay: axes, labels and frame painted over the GPU image
// ---------------------------------------------------------------------------

class TrackOverlay : public QWidget {
public:
  explicit TrackOverlay(TrackView *view) : QWidget(view), view_(view) {
    setAttribute(Qt::WA_TransparentForMouseEvents);
    setAttribute(Qt::WA_NoSystemBackground);
  }

protected:
  void paintEvent(QPaintEvent *) override {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing, false);
    const QColor text = palette().color(QPalette::WindowText);
    QColor faint = text;
    faint.setAlphaF(0.35f);
    QFont small = font();
    small.setPointSizeF(std::max(7.0, small.pointSizeF() * 0.85));
    p.setFont(small);
    const QFontMetricsF fm(small);
    if (view_->render_failed()) {
      p.setPen(text);
      p.drawText(rect(), Qt::AlignCenter | Qt::TextWordWrap,
                 tr("The graphics device could not be initialised, so the paths cannot be drawn.\n"
                    "Set the environment variable AMCSET_RHI to opengl, vulkan, d3d11 or d3d12\n"
                    "to choose another graphics API."));
      return;
    }
    const auto vps = view_->viewports();
    for (const auto &vp : vps) {
      const QRectF r = vp.rect;
      p.setPen(faint);
      if (vps.size() > 1) p.drawRect(r.adjusted(0, 0, -1, -1));
      p.setPen(text);
      QString title;
      switch (vp.projection) {
      case TrackView::Projection::depth_y: title = tr("Depth (X) - Y"); break;
      case TrackView::Projection::depth_z: title = tr("Depth (X) - Z"); break;
      case TrackView::Projection::beam: title = tr("Beam view (Y - Z)"); break;
      case TrackView::Projection::perspective: title = tr("3D  (drag: rotate, right-drag: pan, wheel: zoom)"); break;
      }
      p.drawText(r.adjusted(6, 4, -6, -4), Qt::AlignLeft | Qt::AlignTop, title);
      if (vp.projection == TrackView::Projection::perspective) continue;
      // Ticks along the bottom and the left edge.
      const QRectF w = view_->visible_world(vp);
      const auto hx = nice_ticks(w.left(), w.right(), std::max(2, static_cast<int>(r.width() / 90)));
      const auto vy = nice_ticks(w.top(), w.bottom(), std::max(2, static_cast<int>(r.height() / 60)));
      const auto unit = length_axis_unit(std::max(std::abs(w.left()), std::abs(w.right())));
      const auto unit_v = length_axis_unit(std::max(std::abs(w.top()), std::abs(w.bottom())));
      for (double t : hx) {
        const double sx = r.left() + (t - w.left()) / w.width() * r.width();
        if (sx < r.left() + 30 || sx > r.right() - 20) continue;
        p.setPen(faint);
        p.drawLine(QPointF(sx, r.bottom() - 5), QPointF(sx, r.bottom()));
        p.setPen(text);
        const QString s = format_tick(t / unit.factor);
        p.drawText(QPointF(sx - fm.horizontalAdvance(s) / 2, r.bottom() - 7), s);
      }
      for (double t : vy) {
        const double sy = r.bottom() - (t - w.top()) / w.height() * r.height();
        if (sy > r.bottom() - 30 || sy < r.top() + 3 * fm.height()) continue;
        p.setPen(faint);
        p.drawLine(QPointF(r.left(), sy), QPointF(r.left() + 5, sy));
        p.setPen(text);
        p.drawText(QPointF(r.left() + 7, sy + fm.ascent() / 2 - 1), format_tick(t / unit_v.factor));
      }
      const QString hl = vp.projection == TrackView::Projection::beam ? tr("Y") : tr("Depth");
      const QString vl = vp.projection == TrackView::Projection::depth_y ? tr("Y") : tr("Z");
      const QString label = QStringLiteral("%1 (%2)").arg(hl, unit.name);
      p.drawText(QPointF(r.right() - fm.horizontalAdvance(label) - 6, r.bottom() - 7 - fm.height()), label);
      p.drawText(QPointF(r.left() + 7, r.top() + 6 + 2 * fm.height()), QStringLiteral("%1 (%2)").arg(vl, unit_v.name));
    }
    // Legend.
    const auto *store = view_->store();
    if (store == nullptr || store->size() == 0) {
      p.setPen(faint);
      p.drawText(rect(), Qt::AlignCenter, tr("Particle paths appear here while a calculation runs."));
      return;
    }
    QStringList names;
    QList<QColor> colors;
    if (view_->color_by() == TrackView::ColorBy::element) {
      const auto &sp = store->species();
      const unsigned iz = view_->incident_z();
      names << (iz == 0 ? tr("Incident electrons") : tr("Incident %1").arg(element_symbol(iz)));
      colors << TrackView::incident_color();
      for (unsigned z = 1; z < sp.size(); ++z) {
        if (sp[z] == 0) continue;
        names << tr("%1 recoils").arg(element_symbol(z));
        colors << TrackView::element_color(z);
      }
    } else {
      const QStringList g{tr("Incident"), tr("Primary recoils"), tr("2nd generation"), tr("3rd+")};
      for (int i = 0; i < 4; ++i) {
        names << g[i];
        colors << TrackView::generation_color(i);
      }
    }
    double y = 6;
    const double x = width() - 6;
    for (int i = 0; i < names.size(); ++i) {
      const double tw = fm.horizontalAdvance(names[i]);
      p.fillRect(QRectF(x - tw - 18, y + 3, 10, fm.height() - 6), colors[i]);
      p.setPen(text);
      p.drawText(QPointF(x - tw, y + fm.ascent()), names[i]);
      y += fm.height();
    }
    p.setPen(faint);
    const QString n = tr("%L1 paths, %L2 segments%3")
                          .arg(store->histories())
                          .arg(store->size())
                          .arg(store->full() ? tr(" (limit reached)") : QString());
    p.drawText(QPointF(width() - fm.horizontalAdvance(n) - 6, height() - 6), n);
  }

private:
  TrackView *view_;
};

// ---------------------------------------------------------------------------
// Colours
// ---------------------------------------------------------------------------

QColor TrackView::incident_color() { return QColor(230, 50, 40); }

QColor TrackView::element_color(unsigned z) {
  // Distinct hues from the golden ratio, avoiding the incident particle's red.
  const double h = std::fmod(0.13 + 0.618034 * static_cast<double>(z), 1.0);
  const double hue = 0.1 + 0.8 * h;
  return QColor::fromHsvF(static_cast<float>(hue), 0.7f, 0.92f);
}

QColor TrackView::generation_color(int generation) {
  switch (std::min(generation, 3)) {
  case 0: return QColor(230, 60, 50);
  case 1: return QColor(70, 140, 240);
  case 2: return QColor(60, 190, 110);
  default: return QColor(235, 170, 40);
  }
}

quint32 TrackView::color_of(unsigned z, int generation) const {
  const QColor c = color_by_ == ColorBy::generation ? generation_color(generation)
                   : generation == 0                ? incident_color()
                                                    : element_color(z);
  // RGBA bytes in memory order (UNormByte4).
  return static_cast<quint32>(c.red()) | (static_cast<quint32>(c.green()) << 8) |
         (static_cast<quint32>(c.blue()) << 16) | (255u << 24);
}

// ---------------------------------------------------------------------------
// Construction and settings
// ---------------------------------------------------------------------------

TrackView::TrackView(TrackStore *store, QWidget *parent)
    : QRhiWidget(parent), store_(store), overlay_(new TrackOverlay(this)) {
  const QByteArray api = qgetenv("AMCSET_RHI").toLower();
  if (api == "opengl") setApi(Api::OpenGL);
  else if (api == "vulkan") setApi(Api::Vulkan);
  else if (api == "d3d11") setApi(Api::Direct3D11);
  else if (api == "d3d12") setApi(Api::Direct3D12);
  else if (api == "metal") setApi(Api::Metal);
  else if (api == "null") setApi(Api::Null);
  setSampleCount(4);
  setMouseTracking(true);
  setMinimumSize(200, 160);
  connect(this, &QRhiWidget::renderFailed, this, guarded(this, tr("Path view (renderFailed)"), [this] {
    failed_ = true;
    overlay_->update();
  }));
}

TrackView::~TrackView() = default;

void TrackView::set_arrangement(Arrangement a) {
  arrangement_ = a;
  overlay_->update();
  update();
}

void TrackView::set_single_projection(Projection p) {
  single_ = p;
  overlay_->update();
  update();
}

void TrackView::set_color_by(ColorBy c) {
  if (c == color_by_) return;
  color_by_ = c;
  reupload_ = true;
  overlay_->update();
  update();
}

void TrackView::set_show_incident(bool on) {
  show_incident_ = on;
  update();
}

void TrackView::set_show_recoils(bool on) {
  show_recoils_ = on;
  update();
}

void TrackView::set_max_generation(int g) {
  max_generation_ = std::max(1, g);
  update();
}

void TrackView::set_opacity(double a) {
  opacity_ = static_cast<float>(std::clamp(a, 0.02, 1.0));
  update();
}

void TrackView::set_target(const std::vector<double> &boundaries, unsigned incident_z) {
  if (boundaries.size() >= 2) boundaries_ = boundaries;
  if (incident_z != incident_z_) {
    incident_z_ = incident_z;
    reupload_ = true;
  }
  outline_dirty_ = true;
  for (auto &c : cameras_) {
    if (!c.user_moved) c.scale = -1.0; // refit on the next frame
  }
  overlay_->update();
  update();
}

void TrackView::data_changed() {
  outline_dirty_ = true;
  for (auto &c : cameras_) {
    if (!c.user_moved) c.scale = -1.0;
  }
  overlay_->update();
  update();
}

void TrackView::reset_views() {
  for (auto &c : cameras_) {
    c.user_moved = false;
    c.scale = -1.0;
  }
  outline_dirty_ = true;
  overlay_->update();
  update();
}

void TrackView::set_fit_to_paths(bool on) {
  if (on == fit_to_paths_) return;
  fit_to_paths_ = on;
  reset_views();
}

void TrackView::resizeEvent(QResizeEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    QRhiWidget::resizeEvent(e);
    overlay_->setGeometry(rect());
    for (auto &c : cameras_) {
      if (!c.user_moved) c.scale = -1.0;
    }
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (resizeEvent)"), std::current_exception());
  }
}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

std::vector<TrackView::Viewport> TrackView::viewports() const {
  const QRectF r = QRectF(rect());
  if (arrangement_ == Arrangement::single) return {{single_, r}};
  const double w = r.width() / 2.0;
  const double h = r.height() / 2.0;
  return {{Projection::depth_y, QRectF(0, 0, w, h)},
          {Projection::depth_z, QRectF(w, 0, r.width() - w, h)},
          {Projection::beam, QRectF(0, h, w, r.height() - h)},
          {Projection::perspective, QRectF(w, h, r.width() - w, r.height() - h)}};
}

int TrackView::viewport_at(const QPointF &p) const {
  const auto v = viewports();
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (v[i].rect.contains(p)) return static_cast<int>(i);
  }
  return -1;
}

QRectF TrackView::visible_world(const Viewport &v) const {
  const Camera &c = camera(v.projection);
  const double s = c.scale > 0 ? c.scale : 1.0;
  const double hw = v.rect.width() / (2.0 * s);
  const double hh = v.rect.height() / (2.0 * s);
  // QRectF with top = lowest vertical value (not screen orientation).
  return QRectF(c.center.x() - hw, c.center.y() - hh, 2 * hw, 2 * hh);
}

void TrackView::fit(Projection p, const QRectF &rect) {
  Camera &c = camera(p);
  const double thickness = boundaries_.back();
  const float *lo = store_->lower();
  const float *hi = store_->upper();
  const bool have = store_->size() > 0;
  const double x0 = std::min(0.0, have ? static_cast<double>(lo[0]) : 0.0);
  const double x1 = std::max(thickness, have ? static_cast<double>(hi[0]) : thickness);
  double lat = 0.0;
  if (have) {
    for (int k = 1; k < 3; ++k) lat = std::max({lat, std::abs(static_cast<double>(lo[k])), std::abs(static_cast<double>(hi[k]))});
  }
  lat = std::max(lat, 0.05 * (x1 - x0));
  if (lat <= 0) lat = 1.0;
  const double w = std::max(rect.width(), 10.0);
  const double h = std::max(rect.height(), 10.0);
  if (fit_to_paths_ && have) {
    // The box of the paths (TRIM axes, Angstrom), a little larger, with no
    // side shorter than 2 % of the longest (or 1 A).
    double c3[3], d3[3];
    for (int k = 0; k < 3; ++k) {
      c3[k] = 0.5 * (static_cast<double>(lo[k]) + static_cast<double>(hi[k]));
      d3[k] = static_cast<double>(hi[k]) - static_cast<double>(lo[k]);
    }
    const double longest = std::max({d3[0], d3[1], d3[2], 1.0});
    for (double &d : d3) d = 1.1 * std::max(d, 0.02 * longest);
    switch (p) {
    case Projection::depth_y:
    case Projection::depth_z: {
      const int v = p == Projection::depth_y ? 1 : 2;
      c.center = QPointF(c3[0], c3[v]);
      c.scale = 0.9 * std::min(w / d3[0], h / d3[v]);
      break;
    }
    case Projection::beam:
      c.center = QPointF(c3[1], c3[2]);
      c.scale = 0.9 * std::min(w / d3[1], h / d3[2]);
      break;
    case Projection::perspective:
      c.target = scene(static_cast<float>(c3[0]), static_cast<float>(c3[1]), static_cast<float>(c3[2]));
      c.distance = 1.3f * static_cast<float>(std::max({d3[0], d3[1], d3[2]})) / std::tan(qDegreesToRadians(20.0f));
      c.scale = 1.0;
      break;
    }
    return;
  }
  switch (p) {
  case Projection::depth_y:
  case Projection::depth_z:
    c.center = QPointF(0.5 * (x0 + x1), 0.0);
    c.scale = 0.9 * std::min(w / (x1 - x0), h / (2.0 * lat));
    break;
  case Projection::beam:
    c.center = QPointF(0.0, 0.0);
    c.scale = 0.9 * std::min(w, h) / (2.0 * lat);
    break;
  case Projection::perspective: {
    const float depth = static_cast<float>(x1 - x0);
    c.target = scene(static_cast<float>(0.5 * (x0 + x1)), 0.0f, 0.0f);
    const float extent = std::max(depth, 2.0f * static_cast<float>(lat));
    c.distance = 1.3f * extent / std::tan(qDegreesToRadians(20.0f));
    c.scale = 1.0; // marks the camera as fitted
    break;
  }
  }
}

QMatrix4x4 TrackView::matrix(Projection p, const QRectF &rect, const QSize &) const {
  const Camera &c = camera(p);
  QMatrix4x4 m;
  if (p == Projection::perspective) {
    const float aspect = static_cast<float>(rect.width() / std::max(rect.height(), 1.0));
    QMatrix4x4 proj;
    const float near_plane = std::max(c.distance * 0.001f, 1e-3f);
    proj.perspective(40.0f, aspect, near_plane, c.distance * 20.0f + 1.0f);
    const float yaw = qDegreesToRadians(c.yaw);
    const float pitch = qDegreesToRadians(c.pitch);
    const QVector3D eye = c.target + c.distance * QVector3D(std::cos(pitch) * std::sin(yaw), std::sin(pitch),
                                                            std::cos(pitch) * std::cos(yaw));
    QMatrix4x4 view;
    view.lookAt(eye, c.target, QVector3D(0, 1, 0));
    // TRIM (x depth, y, z) -> scene (y, -x, z).
    const QMatrix4x4 model(0, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
    m = proj * view * model;
  } else {
    const double sx = 2.0 * c.scale / std::max(rect.width(), 1.0);
    const double sy = 2.0 * c.scale / std::max(rect.height(), 1.0);
    // Horizontal and vertical world axes of the projection.
    int h = 0, v = 1;
    if (p == Projection::depth_z) v = 2;
    if (p == Projection::beam) {
      h = 1;
      v = 2;
    }
    float r0[4] = {0, 0, 0, static_cast<float>(-c.center.x() * sx)};
    float r1[4] = {0, 0, 0, static_cast<float>(-c.center.y() * sy)};
    r0[h] = static_cast<float>(sx);
    r1[v] = static_cast<float>(sy);
    m = QMatrix4x4(r0[0], r0[1], r0[2], r0[3], r1[0], r1[1], r1[2], r1[3], 0, 0, 0, 0, 0, 0, 0, 1);
  }
  return rhi_->clipSpaceCorrMatrix() * m;
}

void TrackView::append_vertices(std::size_t first, std::size_t last, std::vector<GpuVertex> &out) const {
  const auto &segs = store_->segments();
  out.reserve(out.size() + 2 * (last - first));
  for (std::size_t i = first; i < last; ++i) {
    const auto &s = segs[i];
    const quint32 col = color_of(s.z, s.generation);
    const float g = static_cast<float>(s.generation);
    out.push_back({s.a[0], s.a[1], s.a[2], col, g});
    out.push_back({s.b[0], s.b[1], s.b[2], col, g});
  }
}

void TrackView::build_outline(std::vector<GpuVertex> &out) const {
  out.clear();
  const QColor fg = palette().color(QPalette::WindowText);
  const quint32 strong = static_cast<quint32>(fg.red()) | (static_cast<quint32>(fg.green()) << 8) |
                         (static_cast<quint32>(fg.blue()) << 16) | (150u << 24);
  const quint32 weak = (strong & 0x00ffffffu) | (70u << 24);
  // 2D (depth projections): long lines at the surface and every interface.
  const float big = 1e9f;
  for (std::size_t i = 0; i < boundaries_.size(); ++i) {
    const float x = static_cast<float>(boundaries_[i]);
    const quint32 c = i == 0 || i + 1 == boundaries_.size() ? strong : weak;
    // A line in the x = const plane spanning y and z (projects to a line in
    // both depth views).
    out.push_back({x, -big, -big, c, -1.0f});
    out.push_back({x, big, big, c, -1.0f});
  }
  const auto n2d = out.size();
  // 3D: rectangles at the surface and the interfaces, lateral size from the data.
  float lat = static_cast<float>(0.25 * boundaries_.back());
  if (store_->size() > 0) {
    for (int k = 1; k < 3; ++k) lat = std::max({lat, std::abs(store_->lower()[k]), std::abs(store_->upper()[k])});
  }
  for (std::size_t i = 0; i < boundaries_.size(); ++i) {
    const float x = static_cast<float>(boundaries_[i]);
    const quint32 c = i == 0 || i + 1 == boundaries_.size() ? strong : weak;
    const float corners[4][2] = {{-lat, -lat}, {lat, -lat}, {lat, lat}, {-lat, lat}};
    for (int k = 0; k < 4; ++k) {
      const auto &a = corners[k];
      const auto &b = corners[(k + 1) % 4];
      out.push_back({x, a[0], a[1], c, -1.0f});
      out.push_back({x, b[0], b[1], c, -1.0f});
    }
  }
  // Vertical edges of the target box.
  const float x0 = 0.0f;
  const float x1 = static_cast<float>(boundaries_.back());
  const float corners[4][2] = {{-lat, -lat}, {lat, -lat}, {lat, lat}, {-lat, lat}};
  for (const auto &cc : corners) {
    out.push_back({x0, cc[0], cc[1], weak, -1.0f});
    out.push_back({x1, cc[0], cc[1], weak, -1.0f});
  }
  // (the 3D part follows the 2D part in the buffer)
  (void)n2d;
}

// ---------------------------------------------------------------------------
// Rendering
// ---------------------------------------------------------------------------

namespace {

QShader load_shader(const QString &name) {
  QFile f(name);
  if (f.open(QIODevice::ReadOnly)) return QShader::fromSerialized(f.readAll());
  return {};
}

} // namespace

void TrackView::initialize(QRhiCommandBuffer *) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (rhi_ != rhi() || (pipeline_ && pipeline_->sampleCount() != sampleCount())) {
      pipeline_.reset();
      vbuf_.reset();
      outline_buf_.reset();
      for (auto &b : ubufs_) b.reset();
      for (auto &s : srbs_) s.reset();
      rhi_ = rhi();
      capacity_ = 0;
      uploaded_ = 0;
      outline_dirty_ = true;
    }
    if (rhi_ == nullptr || pipeline_) return;
    const QRhiDriverInfo info = rhi_->driverInfo();
    emit graphics_ready(QStringLiteral("%1, %2").arg(QString::fromLatin1(rhi_->backendName()),
                                                     QString::fromUtf8(info.deviceName)));
    for (std::size_t i = 0; i < ubufs_.size(); ++i) {
      ubufs_[i].reset(rhi_->newBuffer(QRhiBuffer::Dynamic, QRhiBuffer::UniformBuffer, 80));
      ubufs_[i]->create();
      srbs_[i].reset(rhi_->newShaderResourceBindings());
      srbs_[i]->setBindings({QRhiShaderResourceBinding::uniformBuffer(0, QRhiShaderResourceBinding::VertexStage,
                                                                      ubufs_[i].get())});
      srbs_[i]->create();
    }
    pipeline_.reset(rhi_->newGraphicsPipeline());
    pipeline_->setTopology(QRhiGraphicsPipeline::Lines);
    QRhiGraphicsPipeline::TargetBlend blend;
    blend.enable = true;
    blend.srcColor = QRhiGraphicsPipeline::One;
    blend.dstColor = QRhiGraphicsPipeline::OneMinusSrcAlpha;
    blend.srcAlpha = QRhiGraphicsPipeline::One;
    blend.dstAlpha = QRhiGraphicsPipeline::OneMinusSrcAlpha;
    pipeline_->setTargetBlends({blend});
    pipeline_->setDepthTest(false);
    pipeline_->setDepthWrite(false);
    pipeline_->setSampleCount(sampleCount());
    pipeline_->setShaderStages({{QRhiShaderStage::Vertex, load_shader(QStringLiteral(":/shaders/tracks.vert.qsb"))},
                                {QRhiShaderStage::Fragment, load_shader(QStringLiteral(":/shaders/tracks.frag.qsb"))}});
    QRhiVertexInputLayout layout;
    layout.setBindings({{sizeof(GpuVertex)}});
    layout.setAttributes({{0, 0, QRhiVertexInputAttribute::Float3, offsetof(GpuVertex, x)},
                          {0, 1, QRhiVertexInputAttribute::UNormByte4, offsetof(GpuVertex, rgba)},
                          {0, 2, QRhiVertexInputAttribute::Float, offsetof(GpuVertex, generation)}});
    pipeline_->setVertexInputLayout(layout);
    pipeline_->setShaderResourceBindings(srbs_[0].get());
    pipeline_->setRenderPassDescriptor(renderTarget()->renderPassDescriptor());
    if (!pipeline_->create()) {
      failed_ = true;
      pipeline_.reset();
    }
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (initialize)"), std::current_exception());
  }
}

void TrackView::render(QRhiCommandBuffer *cb) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    QColor bg = palette().color(QPalette::Base);
    if (!rhi_ || !pipeline_) {
      cb->beginPass(renderTarget(), bg, {1.0f, 0});
      cb->endPass();
      return;
    }
    QRhiResourceUpdateBatch *u = rhi_->nextResourceUpdateBatch();

    // Paths: upload what is new.
    if (store_->revision() != uploaded_revision_ || reupload_) {
      uploaded_ = 0;
      uploaded_revision_ = store_->revision();
      reupload_ = false;
    }
    const std::size_t n = store_->size();
    if (n > capacity_) {
      capacity_ = std::max<std::size_t>({n + n / 2, 65536, capacity_ * 2});
      capacity_ = std::min(capacity_, std::max(n, store_->max_segments()));
      vbuf_.reset(rhi_->newBuffer(QRhiBuffer::Static, QRhiBuffer::VertexBuffer,
                                  static_cast<quint32>(capacity_ * 2 * sizeof(GpuVertex))));
      if (!vbuf_->create()) {
        vbuf_.reset();
        capacity_ = 0;
      }
      uploaded_ = 0;
    }
    if (vbuf_ && uploaded_ < n) {
      std::vector<GpuVertex> v;
      append_vertices(uploaded_, n, v);
      u->uploadStaticBuffer(vbuf_.get(), static_cast<quint32>(uploaded_ * 2 * sizeof(GpuVertex)),
                            static_cast<quint32>(v.size() * sizeof(GpuVertex)), v.data());
      uploaded_ = n;
    }
    // Outlines.
    if (outline_dirty_) {
      std::vector<GpuVertex> o;
      build_outline(o);
      outline_2d_ = 2 * boundaries_.size();
      outline_3d_ = o.size() - outline_2d_;
      const auto bytes = static_cast<quint32>(o.size() * sizeof(GpuVertex));
      if (!outline_buf_ || outline_buf_->size() < bytes) {
        outline_buf_.reset(rhi_->newBuffer(QRhiBuffer::Dynamic, QRhiBuffer::VertexBuffer, std::max(bytes, 4096u)));
        outline_buf_->create();
      }
      u->updateDynamicBuffer(outline_buf_.get(), 0, bytes, o.data());
      outline_dirty_ = false;
    }

    // Cameras and uniforms.
    const QSize px = renderTarget()->pixelSize();
    const double ratio = width() > 0 ? static_cast<double>(px.width()) / width() : 1.0;
    const auto vps = viewports();
    for (std::size_t i = 0; i < vps.size(); ++i) {
      Camera &c = camera(vps[i].projection);
      if (c.scale <= 0) fit(vps[i].projection, vps[i].rect);
      const QMatrix4x4 mvp = matrix(vps[i].projection, vps[i].rect, px);
      struct {
        float mvp[16];
        float params[4];
      } ub;
      std::memcpy(ub.mvp, mvp.constData(), sizeof ub.mvp);
      // Recoils first, then the incident particles on top of them.
      ub.params[0] = 0.0f;
      ub.params[1] = show_recoils_ ? 1.0f : 0.0f;
      ub.params[2] = opacity_;
      ub.params[3] = static_cast<float>(max_generation_);
      u->updateDynamicBuffer(ubufs_[i].get(), 0, sizeof ub, &ub);
      ub.params[0] = show_incident_ ? 1.0f : 0.0f;
      ub.params[1] = 0.0f;
      u->updateDynamicBuffer(ubufs_[i + 4].get(), 0, sizeof ub, &ub);
    }

    cb->beginPass(renderTarget(), bg, {1.0f, 0}, u);
    cb->setGraphicsPipeline(pipeline_.get());
    for (std::size_t i = 0; i < vps.size(); ++i) {
      const QRectF r = vps[i].rect;
      const float x = static_cast<float>(r.left() * ratio);
      const float y = static_cast<float>(px.height() - r.bottom() * ratio);
      cb->setViewport(QRhiViewport(x, std::max(0.0f, y), static_cast<float>(r.width() * ratio),
                                   static_cast<float>(r.height() * ratio)));
      cb->setShaderResources(srbs_[i].get());
      const bool is3d = vps[i].projection == Projection::perspective;
      if (vps[i].projection != Projection::beam && outline_buf_) {
        const quint32 offset = is3d ? static_cast<quint32>(outline_2d_ * sizeof(GpuVertex)) : 0;
        const QRhiCommandBuffer::VertexInput in(outline_buf_.get(), offset);
        cb->setVertexInput(0, 1, &in);
        cb->draw(static_cast<quint32>(is3d ? outline_3d_ : outline_2d_));
      }
      if (vbuf_ && uploaded_ > 0) {
        const QRhiCommandBuffer::VertexInput in(vbuf_.get(), 0);
        for (std::size_t pass = 0; pass < 2; ++pass) {
          if (pass == 0 && !show_recoils_) continue;
          if (pass == 1 && !show_incident_) continue;
          cb->setShaderResources(srbs_[i + 4 * pass].get());
          cb->setVertexInput(0, 1, &in);
          cb->draw(static_cast<quint32>(2 * uploaded_));
        }
      }
    }
    cb->endPass();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (render)"), std::current_exception());
  }
}

void TrackView::releaseResources() {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    pipeline_.reset();
    vbuf_.reset();
    outline_buf_.reset();
    for (auto &b : ubufs_) b.reset();
    for (auto &s : srbs_) s.reset();
    rhi_ = nullptr;
    capacity_ = 0;
    uploaded_ = 0;
    outline_dirty_ = true;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (releaseResources)"), std::current_exception());
  }
}

// ---------------------------------------------------------------------------
// Interaction
// ---------------------------------------------------------------------------

void TrackView::mousePressEvent(QMouseEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    drag_view_ = viewport_at(e->position());
    drag_last_ = e->position();
    drag_button_ = e->button();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (mousePressEvent)"), std::current_exception());
  }
}

void TrackView::mouseMoveEvent(QMouseEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    const auto vps = viewports();
    if (drag_view_ >= 0 && drag_view_ < static_cast<int>(vps.size()) && (e->buttons() & drag_button_)) {
      const auto &vp = vps[static_cast<std::size_t>(drag_view_)];
      Camera &c = camera(vp.projection);
      const QPointF d = e->position() - drag_last_;
      drag_last_ = e->position();
      c.user_moved = true;
      if (vp.projection == Projection::perspective) {
        if (drag_button_ == Qt::LeftButton) {
          c.yaw -= static_cast<float>(d.x()) * 0.4f;
          c.pitch = std::clamp(c.pitch + static_cast<float>(d.y()) * 0.4f, -89.0f, 89.0f);
        } else {
          // Pan in the camera plane.
          const float yaw = qDegreesToRadians(c.yaw);
          const QVector3D right(std::cos(yaw), 0, -std::sin(yaw));
          const float k = c.distance * 0.0015f;
          c.target += -right * static_cast<float>(d.x()) * k + QVector3D(0, 1, 0) * static_cast<float>(d.y()) * k;
        }
      } else {
        c.center -= QPointF(d.x() / c.scale, -d.y() / c.scale);
      }
      overlay_->update();
      update();
    }
    // Readout.
    const int i = viewport_at(e->position());
    if (i >= 0 && vps[static_cast<std::size_t>(i)].projection != Projection::perspective) {
      const auto &vp = vps[static_cast<std::size_t>(i)];
      const QRectF w = visible_world(vp);
      const double h = w.left() + (e->position().x() - vp.rect.left()) / vp.rect.width() * w.width();
      const double v = w.top() + (vp.rect.bottom() - e->position().y()) / vp.rect.height() * w.height();
      const auto fmt = [](double a) { return format_length_value(a); };
      QString text;
      switch (vp.projection) {
      case Projection::depth_y: text = tr("depth %1, y %2").arg(fmt(h), fmt(v)); break;
      case Projection::depth_z: text = tr("depth %1, z %2").arg(fmt(h), fmt(v)); break;
      default: text = tr("y %1, z %2").arg(fmt(h), fmt(v)); break;
      }
      emit hover_text(text);
    } else {
      emit hover_text({});
    }
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (mouseMoveEvent)"), std::current_exception());
  }
}

void TrackView::mouseReleaseEvent(QMouseEvent *) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    drag_view_ = -1;
    drag_button_ = Qt::NoButton;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (mouseReleaseEvent)"), std::current_exception());
  }
}

void TrackView::mouseDoubleClickEvent(QMouseEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    const int i = viewport_at(e->position());
    if (i < 0) return;
    Camera &c = camera(viewports()[static_cast<std::size_t>(i)].projection);
    c.user_moved = false;
    c.scale = -1.0;
    c.yaw = -35.0f;
    c.pitch = 25.0f;
    overlay_->update();
    update();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (mouseDoubleClickEvent)"), std::current_exception());
  }
}

void TrackView::wheelEvent(QWheelEvent *e) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    const auto vps = viewports();
    const int i = viewport_at(e->position());
    if (i < 0) return;
    const auto &vp = vps[static_cast<std::size_t>(i)];
    Camera &c = camera(vp.projection);
    const double f = std::pow(1.0015, e->angleDelta().y());
    c.user_moved = true;
    if (vp.projection == Projection::perspective) {
      c.distance = static_cast<float>(c.distance / f);
    } else {
      // Zoom about the cursor.
      const QRectF w = visible_world(vp);
      const double h = w.left() + (e->position().x() - vp.rect.left()) / vp.rect.width() * w.width();
      const double v = w.top() + (vp.rect.bottom() - e->position().y()) / vp.rect.height() * w.height();
      c.scale *= f;
      c.center = QPointF(h - (h - c.center.x()) / f, v - (v - c.center.y()) / f);
    }
    overlay_->update();
    update();
    e->accept();
  } catch (...) {
    report_exception_quietly(QStringLiteral("Path view (wheelEvent)"), std::current_exception());
  }
}

} // namespace gui
} // namespace amcset
