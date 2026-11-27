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

#include "views/sample_scene.h"

#include "crystal/grain_map.h"
#include "crystal/orientation.h"

#include <QCoreApplication>

#include <algorithm>
#include <cmath>
#include <exception>
#include <unordered_map>

namespace amcset {
namespace gui {

namespace {

using crystal::Vec3;

constexpr double pi = 3.14159265358979323846;

std::uint32_t rgba(const QColor &c, int alpha) {
  return static_cast<std::uint32_t>(c.red()) | (static_cast<std::uint32_t>(c.green()) << 8) |
         (static_cast<std::uint32_t>(c.blue()) << 16) | (static_cast<std::uint32_t>(std::clamp(alpha, 0, 255)) << 24);
}

QString tr(const char *text) { return QCoreApplication::translate("SampleScene", text); }

//! One layer as the scene samples it.
struct LayerInfo {
  double top = 0.0, bottom = 0.0;
  bool crystal = false;
  bool cubic = false;
  crystal::Mat3 rotation;
  std::shared_ptr<const crystal::GrainMap> grains;
};

//! A cell's content: the layer (from 0) and the grain (-1: none), or -1 outside the target.
std::int64_t key_of(int layer, int grain) {
  return layer < 0 ? -1 : (static_cast<std::int64_t>(layer) << 32) | static_cast<std::uint32_t>(grain + 1);
}

class Builder {
public:
  Builder(const trim::TrimConfig &c, QColor line, const std::atomic<bool> *cancel, SampleScene &s)
      : c_(c), line_(line), cancel_(cancel), s_(s) {}

  void run() {
    const double tilt = common::in_units(c_.ion.angle, common::radian);
    s_.beam = {std::cos(tilt), std::sin(tilt), 0.0};
    s_.spot = c_.ion.spot;
    double top = 0.0;
    for (std::size_t i = 0; i < c_.layers.size(); ++i) {
      const auto &l = c_.layers[i];
      const double width = common::in_units(l.width, common::angstrom);
      LayerInfo info;
      info.top = top;
      info.bottom = top + width;
      SampleScene::Layer out;
      out.top_A = top;
      out.bottom_A = top + width;
      out.name = QString::fromStdString(l.name);
      if (l.crystal) {
        info.crystal = out.crystal = true;
        info.rotation = out.rotation = l.crystal->rotation();
        info.cubic = crystal::lattice_rotations(l.crystal->cell.lattice).size() == 24;
        out.lattice = l.crystal->cell.lattice;
        if (l.crystal->polycrystal) {
          info.grains = out.grains = crystal::make_grain_map(*l.crystal, width);
          s_.has_grains = true;
          if (!polycrystal_period_ && info.grains) {
            polycrystal_period_ = std::max(info.grains->period_y_A(), info.grains->period_z_A());
          }
        }
      }
      layers_.push_back(info);
      s_.layers.push_back(out);
      top += width;
    }
    s_.depth_A = top;
    s_.ipf = std::any_of(layers_.begin(), layers_.end(), [](const LayerInfo &l) { return l.crystal && l.cubic; });

    // The region: the spot with a margin (a point: one period of the grains).
    const double w = common::in_units(s_.spot.width, common::angstrom);
    const double h = common::in_units(s_.spot.height, common::angstrom);
    double hy = 0.0, hz = 0.0;
    switch (s_.spot.shape) {
    case trim::BeamSpot::Shape::point: hy = hz = 0.5 * (polycrystal_period_ > 0.0 ? polycrystal_period_ : 2000.0); break;
    case trim::BeamSpot::Shape::rectangle: hy = 0.5 * w, hz = 0.5 * h; break;
    case trim::BeamSpot::Shape::circle: hy = hz = 0.5 * w; break;
    case trim::BeamSpot::Shape::gaussian: hy = 1.274 * w, hz = 1.274 * h; break; // 3 sigma
    }
    hy *= 1.06;
    hz *= 1.06;
    s_.y0 = -hy, s_.y1 = hy, s_.z0 = -hz, s_.z1 = hz;
    if (!(s_.depth_A > 0.0) || !(hy > 0.0) || !(hz > 0.0)) {
      s_.message = tr("The target or the beam spot has no extent.");
      return;
    }

    const double depth = s_.depth_A;
    const double dy = s_.y1 - s_.y0, dz = s_.z1 - s_.z0;
    // The surface (beam view; 3D, more opaque), the sections through the
    // spot's centre (depth views) and two sides of the block (3D).
    // (Nearly opaque: the colours are the grains' own; the paths are drawn over them.)
    plane({0, s_.y0, s_.z0}, {0, dy, 0}, {0, 0, dz}, 215, {2});
    plane({0, s_.y0, s_.z0}, {0, dy, 0}, {0, 0, dz}, 240, {3});
    plane({0, s_.y0, 0}, {depth, 0, 0}, {0, dy, 0}, 190, {0});
    plane({0, 0, s_.z0}, {depth, 0, 0}, {0, 0, dz}, 190, {1});
    plane({0, s_.y0, s_.z0}, {depth, 0, 0}, {0, 0, dz}, 215, {3});
    plane({0, s_.y0, s_.z1}, {depth, 0, 0}, {0, dy, 0}, 215, {3});
    block_edges();
    spot_outline();
  }

private:
  int layer_at(double x) const {
    if (x < 0.0 || x > s_.depth_A) return -1;
    for (std::size_t i = 0; i < layers_.size(); ++i) {
      if (x <= layers_[i].bottom) return static_cast<int>(i);
    }
    return static_cast<int>(layers_.size()) - 1;
  }

  std::int64_t key_at(const Vec3 &p) const {
    const int li = layer_at(p.x);
    if (li < 0) return -1;
    const auto &l = layers_[static_cast<std::size_t>(li)];
    int grain = -1;
    if (l.grains) grain = l.grains->at({p.x - l.top, p.y, p.z}).grain;
    return key_of(li, grain);
  }

  std::uint32_t color_of(std::int64_t key, int alpha) {
    if (const auto it = colors_.find(key); it != colors_.end()) {
      return (it->second & 0x00ffffffu) | (static_cast<std::uint32_t>(alpha) << 24);
    }
    QColor c(150, 150, 150); // amorphous
    const int li = static_cast<int>(key >> 32);
    const int grain = static_cast<int>(key & 0xffffffff) - 1;
    const auto &l = layers_[static_cast<std::size_t>(li)];
    if (l.crystal) {
      const crystal::Mat3 &r = grain >= 0 ? l.grains->grain(grain).rotation : l.rotation;
      if (l.cubic) {
        c = ipf_color(r.transposed() * s_.beam);
      } else if (grain >= 0) {
        // Not cubic: a colour per grain (golden-ratio hues).
        const double hue = std::fmod(0.11 + 0.618034 * static_cast<double>(grain), 1.0);
        c = QColor::fromHsvF(static_cast<float>(hue), 0.55f, 0.92f);
      } else {
        c = QColor(214, 178, 110);
      }
    }
    const std::uint32_t v = rgba(c, 255);
    colors_.emplace(key, v);
    return (v & 0x00ffffffu) | (static_cast<std::uint32_t>(alpha) << 24);
  }

  /*!
   * A plane of cells: origin + (i + 1/2) / nu u + (j + 1/2) / nv v, at most
   * some 250 000; runs of one content along u become one quad, the edges
   * between different contents lines.
   */
  void plane(const Vec3 &origin, const Vec3 &u, const Vec3 &v, int alpha, std::initializer_list<int> views) {
    if (cancelled()) return;
    const double lu = crystal::norm(u), lv = crystal::norm(v);
    double cell = std::max(lu, lv) / 420.0;
    if ((lu / cell) * (lv / cell) > 250'000.0) cell = std::sqrt(lu * lv / 250'000.0);
    const int nu = std::max(1, static_cast<int>(std::ceil(lu / cell)));
    const int nv = std::max(1, static_cast<int>(std::ceil(lv / cell)));
    std::vector<std::int64_t> k(static_cast<std::size_t>(nu) * static_cast<std::size_t>(nv));
    for (int j = 0; j < nv; ++j) {
      if (cancelled()) return;
      for (int i = 0; i < nu; ++i) {
        k[static_cast<std::size_t>(j) * nu + i] = key_at(origin + u * ((i + 0.5) / nu) + v * ((j + 0.5) / nv));
      }
    }
    const auto at = [&](int i, int j) { return k[static_cast<std::size_t>(j) * nu + i]; };
    const auto point = [&](double fu, double fv) { return origin + u * fu + v * fv; };
    const auto vertex = [](const Vec3 &p, std::uint32_t col) {
      return SceneVertex{static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z), col};
    };
    std::vector<SceneVertex> tri, lin;
    for (int j = 0; j < nv; ++j) {
      for (int i0 = 0; i0 < nu;) {
        int i1 = i0 + 1;
        while (i1 < nu && at(i1, j) == at(i0, j)) ++i1;
        if (at(i0, j) >= 0) {
          const std::uint32_t col = color_of(at(i0, j), alpha);
          const Vec3 a = point(double(i0) / nu, double(j) / nv), b = point(double(i1) / nu, double(j) / nv);
          const Vec3 cc = point(double(i1) / nu, double(j + 1) / nv), d = point(double(i0) / nu, double(j + 1) / nv);
          for (const Vec3 &p : {a, b, cc, a, cc, d}) tri.push_back(vertex(p, col));
        }
        i0 = i1;
      }
    }
    // Boundaries: between rows, then between columns (runs merged).
    const std::uint32_t lc = rgba(line_, 120);
    for (int j = 0; j + 1 < nv; ++j) {
      for (int i0 = 0; i0 < nu;) {
        if (at(i0, j) == at(i0, j + 1)) {
          ++i0;
          continue;
        }
        int i1 = i0 + 1;
        while (i1 < nu && at(i1, j) != at(i1, j + 1)) ++i1;
        lin.push_back(vertex(point(double(i0) / nu, double(j + 1) / nv), lc));
        lin.push_back(vertex(point(double(i1) / nu, double(j + 1) / nv), lc));
        i0 = i1;
      }
    }
    for (int i = 0; i + 1 < nu; ++i) {
      for (int j0 = 0; j0 < nv;) {
        if (at(i, j0) == at(i + 1, j0)) {
          ++j0;
          continue;
        }
        int j1 = j0 + 1;
        while (j1 < nv && at(i, j1) != at(i + 1, j1)) ++j1;
        lin.push_back(vertex(point(double(i + 1) / nu, double(j0) / nv), lc));
        lin.push_back(vertex(point(double(i + 1) / nu, double(j1) / nv), lc));
        j0 = j1;
      }
    }
    for (const int view : views) {
      auto &t = s_.triangles[static_cast<std::size_t>(view)];
      t.insert(t.end(), tri.begin(), tri.end());
      auto &l = s_.lines[static_cast<std::size_t>(view)];
      l.insert(l.end(), lin.begin(), lin.end());
    }
  }

  void block_edges() {
    const std::uint32_t lc = rgba(line_, 170);
    const double xs[2] = {0.0, s_.depth_A}, ys[2] = {s_.y0, s_.y1}, zs[2] = {s_.z0, s_.z1};
    auto &l = s_.lines[3];
    const auto add = [&](double ax, double ay, double az, double bx, double by, double bz) {
      l.push_back({float(ax), float(ay), float(az), lc});
      l.push_back({float(bx), float(by), float(bz), lc});
    };
    for (double x : xs) {
      add(x, ys[0], zs[0], x, ys[1], zs[0]);
      add(x, ys[1], zs[0], x, ys[1], zs[1]);
      add(x, ys[1], zs[1], x, ys[0], zs[1]);
      add(x, ys[0], zs[1], x, ys[0], zs[0]);
    }
    for (double y : ys) {
      for (double z : zs) add(xs[0], y, z, xs[1], y, z);
    }
  }

  void spot_outline() {
    const double w = common::in_units(s_.spot.width, common::angstrom);
    const double h = common::in_units(s_.spot.height, common::angstrom);
    const std::uint32_t col = rgba(SampleScene::spot_color(), 255);
    auto &o = s_.spot_lines;
    const auto add = [&](double y0, double z0, double y1, double z1) {
      o.push_back({0.0f, float(y0), float(z0), col});
      o.push_back({0.0f, float(y1), float(z1), col});
    };
    switch (s_.spot.shape) {
    case trim::BeamSpot::Shape::point: return;
    case trim::BeamSpot::Shape::rectangle:
      add(-w / 2, -h / 2, w / 2, -h / 2);
      add(w / 2, -h / 2, w / 2, h / 2);
      add(w / 2, h / 2, -w / 2, h / 2);
      add(-w / 2, h / 2, -w / 2, -h / 2);
      break;
    case trim::BeamSpot::Shape::circle:
    case trim::BeamSpot::Shape::gaussian: {
      // A circle of diameter w; the Gaussian's half-maximum ellipse (w x h).
      const double ry = 0.5 * w;
      const double rz = s_.spot.shape == trim::BeamSpot::Shape::circle ? 0.5 * w : 0.5 * h;
      const int n = 160;
      for (int k = 0; k < n; ++k) {
        const double a0 = 2.0 * pi * k / n, a1 = 2.0 * pi * (k + 1) / n;
        add(ry * std::cos(a0), rz * std::sin(a0), ry * std::cos(a1), rz * std::sin(a1));
      }
      break;
    }
    }
    s_.spot_drawn = true;
  }

  bool cancelled() const { return cancel_ != nullptr && cancel_->load(std::memory_order_relaxed); }

  const trim::TrimConfig &c_;
  QColor line_;
  const std::atomic<bool> *cancel_;
  SampleScene &s_;
  std::vector<LayerInfo> layers_;
  double polycrystal_period_ = 0.0;
  std::unordered_map<std::int64_t, std::uint32_t> colors_;
};

} // namespace

QColor ipf_color(const crystal::Vec3 &d) {
  // The standard triangle of a cubic crystal: the direction's components,
  // absolute and sorted (a <= b <= c); [001] = (0 0 1) red, [101] green,
  // [111] blue, mixed linearly between them.
  double a = std::abs(d.x), b = std::abs(d.y), c = std::abs(d.z);
  if (a > b) std::swap(a, b);
  if (b > c) std::swap(b, c);
  if (a > b) std::swap(a, b);
  const double n = std::sqrt(a * a + b * b + c * c);
  if (!(n > 0.0)) return QColor(128, 128, 128);
  a /= n, b /= n, c /= n;
  double r = c - b, g = (b - a) * std::sqrt(2.0), bl = a * std::sqrt(3.0);
  const double m = std::max({r, g, bl, 1e-12});
  r /= m, g /= m, bl /= m;
  // (A square root brightens the mixed colours, as EBSD maps do.)
  return QColor::fromRgbF(static_cast<float>(std::sqrt(r)), static_cast<float>(std::sqrt(g)),
                          static_cast<float>(std::sqrt(bl)));
}

QString SampleScene::describe(double x, double y, double z) const {
  if (x < 0.0 || x > depth_A || layers.empty()) return {};
  const Layer *l = &layers.back();
  std::size_t index = layers.size() - 1;
  for (std::size_t i = 0; i < layers.size(); ++i) {
    if (x <= layers[i].bottom_A) {
      l = &layers[i];
      index = i;
      break;
    }
  }
  QString text = QStringLiteral("%1 (layer %2)").arg(l->name).arg(index + 1);
  if (!l->crystal) return text + tr(", amorphous");
  crystal::Mat3 r = l->rotation;
  if (l->grains) {
    const int g = l->grains->at({x - l->top_A, y, z}).grain;
    if (g < 0) return text + tr(", no grain");
    r = l->grains->grain(g).rotation;
    text += tr(", grain %1").arg(g + 1);
  }
  const auto near = crystal::nearest_direction(l->lattice, r.transposed() * beam, 2);
  return text + tr(": beam %1 deg from [%2 %3 %4]")
                    .arg(common::in_units(near.off, common::degree), 0, 'f', 1)
                    .arg(near.uvw[0])
                    .arg(near.uvw[1])
                    .arg(near.uvw[2]);
}

std::shared_ptr<const SampleScene> build_sample_scene(const trim::TrimConfig &c, QColor line_color,
                                                      const std::atomic<bool> *cancel) {
  auto s = std::make_shared<SampleScene>();
  try {
    Builder(c, line_color, cancel, *s).run();
  } catch (const std::exception &e) {
    *s = SampleScene{};
    s->message = QString::fromUtf8(e.what());
  } catch (...) {
    *s = SampleScene{};
    s->message = QStringLiteral("unknown error");
  }
  return s;
}

} // namespace gui
} // namespace amcset
