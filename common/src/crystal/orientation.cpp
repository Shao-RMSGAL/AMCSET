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

#include "crystal/orientation.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <stdexcept>

namespace amcset {
namespace crystal {

using namespace common;

namespace {

void check_finite(const std::vector<double> &v, const std::string &where) {
  for (double x : v) {
    if (!std::isfinite(x)) throw std::invalid_argument(where + ": indices must be finite numbers");
  }
}

double angle_rad(const Vec3 &a, const Vec3 &b) {
  const double c = dot(a, b) / (norm(a) * norm(b));
  return std::acos(std::clamp(c, -1.0, 1.0));
}

} // namespace

Vec3 direction_indices(const Lattice &lattice, const std::vector<double> &uvw, const std::string &where) {
  check_finite(uvw, where);
  if (uvw.size() == 3) {
    const Vec3 d{uvw[0], uvw[1], uvw[2]};
    if (norm(d) == 0.0) throw std::invalid_argument(where + ": [0 0 0] is not a direction");
    return d;
  }
  if (uvw.size() == 4) {
    if (!lattice.is_hexagonal()) {
      throw std::invalid_argument(where + ": four-index (Miller-Bravais) directions need a hexagonal cell");
    }
    const double u = uvw[0], v = uvw[1], t = uvw[2], w = uvw[3];
    if (std::abs(t + u + v) > 1e-9 * std::max({1.0, std::abs(u), std::abs(v)})) {
      throw std::invalid_argument(where + ": in [u v t w] t must be -(u + v)");
    }
    const Vec3 d{u - t, v - t, w};
    if (norm(d) == 0.0) throw std::invalid_argument(where + ": [0 0 0 0] is not a direction");
    return d;
  }
  throw std::invalid_argument(where + ": expected three indices [u v w] (or four [u v t w] for hexagonal cells)");
}

Vec3 plane_indices(const Lattice &lattice, const std::vector<double> &hkl, const std::string &where) {
  check_finite(hkl, where);
  if (hkl.size() == 3) {
    const Vec3 p{hkl[0], hkl[1], hkl[2]};
    if (norm(p) == 0.0) throw std::invalid_argument(where + ": (0 0 0) is not a plane");
    return p;
  }
  if (hkl.size() == 4) {
    if (!lattice.is_hexagonal()) {
      throw std::invalid_argument(where + ": four-index (Miller-Bravais) planes need a hexagonal cell");
    }
    const double h = hkl[0], k = hkl[1], i = hkl[2], l = hkl[3];
    if (std::abs(i + h + k) > 1e-9 * std::max({1.0, std::abs(h), std::abs(k)})) {
      throw std::invalid_argument(where + ": in (h k i l) i must be -(h + k)");
    }
    const Vec3 p{h, k, l};
    if (norm(p) == 0.0) throw std::invalid_argument(where + ": (0 0 0 0) is not a plane");
    return p;
  }
  throw std::invalid_argument(where + ": expected three indices (h k l) (or four (h k i l) for hexagonal cells)");
}

Vec3 direction_vector(const Lattice &lattice, const Vec3 &uvw) { return lattice.to_cartesian_A(uvw); }

Vec3 plane_normal(const Lattice &lattice, const Vec3 &hkl) {
  // h a* + k b* + l c*: the rows of the reciprocal matrix weighted by hkl.
  return lattice.reciprocal_per_A().transposed() * hkl;
}

Mat3 rotation(const Lattice &lattice, const Orientation &o, const std::string &where) {
  Mat3 r0;
  if (o.kind == OrientationKind::miller) {
    const Vec3 x = normalized(plane_normal(lattice, plane_indices(lattice, o.surface, where + ".surface")));
    Vec3 ref;
    if (!o.reference.empty()) {
      ref = direction_vector(lattice, direction_indices(lattice, o.reference, where + ".reference"));
    } else {
      ref = lattice.vectors_A()[0];
      const double ang = angle_rad(ref, x);
      if (ang < 10.0 * pi_v / 180.0 || ang > 170.0 * pi_v / 180.0) ref = lattice.vectors_A()[1];
    }
    const Vec3 in_plane = ref - x * dot(ref, x);
    if (norm(in_plane) < 1e-6 * norm(ref)) {
      throw std::invalid_argument(where + ".reference: parallel to the surface normal; it must have a component "
                                          "in the surface");
    }
    const Vec3 y = normalized(in_plane);
    r0 = Mat3::from_rows(x, y, cross(x, y));
  } else {
    for (double a : o.euler_deg) {
      if (!std::isfinite(a)) throw std::invalid_argument(where + ".euler: angles must be finite");
    }
    const double d2r = pi_v / 180.0;
    const double c1 = std::cos(o.euler_deg[0] * d2r), s1 = std::sin(o.euler_deg[0] * d2r);
    const double c = std::cos(o.euler_deg[1] * d2r), s = std::sin(o.euler_deg[1] * d2r);
    const double c2 = std::cos(o.euler_deg[2] * d2r), s2 = std::sin(o.euler_deg[2] * d2r);
    // Bunge's g (sample -> crystal: v_c = g v_s).
    const Mat3 g = Mat3::from_rows({c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s},
                                   {-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s},
                                   {s1 * s, -c1 * s, c});
    // Laboratory from sample: x = -zs (into the target), y = xs, z = -ys.
    const Mat3 m = Mat3::from_rows({0, 0, -1}, {1, 0, 0}, {0, -1, 0});
    r0 = m * g.transposed();
  }
  const double tw = in_units(o.twist, radian);
  if (!std::isfinite(tw)) throw std::invalid_argument(where + ".twist: not a finite angle");
  return rotation_about({1, 0, 0}, tw) * r0;
}

TiltTwist tilt_twist_for_beam(const Lattice &lattice, const Orientation &o, const Vec3 &beam,
                              const std::string &where) {
  if (!(norm(beam) > 0.0)) throw std::invalid_argument(where + ": the beam direction is zero");
  Orientation untwisted = o;
  untwisted.twist = 0.0 * degree;
  const Vec3 e = normalized(rotation(lattice, untwisted) * beam);
  if (!(e.x > 1e-9)) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.4g", std::acos(std::clamp(e.x, -1.0, 1.0)) * 180.0 / pi_v);
    throw std::invalid_argument(where + ": the beam is " + buf +
                                " deg from the inward surface normal; it must be below 90 deg to enter the target");
  }
  TiltTwist r;
  r.tilt = std::acos(std::clamp(e.x, -1.0, 1.0)) * radian;
  const double rho = std::hypot(e.y, e.z);
  // (+ 0.0: no negative zero when the beam already lies in the x-y plane.)
  r.twist = rho < 1e-12 ? o.twist : (std::atan2(-e.z, e.y) + 0.0) * radian;
  return r;
}

Vec3 beam_in_crystal(const Lattice &lattice, const Orientation &o, angle_quantity tilt) {
  const double t = in_units(tilt, radian);
  return rotation(lattice, o).transposed() * Vec3{std::cos(t), std::sin(t), 0.0};
}

namespace {

//! nearest_direction() (\p plane false) and nearest_plane() (true).
NearestDirection nearest_indices(const Lattice &lattice, const Vec3 &v, int max_index, bool plane) {
  NearestDirection best;
  double best_angle = 10.0;
  // Larger index sets only win when strictly closer, so the lowest indices
  // of equivalent directions are kept.
  for (int n = 1; n <= max_index; ++n) {
    for (int u = -n; u <= n; ++u) {
      for (int w2 = -n; w2 <= n; ++w2) {
        for (int w = -n; w <= n; ++w) {
          if (std::max({std::abs(u), std::abs(w2), std::abs(w)}) != n) continue;
          if (std::gcd(std::gcd(std::abs(u), std::abs(w2)), std::abs(w)) != 1) continue;
          const Vec3 indices{double(u), double(w2), double(w)};
          const double a = angle_rad(plane ? plane_normal(lattice, indices) : direction_vector(lattice, indices), v);
          if (a < best_angle - 1e-9) {
            best_angle = a;
            best.uvw = {u, w2, w};
          }
        }
      }
    }
  }
  best.off = best_angle * radian;
  return best;
}

//! An angle in degrees in [0, 360).
double wrapped_deg(double rad) {
  double d = std::fmod(rad * 180.0 / pi_v, 360.0);
  if (d < 0.0) d += 360.0;
  return d >= 360.0 ? 0.0 : d + 0.0;
}

} // namespace

NearestDirection nearest_direction(const Lattice &lattice, const Vec3 &v, int max_index) {
  return nearest_indices(lattice, v, max_index, false);
}

NearestDirection nearest_plane(const Lattice &lattice, const Vec3 &v, int max_index) {
  return nearest_indices(lattice, v, max_index, true);
}

Orientation orientation_for_frame(const Lattice &lattice, const Vec3 &normal, const Vec3 &lateral, int max_index,
                                  angle_quantity tolerance, const std::string &where) {
  if (!(norm(normal) > 0.0) || !std::isfinite(norm(normal))) {
    throw std::invalid_argument(where + ": the direction of the surface normal is zero");
  }
  const auto in_surface = [&](const Vec3 &x) {
    const Vec3 p = lateral - x * dot(lateral, x);
    if (!(norm(lateral) > 0.0) || !std::isfinite(norm(lateral)) || norm(p) < 1e-6 * norm(lateral)) {
      throw std::invalid_argument(where + ": the lateral direction is parallel to the surface normal; it must "
                                          "have a component in the surface");
    }
    return normalized(p);
  };
  const double tol = in_units(tolerance, radian);
  Orientation o;
  const auto plane = nearest_plane(lattice, normal, max_index);
  if (in_units(plane.off, radian) <= tol) {
    o.kind = OrientationKind::miller;
    o.surface = {double(plane.uvw[0]), double(plane.uvw[1]), double(plane.uvw[2])};
    const Vec3 x = normalized(plane_normal(lattice, Vec3{o.surface[0], o.surface[1], o.surface[2]}));
    const Vec3 y = in_surface(x);
    const auto ref = nearest_direction(lattice, y, max_index);
    if (in_units(ref.off, radian) <= tol) {
      o.reference = {double(ref.uvw[0]), double(ref.uvw[1]), double(ref.uvw[2])};
    } else {
      // The automatic reference, and the twist from its lateral axes (y0,
      // z0) to y: rotation() turns them by Rx(twist), y = cos y0 - sin z0.
      const Mat3 r0 = rotation(lattice, o, where);
      o.twist = (std::atan2(-dot(y, r0[2]), dot(y, r0[1])) + 0.0) * radian;
    }
    return o;
  }
  // Bunge Euler angles of R = m g^T (rotation()): g = R^T m.
  const Vec3 x = normalized(normal);
  const Vec3 y = in_surface(x);
  const Mat3 r = Mat3::from_rows(x, y, cross(x, y));
  const Mat3 m = Mat3::from_rows({0, 0, -1}, {1, 0, 0}, {0, -1, 0});
  const Mat3 g = r.transposed() * m;
  const double big_phi = std::acos(std::clamp(g[2][2], -1.0, 1.0));
  double phi1 = 0.0;
  double phi2 = 0.0;
  if (std::hypot(g[2][0], g[2][1]) > 1e-9) {
    phi1 = std::atan2(g[2][0], -g[2][1]);
    phi2 = std::atan2(g[0][2], g[1][2]);
  } else {
    // Phi = 0 or 180 deg: only phi1 + phi2 (or phi1 - phi2) is defined.
    phi1 = std::atan2(g[0][1], g[0][0]);
  }
  o.kind = OrientationKind::euler;
  o.euler_deg = {wrapped_deg(phi1), big_phi * 180.0 / pi_v, wrapped_deg(phi2)};
  return o;
}

std::string format_indices(const std::vector<double> &v, char open, char close) {
  std::string s(1, open);
  for (std::size_t i = 0; i < v.size(); ++i) {
    char buf[32];
    if (v[i] == std::round(v[i]) && std::abs(v[i]) < 1e9) {
      std::snprintf(buf, sizeof buf, "%lld", static_cast<long long>(std::llround(v[i])));
    } else {
      std::snprintf(buf, sizeof buf, "%.6g", v[i]);
    }
    s += (i ? " " : "") + std::string(buf);
  }
  return s + close;
}

} // namespace crystal
} // namespace amcset
