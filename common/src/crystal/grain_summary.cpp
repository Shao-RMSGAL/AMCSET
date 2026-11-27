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

#include "crystal/grain_summary.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <limits>
#include <stdexcept>

namespace amcset {
namespace crystal {

using namespace common;

namespace {

double nearest_image(double d, double period) noexcept { return d - period * std::nearbyint(d / period); }

//! One point found on the boundary of grains a < b.
struct Contact {
  int a = 0;
  int b = 0;
  Vec3 point;
};

void add_contact(std::vector<Contact> &out, int ga, int gb, const Vec3 &point) {
  out.push_back({std::min(ga, gb), std::max(ga, gb), point});
}

//! Contacts -> one boundary per pair: the contact nearest the middle of the
//! pair's contacts (lateral distances to the nearest image).
std::vector<GrainBoundary> merge(std::vector<Contact> contacts, double py, double pz) {
  std::sort(contacts.begin(), contacts.end(), [](const Contact &u, const Contact &v) {
    if (u.a != v.a) return u.a < v.a;
    if (u.b != v.b) return u.b < v.b;
    // (A fixed order within a pair: the result does not depend on the order found.)
    if (u.point.x != v.point.x) return u.point.x < v.point.x;
    if (u.point.y != v.point.y) return u.point.y < v.point.y;
    return u.point.z < v.point.z;
  });
  std::vector<GrainBoundary> out;
  for (std::size_t first = 0; first < contacts.size();) {
    std::size_t last = first;
    while (last < contacts.size() && contacts[last].a == contacts[first].a && contacts[last].b == contacts[first].b) {
      ++last;
    }
    const Vec3 &p0 = contacts[first].point;
    const auto offset = [&](const Vec3 &p) {
      return Vec3{p.x - p0.x, nearest_image(p.y - p0.y, py), nearest_image(p.z - p0.z, pz)};
    };
    Vec3 mean;
    for (std::size_t k = first; k < last; ++k) mean += offset(contacts[k].point);
    mean = mean / static_cast<double>(last - first);
    std::size_t best = first;
    double best_d2 = std::numeric_limits<double>::infinity();
    for (std::size_t k = first; k < last; ++k) {
      const Vec3 d = offset(contacts[k].point) - mean;
      const double d2 = dot(d, d);
      if (d2 < best_d2) {
        best_d2 = d2;
        best = k;
      }
    }
    out.push_back({contacts[best].a, contacts[best].b, contacts[best].point, last - first});
    first = last;
  }
  return out;
}

void voxel_contacts(const GrainMap &map, double depth_A, std::vector<Contact> &contacts) {
  const auto n = map.voxel_grid();
  const double dx = map.voxel_depth_A() / n[0];
  const double dy = map.period_y_A() / n[1];
  const double dz = map.period_z_A() / n[2];
  for (int i = 0; i < n[0]; ++i) {
    for (int j = 0; j < n[1]; ++j) {
      for (int k = 0; k < n[2]; ++k) {
        const int g = map.voxel(i, j, k);
        if (g < 0) continue;
        const Vec3 centre{(i + 0.5) * dx, (j + 0.5) * dy, (k + 0.5) * dz};
        // The faces towards +x, +y and +z (laterally periodic).
        if (i + 1 < n[0] && (i + 1) * dx <= depth_A) {
          const int h = map.voxel(i + 1, j, k);
          if (h >= 0 && h != g) add_contact(contacts, g, h, centre + Vec3{0.5 * dx, 0, 0});
        }
        if (centre.x > depth_A) continue;
        if (const int h = map.voxel(i, j + 1, k); h >= 0 && h != g) add_contact(contacts, g, h, centre + Vec3{0, 0.5 * dy, 0});
        if (const int h = map.voxel(i, j, k + 1); h >= 0 && h != g) add_contact(contacts, g, h, centre + Vec3{0, 0, 0.5 * dz});
      }
    }
  }
}

//! (1) of grain_boundaries: points of equal power distance between nearby seeds.
void seed_pair_contacts(const GrainMap &map, double depth_A, std::vector<Contact> &contacts) {
  const auto &seeds = map.seeds();
  const double py = map.period_y_A();
  const double pz = map.period_z_A();
  const bool columnar = map.is_columnar();
  // Buckets of about two seeds; the pairs looked at are those in the same
  // or neighbouring buckets.
  double x_lo = std::numeric_limits<double>::infinity();
  double x_hi = -x_lo;
  for (const auto &s : seeds) {
    x_lo = std::min(x_lo, s.position_A.x);
    x_hi = std::max(x_hi, s.position_A.x);
  }
  const double count = static_cast<double>(seeds.size());
  const double lx = columnar ? 1.0 : std::max(x_hi - x_lo, 1e-6);
  const double h = columnar ? std::sqrt(2.0 * py * pz / count) : std::cbrt(2.0 * lx * py * pz / count);
  const int nb[3] = {columnar ? 1 : std::clamp(static_cast<int>(lx / h), 1, 1024),
                     std::clamp(static_cast<int>(py / h), 1, 4096), std::clamp(static_cast<int>(pz / h), 1, 4096)};
  const double bw[3] = {lx / nb[0], py / nb[1], pz / nb[2]};
  const auto index = [&](int ix, int iy, int iz) {
    return (static_cast<std::size_t>(ix) * static_cast<std::size_t>(nb[1]) + static_cast<std::size_t>(iy)) *
               static_cast<std::size_t>(nb[2]) +
           static_cast<std::size_t>(iz);
  };
  const auto bucket_of = [&](const Vec3 &p, int (&b)[3]) {
    b[0] = columnar ? 0 : std::clamp(static_cast<int>((p.x - x_lo) / bw[0]), 0, nb[0] - 1);
    b[1] = std::clamp(static_cast<int>(p.y / bw[1]), 0, nb[1] - 1);
    b[2] = std::clamp(static_cast<int>(p.z / bw[2]), 0, nb[2] - 1);
  };
  const std::size_t buckets = static_cast<std::size_t>(nb[0]) * nb[1] * nb[2];
  std::vector<std::uint32_t> start(buckets + 1, 0);
  std::vector<std::uint32_t> member(seeds.size());
  for (const auto &s : seeds) {
    int b[3];
    bucket_of(s.position_A, b);
    ++start[index(b[0], b[1], b[2]) + 1];
  }
  for (std::size_t i = 0; i < buckets; ++i) start[i + 1] += start[i];
  {
    std::vector<std::uint32_t> next(start.begin(), start.end() - 1);
    for (std::size_t i = 0; i < seeds.size(); ++i) {
      int b[3];
      bucket_of(seeds[i].position_A, b);
      member[next[index(b[0], b[1], b[2])]++] = static_cast<std::uint32_t>(i);
    }
  }
  std::vector<std::uint32_t> near;
  for (std::size_t i = 0; i < seeds.size(); ++i) {
    const auto &si = seeds[i];
    if (si.grain < 0) continue;
    int b[3];
    bucket_of(si.position_A, b);
    near.clear();
    for (int ox = -1; ox <= 1; ++ox) {
      const int ix = b[0] + ox;
      if (ix < 0 || ix >= nb[0]) continue;
      for (int oy = -1; oy <= 1; ++oy) {
        const int iy = ((b[1] + oy) % nb[1] + nb[1]) % nb[1];
        for (int oz = -1; oz <= 1; ++oz) {
          const int iz = ((b[2] + oz) % nb[2] + nb[2]) % nb[2];
          const std::size_t bucket = index(ix, iy, iz);
          for (std::uint32_t k = start[bucket]; k < start[bucket + 1]; ++k) {
            if (member[k] > i) near.push_back(member[k]);
          }
        }
      }
    }
    // (With fewer than three buckets along an axis a bucket is reached twice.)
    std::sort(near.begin(), near.end());
    near.erase(std::unique(near.begin(), near.end()), near.end());
    for (const std::uint32_t j : near) {
      const auto &sj = seeds[j];
      if (sj.grain < 0 || sj.grain == si.grain) continue;
      const Vec3 d{columnar ? 0.0 : sj.position_A.x - si.position_A.x,
                   nearest_image(sj.position_A.y - si.position_A.y, py),
                   nearest_image(sj.position_A.z - si.position_A.z, pz)};
      const double d2 = dot(d, d);
      if (!(d2 > 0.0)) continue;
      // Equal power distance |r - p|^2 - w on the line from seed i to seed j.
      const double t = 0.5 + (si.weight_A2 - sj.weight_A2) / (2.0 * d2);
      if (!(t > 0.0) || !(t < 1.0)) continue;
      const Vec3 point = si.position_A + d * t;
      if (!columnar && (point.x < 0.0 || point.x > depth_A)) continue;
      // Just on either side, the two grains: no other seed claims the point.
      const Vec3 step = d * 1e-6;
      if (map.at(point - step).grain != si.grain || map.at(point + step).grain != sj.grain) continue;
      add_contact(contacts, si.grain, sj.grain, point);
    }
  }
}

//! (2) of grain_boundaries: neighbouring points of a grid in different grains.
void grid_contacts(const GrainMap &map, double depth_A, std::size_t points, std::vector<Contact> &contacts) {
  const bool columnar = map.is_columnar();
  const double py = map.period_y_A();
  const double pz = map.period_z_A();
  const double budget = static_cast<double>(std::max<std::size_t>(points, 1000));
  int n[3] = {1, 1, 1};
  if (columnar) {
    const double h = std::sqrt(py * pz / budget);
    n[1] = std::max(1, static_cast<int>(std::ceil(py / h)));
    n[2] = std::max(1, static_cast<int>(std::ceil(pz / h)));
  } else {
    const double h = std::cbrt(depth_A * py * pz / budget);
    n[0] = std::max(1, static_cast<int>(std::ceil(depth_A / h)));
    n[1] = std::max(1, static_cast<int>(std::ceil(py / h)));
    n[2] = std::max(1, static_cast<int>(std::ceil(pz / h)));
  }
  const double hx = depth_A / n[0];
  const double hy = py / n[1];
  const double hz = pz / n[2];
  const auto point = [&](int i, int j, int k) {
    return Vec3{columnar ? 0.0 : (i + 0.5) * hx, (j + 0.5) * hy, (k + 0.5) * hz};
  };
  // Two slabs of grain numbers (this depth and the next).
  const std::size_t slab = static_cast<std::size_t>(n[1]) * static_cast<std::size_t>(n[2]);
  std::vector<int> here(slab), below(slab);
  const auto fill = [&](int i, std::vector<int> &g) {
    for (int j = 0; j < n[1]; ++j) {
      for (int k = 0; k < n[2]; ++k) g[static_cast<std::size_t>(j) * n[2] + k] = map.at(point(i, j, k)).grain;
    }
  };
  fill(0, here);
  for (int i = 0; i < n[0]; ++i) {
    if (i + 1 < n[0]) fill(i + 1, below);
    for (int j = 0; j < n[1]; ++j) {
      for (int k = 0; k < n[2]; ++k) {
        const int g = here[static_cast<std::size_t>(j) * n[2] + k];
        if (g < 0) continue;
        const Vec3 p = point(i, j, k);
        const int gy = here[static_cast<std::size_t>((j + 1) % n[1]) * n[2] + k];
        if (gy >= 0 && gy != g) add_contact(contacts, g, gy, p + Vec3{0, 0.5 * hy, 0});
        const int gz = here[static_cast<std::size_t>(j) * n[2] + (k + 1) % n[2]];
        if (gz >= 0 && gz != g) add_contact(contacts, g, gz, p + Vec3{0, 0, 0.5 * hz});
        if (i + 1 < n[0]) {
          const int gx = below[static_cast<std::size_t>(j) * n[2] + k];
          if (gx >= 0 && gx != g) add_contact(contacts, g, gx, p + Vec3{0.5 * hx, 0, 0});
        }
      }
    }
    std::swap(here, below);
  }
}

double percentile(const std::vector<double> &sorted, double q) {
  if (sorted.empty()) return 0.0;
  const double at = q * static_cast<double>(sorted.size() - 1);
  const auto lo = static_cast<std::size_t>(std::floor(at));
  const auto hi = std::min(lo + 1, sorted.size() - 1);
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (at - static_cast<double>(lo));
}

//! b S with the lattice rotation S that maximises trace(a^T b S): the
//! orientation equivalent to b nearest a.
Mat3 nearest_equivalent(const Mat3 &a, const Mat3 &b, const std::vector<Mat3> &symmetry) {
  const Mat3 d = a.transposed() * b;
  Mat3 best = Mat3::identity();
  double best_trace = -std::numeric_limits<double>::infinity();
  for (const auto &s : symmetry) {
    double t = 0.0;
    for (int i = 0; i < 3; ++i) t += dot(d[i], s.column(i));
    if (t > best_trace + 1e-12) {
      best_trace = t;
      best = s;
    }
  }
  return b * best;
}

} // namespace

std::vector<GrainBoundary> grain_boundaries(const GrainMap &map, double depth_A, std::size_t grid_points) {
  if (!(depth_A > 0.0)) throw std::invalid_argument("grain boundaries: the layer's depth must be positive");
  std::vector<Contact> contacts;
  if (map.is_voxels()) {
    voxel_contacts(map, depth_A, contacts);
  } else if (!map.seeds().empty()) {
    seed_pair_contacts(map, depth_A, contacts);
    if (grid_points == 0) {
      // Automatic: about 500 points per grain (some 8 across one), from
      // 50 000 to 2 000 000 points.
      grid_points = std::clamp<std::size_t>(500 * map.grain_count(), 50'000, 2'000'000);
    }
    grid_contacts(map, depth_A, grid_points, contacts);
  }
  return merge(std::move(contacts), map.period_y_A(), map.period_z_A());
}

GrainSummary summarize_grains(const GrainMap &map, const Lattice &lattice, const Mat3 &reference, double depth_A,
                              std::size_t max_samples) {
  if (!(depth_A > 0.0)) throw std::invalid_argument("grain summary: the layer's depth must be positive");
  GrainSummary s;
  s.grains = map.grain_count();
  s.columnar = map.is_columnar();
  s.period_y_A = map.period_y_A();
  s.period_z_A = map.period_z_A();
  const double py = s.period_y_A;
  const double pz = s.period_z_A;
  const double budget = static_cast<double>(
      std::max<std::size_t>(std::min(max_samples, std::max<std::size_t>(64 * map.grain_count(), 100'000)), 1000));

  // A grid of sample points (cell centres): over the surface for columns,
  // over the layer otherwise.
  int n[3] = {1, 1, 1};
  if (s.columnar) {
    const double h = std::sqrt(py * pz / budget);
    n[1] = std::max(1, static_cast<int>(std::ceil(py / h)));
    n[2] = std::max(1, static_cast<int>(std::ceil(pz / h)));
  } else {
    const double h = std::cbrt(depth_A * py * pz / budget);
    n[0] = std::max(1, static_cast<int>(std::ceil(depth_A / h)));
    n[1] = std::max(1, static_cast<int>(std::ceil(py / h)));
    n[2] = std::max(1, static_cast<int>(std::ceil(pz / h)));
  }
  const double hx = depth_A / n[0];
  const double hy = py / n[1];
  const double hz = pz / n[2];
  s.sample_spacing_A = s.columnar ? std::max(hy, hz) : std::max({hx, hy, hz});
  const double cell = s.columnar ? hy * hz : hx * hy * hz;
  std::vector<std::size_t> hits(map.grain_count(), 0);
  std::size_t none = 0;
  for (int i = 0; i < n[0]; ++i) {
    for (int j = 0; j < n[1]; ++j) {
      for (int k = 0; k < n[2]; ++k) {
        const int g = map.at({s.columnar ? 0.0 : (i + 0.5) * hx, (j + 0.5) * hy, (k + 0.5) * hz}).grain;
        if (g < 0) ++none;
        else ++hits[static_cast<std::size_t>(g)];
        ++s.samples;
      }
    }
  }
  s.ungrained_fraction = s.samples > 0 ? static_cast<double>(none) / static_cast<double>(s.samples) : 0.0;
  s.coarse = s.samples < 32 * map.grain_count();

  // Sizes and texture of the grains met.
  const auto symmetry = lattice_rotations(lattice);
  const Vec3 normal_ref = reference[0]; // the crystal direction along the depth (row 0 of crystal -> laboratory)
  std::vector<double> diameters;
  double weight = 0.0, mis_sum = 0.0, mis2_sum = 0.0, tilt_sum = 0.0;
  for (std::size_t g = 0; g < hits.size(); ++g) {
    if (hits[g] == 0) continue;
    ++s.grains_in_layer;
    const double size = static_cast<double>(hits[g]) * cell;
    diameters.push_back(s.columnar ? std::sqrt(4.0 * size / pi_v) : std::cbrt(6.0 * size / pi_v));
    const Mat3 &r = map.grain(static_cast<int>(g)).rotation;
    const double mis = misorientation_rad(r, reference, symmetry);
    double tilt = pi_v;
    for (const auto &sym : symmetry) {
      tilt = std::min(tilt, std::acos(std::clamp(std::abs(dot(normalized(r[0]), normalized(sym * normal_ref))), 0.0, 1.0)));
    }
    const double w = static_cast<double>(hits[g]);
    weight += w;
    mis_sum += w * mis;
    mis2_sum += w * mis * mis;
    tilt_sum += w * tilt;
  }
  std::sort(diameters.begin(), diameters.end());
  if (!diameters.empty()) {
    double sum = 0.0;
    for (double d : diameters) sum += d;
    s.diameter_mean_A = sum / static_cast<double>(diameters.size());
    s.diameter_median_A = percentile(diameters, 0.5);
    s.diameter_p10_A = percentile(diameters, 0.1);
    s.diameter_p90_A = percentile(diameters, 0.9);
  }
  if (weight > 0.0) {
    const double to_deg = 180.0 / pi_v;
    s.misorientation_mean_deg = mis_sum / weight * to_deg;
    s.misorientation_rms_deg = std::sqrt(mis2_sum / weight) * to_deg;
    s.normal_tilt_mean_deg = tilt_sum / weight * to_deg;
  }

  // Boundaries (not for the largest maps: tens of seconds).
  if (map.grain_count() > max_grains_for_boundaries) return s;
  s.boundaries_counted = true;
  const auto boundaries = grain_boundaries(map, depth_A);
  s.boundaries = boundaries.size();
  std::size_t low = 0;
  double sum = 0.0;
  for (const auto &b : boundaries) {
    const double mis = misorientation_rad(map.grain(b.a).rotation, map.grain(b.b).rotation, symmetry) * 180.0 / pi_v;
    sum += mis;
    if (mis < 15.0) ++low;
  }
  if (!boundaries.empty()) {
    s.boundary_misorientation_mean_deg = sum / static_cast<double>(boundaries.size());
    s.low_angle_fraction = static_cast<double>(low) / static_cast<double>(boundaries.size());
  }
  return s;
}

std::string to_json(const GrainSummary &s) {
  const auto num = [](double v) {
    if (!std::isfinite(v)) return std::string("null");
    char buf[48];
    std::snprintf(buf, sizeof buf, "%.6g", v);
    return std::string(buf);
  };
  std::string j = "{\"grains\": " + std::to_string(s.grains) + ", \"grains_in_layer\": " +
                  std::to_string(s.grains_in_layer) + ", \"columnar\": " + (s.columnar ? "true" : "false") +
                  ", \"period_A\": [" + num(s.period_y_A) + ", " + num(s.period_z_A) + "], \"samples\": " +
                  std::to_string(s.samples) + ", \"sample_spacing_A\": " + num(s.sample_spacing_A) +
                  ", \"coarse\": " + (s.coarse ? "true" : "false") +
                  ", \"ungrained_fraction\": " + num(s.ungrained_fraction) + ", \"diameter_A\": {\"mean\": " +
                  num(s.diameter_mean_A) + ", \"median\": " + num(s.diameter_median_A) + ", \"p10\": " +
                  num(s.diameter_p10_A) + ", \"p90\": " + num(s.diameter_p90_A) +
                  "}, \"misorientation_from_reference_deg\": {\"mean\": " + num(s.misorientation_mean_deg) +
                  ", \"rms\": " + num(s.misorientation_rms_deg) + "}, \"normal_tilt_mean_deg\": " +
                  num(s.normal_tilt_mean_deg);
  if (s.boundaries_counted) {
    j += ", \"boundaries\": " + std::to_string(s.boundaries) + ", \"boundary_misorientation_mean_deg\": " +
         num(s.boundary_misorientation_mean_deg) + ", \"low_angle_fraction\": " + num(s.low_angle_fraction);
  }
  j += ", \"registry_below_deg\": " + num(s.registry_below_deg) + ", \"registry_groups\": " +
       std::to_string(s.registry_groups) + ", \"registry_grains\": " + std::to_string(s.registry_grains) + "}";
  return j;
}

GrainMap put_in_registry(const GrainMap &map, const Lattice &lattice, double below_rad, double depth_A,
                         RegistryResult *result) {
  if (result) *result = {};
  if (!(below_rad > 0.0)) return map;
  const auto symmetry = lattice_rotations(lattice);
  const auto boundaries = grain_boundaries(map, depth_A);
  std::vector<Grain> grains = map.grains();
  const std::size_t n = grains.size();

  // The low-angle neighbours of every grain (by the original orientations;
  // the misorientation does not change when a grain takes an equivalent one).
  struct Edge {
    int other;
    std::size_t boundary;
  };
  std::vector<std::vector<Edge>> edges(n);
  for (std::size_t k = 0; k < boundaries.size(); ++k) {
    const auto &b = boundaries[k];
    if (misorientation_rad(grains[static_cast<std::size_t>(b.a)].rotation, grains[static_cast<std::size_t>(b.b)].rotation,
                           symmetry) < below_rad) {
      edges[static_cast<std::size_t>(b.a)].push_back({b.b, k});
      edges[static_cast<std::size_t>(b.b)].push_back({b.a, k});
    }
  }
  for (auto &e : edges) std::sort(e.begin(), e.end(), [](const Edge &u, const Edge &v) { return u.other < v.other; });

  const Mat3 cartesian = lattice.vectors_A().transposed(); // fractional -> crystal frame
  const Mat3 fractional = lattice.reciprocal_per_A();      // crystal frame -> fractional
  std::vector<char> done(n, 0);
  std::deque<int> queue;
  for (std::size_t root = 0; root < n; ++root) {
    if (done[root] || edges[root].empty()) continue;
    done[root] = 1;
    std::size_t members = 1;
    queue.assign(1, static_cast<int>(root));
    while (!queue.empty()) {
      const int p = queue.front();
      queue.pop_front();
      const Grain &parent = grains[static_cast<std::size_t>(p)];
      for (const auto &e : edges[static_cast<std::size_t>(p)]) {
        if (done[static_cast<std::size_t>(e.other)]) continue;
        auto &child = grains[static_cast<std::size_t>(e.other)];
        const Vec3 &point = boundaries[e.boundary].point_A;
        child.rotation = nearest_equivalent(parent.rotation, child.rotation, symmetry);
        // The parent's lattice point nearest the boundary point: the child's
        // lattice passes through it.
        Vec3 f = fractional * (parent.rotation.transposed() * (point - parent.origin_A));
        f = {std::nearbyint(f.x), std::nearbyint(f.y), std::nearbyint(f.z)};
        child.origin_A = parent.origin_A + parent.rotation * (cartesian * f);
        done[static_cast<std::size_t>(e.other)] = 1;
        queue.push_back(e.other);
        ++members;
      }
    }
    if (result && members > 1) {
      ++result->groups;
      result->grains += members;
    }
  }
  return map.with_grains(std::move(grains));
}

} // namespace crystal
} // namespace amcset
