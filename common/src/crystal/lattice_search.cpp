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

#include "crystal/lattice_search.h"

#include "amcset_units.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace amcset {
namespace crystal {

namespace {

//! SplitMix64's output function (Steele, Lea and Flood 2014): a bijective mixer.
constexpr std::uint64_t mix64(std::uint64_t z) noexcept {
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

//! Next number of the SplitMix64 stream \p state, uniform in [0, 1).
inline double next_uniform(std::uint64_t &state) noexcept {
  state += 0x9E3779B97F4A7C15ull;
  return static_cast<double>(mix64(state) >> 11) * (1.0 / 9007199254740992.0); // 2^-53
}

//! Three standard normal deviates (Box-Muller), each cut at +-cut.
inline Vec3 next_gaussians(std::uint64_t &state, double cut) noexcept {
  const double two_pi = 2.0 * common::pi_v;
  // 1 - u is in (0, 1]: the logarithm is finite.
  const double r1 = std::sqrt(-2.0 * std::log(1.0 - next_uniform(state)));
  const double a1 = two_pi * next_uniform(state);
  const double r2 = std::sqrt(-2.0 * std::log(1.0 - next_uniform(state)));
  const double a2 = two_pi * next_uniform(state);
  return {std::clamp(r1 * std::cos(a1), -cut, cut), std::clamp(r1 * std::sin(a1), -cut, cut),
          std::clamp(r2 * std::cos(a2), -cut, cut)};
}

} // namespace

std::uint64_t site_hash(const SiteId &id, std::uint64_t seed) noexcept {
  std::uint64_t h = mix64(seed + 0x9E3779B97F4A7C15ull);
  for (const std::uint64_t v : {static_cast<std::uint64_t>(static_cast<std::uint32_t>(id.i)),
                                static_cast<std::uint64_t>(static_cast<std::uint32_t>(id.j)),
                                static_cast<std::uint64_t>(static_cast<std::uint32_t>(id.k)),
                                static_cast<std::uint64_t>(id.site)}) {
    h = mix64(h ^ (v + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2)));
  }
  return h;
}

LatticeSearch::LatticeSearch(const UnitCell &cell, const Mat3 &rotation, const std::map<unsigned, double> &rms_A,
                             double max_impact_A, double reach_A)
    : max_impact_(max_impact_A), reach_(reach_A) {
  if (!(max_impact_A > 0.0) || !std::isfinite(max_impact_A)) {
    throw std::invalid_argument("lattice search: the maximum impact parameter must be positive");
  }
  if (!(reach_A > 0.0) || !std::isfinite(reach_A)) {
    throw std::invalid_argument("lattice search: the reach must be positive");
  }
  // Lattice and reciprocal vectors in the search frame (a rotation keeps
  // their scalar products).
  for (int i = 0; i < 3; ++i) {
    vectors_[i] = rotation * cell.lattice.vectors_A()[i];
    reciprocal_[i] = rotation * cell.lattice.reciprocal_per_A()[i];
  }
  double sigma_max = 0.0;
  for (const auto &s : cell.sites) {
    SiteData sd;
    sd.fractional = s.position;
    double sum = 0.0;
    for (const auto &o : s.occupants) {
      sum += o.fraction;
      sd.z.push_back(o.z);
      sd.cumulative.push_back(sum);
      const auto it = rms_A.find(o.z);
      const double sigma = it != rms_A.end() ? it->second : 0.0;
      if (!(sigma >= 0.0) || !std::isfinite(sigma)) {
        throw std::invalid_argument("lattice search: thermal amplitudes must not be negative");
      }
      sd.sigma.push_back(sigma);
      sigma_max = std::max(sigma_max, sigma);
    }
    sites_.push_back(std::move(sd));
  }
  margin_ = displacement_cut * std::sqrt(3.0) * sigma_max;
  // A site within D of a point has fractional coordinates within D |a*_i|
  // of the point's along each axis (|a*_i| is the reciprocal of the spacing
  // of the lattice planes normal to a*_i): the cells a search must visit.
  const double d = std::hypot(reach_ + margin_, max_impact_ + margin_);
  for (int i = 0; i < 3; ++i) cell_reach_[i] = d * norm(reciprocal_[i]);
  const Mat3 to_cartesian = vectors_.transposed();
  for (std::size_t s = 0; s < sites_.size(); ++s) {
    if (sites_[s].z.empty()) continue; // an empty site never holds an atom
    occupied_.push_back(static_cast<std::uint32_t>(s));
    offset_.push_back(to_cartesian * sites_[s].fractional);
  }
}

std::size_t LatticeSearch::sites_per_search() const noexcept {
  std::size_t cells = 1;
  for (int i = 0; i < 3; ++i) cells *= static_cast<std::size_t>(std::floor(2.0 * cell_reach_[i])) + 2;
  return cells * occupied_.size();
}

bool LatticeSearch::occupant(const SiteData &sd, std::uint64_t &state, unsigned &z, double &sigma) const noexcept {
  // One number per site, drawn also for fully occupied sites: the stream
  // position of the thermal displacement does not depend on the occupancy.
  const double u = next_uniform(state);
  for (std::size_t k = 0; k < sd.z.size(); ++k) {
    if (u < sd.cumulative[k]) {
      z = sd.z[k];
      sigma = sd.sigma[k];
      return true;
    }
  }
  return false; // the rest of the site's probability: a vacancy
}

Vec3 LatticeSearch::site_position(const SiteId &id) const noexcept {
  const Vec3 f = sites_[id.site].fractional + Vec3{double(id.i), double(id.j), double(id.k)};
  return vectors_.transposed() * f;
}

bool LatticeSearch::atom_at(const SiteId &id, std::uint64_t seed, unsigned &z, Vec3 &position) const noexcept {
  std::uint64_t state = site_hash(id, seed);
  double sigma = 0.0;
  if (!occupant(sites_[id.site], state, z, sigma)) return false;
  position = site_position(id);
  if (sigma > 0.0) position += next_gaussians(state, displacement_cut) * sigma;
  return true;
}

bool LatticeSearch::nearest_atom(const Vec3 &r, std::uint64_t seed, unsigned z_wanted, double x_lo, double x_hi,
                                 SiteId &id, Vec3 &position) const noexcept {
  const Vec3 f = reciprocal_ * r;
  const double fi = std::floor(f.x), fj = std::floor(f.y), fk = std::floor(f.z);
  const auto ci = static_cast<std::int32_t>(fi);
  const auto cj = static_cast<std::int32_t>(fj);
  const auto ck = static_cast<std::int32_t>(fk);
  const Vec3 local = r - vectors_.transposed() * Vec3{fi, fj, fk};
  const double eps = 1e-7; // (as candidates())
  double best = std::numeric_limits<double>::infinity();
  bool found = false;
  for (int di = -2; di <= 2; ++di) {
    for (int dj = -2; dj <= 2; ++dj) {
      for (int dk = -2; dk <= 2; ++dk) {
        const Vec3 cell = vectors_[0] * double(di) + vectors_[1] * double(dj) + vectors_[2] * double(dk) - local;
        for (std::size_t s = 0; s < occupied_.size(); ++s) {
          const Vec3 d = cell + offset_[s];
          const double site_x = r.x + d.x;
          if (site_x < x_lo - eps || site_x >= x_hi - eps) continue;
          const SiteId sid{ci + di, cj + dj, ck + dk, occupied_[s]};
          std::uint64_t state = site_hash(sid, seed);
          unsigned z = 0;
          double sigma = 0.0;
          if (!occupant(sites_[occupied_[s]], state, z, sigma)) continue;
          if (z_wanted != 0 && z != z_wanted) continue;
          Vec3 dd = d;
          if (sigma > 0.0) dd += next_gaussians(state, displacement_cut) * sigma;
          const double d2 = dot(dd, dd);
          if (d2 < best) {
            best = d2;
            id = sid;
            position = r + dd;
            found = true;
          }
        }
      }
    }
  }
  return found;
}

void LatticeSearch::candidates(const Vec3 &r, const Vec3 &dir, std::uint64_t seed, double x_lo, double x_hi,
                               std::vector<Candidate> &out) const {
  out.clear();
  const Vec3 f = reciprocal_ * r;
  const double fi = std::floor(f.x), fj = std::floor(f.y), fk = std::floor(f.z);
  const auto ci = static_cast<std::int32_t>(fi);
  const auto cj = static_cast<std::int32_t>(fj);
  const auto ck = static_cast<std::int32_t>(fk);
  const Vec3 local = r - vectors_.transposed() * Vec3{fi, fj, fk}; // r relative to its cell's origin
  const double xi_hi = reach_ + margin_;
  const double coarse2 = (max_impact_ + margin_) * (max_impact_ + margin_);
  const double impact2 = max_impact_ * max_impact_;
  const double eps = 1e-7; // Angstrom: sites on a boundary plane belong to the layer below it
  // The cells that can hold a site within reach (a site of cell d has
  // fractional coordinates in [d, d + 1) relative to r's cell).
  const double frac[3] = {f.x - fi, f.y - fj, f.z - fk};
  int lo[3], hi[3];
  for (int i = 0; i < 3; ++i) {
    lo[i] = static_cast<int>(std::floor(frac[i] - cell_reach_[i]));
    hi[i] = static_cast<int>(std::floor(frac[i] + cell_reach_[i]));
  }
  const std::size_t n_sites = occupied_.size();
  for (int di = lo[0]; di <= hi[0]; ++di) {
    for (int dj = lo[1]; dj <= hi[1]; ++dj) {
      for (int dk = lo[2]; dk <= hi[2]; ++dk) {
        // The cell's origin relative to r.
        const Vec3 cell = vectors_[0] * double(di) + vectors_[1] * double(dj) + vectors_[2] * double(dk) - local;
        for (std::size_t s = 0; s < n_sites; ++s) {
          const Vec3 d = cell + offset_[s]; // the lattice site relative to r
          const double xi0 = dot(d, dir);
          if (xi0 < -margin_ || xi0 > xi_hi) continue;
          if (dot(d, d) - xi0 * xi0 > coarse2) continue;
          const double site_x = r.x + d.x;
          if (site_x < x_lo - eps || site_x >= x_hi - eps) continue;
          const SiteId id{ci + di, cj + dj, ck + dk, occupied_[s]};
          std::uint64_t state = site_hash(id, seed);
          unsigned z = 0;
          double sigma = 0.0;
          if (!occupant(sites_[occupied_[s]], state, z, sigma)) continue;
          Vec3 dd = d;
          if (sigma > 0.0) dd += next_gaussians(state, displacement_cut) * sigma;
          const double xi = dot(dd, dir);
          if (!(xi > 0.0) || xi > reach_) continue;
          const double s2 = dot(dd, dd) - xi * xi;
          if (s2 > impact2) continue;
          out.push_back({id, z, r + dd, xi, std::sqrt(std::max(0.0, s2))});
        }
      }
    }
  }
  std::sort(out.begin(), out.end(), [](const Candidate &a, const Candidate &b) {
    if (a.xi != b.xi) return a.xi < b.xi;
    // A total order (symmetric sites have equal xi): results must not
    // depend on the order of the template.
    if (a.id.i != b.id.i) return a.id.i < b.id.i;
    if (a.id.j != b.id.j) return a.id.j < b.id.j;
    if (a.id.k != b.id.k) return a.id.k < b.id.k;
    return a.id.site < b.id.site;
  });
}

double largest_hole_radius_A(const UnitCell &cell) {
  std::vector<Vec3> atoms;
  for (const auto &s : cell.sites) {
    if (s.kind == SiteKind::lattice && !s.occupants.empty()) atoms.push_back(s.position);
  }
  if (atoms.empty()) throw std::invalid_argument("largest hole: the cell has no occupied lattice site");
  // The nearest image of an atom to a point of the cell lies in the
  // neighbouring cells; strongly skewed cells get one more shell.
  const double deg = common::pi_v / 180.0;
  bool skewed = false;
  for (const auto &angle : {cell.lattice.alpha(), cell.lattice.beta(), cell.lattice.gamma()}) {
    const double a = common::in_units(angle, common::radian);
    skewed = skewed || a < 60.0 * deg - 1e-9 || a > 120.0 * deg + 1e-9;
  }
  const int n[3] = {skewed ? 2 : 1, skewed ? 2 : 1, skewed ? 2 : 1};
  const auto nearest = [&](const Vec3 &f) {
    double best = std::numeric_limits<double>::max();
    for (const auto &a : atoms) {
      for (int i = -n[0]; i <= n[0]; ++i) {
        for (int j = -n[1]; j <= n[1]; ++j) {
          for (int k = -n[2]; k <= n[2]; ++k) {
            const Vec3 d = cell.lattice.to_cartesian_A(a + Vec3{double(i), double(j), double(k)} - f);
            best = std::min(best, dot(d, d));
          }
        }
      }
    }
    return std::sqrt(best);
  };
  // Coarse grid over the cell, then a shrinking grid around the best point.
  const int grid = 16;
  Vec3 best_f;
  double best = -1.0;
  for (int i = 0; i < grid; ++i) {
    for (int j = 0; j < grid; ++j) {
      for (int k = 0; k < grid; ++k) {
        const Vec3 f{(i + 0.5) / grid, (j + 0.5) / grid, (k + 0.5) / grid};
        const double d = nearest(f);
        if (d > best) {
          best = d;
          best_f = f;
        }
      }
    }
  }
  double half = 1.0 / grid;
  for (int level = 0; level < 12; ++level) {
    Vec3 centre = best_f;
    for (int i = -2; i <= 2; ++i) {
      for (int j = -2; j <= 2; ++j) {
        for (int k = -2; k <= 2; ++k) {
          const Vec3 f = centre + Vec3{i * half / 2, j * half / 2, k * half / 2};
          const double d = nearest(f);
          if (d > best) {
            best = d;
            best_f = f;
          }
        }
      }
    }
    half *= 0.5;
  }
  return best;
}

} // namespace crystal
} // namespace amcset
