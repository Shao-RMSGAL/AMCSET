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

#include "crystal/grain_io.h"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <deque>
#include <fstream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace amcset {
namespace crystal {

using namespace common;

namespace {

constexpr double deg = pi_v / 180.0;
constexpr double micrometre_A = 1.0e4;

struct Token {
  std::string_view text;
  std::size_t line = 0;
};

std::string at_line(const std::string &source, std::size_t line) { return source + ", line " + std::to_string(line); }

bool separator(char c) { return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r'; }

//! The lines of \p text (without their ends), numbered from 1.
std::vector<Token> lines_of(std::string_view text) {
  std::vector<Token> lines;
  std::size_t start = 0;
  std::size_t number = 1;
  while (start <= text.size()) {
    std::size_t end = text.find('\n', start);
    if (end == std::string_view::npos) end = text.size();
    lines.push_back({text.substr(start, end - start), number++});
    start = end + 1;
  }
  return lines;
}

//! The words of a line, up to a '#'.
void words_of(std::string_view line, std::vector<std::string_view> &out, bool comments = true) {
  out.clear();
  std::size_t i = 0;
  while (i < line.size()) {
    while (i < line.size() && separator(line[i])) ++i;
    if (i >= line.size() || (comments && line[i] == '#')) break;
    const std::size_t s = i;
    while (i < line.size() && !separator(line[i]) && !(comments && line[i] == '#')) ++i;
    out.push_back(line.substr(s, i - s));
  }
}

bool to_double(std::string_view t, double &v) {
  if (!t.empty() && t.front() == '+') t.remove_prefix(1);
  if (t.empty()) return false;
  const auto r = std::from_chars(t.data(), t.data() + t.size(), v);
  return r.ec == std::errc() && r.ptr == t.data() + t.size() && std::isfinite(v);
}

double number(std::string_view t, const std::string &where) {
  double v = 0.0;
  if (!to_double(t, v)) throw std::invalid_argument(where + ": '" + std::string(t) + "' is not a number");
  return v;
}

std::string lower(std::string_view s) {
  std::string out(s);
  for (auto &c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return out;
}

double unit_A(std::string_view unit, const std::string &where) {
  const std::string u = lower(unit);
  if (u == "a" || u == "angstrom") return 1.0;
  if (u == "nm") return 10.0;
  if (u == "um" || u == "micron" || u == "micrometer" || u == "micrometre") return micrometre_A;
  throw std::invalid_argument(where + ": unknown unit '" + std::string(unit) + "' (A, nm, um)");
}

//! Bunge angles (radians) of a crystal relative to the sample -> its rotation crystal -> laboratory.
Mat3 bunge_rotation(double phi1, double Phi, double phi2) {
  const double c1 = std::cos(phi1), s1 = std::sin(phi1);
  const double c = std::cos(Phi), s = std::sin(Phi);
  const double c2 = std::cos(phi2), s2 = std::sin(phi2);
  // Bunge's g (sample -> crystal: v_c = g v_s), as in crystal::rotation().
  const Mat3 g = Mat3::from_rows({c1 * c2 - s1 * s2 * c, s1 * c2 + c1 * s2 * c, s2 * s},
                                 {-c1 * s2 - s1 * c2 * c, -s1 * s2 + c1 * c2 * c, c2 * s},
                                 {s1 * s, -c1 * s, c});
  return sample_to_laboratory() * g.transposed();
}

//! A point of the unit cell that follows from the grain's number: the grains' lattices are not in registry.
Vec3 lattice_offset(std::size_t grain, const UnitCell &cell) {
  std::uint64_t z = (static_cast<std::uint64_t>(grain) + 1) * 0x9E3779B97F4A7C15ull;
  const auto next = [&z] {
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) * (1.0 / 9007199254740992.0);
  };
  const double u = next();
  const double v = next();
  const double w = next();
  return cell.lattice.to_cartesian_A({u, v, w});
}

double period_or(length_quantity given, double from_file) {
  const double g = in_units(given, angstrom);
  return g > 0.0 ? g : from_file;
}

std::string read_text(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read the grain file " + path.string());
  std::ostringstream out;
  out << in.rdbuf();
  return out.str();
}

} // namespace

Polycrystal::Source grain_source_from_path(const std::filesystem::path &path) {
  const std::string ext = lower(path.extension().string());
  if (ext == ".tess") return Polycrystal::Source::neper;
  if (ext == ".ang" || ext == ".ctf") return Polycrystal::Source::ebsd;
  return Polycrystal::Source::grain_file;
}

std::shared_ptr<const GrainMap> read_grain_map(const Polycrystal &p, const std::filesystem::path &path,
                                               const UnitCell &cell) {
  const std::string text = read_text(path);
  const std::string name = path.filename().string();
  switch (p.source) {
  case Polycrystal::Source::grain_file: return std::make_shared<const GrainMap>(parse_grain_file(text, p, cell, name));
  case Polycrystal::Source::neper: return std::make_shared<const GrainMap>(parse_neper_tess(text, p, cell, name));
  case Polycrystal::Source::ebsd:
    return std::make_shared<const GrainMap>(parse_ebsd(text, lower(path.extension().string()) == ".ctf", p, cell, name));
  case Polycrystal::Source::voronoi: break;
  }
  throw std::invalid_argument("polycrystal: generated grains are not read from a file");
}

// ---------------------------------------------------------------------------
// AMCSET's grain file
// ---------------------------------------------------------------------------

GrainMap parse_grain_file(std::string_view text, const Polycrystal &p, const UnitCell &cell, const std::string &source) {
  double unit = 1.0;
  double period[2] = {0.0, 0.0};
  double depth = 0.0;
  bool columnar = false;
  std::array<int, 3> n{0, 0, 0};
  std::vector<Grain> grains;
  std::vector<GrainMap::Seed> seeds;
  std::vector<int> ids;
  bool in_ids = false;
  std::vector<std::string_view> w;
  for (const auto &line : lines_of(text)) {
    words_of(line.text, w);
    if (w.empty()) continue;
    const std::string where = at_line(source, line.line);
    if (in_ids) {
      for (const auto t : w) {
        const double v = number(t, where);
        if (v != std::floor(v)) throw std::invalid_argument(where + ": grain numbers are whole numbers");
        ids.push_back(static_cast<int>(v));
      }
      continue;
    }
    const std::string key = lower(w[0]);
    const auto need = [&](std::size_t count) {
      if (w.size() < count + 1) {
        throw std::invalid_argument(where + ": '" + key + "' needs " + std::to_string(count) + " values");
      }
    };
    if (key == "unit") {
      need(1);
      unit = unit_A(w[1], where);
    } else if (key == "period") {
      need(1);
      period[0] = number(w[1], where) * unit;
      period[1] = w.size() > 2 ? number(w[2], where) * unit : period[0];
    } else if (key == "depth") {
      need(1);
      depth = number(w[1], where) * unit;
    } else if (key == "columnar") {
      columnar = true;
    } else if (key == "grain") {
      need(6);
      GrainMap::Seed s;
      s.position_A = {number(w[1], where) * unit, number(w[2], where) * unit, number(w[3], where) * unit};
      if (w.size() > 7) {
        const double radius = number(w[7], where) * unit;
        s.weight_A2 = radius * radius;
      }
      s.grain = static_cast<int>(grains.size());
      Grain g;
      g.rotation = bunge_rotation(number(w[4], where) * deg, number(w[5], where) * deg, number(w[6], where) * deg);
      g.origin_A = s.position_A;
      grains.push_back(g);
      seeds.push_back(s);
    } else if (key == "voxels") {
      need(3);
      for (int i = 0; i < 3; ++i) {
        const double v = number(w[static_cast<std::size_t>(i) + 1], where);
        if (!(v >= 1.0) || v != std::floor(v) || v > 4096.0) {
          throw std::invalid_argument(where + ": voxels: three whole numbers from 1 to 4096");
        }
        n[static_cast<std::size_t>(i)] = static_cast<int>(v);
      }
    } else if (key == "orientation") {
      need(3);
      Grain g;
      g.rotation = bunge_rotation(number(w[1], where) * deg, number(w[2], where) * deg, number(w[3], where) * deg);
      g.origin_A = lattice_offset(grains.size(), cell);
      grains.push_back(g);
    } else if (key == "ids") {
      in_ids = true;
      for (std::size_t i = 1; i < w.size(); ++i) ids.push_back(static_cast<int>(number(w[i], where)));
    } else {
      throw std::invalid_argument(where + ": unknown entry '" + std::string(w[0]) +
                                  "' (unit, period, depth, columnar, grain, voxels, orientation, ids)");
    }
  }
  const double py = period_or(p.period_y, period[0]);
  const double pz = period_or(p.period_z, period[1]);
  if (!(py > 0.0) || !(pz > 0.0)) {
    throw std::invalid_argument(source + ": the lateral periods are missing (a line \"period y z\", or the "
                                         "polycrystal's period)");
  }
  if (n[0] > 0) {
    if (!seeds.empty()) throw std::invalid_argument(source + ": give seeds (grain ...) or voxels, not both");
    if (!(depth > 0.0)) throw std::invalid_argument(source + ": voxels need a line \"depth ...\"");
    return GrainMap::from_voxels(std::move(grains), std::move(ids), n, depth, py, pz);
  }
  if (seeds.empty()) throw std::invalid_argument(source + ": no grains (lines \"grain x y z phi1 Phi phi2\")");
  return GrainMap::from_seeds(std::move(grains), std::move(seeds), py, pz, columnar);
}

// ---------------------------------------------------------------------------
// Neper .tess
// ---------------------------------------------------------------------------

GrainMap parse_neper_tess(std::string_view text, const Polycrystal &p, const UnitCell &, const std::string &source) {
  std::vector<Token> tokens;
  std::vector<std::string_view> w;
  for (const auto &line : lines_of(text)) {
    words_of(line.text, w, false);
    for (const auto t : w) tokens.push_back({t, line.line});
  }
  const auto find = [&](std::string_view name, std::size_t from) {
    for (std::size_t i = from; i < tokens.size(); ++i) {
      if (tokens[i].text == name) return i;
    }
    return tokens.size();
  };
  const auto value = [&](std::size_t i, const char *what) {
    if (i >= tokens.size()) throw std::invalid_argument(source + ": the file ends in " + what);
    return number(tokens[i].text, at_line(source, tokens[i].line) + " (" + what + ")");
  };
  if (tokens.empty() || tokens.front().text != "***tess") {
    throw std::invalid_argument(source + ": not a Neper tessellation (it does not start with ***tess)");
  }
  const std::size_t cell_at = find("**cell", 0);
  if (cell_at == tokens.size()) throw std::invalid_argument(source + ": no **cell section");
  const double count = value(cell_at + 1, "the number of cells");
  if (!(count >= 1.0) || count != std::floor(count) || count > static_cast<double>(max_grains)) {
    throw std::invalid_argument(source + ": " + std::string(tokens[cell_at + 1].text) + " cells");
  }
  const std::size_t n = static_cast<std::size_t>(count);
  const std::size_t seed_at = find("*seed", cell_at);
  if (seed_at == tokens.size()) throw std::invalid_argument(source + ": no *seed section (the cells' seeds)");
  const std::size_t ori_at = find("*ori", cell_at);
  if (ori_at == tokens.size()) {
    throw std::invalid_argument(source + ": no *ori section (the cells' orientations)");
  }
  if (ori_at + 1 >= tokens.size()) throw std::invalid_argument(source + ": the file ends in *ori");

  // The orientation descriptor: name[:active|:passive].
  std::string descriptor = lower(tokens[ori_at + 1].text);
  bool passive = false;
  if (const auto colon = descriptor.find(':'); colon != std::string::npos) {
    const std::string convention = descriptor.substr(colon + 1);
    if (convention == "passive") passive = true;
    else if (convention != "active") {
      throw std::invalid_argument(at_line(source, tokens[ori_at + 1].line) + ": unknown convention '" + convention + "'");
    }
    descriptor.resize(colon);
  }
  enum class Kind { bunge, rodrigues, quaternion, axis_angle, matrix };
  Kind kind = Kind::bunge;
  std::size_t width = 3;
  if (descriptor == "euler-bunge" || descriptor == "e") {
    kind = Kind::bunge;
  } else if (descriptor == "rodrigues" || descriptor == "r") {
    kind = Kind::rodrigues;
  } else if (descriptor == "quaternion" || descriptor == "q") {
    kind = Kind::quaternion;
    width = 4;
  } else if (descriptor == "axis-angle" || descriptor == "rtheta") {
    kind = Kind::axis_angle;
    width = 4;
  } else if (descriptor == "rotmat" || descriptor == "g") {
    kind = Kind::matrix;
    width = 9;
  } else {
    throw std::invalid_argument(at_line(source, tokens[ori_at + 1].line) + ": orientations as '" + descriptor +
                                "' are not read (euler-bunge, rodrigues, quaternion, axis-angle, rotmat)");
  }

  // The box: that of the vertices (of the seeds, if there are none).
  Vec3 lo{1e300, 1e300, 1e300};
  Vec3 hi{-1e300, -1e300, -1e300};
  const auto extend = [&](const Vec3 &v) {
    for (int k = 0; k < 3; ++k) {
      lo[k] = std::min(lo[k], v[k]);
      hi[k] = std::max(hi[k], v[k]);
    }
  };
  std::vector<Vec3> position(n);
  std::vector<double> weight(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t at = seed_at + 1 + 5 * i;
    position[i] = {value(at + 1, "a seed"), value(at + 2, "a seed"), value(at + 3, "a seed")};
    weight[i] = value(at + 4, "a seed's weight");
  }
  if (const std::size_t vertex_at = find("**vertex", 0); vertex_at != tokens.size()) {
    const double vertices = value(vertex_at + 1, "the number of vertices");
    for (std::size_t i = 0; i < static_cast<std::size_t>(std::max(vertices, 0.0)); ++i) {
      const std::size_t at = vertex_at + 2 + 5 * i;
      extend({value(at + 1, "a vertex"), value(at + 2, "a vertex"), value(at + 3, "a vertex")});
    }
  }
  if (!(hi.x > lo.x) || !(hi.y > lo.y) || !(hi.z > lo.z)) {
    for (const auto &v : position) extend(v);
  }
  if (!(hi.x > lo.x) || !(hi.y > lo.y)) throw std::invalid_argument(source + ": the tessellation has no extent");

  const double scale = in_units(p.scale, angstrom);
  std::vector<Grain> grains(n);
  std::vector<GrainMap::Seed> seeds(n);
  for (std::size_t i = 0; i < n; ++i) {
    const std::size_t at = ori_at + 2 + width * i;
    double v[9] = {0};
    for (std::size_t k = 0; k < width; ++k) v[k] = value(at + k, "an orientation");
    Mat3 active; // crystal -> sample
    switch (kind) {
    case Kind::bunge: active = sample_to_laboratory().transposed() * bunge_rotation(v[0] * deg, v[1] * deg, v[2] * deg); break;
    case Kind::rodrigues: {
      const Vec3 r{v[0], v[1], v[2]};
      active = norm(r) > 0.0 ? rotation_about(r, 2.0 * std::atan(norm(r))) : Mat3::identity();
      break;
    }
    case Kind::quaternion: {
      const double q0 = v[0];
      const Vec3 q{v[1], v[2], v[3]};
      const double s = std::sqrt(q0 * q0 + dot(q, q));
      active = norm(q) > 0.0 && s > 0.0 ? rotation_about(q, 2.0 * std::atan2(norm(q), q0)) : Mat3::identity();
      break;
    }
    case Kind::axis_angle: {
      const Vec3 axis{v[0], v[1], v[2]};
      active = norm(axis) > 0.0 ? rotation_about(axis, v[3] * deg) : Mat3::identity();
      break;
    }
    case Kind::matrix: active = Mat3::from_rows({v[0], v[1], v[2]}, {v[3], v[4], v[5]}, {v[6], v[7], v[8]}); break;
    }
    if (passive) active = active.transposed();
    if (rotation_error(active) > 1e-3) {
      throw std::invalid_argument(at_line(source, tokens[std::min(at, tokens.size() - 1)].line) +
                                  ": the orientation of cell " + std::to_string(i + 1) + " is not a rotation");
    }
    grains[i].rotation = sample_to_laboratory() * active;
    // Sample (x, y, z) -> layer (z_top - z, x, -y), from the box's corner.
    auto &s = seeds[i];
    s.position_A = {(hi.z - position[i].z) * scale, (position[i].x - lo.x) * scale, (hi.y - position[i].y) * scale};
    s.weight_A2 = weight[i] * weight[i] * scale * scale;
    s.grain = static_cast<int>(i);
    grains[i].origin_A = s.position_A;
  }
  return GrainMap::from_seeds(std::move(grains), std::move(seeds), period_or(p.period_y, (hi.x - lo.x) * scale),
                              period_or(p.period_z, (hi.y - lo.y) * scale), false);
}

// ---------------------------------------------------------------------------
// EBSD maps
// ---------------------------------------------------------------------------

GrainMap parse_ebsd(std::string_view text, bool ctf, const Polycrystal &p, const UnitCell &cell, const std::string &source) {
  struct Point {
    double x = 0.0, y = 0.0; // micrometres
    Mat3 rotation;
    bool valid = false;
  };
  std::vector<Point> points;
  double step = 0.0;
  std::vector<std::string_view> w;
  int columns[6] = {-1, -1, -1, -1, -1, -1}; // ctf: Phase, X, Y, Euler1..3
  bool in_data = !ctf;
  for (const auto &line : lines_of(text)) {
    const std::string where = at_line(source, line.line);
    if (!ctf) {
      // TSL .ang: "# XSTEP: 0.5" in the header; phi1 PHI phi2 x y IQ CI Phase ...
      if (!line.text.empty() && line.text.front() == '#') {
        words_of(line.text.substr(1), w, false);
        if (w.size() >= 2 && lower(w[0]) == "xstep:") step = number(w[1], where);
        continue;
      }
      words_of(line.text, w, false);
      if (w.empty()) continue;
      if (w.size() < 5) throw std::invalid_argument(where + ": expected phi1 PHI phi2 x y ...");
      Point pt;
      const double e1 = number(w[0], where), e2 = number(w[1], where), e3 = number(w[2], where);
      pt.x = number(w[3], where);
      pt.y = number(w[4], where);
      const double ci = w.size() > 6 ? number(w[6], where) : 1.0;
      const int phase = w.size() > 7 ? static_cast<int>(number(w[7], where)) : 0;
      // Points that are not indexed: a negative confidence index, or the
      // angles 4 pi that TSL writes for them.
      const bool indexed = ci >= 0.0 && !(e1 > 12.5 && e2 > 12.5 && e3 > 12.5);
      pt.valid = indexed && (p.phase == 0 || phase == p.phase);
      if (pt.valid) pt.rotation = bunge_rotation(e1, e2, e3);
      points.push_back(pt);
      continue;
    }
    // Oxford .ctf: a header, then a line of column names, then the points.
    words_of(line.text, w, false);
    if (w.empty()) continue;
    if (!in_data) {
      if (lower(w[0]) == "xstep" && w.size() >= 2) step = number(w[1], where);
      if (lower(w[0]) == "phase" && w.size() >= 6) {
        static const char *const names[6] = {"phase", "x", "y", "euler1", "euler2", "euler3"};
        for (std::size_t c = 0; c < w.size(); ++c) {
          const std::string name = lower(w[c]);
          for (int k = 0; k < 6; ++k) {
            if (name == names[k]) columns[k] = static_cast<int>(c);
          }
        }
        for (int k = 0; k < 6; ++k) {
          if (columns[k] < 0) throw std::invalid_argument(where + ": no column " + std::string(names[k]));
        }
        in_data = true;
      }
      continue;
    }
    const auto column = [&](int k) {
      if (static_cast<std::size_t>(columns[k]) >= w.size()) throw std::invalid_argument(where + ": too few columns");
      return number(w[static_cast<std::size_t>(columns[k])], where);
    };
    Point pt;
    const int phase = static_cast<int>(column(0));
    pt.x = column(1);
    pt.y = column(2);
    pt.valid = phase != 0 && (p.phase == 0 || phase == p.phase);
    if (pt.valid) pt.rotation = bunge_rotation(column(3) * deg, column(4) * deg, column(5) * deg);
    points.push_back(pt);
  }
  if (ctf && !in_data) throw std::invalid_argument(source + ": no line of column names (Phase X Y ... Euler1 Euler2 Euler3)");
  if (points.empty()) throw std::invalid_argument(source + ": no points");
  if (std::none_of(points.begin(), points.end(), [](const Point &q) { return q.valid; })) {
    throw std::invalid_argument(source + ": no indexed point" + (p.phase != 0 ? " of phase " + std::to_string(p.phase) : ""));
  }

  double x_lo = 1e300, x_hi = -1e300, y_lo = 1e300, y_hi = -1e300;
  for (const auto &q : points) {
    x_lo = std::min(x_lo, q.x);
    x_hi = std::max(x_hi, q.x);
    y_lo = std::min(y_lo, q.y);
    y_hi = std::max(y_hi, q.y);
  }
  // The step: from the header, else from the density of the points.
  if (!(step > 0.0)) {
    const double area = std::max(x_hi - x_lo, 1e-12) * std::max(y_hi - y_lo, 1e-12);
    step = points.size() > 1 ? std::sqrt(area / static_cast<double>(points.size())) : 1.0;
  }

  // Neighbours: the points within 1.5 steps (4 or 8 on a square grid, 6 on
  // a hexagonal one), found through a grid of buckets one and a half steps wide.
  const double cell_size = 1.5 * step;
  const auto key = [&](double x, double y) {
    const auto ix = static_cast<std::int64_t>(std::floor((x - x_lo) / cell_size));
    const auto iy = static_cast<std::int64_t>(std::floor((y - y_lo) / cell_size));
    return (static_cast<std::uint64_t>(ix) << 32) ^ static_cast<std::uint64_t>(static_cast<std::uint32_t>(iy));
  };
  std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> buckets;
  buckets.reserve(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) buckets[key(points[i].x, points[i].y)].push_back(static_cast<std::uint32_t>(i));
  const auto neighbours = [&](std::size_t i, auto &&f) {
    for (int a = -1; a <= 1; ++a) {
      for (int b = -1; b <= 1; ++b) {
        const auto it = buckets.find(key(points[i].x + a * cell_size, points[i].y + b * cell_size));
        if (it == buckets.end()) continue;
        for (const std::uint32_t j : it->second) {
          if (j == i) continue;
          const double dx = points[j].x - points[i].x;
          const double dy = points[j].y - points[i].y;
          if (dx * dx + dy * dy <= cell_size * cell_size) f(static_cast<std::size_t>(j));
        }
      }
    }
  };

  // Grains: neighbouring points closer in orientation than the threshold.
  const auto symmetry = lattice_rotations(cell.lattice);
  const double threshold = in_units(p.misorientation, radian);
  std::vector<std::uint32_t> parent(points.size());
  std::iota(parent.begin(), parent.end(), 0u);
  const auto root = [&parent](std::uint32_t i) {
    while (parent[i] != i) {
      parent[i] = parent[parent[i]];
      i = parent[i];
    }
    return i;
  };
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!points[i].valid) continue;
    neighbours(i, [&](std::size_t j) {
      if (j < i || !points[j].valid) return;
      const std::uint32_t a = root(static_cast<std::uint32_t>(i));
      const std::uint32_t b = root(static_cast<std::uint32_t>(j));
      if (a == b) return;
      if (misorientation_rad(points[i].rotation, points[j].rotation, symmetry) <= threshold) parent[std::max(a, b)] = std::min(a, b);
    });
  }
  std::vector<int> grain_of(points.size(), -1);
  std::vector<Grain> grains;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (!points[i].valid) continue;
    const std::uint32_t r = root(static_cast<std::uint32_t>(i));
    if (grain_of[r] < 0) {
      // (The root is the grain's first point: its orientation is the grain's.)
      grain_of[r] = static_cast<int>(grains.size());
      Grain g;
      g.rotation = points[r].rotation;
      g.origin_A = lattice_offset(grains.size(), cell);
      grains.push_back(g);
    }
    grain_of[i] = grain_of[r];
  }
  // Points without an orientation take the grain of the nearest indexed point.
  std::deque<std::size_t> queue;
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (grain_of[i] >= 0) queue.push_back(i);
  }
  while (!queue.empty()) {
    const std::size_t i = queue.front();
    queue.pop_front();
    neighbours(i, [&](std::size_t j) {
      if (grain_of[j] >= 0) return;
      grain_of[j] = grain_of[i];
      queue.push_back(j);
    });
  }

  std::vector<GrainMap::Seed> seeds(points.size());
  for (std::size_t i = 0; i < points.size(); ++i) {
    // Map (x, y) -> layer (y, z) = (x, -y), from the map's corner.
    seeds[i].position_A = {0.0, (points[i].x - x_lo) * micrometre_A, (y_hi - points[i].y) * micrometre_A};
    seeds[i].grain = grain_of[i];
  }
  return GrainMap::from_seeds(std::move(grains), std::move(seeds),
                              period_or(p.period_y, (x_hi - x_lo + step) * micrometre_A),
                              period_or(p.period_z, (y_hi - y_lo + step) * micrometre_A), true);
}

} // namespace crystal
} // namespace amcset
