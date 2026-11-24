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

/*! \file trim_dat.cpp
 *  \brief TRIM.DAT reader and writer.
 */

#include "trim/trim_dat.h"

#include <cctype>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace trim {

using namespace common;

namespace {

// Split on blanks and commas (TRIM accepts both separators).
std::vector<std::string> tokens(const std::string &line) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : line) {
    if (c == ',' || std::isspace(static_cast<unsigned char>(c)) || c == '\x1a') {
      if (!cur.empty()) out.push_back(std::move(cur));
      cur.clear();
    } else {
      cur.push_back(c);
    }
  }
  if (!cur.empty()) out.push_back(std::move(cur));
  return out;
}

bool to_number(const std::string &s, double &v) {
  if (s.empty()) return false;
  // VB's Val/Input accept a leading '.', e.g. ".1016364E+02" or "-.7111107".
  char *end = nullptr;
  errno = 0;
  v = std::strtod(s.c_str(), &end);
  return end != s.c_str() && *end == '\0' && errno != ERANGE;
}

void strip_line(std::string &line) {
  if (!line.empty() && line.back() == '\r') line.pop_back();
}

bool blank_or_eof(const std::string &line) {
  for (char c : line) {
    if (c != '\x1a' && !std::isspace(static_cast<unsigned char>(c))) return false;
  }
  return true;
}

/*!
 * A particle line that can be transported: whole Z in 1..118 (checked
 * before any integer conversion: NaN or 1e30 in an integer cast is
 * undefined behaviour), finite energy >= 0, finite position and a finite,
 * non-zero direction (a NaN direction would never let the particle stop).
 */
bool plausible_particle(double z, double energy, const double *xyz, const double *dir) {
  if (!(z >= 1.0 && z <= 118.0) || z != std::floor(z)) return false;
  if (!(energy >= 0.0) || !std::isfinite(energy)) return false;
  for (int k = 0; k < 3; ++k) {
    if (!std::isfinite(xyz[k]) || !std::isfinite(dir[k])) return false;
  }
  return dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2] > 0.0;
}

} // namespace

ParticleList parse_trim_dat(std::istream &in) {
  ParticleList list;
  std::string line;
  std::size_t number = 0;
  for (int i = 0; i < 10; ++i) {
    if (!std::getline(in, line)) {
      throw std::runtime_error("TRIM.DAT: the file must start with ten comment lines (line 8 "
                               "is the title); it ends after line " +
                               std::to_string(number));
    }
    ++number;
    strip_line(line);
    if (number == 8) {
      const auto b = line.find_first_not_of(" \t");
      const auto e = line.find_last_not_of(" \t");
      list.title = b == std::string::npos ? std::string() : line.substr(b, e - b + 1);
    }
  }
  while (std::getline(in, line)) {
    ++number;
    strip_line(line);
    if (blank_or_eof(line)) continue;
    const auto tok = tokens(line);
    // Name (one token) + 8 numbers. A name with blanks inside (allowed: any
    // 5 characters) leaves extra tokens at the front.
    double v[8];
    bool ok = tok.size() >= 9;
    if (ok) {
      const std::size_t first = tok.size() - 8;
      for (std::size_t k = 0; k < 8 && ok; ++k) ok = to_number(tok[first + k], v[k]);
    }
    if (!ok || !plausible_particle(v[0], v[1], v + 2, v + 5)) {
      list.warnings.push_back("TRIM.DAT line " + std::to_string(number) +
                              ": illegal data line skipped: " + line);
      continue;
    }
    SourceParticle p;
    {
      const std::size_t first = tok.size() - 8;
      std::string name;
      for (std::size_t k = 0; k < first; ++k) {
        if (k) name += ' ';
        name += tok[k];
      }
      p.name = name;
    }
    p.z = static_cast<unsigned>(v[0]);
    p.energy = v[1] * electron_volt;
    p.x = v[2] * angstrom;
    p.y = v[3] * angstrom;
    p.z_position = v[4] * angstrom;
    p.cos_x = v[5];
    p.cos_y = v[6];
    p.cos_z = v[7];
    list.particles.push_back(std::move(p));
  }
  return list;
}

ParticleList read_trim_dat(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open " + path.string());
  try {
    return parse_trim_dat(in);
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(path.string() + ": " + e.what());
  }
}

ParticleList read_exit_file_as_particles(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open " + path.string());
  ParticleList list;
  list.title = path.filename().string();
  std::string line;
  std::size_t number = 0;
  while (std::getline(in, line)) {
    ++number;
    strip_line(line);
    const auto tok = tokens(line);
    if (tok.size() < 10 || tok[0].size() != 1 || std::string("BTS").find(tok[0][0]) == std::string::npos) {
      continue;
    }
    double v[9];
    bool ok = true;
    for (std::size_t k = 0; k < 9 && ok; ++k) ok = to_number(tok[1 + k], v[k]);
    if (!ok || !(v[0] >= 0.0 && v[0] < 4294967296.0) || !plausible_particle(v[1], v[2], v + 3, v + 6)) {
      continue;
    }
    SourceParticle p;
    char buf[16];
    std::snprintf(buf, sizeof buf, "%c%04u", tok[0][0], static_cast<unsigned>(v[0]) % 10000u);
    p.name = buf;
    p.z = static_cast<unsigned>(v[1]);
    p.energy = v[2] * electron_volt;
    p.x = v[3] * angstrom;
    p.y = v[4] * angstrom;
    p.z_position = v[5] * angstrom;
    p.cos_x = v[6];
    p.cos_y = v[7];
    p.cos_z = v[8];
    list.particles.push_back(std::move(p));
  }
  if (list.particles.empty()) {
    throw std::runtime_error(path.string() + ": no B/T/S exit records found");
  }
  return list;
}

void write_trim_dat(std::ostream &out, const ParticleList &list) {
  const char *nl = "\r\n";
  out << "TRIM.DAT written by AMCSET: one particle per line after the 10 header lines" << nl
      << "Comment: Line #8 is the title written into the output files." << nl
      << "Comment: Data lines: EventName(5 char.) Z Energy(eV) X Y Z (Angstrom) cos(X) cos(Y) cos(Z)."
      << nl << "Comment: cos(X) = 1 for normal incidence; cos(X) = -1 towards the target surface." << nl
      << "Comment:" << nl << "Comment:" << nl << "Comment:" << nl << list.title << nl
      << "Event  Atom   Energy      Depth    -- Lateral Position --  ---- Atom Direction ----" << nl
      << "Name   Numb    (eV)       X (A)      Y (A)      Z (A)      cos(X)    cos(Y)    cos(Z)" << nl;
  char buf[256];
  for (const auto &p : list.particles) {
    std::string name = p.name.empty() ? std::string("-") : p.name;
    for (char &c : name) {
      if (std::isspace(static_cast<unsigned char>(c)) || c == ',') c = '_';
    }
    std::snprintf(buf, sizeof buf, "%-6s %3u %.7E %.6E %.6E %.6E %10.7f %10.7f %10.7f", name.c_str(),
                  p.z, in_units(p.energy, electron_volt), in_units(p.x, angstrom),
                  in_units(p.y, angstrom), in_units(p.z_position, angstrom), p.cos_x, p.cos_y,
                  p.cos_z);
    out << buf << nl;
  }
  out << '\x1a';
}

} // namespace trim
} // namespace amcset
