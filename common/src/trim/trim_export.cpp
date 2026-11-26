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

/*! \file trim_export.cpp
 *  \brief The trajectory and MD event files: their columns (documented in
 *  docs/OUTPUT.md, which must be kept in step with this file) and rows.
 */

#include "trim_export.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <unordered_map>

namespace amcset {
namespace trim {

using common::CollisionPoint;
using output::Column;
using output::Dimension;
using output::Table;

namespace {

double parent_value(std::uint32_t v) { return v == CollisionPoint::no_parent ? -1.0 : static_cast<double>(v); }

const char *kind_text(CollisionPoint::Kind k) {
  switch (k) {
  case CollisionPoint::Kind::collision: return "collision";
  case CollisionPoint::Kind::birth: return "birth";
  case CollisionPoint::Kind::stop: return "stop";
  case CollisionPoint::Kind::exit: return "exit";
  }
  return "collision";
}

//! Crystal transport: where the crystal lattice is for this ion.
void add_lattice_offset_columns(Table &t) {
  const char *what = "Crystal transport: the lateral offset of the crystal lattice for this ion (its entry point "
                     "on the beam spot). A point (x, y, z) of this file is at (x + lattice_dx_A, y + lattice_dy_A, "
                     "z - top of the layer) in the frame of the layer's lattice, whose origin is a lattice point; "
                     "the crystal's orientation is in the file's metadata and in summary.json";
  t.columns.push_back(Column::real("lattice_dx_A", "Lattice offset x", Dimension::length, "A", {}, what));
  t.columns.push_back(Column::real("lattice_dy_A", "Lattice offset y", Dimension::length, "A", {}, what));
}

//! The trajectory file's columns (keys carry the units; docs/OUTPUT.md).
Table trajectory_schema(bool lattice_offsets) {
  Table t;
  t.key = "trajectories";
  t.title = "Trajectory points (SRIM mode)";
  t.description = "One row per recorded point. Axes: z is the depth below the surface, x and y lateral "
                  "(TRIM's y and z). Energies and directions are those after the collision at the point.";
  const auto i = [&](const char *k, const char *label, const char *d) { t.columns.push_back(Column::integer(k, label, {}, d)); };
  const auto r = [&](const char *k, const char *label, Dimension dim, const char *u, const char *d) {
    t.columns.push_back(Column::real(k, label, dim, u, {}, d));
  };
  i("ion", "Ion", "1-based number of the incident particle");
  i("point", "Point", "Index of the point in its ion's trajectory (0 = the ion's entry)");
  i("track", "Track", "Path within the ion's trajectory: 0 = the ion, one number per recoil");
  i("parent_track", "Parent track", "Track of the atom whose collision set this one in motion (-1: none)");
  i("parent_point", "Parent point", "Point (index) of that collision (-1: none)");
  t.columns.push_back(Column::strings("kind", "Kind", {}, "birth, collision, stop or exit"));
  i("gen", "Generation", "0 = the ion, 1 = atoms it struck, 2 = atoms they struck, ...");
  i("Z", "Z", "Atomic number");
  r("mass_amu", "Mass", Dimension::dimensionless, "1", "Atomic mass in amu");
  t.columns.back().unit = "amu"; // (no mass dimension in output/units.h; metadata only)
  r("x_A", "x", Dimension::length, "A", "Lateral position");
  r("y_A", "y", Dimension::length, "A", "Lateral position");
  r("z_A", "z (depth)", Dimension::length, "A", "Depth below the surface");
  r("ux", "ux", Dimension::dimensionless, "1", "Direction cosine (0 0 0 if at rest)");
  r("uy", "uy", Dimension::dimensionless, "1", "Direction cosine");
  r("uz", "uz", Dimension::dimensionless, "1", "Direction cosine along the depth");
  r("E_eV", "Energy", Dimension::energy, "eV", "Kinetic energy (after the collision at a collision point)");
  r("E_before_eV", "Energy before", Dimension::energy, "eV", "Collision points: energy just before the collision");
  r("T_eV", "Transferred", Dimension::energy, "eV", "Collision points: energy given to the struck atom");
  r("dE_e_eV", "Electronic loss", Dimension::energy, "eV",
    "Collision points: electronic energy loss over the flight that led to the point");
  r("t_ps", "Time", Dimension::time, "ps", "Time since the ion entered (flight length / speed, summed)");
  if (lattice_offsets) add_lattice_offset_columns(t);
  return t;
}

//! The MD event file's columns.
Table event_schema(bool lattice_offsets) {
  Table t;
  t.key = "md_events";
  t.title = "Recoil events for molecular dynamics (hybrid BCA-MD)";
  t.description = "One row per event and transition energy E_c. Axes as in the trajectory file: z is the depth.";
  const auto i = [&](const char *k, const char *label, const char *d) { t.columns.push_back(Column::integer(k, label, {}, d)); };
  const auto r = [&](const char *k, const char *label, Dimension dim, const char *u, const char *d) {
    t.columns.push_back(Column::real(k, label, dim, u, {}, d));
  };
  i("ion_id", "Ion", "1-based number of the incident particle");
  i("event_id", "Event", "1-based within the ion and E_c, in time order");
  r("E_c_eV", "E_c", Dimension::energy, "eV", "The transition energy the event belongs to");
  t.columns.push_back(Column::strings("kind", "Kind", {},
                                      "birth: a recoil set in motion below E_c (E = the energy transferred); "
                                      "crossing: a moving atom whose energy fell below E_c at a collision"));
  r("t_ps", "Time", Dimension::time, "ps", "Since the ion entered the target");
  r("x_A", "x", Dimension::length, "A", "Lateral position");
  r("y_A", "y", Dimension::length, "A", "Lateral position");
  r("z_A", "z (depth)", Dimension::length, "A", "Depth below the surface");
  r("E_eV", "Energy", Dimension::energy, "eV", "Kinetic energy handed to MD");
  r("ux", "ux", Dimension::dimensionless, "1", "Direction cosine");
  r("uy", "uy", Dimension::dimensionless, "1", "Direction cosine");
  r("uz", "uz", Dimension::dimensionless, "1", "Direction cosine along the depth");
  i("gen", "Generation", "0 = the ion, 1 = its primaries, ...");
  i("Z", "Z", "Atomic number");
  r("mass_amu", "Mass", Dimension::dimensionless, "1", "Atomic mass in amu");
  t.columns.back().unit = "amu"; // (no mass dimension in output/units.h; metadata only)
  i("track", "Track", "The atom's track (as in the trajectory file)");
  i("parent_track", "Parent track", "Track of the atom that set it in motion (-1: none)");
  i("parent_event", "Parent event",
    "Always -1: an event's atom is not followed further for its E_c, so no event descends from another");
  if (lattice_offsets) add_lattice_offset_columns(t);
  return t;
}

void push(Table &t, std::size_t k, double v) { t.columns[k].values.push_back(v); }

} // namespace

std::filesystem::path export_path(const std::string &file, const std::string &output_directory) {
  const std::filesystem::path f(file);
  if (f.is_absolute() || output_directory.empty()) return f;
  return std::filesystem::path(output_directory) / f;
}

bool export_keeps(const TrimConfig::TrajectoryExport &e, const CollisionPoint &p) {
  if (e.births_only && p.kind != CollisionPoint::Kind::birth) return false;
  return p.energy_eV >= common::in_units(e.min_energy, common::electron_volt);
}

void number_events(std::vector<MdEvent> &events, std::size_t first) {
  if (first >= events.size()) return;
  std::stable_sort(events.begin() + static_cast<std::ptrdiff_t>(first), events.end(),
                   [](const MdEvent &a, const MdEvent &b) {
                     if (a.transition_energy_eV != b.transition_energy_eV) return a.transition_energy_eV < b.transition_energy_eV;
                     return a.time_s < b.time_s;
                   });
  std::size_t n = 0;
  double ec = -1.0;
  for (std::size_t i = first; i < events.size(); ++i) {
    if (events[i].transition_energy_eV != ec) {
      ec = events[i].transition_energy_eV;
      n = 0;
    }
    events[i].event = ++n;
  }
}

bool exports_lattice_offsets(const TrimConfig &config) {
  if (config.transport != TrimConfig::Transport::crystal) return false;
  return std::any_of(config.layers.begin(), config.layers.end(), [](const TargetLayer &l) { return l.crystal.has_value(); });
}

std::vector<std::pair<std::string, std::string>> crystal_metadata(const TrimConfig &config) {
  std::vector<std::pair<std::string, std::string>> m;
  const auto num = [](double v) {
    char buf[40];
    std::snprintf(buf, sizeof buf, "%.10g", v);
    return std::string(buf);
  };
  const auto matrix = [&](const crystal::Mat3 &r) {
    std::string s;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) s += (i == 0 && j == 0 ? "" : " ") + num(r[i][j]);
    }
    return s;
  };
  for (const auto &r : crystal_reports(config)) {
    const std::string k = "crystal.layer" + std::to_string(r.layer + 1) + ".";
    m.emplace_back(k + "name", r.cell_name);
    m.emplace_back(k + "followed", r.followed ? "yes (crystal transport)" : "no (treated as amorphous)");
    m.emplace_back(k + "top_A", num(r.top_A));
    m.emplace_back(k + "lattice_vectors_A", matrix(r.lattice_A));
    m.emplace_back(k + "rotation_crystal_to_file_axes", matrix(r.to_output));
    m.emplace_back(k + "euler_bunge_deg", num(r.euler_deg[0]) + " " + num(r.euler_deg[1]) + " " + num(r.euler_deg[2]));
    m.emplace_back(k + "depth_direction", uvw_text(r.depth_uvw) + " + " + num(r.depth_off_deg) + " deg");
    m.emplace_back(k + "beam_direction", uvw_text(r.beam_uvw) + " + " + num(r.beam_off_deg) + " deg");
    m.emplace_back(k + "tilt_twist_deg", num(r.tilt_deg) + " " + num(r.twist_deg));
    m.emplace_back(k + "temperature_K", num(r.temperature_K));
  }
  if (!m.empty()) {
    m.emplace_back("crystal.convention",
                   "rotation_crystal_to_file_axes maps the crystal's Cartesian frame (lattice_vectors_A, rows) to "
                   "the file's axes (x, y lateral, z depth); the lattice origin of a layer is at its top, shifted "
                   "laterally per ion by the columns lattice_dx_A, lattice_dy_A");
  }
  return m;
}

TrimSimulation::ExportWriters::ExportWriters(const TrimConfig &c,
                                             const std::vector<std::pair<std::string, std::string>> &given) {
  lattice_offsets = exports_lattice_offsets(c);
  auto metadata = given;
  for (auto &kv : crystal_metadata(c)) metadata.push_back(std::move(kv));
  if (!c.trajectory_export.file.empty()) {
    auto m = metadata;
    m.emplace_back("points", c.trajectory_export.births_only ? "births" : "all");
    m.emplace_back("min_energy_eV", std::to_string(common::in_units(c.trajectory_export.min_energy, common::electron_volt)));
    trajectories = std::make_unique<output::StreamTableWriter>(export_path(c.trajectory_export.file, c.output_directory),
                                                               trajectory_schema(lattice_offsets), m);
  }
  if (c.md_events.enabled()) {
    auto m = metadata;
    std::string ec;
    for (const auto &e : c.md_events.transition_energies) {
      ec += (ec.empty() ? "" : ",") + std::to_string(common::in_units(e, common::electron_volt));
    }
    m.emplace_back("E_c_eV", ec);
    m.emplace_back("E_min_eV", c.md_events.min_energy > 0.0 * common::electron_volt
                                   ? std::to_string(common::in_units(c.md_events.min_energy, common::electron_volt))
                                   : std::string("Ed of the struck element"));
    m.emplace_back("mode", c.md_events.handoff ? "handoff" : "record");
    events = std::make_unique<output::StreamTableWriter>(export_path(c.md_events.file, c.output_directory),
                                                         event_schema(lattice_offsets), m);
  }
}

void TrimSimulation::ExportWriters::write(const ExportBlock &b) {
  std::unordered_map<std::size_t, std::array<double, 2>> offset;
  for (const auto &[ion, o] : b.lattice_offsets) offset[ion] = o;
  const auto offsets = [&](Table &t, std::size_t &k, std::size_t ion) {
    if (!lattice_offsets) return;
    const auto it = offset.find(ion);
    push(t, k++, it != offset.end() ? it->second[0] : 0.0);
    push(t, k++, it != offset.end() ? it->second[1] : 0.0);
  };
  if (trajectories && !b.points.empty()) {
    Table t = trajectory_schema(lattice_offsets);
    for (auto &col : t.columns) {
      if (col.type == output::ColumnType::text) col.text.reserve(b.points.size());
      else col.values.reserve(b.points.size());
    }
    const double m_to_A = 1e10; // CollisionPoint positions are SI (metres)
    for (const auto &p : b.points) {
      std::size_t k = 0;
      push(t, k++, static_cast<double>(p.ion));
      push(t, k++, static_cast<double>(p.point_index));
      push(t, k++, static_cast<double>(p.track));
      push(t, k++, parent_value(p.parent_track));
      push(t, k++, parent_value(p.parent_point));
      t.columns[k++].text.emplace_back(kind_text(p.kind));
      push(t, k++, static_cast<double>(p.generation));
      push(t, k++, static_cast<double>(p.z_number));
      push(t, k++, p.atomic_mass_number);
      push(t, k++, p.coordinate.x_.value() * m_to_A);
      push(t, k++, p.coordinate.y_.value() * m_to_A);
      push(t, k++, p.coordinate.z_.value() * m_to_A);
      const double v = std::sqrt(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
      push(t, k++, v > 0.0 ? p.vx / v : 0.0);
      push(t, k++, v > 0.0 ? p.vy / v : 0.0);
      push(t, k++, v > 0.0 ? p.vz / v : 0.0);
      push(t, k++, p.energy_eV);
      push(t, k++, p.energy_before_eV);
      push(t, k++, p.energy_transferred_eV);
      push(t, k++, p.electronic_loss_eV);
      push(t, k++, p.time * 1e12);
      offsets(t, k, static_cast<std::size_t>(p.ion));
    }
    trajectories->append(t);
  }
  if (events && !b.events.empty()) {
    Table t = event_schema(lattice_offsets);
    for (const auto &e : b.events) {
      std::size_t k = 0;
      push(t, k++, static_cast<double>(e.ion));
      push(t, k++, static_cast<double>(e.event));
      push(t, k++, e.transition_energy_eV);
      t.columns[k++].text.emplace_back(e.kind == MdEvent::Kind::birth ? "birth" : "crossing");
      push(t, k++, e.time_s * 1e12);
      push(t, k++, e.x_A);
      push(t, k++, e.y_A);
      push(t, k++, e.z_A);
      push(t, k++, e.energy_eV);
      push(t, k++, e.ux);
      push(t, k++, e.uy);
      push(t, k++, e.uz);
      push(t, k++, static_cast<double>(e.generation));
      push(t, k++, static_cast<double>(e.z));
      push(t, k++, e.mass_amu);
      push(t, k++, static_cast<double>(e.track));
      push(t, k++, parent_value(e.parent_track));
      push(t, k++, -1.0);
      offsets(t, k, e.ion);
    }
    events->append(t);
  }
}

void TrimSimulation::ExportWriters::finish() {
  if (trajectories) trajectories->finish();
  if (events) events->finish();
}

} // namespace trim
} // namespace amcset
