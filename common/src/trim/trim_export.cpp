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

//! The trajectory file's columns (keys carry the units; docs/OUTPUT.md).
Table trajectory_schema() {
  Table t;
  t.key = "trajectories";
  t.title = "Trajectory points (SRIM mode, trim2013 transport)";
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
  return t;
}

//! The MD event file's columns.
Table event_schema() {
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

TrimSimulation::ExportWriters::ExportWriters(const TrimConfig &c,
                                             const std::vector<std::pair<std::string, std::string>> &metadata) {
  if (!c.trajectory_export.file.empty()) {
    auto m = metadata;
    m.emplace_back("points", c.trajectory_export.births_only ? "births" : "all");
    m.emplace_back("min_energy_eV", std::to_string(common::in_units(c.trajectory_export.min_energy, common::electron_volt)));
    trajectories = std::make_unique<output::StreamTableWriter>(export_path(c.trajectory_export.file, c.output_directory),
                                                               trajectory_schema(), m);
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
                                                         event_schema(), m);
  }
}

void TrimSimulation::ExportWriters::write(const ExportBlock &b) {
  if (trajectories && !b.points.empty()) {
    Table t = trajectory_schema();
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
    }
    trajectories->append(t);
  }
  if (events && !b.events.empty()) {
    Table t = event_schema();
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
