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

/*! \file trim_export.h
 *  \brief Internal to SRIM mode: the trajectory and MD event files
 *  (TrimConfig::trajectory_export, TrimConfig::md_events).
 *
 * Every block of 16 ions collects its points and events in an ExportBlock;
 * run() hands the blocks to ExportWriters in block order (the order the
 * tallies are merged in), so the files list the ions in order for any
 * number of threads. Nothing is kept after a block is written.
 */

#pragma once

#include "amcset_common.h"
#include "output/stream.h"
#include "trim/trim.h"

#include <filesystem>
#include <memory>
#include <vector>

namespace amcset {
namespace trim {

struct TrimSimulation::ExportBlock {
  std::vector<common::CollisionPoint> points; //!< Kept points (the export's filter applied)
  std::vector<MdEvent> events;
};

struct TrimSimulation::ExportWriters {
  //! Open the files the configuration asks for (relative names: in the
  //! output directory). \throws std::runtime_error, std::invalid_argument
  ExportWriters(const TrimConfig &config, const std::vector<std::pair<std::string, std::string>> &metadata);
  void write(const ExportBlock &block);
  //! Complete the files (Parquet footers). \throws std::runtime_error
  void finish();

  std::unique_ptr<output::StreamTableWriter> trajectories;
  std::unique_ptr<output::StreamTableWriter> events;
};

//! Where an export file goes: \p file itself if absolute or if there is no
//! output directory, else inside \p output_directory.
std::filesystem::path export_path(const std::string &file, const std::string &output_directory);

//! Whether the trajectory export keeps \p p.
bool export_keeps(const TrimConfig::TrajectoryExport &e, const common::CollisionPoint &p);

//! Sort the events [first, end) of one ion by E_c and time and number them
//! (MdEvent::event, from 1 for every E_c).
void number_events(std::vector<MdEvent> &events, std::size_t first);

} // namespace trim
} // namespace amcset
