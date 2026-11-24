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

/*! \file simulation_engine.h
 *  \brief Front-end-neutral interface implemented by every simulation engine
 *  (the native AMCSET Simulation and the SRIM-mode trim::TrimSimulation).
 *
 * The GUI and the CLI drive simulations only through this interface: run on
 * a worker thread, poll progress, live statistics and trajectories from
 * another thread, and cancel.
 */

#pragma once

#include "amcset_units.h"
#include "tallies.h"

#include <cstddef>
#include <string>
#include <vector>

namespace amcset {
namespace common {

struct CollisionPoint;

/*!
 * \brief Abstract simulation engine.
 *
 * Thread-safety: run() is called once, on one thread. Every other member
 * may be called concurrently from any thread while run() executes.
 */
class SimulationEngine {
public:
  virtual ~SimulationEngine() = default;

  /*!
   * \brief Run the whole simulation (blocking) and write its output files.
   * \return 0 on success (including cancellation).
   */
  virtual int run() = 0;

  //! Ask the engine to stop as soon as possible (pending histories are skipped).
  virtual void request_cancel() noexcept = 0;
  //! Whether request_cancel() was called.
  virtual bool is_cancelled() const noexcept = 0;
  //! Whether set_paused() has an effect.
  virtual bool supports_pause() const noexcept { return false; }
  /*!
   * \brief Pause (true) or resume (false) the calculation: the workers stop
   * between groups of incident particles. Cancelling also resumes.
   */
  virtual void set_paused(bool /*paused*/) {}
  //! Whether the calculation is paused.
  virtual bool is_paused() const noexcept { return false; }

  //! Number of incident particles finished so far.
  virtual std::size_t completed_count() const noexcept = 0;
  //! Total number of incident particles.
  virtual std::size_t total_count() const noexcept = 0;

  //! Snapshot of the live statistics (depths in metres, energies in joules).
  virtual Tallies snapshot_tallies() const = 0;

  /*!
   * \brief Pop the trajectory (collision points) of one finished incident
   * particle, if any is queued. Coordinates follow the AMCSET convention
   * (+z = depth). \return false if the queue is empty.
   */
  virtual bool try_pop_trajectory(std::vector<CollisionPoint> &out) = 0;
  //! Drop all queued trajectories.
  virtual void clear_trajectories() = 0;

  //! Depth extent of the plots and histograms.
  virtual length_quantity depth_extent() const = 0;

  //! Human-readable description of the configuration (for logs).
  virtual std::string describe() const = 0;
};

} // namespace common
} // namespace amcset
