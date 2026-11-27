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

/*! \file track_store.h
 *  \brief The particle paths shown by the trajectory views.
 *
 * Trajectories arrive from the engine as lists of collision points; the store
 * turns them into straight segments (consecutive points of the same path)
 * kept in one append-only array, so the GPU views only upload what was added
 * since their last frame. Coordinates follow TRIM: x = depth, y and z
 * lateral, in Angstrom.
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace amcset {
namespace common {
struct CollisionPoint;
}
namespace gui {

//! One straight piece of a particle path.
struct TrackSegment {
  float a[3];               //!< Start (depth x, lateral y, lateral z), Angstrom
  float b[3];               //!< End
  std::uint8_t z;           //!< Atomic number of the moving particle (0 = electron)
  std::uint8_t generation;  //!< 0 = incident particle, 1 = primary recoil, ...
  std::uint16_t flags;      //!< Reserved
  std::uint32_t history;    //!< Incident particle number (0-based, in arrival order)
};

//! Append-only segment store with a size cap.
class TrackStore {
public:
  explicit TrackStore(std::size_t max_segments = 4'000'000);

  /*!
   * \brief Add the trajectory of one incident particle (points in AMCSET's
   * convention: coordinate z = depth). Points with the same track number are
   * joined in order. Returns false (and adds nothing) once the cap is reached.
   */
  bool add_history(const std::vector<common::CollisionPoint> &points);
  //! Remove everything (views re-upload: revision() changes).
  void clear();

  const std::vector<TrackSegment> &segments() const noexcept { return segments_; }
  std::size_t size() const noexcept { return segments_.size(); }
  std::size_t histories() const noexcept { return histories_; }
  //! Changes when existing segments are removed (appends keep it).
  std::uint64_t revision() const noexcept { return revision_; }
  bool full() const noexcept { return segments_.size() >= max_segments_; }
  std::size_t max_segments() const noexcept { return max_segments_; }
  void set_max_segments(std::size_t n) { max_segments_ = n; }

  //! Bounding box of everything stored (empty store: zero box).
  const float *lower() const noexcept { return lo_; }
  const float *upper() const noexcept { return hi_; }
  /*!
   * The entry point of each incident particle on the beam spot (TRIM's
   * lateral y, z; Angstrom; index = TrackSegment::history): the segments are
   * relative to it, their true place on the sample is a + (0, y, z). Zero
   * outside the crystal transport.
   */
  struct Entry {
    float y = 0.0f;
    float z = 0.0f;
  };
  const std::vector<Entry> &entries() const noexcept { return entries_; }
  //! Bounding box of the segments at their true places (with the entry points).
  const float *true_lower() const noexcept { return true_lo_; }
  const float *true_upper() const noexcept { return true_hi_; }
  //! Atomic numbers of the recoils seen so far (index = Z, value = segment count).
  const std::vector<std::size_t> &species() const noexcept { return species_; }

private:
  std::vector<TrackSegment> segments_;
  std::size_t max_segments_;
  std::size_t histories_ = 0;
  std::uint64_t revision_ = 1;
  float lo_[3] = {0, 0, 0};
  float hi_[3] = {0, 0, 0};
  float true_lo_[3] = {0, 0, 0};
  float true_hi_[3] = {0, 0, 0};
  std::vector<std::size_t> species_;
  std::vector<Entry> entries_;
};

} // namespace gui
} // namespace amcset
