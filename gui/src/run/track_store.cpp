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

#include "run/track_store.h"

#include "amcset_common.h"

#include <algorithm>
#include <unordered_map>

namespace amcset {
namespace gui {

TrackStore::TrackStore(std::size_t max_segments) : max_segments_(max_segments) {}

void TrackStore::clear() {
  segments_.clear();
  segments_.shrink_to_fit();
  histories_ = 0;
  species_.clear();
  entries_.clear();
  entries_.shrink_to_fit();
  for (int k = 0; k < 3; ++k) lo_[k] = hi_[k] = true_lo_[k] = true_hi_[k] = 0.0f;
  ++revision_;
}

bool TrackStore::add_history(const std::vector<common::CollisionPoint> &points) {
  if (full()) return false;
  const auto history = static_cast<std::uint32_t>(histories_++);
  // The entry point on the beam spot (AMCSET x, y = TRIM y, z).
  Entry entry;
  if (!points.empty()) {
    entry.y = static_cast<float>(points.front().entry_x_A);
    entry.z = static_cast<float>(points.front().entry_y_A);
  }
  entries_.push_back(entry);
  const float shift[3] = {0.0f, entry.y, entry.z};
  // Last point of every open path.
  struct Last {
    float p[3];
  };
  std::unordered_map<std::uint32_t, Last> last;
  last.reserve(64);
  const bool first_history = segments_.empty();
  for (const auto &cp : points) {
    // AMCSET (x, y lateral, z depth) -> TRIM (x depth, y, z lateral), Angstrom.
    const float p[3] = {
        static_cast<float>(common::in_units(cp.coordinate.z_, common::angstrom)),
        static_cast<float>(common::in_units(cp.coordinate.x_, common::angstrom)),
        static_cast<float>(common::in_units(cp.coordinate.y_, common::angstrom))};
    auto it = last.find(cp.track);
    if (it != last.end()) {
      if (segments_.size() >= max_segments_) return false;
      TrackSegment s{};
      for (int k = 0; k < 3; ++k) {
        s.a[k] = it->second.p[k];
        s.b[k] = p[k];
      }
      s.z = static_cast<std::uint8_t>(std::min<std::size_t>(cp.z_number, 255));
      s.generation = static_cast<std::uint8_t>(std::min<std::size_t>(cp.generation, 255));
      s.flags = 0;
      s.history = history;
      segments_.push_back(s);
      if (s.generation > 0) {
        if (species_.size() <= s.z) species_.resize(s.z + 1, 0);
        ++species_[s.z];
      }
      for (int k = 0; k < 3; ++k) {
        const float ta = s.a[k] + shift[k];
        const float tb = s.b[k] + shift[k];
        if (first_history && segments_.size() == 1) {
          lo_[k] = std::min(s.a[k], s.b[k]);
          hi_[k] = std::max(s.a[k], s.b[k]);
          true_lo_[k] = std::min(ta, tb);
          true_hi_[k] = std::max(ta, tb);
        } else {
          lo_[k] = std::min(lo_[k], std::min(s.a[k], s.b[k]));
          hi_[k] = std::max(hi_[k], std::max(s.a[k], s.b[k]));
          true_lo_[k] = std::min(true_lo_[k], std::min(ta, tb));
          true_hi_[k] = std::max(true_hi_[k], std::max(ta, tb));
        }
      }
      for (int k = 0; k < 3; ++k) it->second.p[k] = p[k];
    } else {
      Last l;
      for (int k = 0; k < 3; ++k) l.p[k] = p[k];
      last.emplace(cp.track, l);
    }
  }
  return true;
}

} // namespace gui
} // namespace amcset
