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

/*! \file volume.cpp
 *  \brief Volume: ordered layer stack and depth -> layer lookup.
 */

#include "amcset_common.h"
#include <absl/log/log.h>

namespace amcset {
namespace common {

Volume::Volume(std::vector<Layer> &&layers) try : layers_(std::move(layers)) {
  VLOG(1) << "Created volume with " << layers_.size() << " layers.";
  if (layers_.size() == 0) {
    throw std::invalid_argument(
        "Vector cannot be empty to Volume constructor.");
  }

  // Layer depths are cumulative (bottom boundary of each layer), so they must
  // be non-decreasing from the surface down.
  auto prev_depth = layers_.front().get_depth();
  for (const auto &layer : layers_) {
    if (layer.get_depth() < prev_depth) {
      throw std::invalid_argument(EXCEPTION_MESSAGE(
          "Provided vector of layers contains a layer which has a smaller "
          "depth than the previous layer. Depth: " +
          to_string(layer.get_depth()) +
          ". Previous depth: " + to_string(prev_depth)));
    }
    prev_depth = layer.get_depth();
  }
} catch (...) {
  rethrow();
}

const Layer &Volume::get_layer(size_t index) const try {
  if (index >= layers_.size()) {
    throw std::out_of_range(EXCEPTION_MESSAGE(
        "Index out of range. Index: " + std::to_string(index) +
        ". Vector size: " + std::to_string(layers_.size())));
  }
  return layers_[index];
} catch (...) {
  rethrow();
}

length_quantity Volume::get_depth() const { return layers_.back().get_depth(); }

const std::pair<const Layer *, size_t>
Volume::get_layer(length_quantity depth) const try {
  VLOG(2) << "Getting layer at depth " << depth;
  if (depth < length_quantity(0.0 * angstrom)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Depth: " + to_string(depth) + " cannot be less than zero."));
  }

  if (depth > layers_.back().get_depth()) {
    throw std::out_of_range(EXCEPTION_MESSAGE(
        "Depth: " + to_string(depth) + " is deeper than the deepest layer."));
  }

  // Rely on the ordering of the layers (depths are non-decreasing). A depth
  // exactly on a boundary belongs to the upper layer.
  for (size_t i = 0; i < layers_.size(); i++) {
    VLOG(3) << "Checking layer " << i << " with depth "
            << layers_[i].get_depth();
    if (depth <= layers_[i].get_depth()) {
      VLOG(2) << "Found match at " << i << " with depth "
              << layers_[i].get_depth() << " greater than " << depth;
      return std::make_pair(&layers_[i], i);
    }
  }

  // There is a logic error if this is reached
  throw std::logic_error(
      EXCEPTION_MESSAGE("End of layer vector reached. Maximum depth: " +
                        to_string(layers_.back().get_depth()) +
                        ". Provided depth: " + to_string(depth)));
} catch (...) {
  rethrow();
}

} // namespace common
} // namespace amcset
