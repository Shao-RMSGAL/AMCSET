// Copyright 2024, Texas A&M University
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

#include "amcset_common.h"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <stdexcept>

namespace amcset {
namespace common {
const std::vector<double> Volume::Layer::get_relative_compositions() const {
  std::vector<double> doubles;
  doubles.reserve(material_.size());
  std::transform(material_.begin(), material_.end(),
                 std::back_inserter(doubles),
                 [](const std::pair<double, Particle::Properties>& pair) {
                   return pair.first;
                 });
  return doubles;
}

Volume::Layer::Layer(const material_vector&& material, length_quantity depth)
    : material_(std::move(material)), depth_(depth) {
  if (depth < length_quantity(0)) {
    throw std::invalid_argument("Depth: " + std::to_string(depth.value()) +
                                " is less than zero.");
  }

  if (material.size() == 0) {
    throw std::invalid_argument("Cannot pass a vector of size 0");
  }

  auto relative_compositions = get_relative_compositions();
  double sum = std::accumulate(relative_compositions.begin(),
                               relative_compositions.end(), 0);

  for (auto& arg : material_) {
    if (arg.first < 0) {
      throw std::invalid_argument(
          "Relative composition: " + std::to_string(arg.first) +
          " is less than zero.");
    }
    arg.first /= sum;
  }
}
}  // namespace common
}  // namespace amcset
