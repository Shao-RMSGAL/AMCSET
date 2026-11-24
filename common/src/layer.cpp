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

/*! \file layer.cpp
 *  \brief Layer construction: composition normalization, number density and
 *  the isotope-selection distribution.
 */

#include "amcset_common.h"
#include <absl/log/log.h>

#include <numeric>

namespace amcset {
namespace common {

Layer::Layer(material_vector &&material, length_quantity depth,
             mass_density_quantity mass_density) try
    : material_([&material, depth, mass_density]() {
        VLOG(1) << "Creating Layer material with vector of length "
                << material.size() << ", depth " << depth << ", and density "
                << mass_density;
        auto relative_compositions =
            material | std::views::transform(
                           &std::pair<double, Particle::Properties>::first);
        double sum = std::accumulate(relative_compositions.begin(),
                                     relative_compositions.end(), 0.0);
        VLOG(2) << "Sum of compositions is " << sum;
        // Guard against NaN fractions (0/0) when the material is non-empty
        // but every fraction is zero. An empty material is rejected below.
        if (!material.empty() && !(sum > 0.0)) {
          throw std::invalid_argument(EXCEPTION_MESSAGE(
              "Relative compositions must sum to a positive value. Sum: " +
              std::to_string(sum)));
        }

        for (auto &pair : material) {
          if (pair.first < 0) {
            throw std::invalid_argument(EXCEPTION_MESSAGE(
                "Relative composition: " + std::to_string(pair.first) +
                " is less than zero."));
          }
          VLOG(3) << "Initial composition is " << pair.first;
          pair.first /= sum;
          VLOG(3) << "Final composition is " << pair.first;
        }

        return std::move(material);
      }()),
      depth_(depth), mass_density_(mass_density), number_density_([&]() {
        // Atoms per volume: mass density / mean atomic mass.
        mass_quantity mean_mass = 0.0 * si::kilogram;
        for (const auto &component : material_) {
          mean_mass += component.first * component.second.mass_;
        }
        return mass_density / mean_mass;
      }()) {
  VLOG(1) << "Created Layer with " << material_.size()
          << " components, depth of " << depth << ", density of "
          << mass_density << ", and number density of " << number_density_;

  if (depth < length_quantity(0)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE("Depth: " + to_string(depth) +
                                                  " is less than zero."));
  }

  if (material_.size() == 0) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Cannot pass a vector of size 0 to Layer constructor"));
  }

  VLOG(2) << "Passing relative compositions to discrete distribution.";
  auto composition_view = get_relative_compositions();
  discrete_distribution_ = boost::random::discrete_distribution<size_t, double>(
      composition_view.begin(), composition_view.end());

} catch (...) {
  rethrow();
}

void Layer::set_displacement_energies(std::vector<energy_quantity> energies) {
  if (!energies.empty() && energies.size() != material_.size()) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Expected " + std::to_string(material_.size()) +
        " displacement energies, got " + std::to_string(energies.size())));
  }
  displacement_energies_ = std::move(energies);
}

auto Layer::get_relative_compositions() const
    -> decltype(std::declval<const material_vector &>() |
                std::views::keys) try {
  VLOG(2) << "Getting relative compositions for layer with " << material_.size()
          << " components";
  return material_ | std::views::keys;
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

} // namespace common
} // namespace amcset
