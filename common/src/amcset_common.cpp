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

#include <numeric>
#include <stdexcept>

#include "amcset_utilities.h"

namespace amcset {
namespace common {

Volume::Layer::Layer(material_vector&& material, length_quantity depth) try
    : material_(std::move(material)), depth_(depth) {
  if (depth < length_quantity(0)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Depth: " + to_string_with_unit(depth) + " is less than zero."));
  }

  if (material_.size() == 0) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Cannot pass a vector of size 0 to Volume::Layer constructor"));
  }

  auto relative_compositions = get_relative_compositions();
  double sum = std::accumulate(relative_compositions.begin(),
                               relative_compositions.end(), 0);

  for (auto& arg : material_) {
    if (arg.first < 0) {
      throw std::invalid_argument(EXCEPTION_MESSAGE(
          "Relative composition: " + std::to_string(arg.first) +
          " is less than zero."));
    }
    arg.first /= sum;
  }
} catch (...) {
  rethrow();
}

const Volume::Layer& Volume::get_layer(length_quantity depth) const try {
  if (depth < length_quantity(0.0 * angstrom)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Depth: " + to_string_with_unit(depth) + " cannot be less than zero."));
  }

  if (depth > layers_.back().get_depth()) {
    throw std::out_of_range(
        EXCEPTION_MESSAGE("Depth: " + to_string_with_unit(depth) +
                          " is deeper than the deepest layer."));
  }

  // Rely on the ordering of the layers (depths are strictly increasing)
  for (const auto& layer : layers_) {
    if (depth <= layer.get_depth()) {
      return layer;
    }
  }

  // There is a logic error if this is reached
  throw std::logic_error(
      EXCEPTION_MESSAGE("End of layer vector vector reached. Maximum depth: " +
                        to_string_with_unit(layers_.back().get_depth()) +
                        ". Provided depth: " + to_string_with_unit(depth)));
} catch (...) {
  rethrow();
}

auto Volume::Layer::get_relative_compositions() const
    -> decltype(std::declval<const material_vector&>() |
                std::views::transform(
                    &std::pair<double, Particle::Properties>::first)) try {
  return material_ |
         std::views::transform(&std::pair<double, Particle::Properties>::first);
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

Volume::Volume(std::vector<Layer>&& layers) try : layers_(std::move(layers)) {
  if (layers_.size() == 0) {
    throw std::invalid_argument(
        "Vector cannot be empty to Volume constructor.");
  }

  auto prev_depth = layers_.at(0).get_depth();
  for (const auto& layer : layers_) {
    if (layer.get_depth() < prev_depth) {
      throw std::invalid_argument(EXCEPTION_MESSAGE(
          "Provided vector of layers contains a layer which has a smaller "
          "depth than the previous layer. Depth: " +
          to_string_with_unit(layer.get_depth()) +
          ". Previous depth: " + to_string_with_unit(prev_depth)));
    }
  }
} catch (...) {
  rethrow();
}

Simulation::Settings::Settings(
    energy_quantity electron_stopping_energy,
    Particle::Properties incident_particle_properties,
    bool enable_damage_cascade, energy_quantity ion_stopping_energy,
    energy_quantity ion_displacement_energy, bool log_single_displacement,
    size_t divisor_angle_number, size_t flying_distance_number,
    length_quantity range, size_t bombardment_count, bool is_electron,
    energy_quantity incident_energy) try
    : electron_stopping_energy_(electron_stopping_energy),
      incident_particle_properties_(incident_particle_properties),
      enable_damage_cascade_(enable_damage_cascade),
      ion_stopping_energy_(ion_stopping_energy),
      ion_displacement_energy_(ion_displacement_energy),
      log_single_displacement_(log_single_displacement),
      divisor_angle_number_(divisor_angle_number),
      flying_distance_number_(flying_distance_number),
      range_(range),
      bombardment_count_(bombardment_count),
      is_electron_(is_electron),
      incident_energy_(incident_energy) {
  if (incident_energy <= energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Incident energy should be greater than 0. Value: " +
                          to_string_with_unit(incident_energy)));
  }
  if (bombardment_count == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Bombardment count must be greater than 0"));
  }
  if (divisor_angle_number == 0) {
    EXCEPTION_MESSAGE("Divisor angle number must be greater than 0");
  }
  if (ion_displacement_energy < energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Ion displacement energy cannot be less than 0. Value: " +
        to_string_with_unit(ion_displacement_energy)));
  }
  if (range < length_quantity(0.0 * angstrom)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Range cannot be less than 0. Value: " + to_string_with_unit(range)));
  }
  if (flying_distance_number == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Flying distance number must be greater than 0"));
  }
  if (electron_stopping_energy < energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Electron stopping energy cannot be less than 0. Value: " +
        to_string_with_unit(electron_stopping_energy)));
  }
} catch (...) {
  rethrow();
}

void Simulation::run() {}

}  // namespace common
}  // namespace amcset
