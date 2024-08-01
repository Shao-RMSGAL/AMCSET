// Copyright 2024, Texas A&M University
//
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

#include <glog/logging.h>

#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "amcset_utilities.h"

namespace amcset {
namespace common {

Layer::Layer(material_vector &&material, length_quantity depth,
             mass_density_quantity mass_density) try
    : material_([&material]() {
          auto relative_compositions =
              material | std::views::transform(
                             &std::pair<double, Particle::Properties>::first);
          double sum = std::accumulate(relative_compositions.begin(),
                                       relative_compositions.end(), 0.0);

          for (auto &pair : material) {
              if (pair.first < 0) {
                  throw std::invalid_argument(EXCEPTION_MESSAGE(
                      "Relative composition: " + std::to_string(pair.first) +
                      " is less than zero."));
              }
              pair.first /= sum;
          }
          return std::move(material);
      }()),
      depth_(depth),
      mass_density_(mass_density),
      number_density_([&]() {
          return mass_density /
                 (std::accumulate(
                      material_.begin(), material_.end(), mass_quantity(0),
                      [](const mass_quantity &sum, const auto &mat) {
                          return sum + mat.first * mat.second.mass_;
                      }) *
                  constants::N_A);
      }()) {
    if (depth < length_quantity(0)) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Depth: " + to_string(depth) + " is less than zero."));
    }

    if (material_.size() == 0) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Cannot pass a vector of size 0 to Layer constructor"));
    }

    discrete_distribution_ =
        boost::random::discrete_distribution<size_t, double>(
            get_relative_compositions().begin(),
            get_relative_compositions().end());
} catch (...) {
    rethrow();
}

std::pair<const Layer *, size_t> Volume::get_layer(length_quantity depth) const
    try {
    if (depth < length_quantity(0.0 * angstrom)) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Depth: " + to_string(depth) + " cannot be less than zero."));
    }

    if (depth > layers_.back().get_depth()) {
        throw std::out_of_range(
            EXCEPTION_MESSAGE("Depth: " + to_string(depth) +
                              " is deeper than the deepest layer."));
    }

    // Rely on the ordering of the layers (depths are strictly increasing)
    for (size_t i = 0; i < layers_.size(); i++) {
        if (depth <= layers_[i].get_depth()) {
            return std::make_pair(&layers_[i], i);
        }
    }

    // There is a logic error if this is reached
    throw std::logic_error(EXCEPTION_MESSAGE(
        "End of layer vector vector reached. Maximum depth: " +
        to_string(layers_.back().get_depth()) +
        ". Provided depth: " + to_string(depth)));
} catch (...) {
    rethrow();
}

const Layer &Volume::get_layer(size_t index) const try {
    if (index > layers_.size()) {
        throw std::out_of_range(EXCEPTION_MESSAGE(
            "Index out of range. Index: " + std::to_string(index) +
            ". Vector size: " + std::to_string(layers_.size())));
    }
    return layers_[index];
} catch (...) {
    rethrow();
}
auto Layer::get_relative_compositions() const
    -> decltype(std::declval<const material_vector &>() |
                std::views::transform(
                    &std::pair<double, Particle::Properties>::first)) try {
    return material_ | std::views::transform(
                           &std::pair<double, Particle::Properties>::first);
} catch (...) {  //!<
    rethrow(EXCEPTION_MESSAGE(""));
}

Volume::Volume(std::vector<Layer> &&layers) try : layers_(std::move(layers)) {
    if (layers_.size() == 0) {
        throw std::invalid_argument(
            "Vector cannot be empty to Volume constructor.");
    }

    auto prev_depth = layers_.at(0).get_depth();
    for (const auto &layer : layers_) {
        if (layer.get_depth() < prev_depth) {
            throw std::invalid_argument(EXCEPTION_MESSAGE(
                "Provided vector of layers contains a layer which has a "
                "smaller "
                "depth than the previous layer. Depth: " +
                to_string(layer.get_depth()) +
                ". Previous depth: " + to_string(prev_depth)));
        }
    }
} catch (...) {
    rethrow();
}

Simulation::Settings::Settings(
    energy_quantity electron_stopping_energy, size_t z_number,
    size_t mass_number, bool enable_damage_cascade,
    energy_quantity ion_stopping_energy,
    energy_quantity ion_displacement_energy, bool log_single_displacement,
    size_t divisor_angle_number, size_t flying_distance_number,
    length_quantity range, size_t bombardment_count, bool is_electron,
    energy_quantity incident_energy, size_t thread_count) try
    : electron_stopping_energy_(electron_stopping_energy),
      z_number_(z_number),
      mass_number_(mass_number),
      incident_particle_properties_(
          Particle::Properties(z_number, mass_number)),
      enable_damage_cascade_(enable_damage_cascade),
      ion_stopping_energy_(ion_stopping_energy),
      ion_displacement_energy_(ion_displacement_energy),
      log_single_displacement_(log_single_displacement),
      divisor_angle_number_(divisor_angle_number),
      flying_distance_number_(flying_distance_number),
      range_(range),
      bombardment_count_(bombardment_count),
      is_electron_(is_electron),
      incident_energy_(incident_energy),
      thread_count_(thread_count) {
    if (incident_energy <= energy_quantity(0.0 * kilo_electron_volt)) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Incident energy should be greater than 0. Value: " +
            to_string(incident_energy)));
    }
    if (bombardment_count == 0) {
        throw std::invalid_argument(
            EXCEPTION_MESSAGE("Bombardment count must be greater than 0"));
    }
    if (divisor_angle_number == 0) {
        throw std::invalid_argument(
            EXCEPTION_MESSAGE("Divisor angle number must be greater than 0"));
    }
    if (ion_displacement_energy < energy_quantity(0.0 * kilo_electron_volt)) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Ion displacement energy cannot be less than 0. Value: " +
            to_string(ion_displacement_energy)));
    }
    if (range < length_quantity(0.0 * angstrom)) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Range cannot be less than 0. Value: " + to_string(range)));
    }
    if (flying_distance_number == 0) {
        throw std::invalid_argument(
            EXCEPTION_MESSAGE("Flying distance number must be greater than 0"));
    }
    if (electron_stopping_energy < energy_quantity(0.0 * kilo_electron_volt)) {
        throw std::invalid_argument(EXCEPTION_MESSAGE(
            "Electron stopping energy cannot be less than 0. Value: " +
            to_string(electron_stopping_energy)));
    }
    if (thread_count == 0) {
        throw std::invalid_argument(
            EXCEPTION_MESSAGE("Thread count must be greater than 0"));
    }
} catch (...) {
    rethrow();
}

Simulation::Simulation(const Settings settings, Volume &&volume) try
    : settings_(settings),
      volume_(std::move(volume)),
      thread_pool_(settings.thread_count_) {
    std::vector<Bombardment> bombardments;
    bombardments.reserve(settings.bombardment_count_);
    for (size_t i = 0; i < settings.bombardment_count_; i++) {
        bombardments.push_back(Bombardment(*this, i));
    }
    bombardments_ = std::move(bombardments);

} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

void Simulation::run_simulation() try {
    LOG(INFO) << "Starting Simulation...\n"
              << "Settings:\n"
              << Simulation::print_settings();
    for (auto &bombardment : bombardments_) {
        // boost::asio::post(thread_pool_, [&bombardment]() {
        LOG(INFO) << "Starting Bombardment...";
        bombardment.run_bombardment();
        LOG(INFO) << "...Bombardment Complete!";
        // });
    }
    thread_pool_.join();
    LOG(INFO) << "...Simulation complete";
} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

std::string Simulation::print_settings() const try {
    std::ostringstream ss;
    ss << "Electron stopping energy: " << settings_.electron_stopping_energy_
       << std::endl;
    ss << "Incident particle z number: " << settings_.z_number_ << std::endl;
    ss << "Incident particle mass number: " << settings_.mass_number_
       << std::endl;
    ss << "Incident particle properties: " << std::endl;
    ss << "\tMass: " << settings_.incident_particle_properties_.mass_
       << std::endl;
    ss << "\tCharge: " << settings_.incident_particle_properties_.charge_
       << std::endl;
    ss << "Enable damage cascade? "
       << (settings_.enable_damage_cascade_ ? "Yes" : "No") << std::endl;
    ss << "Ion stopping energy: " << settings_.ion_stopping_energy_
       << std::endl;
    ss << "Ion displacement energy: " << settings_.ion_displacement_energy_
       << std::endl;
    ss << "Log single displacement? "
       << (settings_.log_single_displacement_ ? "Yes" : "No") << std::endl;
    ss << "Divisor angle number: " << settings_.divisor_angle_number_
       << std::endl;
    ss << "Flying distance number: " << settings_.flying_distance_number_
       << std::endl;
    ss << "Range: " << settings_.range_ << std::endl;
    ss << "Bombardment count: " << settings_.bombardment_count_ << std::endl;
    ss << "Is electron? " << (settings_.is_electron_ ? "Yes" : "No")
       << std::endl;
    ss << "Incident energy: " << settings_.incident_energy_ << std::endl;
    ss << "Thread count: " << settings_.thread_count_ << std::endl;
    return ss.str();
} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

Bombardment::Bombardment(const Simulation &simulation, size_t id) try
    : simulation_(simulation),
      random_number_generator_(),
      uniform_distribution_(),
      id_(id) {
    std::cout << "Attempting random number: "
              << uniform_distribution_(random_number_generator_) << std::endl;

} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

void Bombardment::run_bombardment() try {
    if (simulation_.get_settings(&Simulation::Settings::is_electron_)) {
        incident_particle_ = std::make_unique<Electron>(
            simulation_.get_settings(&Simulation::Settings::z_number_),
            simulation_.get_settings(&Simulation::Settings::mass_number_),
            simulation_, uniform_distribution_, random_number_generator_);
    } else {
        incident_particle_ = std::make_unique<Ion>(
            simulation_.get_settings(&Simulation::Settings::z_number_),
            simulation_.get_settings(&Simulation::Settings::mass_number_),
            simulation_, uniform_distribution_, random_number_generator_);
    }

    incident_particle_->fire();

} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

Particle::Particle(
    size_t z_number, size_t mass_number, const Simulation &simulation,
    const boost::random::uniform_01<double> &uniform_distribution,
    boost::random::mt19937 &random_number_generator) try
    : properties_(Particle::Properties(z_number, mass_number)),
      particles_(std::vector<std::unique_ptr<Particle>>()),
      coordinates_(std::vector<Coordinate>()),
      simulation_(simulation),
      velocity_(Velocity(
          0, 0,
          simulation.get_settings(&Simulation::Settings::incident_energy_))),
      coordinate_(Coordinate(0, 0, 0)),
      uniform_distribution_(uniform_distribution),
      random_number_generator_(random_number_generator) {
} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

void Electron::fire() try {
    // TODO: Implement Electron simulation
} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

void Ion::fire() try {
    const auto ion_stopping_energy =
        simulation_.get_settings(&Simulation::Settings::ion_stopping_energy_);
    auto current_layer_pair =
        simulation_.get_volume().get_layer(coordinate_.z_);
    const Layer &current_layer = *(current_layer_pair.first);

    auto comp_view = current_layer.get_relative_compositions().begin();

    // TODO: Remove
    for (size_t i = 0; i < current_layer.get_number_of_components(); i++) {
        std::cout << "Index: " << i << std::endl;
        std::cout << "Relative (direct): "
                  << current_layer.get_relative_composition(i) << std::endl;
        std::cout << "Relative (view): " << *comp_view << std::endl;
        comp_view++;
    }

    while (velocity_.energy_ > ion_stopping_energy) {
        // TODO: Implement simulation calculations
        //
        // 1. Subtract electronic stopping energy
        velocity_.energy_ =
            electronic_stopping_energy(properties_.charge_, properties_.mass_,
                                       velocity_.energy_, current_layer);
        // 2. Calculate recoil energy and velocity
        // 3. Calculate new coordinate (push old coordinate to vector)
        // 4. Determine if sputtering occured
        // 5. Check if damage cascade occurs
        //    a. If so, if energy difference is larger than displacement energy,
        //    create a cascade particle and push it to the particles_ vector
    }
} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

length_quantity Ion::interatomic_spacing(
    number_density_quantity number_density) const {
    return root<-3>(number_density * constants::N_A);
};

energy_quantity Ion::electronic_stopping_energy(charge_quantity charge,
                                                mass_quantity mass,
                                                energy_quantity energy,
                                                const Layer &layer) const try {
    const size_t atom_index =
        layer.get_discrete_distribution()(random_number_generator_);
    std::cout << "Index: " << atom_index << std::endl;
    const auto properties = layer.get_property(atom_index);
    const auto c_1 = 1.212 * root<2>(constants::m_u) * root<2>(electron_volt) /
                     angstrom / pow<static_rational<7, 6>>(constants::e);
    const auto k_l = c_1 * pow<static_rational<7, 6>>(charge) *
                     properties.charge_ /
                     (pow<static_rational<3, 2>>(
                          pow<static_rational<2, 3>>(charge) +
                          pow<static_rational<2, 3>>(properties.charge_)) *
                      root<2>(mass));
    const auto stopping_energy = k_l * root<2>(energy);
    const auto c_2 = 1.59 * pow<3>(angstrom);
    const auto travel_length = interatomic_spacing(layer.get_number_density());
    const auto energy_loss = c_2 * travel_length * layer.get_number_density() *
                             constants::N_A * stopping_energy;
    return energy - energy_loss;
} catch (...) {
    rethrow(EXCEPTION_MESSAGE(""));
}

}  // namespace common
}  // namespace amcset
