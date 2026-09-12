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

// Includes{{{

#include "amcset_common.h"

#include <glog/logging.h>

#include <boost/math/constants/constants.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/systems/si/amount.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/systems/si/energy.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass.hpp>

#include <algorithm>
// #include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
// #include <limits>
#include <numeric>
#include <random>
#include <ranges>
#include <sstream>
#include <stdexcept>

#include "amcset_utilities.h" // }}}

namespace amcset {
namespace common {

// Layer function implementations{{{
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

        for (auto &pair : material) {
          if (pair.first < 0) {
            throw std::invalid_argument(EXCEPTION_MESSAGE(
                "Relative composition: " + std::to_string(pair.first) +
                " is less than zero."));
          }
          VLOG(3) << "Initial composition is " << pair.first;
          pair.first /= sum;
          VLOG(3) << "Final composition is "
                  << pair.first;        }

        return std::move(material);
      }()),
      depth_(depth), mass_density_(mass_density), number_density_([&]() {
        VLOG(2) << "Calculating mass density using " << mass_density;
        return mass_density /
               (std::accumulate(material_.begin(), material_.end(),
                                mass_quantity(0),
                                [](const mass_quantity &sum, const auto &mat) {
                                  VLOG(3)
                                      << "Number density calculation. Sum of "
                                      << sum << ", fraction " << mat.first
                                      << ", and mass of " << mat.second.mass_;
                                  return sum + mat.first * mat.second.mass_;
                                }) *
                constants::N_A);
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

auto Layer::get_relative_compositions() const
    -> decltype(std::declval<const material_vector &>() |
                std::views::keys) try {
  VLOG(2) << "Getting relative compositions for layer with " << material_.size()
          << " components";
  // return material_ |
  //        std::views::transform(&std::pair<double,
  //        Particle::Properties>::first);
  return material_ | std::views::keys;
} catch (...) { //!<
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Volume function implementations{{{
Volume::Volume(std::vector<Layer> &&layers) try : layers_(std::move(layers)) {
  VLOG(1) << "Created volume with " << layers_.size() << " layers.";
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

  // Rely on the ordering of the layers (depths are strictly increasing)
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
      EXCEPTION_MESSAGE("End of layer vector vector reached. Maximum depth: " +
                        to_string(layers_.back().get_depth()) +
                        ". Provided depth: " + to_string(depth)));
} catch (...) {
  rethrow();
} // }}}

// Simulation function implementations{{{
Simulation::Settings::Settings(
    energy_quantity electron_stopping_energy, size_t z_number,
    size_t mass_number, bool enable_damage_cascade,
    energy_quantity ion_stopping_energy,
    energy_quantity ion_displacement_energy, bool log_single_displacement,
    size_t divisor_angle_number, size_t flying_distance_number,
    length_quantity range, size_t bombardment_count, bool is_electron,
    energy_quantity incident_energy, size_t thread_count,
    angle_quantity incident_angle, size_t displacement_bin_count,
    size_t electronic_energy_loss_bin_count,
    size_t nuclear_energy_loss_bin_count,
    size_t total_energy_loss_bin_count,
    std::string displacement_histogram_output_file,
    std::string electronic_energy_loss_histogram_output_file,
    std::string nuclear_energy_loss_histogram_output_file,
    std::string total_energy_loss_histogram_output_file,
    bool log_particle_tree, std::string particle_tree_output_file) try
    : electron_stopping_energy_(electron_stopping_energy), z_number_(z_number),
      mass_number_(mass_number),
      incident_particle_properties_(
          Particle::Properties(z_number, mass_number)),
      enable_damage_cascade_(enable_damage_cascade),
      ion_stopping_energy_(ion_stopping_energy),
      ion_displacement_energy_(ion_displacement_energy),
      log_single_displacement_(log_single_displacement),
      divisor_angle_number_(divisor_angle_number),
      flying_distance_number_(flying_distance_number), range_(range),
      bombardment_count_(bombardment_count), is_electron_(is_electron),
      incident_energy_(incident_energy), thread_count_(thread_count),
      incident_angle_(incident_angle),
      displacement_bin_count_(displacement_bin_count),
      electronic_energy_loss_bin_count_(electronic_energy_loss_bin_count),
      nuclear_energy_loss_bin_count_(nuclear_energy_loss_bin_count),
      total_energy_loss_bin_count_(total_energy_loss_bin_count),
      displacement_histogram_output_file_(std::move(displacement_histogram_output_file)),
      electronic_energy_loss_histogram_output_file_(
          std::move(electronic_energy_loss_histogram_output_file)),
      nuclear_energy_loss_histogram_output_file_(
          std::move(nuclear_energy_loss_histogram_output_file)),
      total_energy_loss_histogram_output_file_(
          std::move(total_energy_loss_histogram_output_file)),
      log_particle_tree_(log_particle_tree),
      particle_tree_output_file_(std::move(particle_tree_output_file)) {
  if (incident_energy <= energy_quantity(0.0 * kilo_electron_volt)) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Incident energy should be greater than 0. Value: " +
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
  if (displacement_bin_count == 0 || electronic_energy_loss_bin_count == 0 ||
      nuclear_energy_loss_bin_count == 0 || total_energy_loss_bin_count == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Histogram bin counts must be greater than 0"));
  }
} catch (...) {
  rethrow();
}

Simulation::Simulation(const Settings settings, Volume &&volume,
                       DrawFunction draw_function) try
    : settings_(settings), volume_(std::move(volume)),
      tallies_(100, settings.range_.value(), settings.displacement_bin_count_,
               settings.electronic_energy_loss_bin_count_,
               settings.nuclear_energy_loss_bin_count_,
               settings.total_energy_loss_bin_count_),
      thread_pool_(settings.thread_count_), draw_function_(draw_function) {
  LOG(INFO) << "Creating simulation with settings:\n"
            << Simulation::print_settings(); // TODO: Print volume
                                             // information as well.
  std::vector<Bombardment> bombardments;
  bombardments.reserve(settings.bombardment_count_);

  for (size_t i = 0; i < settings.bombardment_count_; i++) {
    VLOG(1) << "Creating bombardment " << i;
    bombardments.push_back(Bombardment(*this, i));
  }
  bombardments_ = std::move(bombardments);
  LOG(INFO) << settings.bombardment_count_ << " Bombardments Created";

  LOG(INFO) << "Tallies created";

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

int Simulation::run_simulation() try {
  LOG(INFO) << "Starting Simulation...";
  for (auto &bombardment : bombardments_) {
    // TODO: Enable multithreading
    boost::asio::post(thread_pool_, [this, &bombardment]() {
      VLOG(1) << "Starting Bombardment...";
      bombardment.run_bombardment();
      increment_completed_bombardment_count();
      VLOG(1) << "...Bombardment Complete! Completed: "
              << get_completed_bombardment_count() << "/"
              << settings_.bombardment_count_;
    });
  }
  thread_pool_.join();
  LOG(INFO) << "...Simulation complete";
  // TODO: Come up with a more elegant solution to this.
  LOG(INFO) << "Saving simulation data...";
  try {
    const std::filesystem::path output_root = "output";
    std::filesystem::create_directories(output_root);
    std::filesystem::remove(output_root / "collisions.csv");

    auto write_histogram_with_log = [this](const std::string &label,
                                          const std::string &path,
                                          const Histogram &histogram,
                                          const std::string &value_label) {
      const std::filesystem::path output_path(path);
      if (!output_path.parent_path().empty()) {
        std::filesystem::create_directories(output_path.parent_path());
      }
      LOG(INFO) << "Writing " << label << " histogram to: "
                << output_path.string();
      write_histogram_csv(path, histogram, value_label);
    };

    write_histogram_with_log("displacement",
                             settings_.displacement_histogram_output_file_,
                             get_tallies().displacement_histogram, "count");
    write_histogram_with_log("electronic loss",
                             settings_.electronic_energy_loss_histogram_output_file_,
                             get_tallies().electronic_energy_loss_histogram,
                             "energy_loss_eV");
    write_histogram_with_log("nuclear loss",
                             settings_.nuclear_energy_loss_histogram_output_file_,
                             get_tallies().nuclear_energy_loss_histogram,
                             "energy_loss_eV");
    write_histogram_with_log("total loss",
                             settings_.total_energy_loss_histogram_output_file_,
                             get_tallies().total_energy_loss_histogram,
                             "energy_loss_eV");

    if (settings_.log_particle_tree_) {
      write_particle_tree_log();
    }

    const std::filesystem::path tallies_path = output_root / "tallies.json";
    LOG(INFO) << "Writing tallies summary to: " << tallies_path.string();
    std::ofstream tally_stream(tallies_path);
    if (!tally_stream.is_open()) {
      throw std::runtime_error("Failed to open tallies output: " + tallies_path.string());
    }
    tally_stream << get_tallies().to_json();
  } catch (const std::exception &e) {
    VLOG(1) << "Failed to write simulation output: " << e.what();
  } catch (...) {
    VLOG(1) << "Failed to write simulation output";
  }

  return 0;
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
  ss << "Ion stopping energy: " << settings_.ion_stopping_energy_ << std::endl;
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
  ss << "Is electron? " << (settings_.is_electron_ ? "Yes" : "No") << std::endl;
  ss << "Incident energy: " << settings_.incident_energy_ << std::endl;
  ss << "Thread count: " << settings_.thread_count_ << std::endl;
  ss << "Displacement histogram bins: "
     << settings_.displacement_bin_count_ << std::endl;
  ss << "Electronic loss histogram bins: "
     << settings_.electronic_energy_loss_bin_count_ << std::endl;
  ss << "Nuclear loss histogram bins: "
     << settings_.nuclear_energy_loss_bin_count_ << std::endl;
  ss << "Total loss histogram bins: "
     << settings_.total_energy_loss_bin_count_ << std::endl;
  ss << "Log particle tree? "
     << (settings_.log_particle_tree_ ? "Yes" : "No") << std::endl;
  return ss.str();
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

void Simulation::merge_tallies(const Tallies &tallies) {
  std::lock_guard<std::mutex> lk(tallies_mutex_);
  tallies_.merge(tallies);
}

const Tallies &Simulation::get_tallies() const noexcept { return tallies_; }

void Simulation::write_histogram_csv(const std::string &output_file,
                                    const Histogram &histogram,
                                    const std::string &value_label) const {
  if (output_file.empty()) {
    throw std::invalid_argument("Histogram output file path cannot be empty.");
  }

  const std::filesystem::path output_path(output_file);
  if (!output_path.parent_path().empty()) {
    std::filesystem::create_directories(output_path.parent_path());
  }

  std::ofstream output(output_path);
  if (!output.is_open()) {
    throw std::runtime_error("Failed to open histogram output: " + output_file);
  }

  const auto left_to_angstrom = [](double meters) {
    return (length_quantity(meters * si::meter) / angstrom).value();
  };

  const auto value_name = value_label.empty() ? std::string("count") : value_label;
  output << "depth_min_angstrom,depth_max_angstrom," << value_name << "\n";
  const auto width = histogram.max - histogram.min;
  if (histogram.bins.empty()) {
    output << "0,0,0\n";
    return;
  }
  for (size_t i = 0; i < histogram.bins.size(); ++i) {
    const auto left = histogram.min + width * double(i) / double(histogram.bins.size());
    const auto right = histogram.min + width * double(i + 1) / double(histogram.bins.size());
    const auto value = value_name == "energy_loss_eV" ? histogram.bins[i] : histogram.bins[i];
    output << left_to_angstrom(left) << "," << left_to_angstrom(right) << "," << value << "\n";
  }
  output.flush();
  LOG(INFO) << "Histogram written to " << output_path.string() << " with "
            << histogram.bins.size() << " bins";
}

void Simulation::write_particle_tree_log() const {
  if (!settings_.log_particle_tree_) return;
  std::filesystem::create_directories(
      std::filesystem::path(settings_.particle_tree_output_file_).parent_path());
  std::ofstream output(settings_.particle_tree_output_file_);
  output << "{\n  \"bombardments\": [\n";
  for (size_t i = 0; i < bombardments_.size(); ++i) {
    if (i != 0) output << ",\n";
    output << "    {\n";
    output << "      \"bombardment_id\": " << i << ",\n";
    output << "      \"tree\": " << bombardments_[i].get_data_leaf().to_json(8)
           << "\n    }";
  }
  output << "\n  ]\n}\n";
}

// Bombardment tally helpers
void Bombardment::record_incident(double range_meters) {
  local_tallies_.incidents += 1;
  local_tallies_.range_stats.add((long double)range_meters);
  local_tallies_.depth_histogram.add(range_meters);
}

void Bombardment::record_displacement(double depth_meters) {
  local_tallies_.total_displacements += 1;
  local_tallies_.depth_histogram.add(depth_meters);
  local_tallies_.displacement_histogram.add(depth_meters);
}

void Bombardment::record_sputtered() { local_tallies_.sputtered += 1; }

void Bombardment::record_energy_loss(long double energy_joules, double depth_meters,
                                     EnergyLossType type) {
  local_tallies_.energy_loss_stats.add(energy_joules);
  const auto energy_eV = (energy_quantity(static_cast<double>(energy_joules) * si::joule) / electron_volt).value();
  if (type == EnergyLossType::electronic) {
    local_tallies_.electronic_energy_loss_histogram.add_value(depth_meters, energy_eV);
  } else if (type == EnergyLossType::nuclear) {
    local_tallies_.nuclear_energy_loss_histogram.add_value(depth_meters, energy_eV);
  } else {
    local_tallies_.total_energy_loss_histogram.add_value(depth_meters, energy_eV);
  }
}

// Bombardment function implementations{{{
Bombardment::Bombardment(Simulation &simulation, size_t id) try
    : simulation_(simulation),
      local_tallies_(
          100, simulation.get_settings(&Simulation::Settings::range_).value()),
      random_number_generator_(EntropySeeder::generate()),
      uniform_distribution_(), id_(id),
      data_leaf_() { // TODO: Provide the option to set a
  // seed instead of random generation.
  VLOG(1) << "Attempting random number: "
          << uniform_distribution_(random_number_generator_);

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

void Bombardment::run_bombardment() try {
  if (simulation_.get_settings(&Simulation::Settings::is_electron_)) {
    VLOG(1) << "Simulating electron bombardment...";
    incident_particle_ = std::make_unique<Electron>(
        simulation_.get_settings(&Simulation::Settings::z_number_),
        simulation_.get_settings(&Simulation::Settings::mass_number_),
        simulation_, uniform_distribution_, random_number_generator_,
        get_data_leaf(), *(this));
  } else {
    VLOG(1) << "Simulating ion bombardment";
    incident_particle_ = std::make_unique<Ion>(
        simulation_.get_settings(&Simulation::Settings::z_number_),
        simulation_.get_settings(&Simulation::Settings::mass_number_),
        simulation_, uniform_distribution_, random_number_generator_,
        get_data_leaf(), *(this));
  }

  VLOG(1) << "Firing incident particle with velocity quaternion (" << incident_particle_->get_velocity().quaternion_.r_ << ", "
            << incident_particle_->get_velocity().quaternion_.x_ << ", "
            << incident_particle_->get_velocity().quaternion_.y_ << ", "
            << incident_particle_->get_velocity().quaternion_.z_ << ")."
            << " and energy " << incident_particle_->get_velocity().energy_;
  incident_particle_->fire();
  std::this_thread::yield(); // TODO: Remove
  VLOG(1) << "...Bombardment complete.";

  // Simple local tally collection: record final range and counts
  try {
    auto &coords = incident_particle_->get_coordinates();
    if (!coords.empty()) {
      auto last_z = coords.back().z_.value();
      // record via Bombardment helper
      record_incident(last_z);
    }
    // Note: individual displacements and sputters are recorded during
    // the particle history via Bombardment::record_displacement/_sputtered
  } catch (...) {
    VLOG(1) << "Failed to collect local tallies";
  }

  // Merge local tallies into simulation-level accumulator
  simulation_.merge_tallies(local_tallies_);

  // Post-bombardment processing.
  // This goes through the incident_particle_ data_leaf_ and
  // looks for knock-ons through the list of pointers
  // contained within the ParticleDataLeaf class. It is a tree-structure,
  // hence why it is called a "leaf".
  // TODO:Remove this comment


} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

Particle::Particle(size_t z_number, size_t mass_number,
                   const Simulation &simulation,
                   boost::random::uniform_01<double> &uniform_distribution,
                   boost::random::mt19937 &random_number_generator,
                   ParticleDataLeaf &data_leaf, Bombardment &bombardment) try
    : properties_(Particle::Properties(z_number, mass_number)),
      coordinates_(std::vector<Coordinate>()), simulation_(simulation),
      velocity_(Velocity(
          Quaternion(
              0.0,
              1.0 * sin(simulation_.get_settings(
                        &Simulation::Settings::incident_angle_)),
              0.0,
              1.0 * cos(simulation_.get_settings(
                        &Simulation::Settings::
                            incident_angle_))), // Point straight towards
                                                // z-axis. TODO: Allow user to
                                                // specify an incident angle
          simulation.get_settings(&Simulation::Settings::incident_energy_))),
      coordinate_(Coordinate(angstrom * 0.0, angstrom * 0.0, angstrom * 0.0)),
      uniform_distribution_(uniform_distribution),
      random_number_generator_(random_number_generator), data_leaf_(data_leaf),
      bombardment_(bombardment) {
  VLOG(1) << "Creating particle with velocity (" << velocity_.quaternion_.r_
          << ", " << velocity_.quaternion_.z_ << ", "
          << velocity_.quaternion_.y_ << ", " << velocity_.quaternion_.z_
          << ").";
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

velocity_quantity Particle::speed_from_energy(energy_quantity energy,
                                              mass_quantity mass) const {
  return boost::units::sqrt(2.0 * energy / mass);
}

mass_quantity Particle::reduced_mass(mass_quantity mass_1,
                                     mass_quantity mass_2) const {
  return mass_1 * mass_2 / (mass_1 + mass_2);
}

energy_quantity Particle::cm_energy(mass_quantity reduced_mass,
                                    velocity_quantity velocity) const {
  return reduced_mass * boost::units::pow<2>(velocity) / 2.0;
}

// Ensure this does not fail with quaternions aligned with axes, nan components,
// and zero-length vector.
Particle &Particle::rotate_inplace(angle_quantity scattering_angle) {
  VLOG(2) << "Beginning rotation with scattering_angle " << scattering_angle
          << " velocity (" << velocity_.quaternion_.r_ << ", "
          << velocity_.quaternion_.x_ << ", " << velocity_.quaternion_.y_
          << ", " << velocity_.quaternion_.z_ << ")";

  auto x = velocity_.quaternion_.x_;
  auto y = velocity_.quaternion_.y_;
  auto z = velocity_.quaternion_.z_;

  if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
    // throw std::invalid_argument("Input vector contains NaNs");
    // HACK: Sometimes nans are in the array. Redirect in the z direction
    // if so.
    LOG(ERROR) << "Nans detected in bombardment. Position (" << x << ", " << y
               << ", " << z
               << "). Energy: " << velocity_.energy_ / electron_volt << " eV"
               << ". Bombardment ID: " << bombardment_.get_id();
    x = 0.0;
    y = 0.0;
    z = 1.0;
  }

  auto length = std::sqrt(x * x + y * y + z * z);
  if (length == 0.0) {
    throw std::invalid_argument("Input vector has zero length");
  }

  Quaternion rotation_quaternion;

  // Special cases: TODO: Add [[likely]] and [[unlikely]] optimizations
  if (x == 0 && y == 0) { // z-axis aligned
    rotation_quaternion = Quaternion(0.0, 0.0, 1.0, 0.0);
  } else if (x == 0 && z == 0) { // y-axis aligned
    rotation_quaternion = Quaternion(0.0, 1.0, 0.0, 0.0);
  } else if (y == 0.0 && z == 0.0) { // x-axis aligned
    rotation_quaternion = Quaternion(0.0, 0.0, 0.0, 1.0);
  } else { // Non-axis aligned, manual alignment instead.
    rotation_quaternion = Quaternion(0.0, y, x,
                                     0.0); // Direction is not important,
                                           // rotation has random alpha tilt.
  }

  const auto alpha_angle =
      radian * 2.0 * pi *
      (0.5 - uniform_distribution_(random_number_generator_));

  auto alpha_scalar = cos(alpha_angle / 2.0);
  auto alpha_pos = sin(alpha_angle / 2.0);
  auto quaternion_alpha =
      Quaternion(alpha_scalar, x * alpha_pos, y * alpha_pos, z * alpha_pos);

  // TODO: Try pythagorean identity to make more efficient.
  auto beta_angle = scattering_angle;

  auto beta_scalar = cos(beta_angle / 2.0);
  auto beta_pos = sin(beta_angle / 2.0);
  auto quaternion_beta = Quaternion(
      beta_scalar, rotation_quaternion.x_ * beta_pos,
      rotation_quaternion.y_ * beta_pos, rotation_quaternion.z_ * beta_pos);

  velocity_.quaternion_ =
      (quaternion_alpha * (quaternion_beta * velocity_.quaternion_) *
       quaternion_beta.conjugate_copy()) *
      quaternion_alpha.conjugate_copy();
  velocity_.quaternion_.r_ = 0.0;
  velocity_.quaternion_.normalize_inplace();

  VLOG(2) << "Post-rotation position-only normalized result quaternion ("
          << velocity_.quaternion_.r_ << ", " << velocity_.quaternion_.x_
          << ", " << velocity_.quaternion_.y_ << ", "
          << velocity_.quaternion_.z_ << ")";

  return *this;
}

dimensionless_quantity
Ion::screening_function(dimensionless_quantity reduced_radius) const {
  const auto result = 0.1818 * exp(-3.2 * reduced_radius) +
                      0.5099 * exp(-0.9423 * reduced_radius) +
                      0.2802 * exp(-0.4028 * reduced_radius) +
                      0.2817 * exp(-0.2016 * reduced_radius);
  VLOG(2) << "Screening function with reduced radius " << reduced_radius
          << " is " << result;
  return result;
}

// TODO: Find a way to cache the results of this function for speedup. Consider
// calculating ahead of time and storing the results in a data structure
// constructed at runtime for each Layer material.
length_quantity Ion::screening_length(dimensionless_quantity z_1,
                                      dimensionless_quantity z_2) const {
  const auto result = 0.25 * cbrt(9.0 * pi * pi / 2.0) * bohr_radius /
                      (pow(z_1, 0.23) + pow(z_2, 0.23));
  VLOG(2) << "Screening length with z_1 " << z_1 << " and z_2 " << z_2 << " is "
          << result;
  return result;
}

energy_quantity Ion::screening_potential(length_quantity radius,
                                         dimensionless_quantity z_1,
                                         dimensionless_quantity z_2) const {
  auto reduced_length = radius / screening_length(z_1, z_2);
  auto result = screening_function(reduced_length) * z_1 * z_2 * e_statcoulomb *
                e_statcoulomb / radius;
  VLOG(2) << "Screening potential with radius " << radius << ", z_1 " << z_1
          << ", and z_2 " << z_2 << " is " << result;
  return result;
}

voltage_quantity
Ion::screening_potential_derivative(length_quantity radius,
                                    dimensionless_quantity z_1,
                                    dimensionless_quantity z_2) const {
  // TODO: Consolidate calls to screening_length to above-mentioned method.
  // Optimize this function as well.
  auto s_length = screening_length(z_1, z_2);
  auto reduced_radius = radius / s_length;
  auto result =
      (0.1818 * exp(-3.2 * reduced_radius) * (-3.2 / s_length) +
       0.5099 * exp(-0.9423 * reduced_radius) * (-0.9423 / s_length) +
       0.2802 * exp(-0.4028 * reduced_radius) * (-0.4028 / s_length) +
       0.2817 * exp(-0.2016 * reduced_radius) * (-0.2016 / s_length)) *
          z_1 * z_2 * e_statcoulomb * e_statcoulomb / radius -
      (0.1818 * exp(-3.2 * reduced_radius) +
       0.5099 * exp(-0.9423 * reduced_radius) +
       0.2802 * exp(-0.4028 * reduced_radius) +
       0.2817 * exp(-0.2016 * reduced_radius)) *
          z_1 * z_2 * e_statcoulomb * e_statcoulomb / (radius * radius);
  VLOG(2) << "Screening potential derivative with radius " << radius << ", z_1 "
          << z_1 << ", and z_2 " << z_2 << " is " << result;
  return result;
}

dimensionless_quantity Ion::reduced_energy(length_quantity screening_length,
                                           energy_quantity cm_energy,
                                           dimensionless_quantity z_1,
                                           dimensionless_quantity z_2) const {
  auto result = screening_length * cm_energy /
                (z_1 * z_2 * e_statcoulomb * e_statcoulomb);
  VLOG(2) << "Reduced energy with screening length " << screening_length
          << ", center of mass energy " << cm_energy << ", z_1 " << z_1
          << ", z_2 " << z_2 << ", result is " << result;
  return result;
}

// TODO: Verify that dividing by Avogadro's number is correct.
length_quantity
Ion::free_flying_path_length(mass_quantity m_1, mass_quantity m_2,
                             dimensionless_quantity reduced_energy,
                             length_quantity screening_length,
                             number_density_quantity number_density) const {
  length_quantity L;
  // if (reduced_energy > 100) { // TODO: Figure out this (odd behavior at high
  // E)
  //   L = (0.02 * pow(1.0 + (m_1 + m_2) / atomic_mass_unit, 2) *
  //   reduced_energy *
  //            reduced_energy +
  //        0.1 * pow(reduced_energy, 1.38)) /
  //       (4.0 * pi * screening_length * screening_length * number_density *
  //        log(1.0 + reduced_energy)) /
  //       si::constants::codata::N_A;
  //   VLOG(3) << "High energy free flying path length with m_1 " << m_1 << ",
  //   m_2"
  //           << m_2 << ", reduced_energy " << reduced_energy
  //           << ", screening_length " << screening_length
  //           << ", and number_density " << number_density << " with result "
  //           << L;
  // } else {
  L = pow<static_rational<1, 3>>(1.0 / number_density /
                                 si::constants::codata::N_A);
  VLOG(3) << "Low energy free flying path length with m_1 " << m_1 << ", m_2 "
          << m_2 << ", reduced_energy " << reduced_energy
          << ", screening_length " << screening_length
          << ", and number_density " << number_density << " with result " << L;
  // }
  return L;
}

length_quantity Ion::collision_diameter(dimensionless_quantity z_1,
                                        dimensionless_quantity z_2,
                                        mass_quantity rd_mass,
                                        velocity_quantity speed) const {
  auto result =
      4 * z_1 * z_2 * e_statcoulomb * e_statcoulomb / (rd_mass * speed * speed);
  VLOG(2) << "Collision diameter given z_1 " << z_1 << ", z_2 " << z_2
          << " reduced mass " << rd_mass << ", and speed " << speed << " is "
          << result;
  return result;
}

// TODO: Fix magnitude error. This code reports closest approach on the
// order of meters. Find out why.
length_quantity Ion::closest_approach(dimensionless_quantity z_1,
                                      dimensionless_quantity z_2,
                                      energy_quantity cm_energy,
                                      length_quantity impact_param) const {
  // Potential function
  auto F = [&](length_quantity closest_approach) {
    return 1 - screening_potential(closest_approach, z_1, z_2) / cm_energy -
           pow<2>(impact_param / closest_approach);
  };
  // Derivative of potential function
  auto dF = [&](length_quantity closest_approach) {
    return -screening_potential_derivative(closest_approach, z_1, z_2) /
               cm_energy +
           2.0 * impact_param * impact_param /
               pow<static_rational<3>>(closest_approach);
  };
  auto old_result = impact_param; // Initial guess, adjust if encountering
                                  // convergence issues
  VLOG(3) << "Closest approach initial guess " << old_result;
  auto old_pot = F(old_result);
  auto result = old_result - old_pot / dF(old_result); // First iteration
  auto pot = F(result);
  // Newtonian iteration
  size_t iterations = 0;
  constexpr size_t max_iterations = 100; // If encountering convergence errors,
                                         // increase. Will decrease performance.
  constexpr double rel_tolerance =
      1e-3; // Threshold for relative tolerance. Increase if encountering
            // convergence issues.
  while (iterations < max_iterations) {
    auto rel_diff = abs(result - old_result) / result;

    if (rel_diff < rel_tolerance) {
      VLOG(2) << "Closest approach found after " << iterations
              << " iterations. Result " << result;
      return result;
    }

    auto d_pot = dF(old_result);
    old_pot = pot;
    pot = F(result);

    // Based on verbose log inspection, this is unlikely to occur.
    if ((old_pot < 0.0 && pot > 0.0) || (old_pot > 0.0 && pot < 0.0)) {
      // Go into binary search (Solves oscillation problem)
      size_t binary_search_iterations = 0;

      while (true) {
        if (iterations > max_iterations) {
          break;
        }
        if (rel_diff < rel_tolerance) {
          VLOG(2) << "Closest approach found after " << iterations << " and "
                  << binary_search_iterations << " binary search iterations."
                  << " iterations. Result " << result;
          rel_diff = abs(result - old_result) / result;
          return result;
        }
        auto mid_result = old_result + (result - old_result) / 2.0;
        auto mid_pot = F(mid_result);

        VLOG(3) << "Binary search old_result " << old_result << " old_pot "
                << old_pot << " mid_result " << mid_result << " mid_pot "
                << mid_pot << " result " << result << " pot " << pot;
        if (copysign(old_pot, mid_pot) == old_pot) {
          // If old_pot and mid_pot have the same sign
          old_result = mid_result; // Search the "right" side
          old_pot = mid_pot;
          VLOG(3) << "Searching pot side";
        } else {
          result = mid_result; // Search the "left" side
          pot = mid_pot;
          VLOG(3) << "Searching old_pot side";
        }
        binary_search_iterations++;
      }
    }

    old_result = result;
    result = old_result - pot / d_pot;
    iterations++;

    VLOG(3) << "Zeroing function value " << pot << ", derivative " << d_pot;
    VLOG(3) << "Newtonian iteration " << iterations << ". Previous result "
            << old_result << ", current result " << result << ". Tolerance "
            << rel_tolerance << ", relative difference " << rel_diff;
  }

  // TODO: Write more graceful code to handle this. And throw an exception to
  // be handled by the GUI if convergence doesn't occur.
  LOG(ERROR) << "Convergence error. Please adjust closest approach parameters. "
                "Previous result "
             << old_result << ", current result " << result;

  return result;
}

length_quantity Ion::radius_of_curvature(length_quantity radius,
                                         energy_quantity cm_energy,
                                         dimensionless_quantity z_1,
                                         dimensionless_quantity z_2) const {
  auto result = 2.0 * abs(cm_energy - screening_potential(radius, z_1, z_2)) /
                (-screening_potential_derivative(radius, z_1, z_2));
  VLOG(2) << "Radius of curvature for radius " << radius << " cm_energy "
          << cm_energy << " z_1 " << z_1 << " z_2 " << z_2 << " has result "
          << result;
  return result;
}

dimensionless_quantity Ion::magic_formula_correction_parameter(
    dimensionless_quantity reduced_energy, length_quantity impact_parameter,
    length_quantity screening_length, length_quantity closest_approach) const {
  auto C_1 = 0.99229;
  auto C_2 = 0.011615;
  auto C_3 = 0.007122;
  auto C_4 = 9.3066;
  auto C_5 = 14.813;

  auto alpha = 1 + C_1 / sqrt(reduced_energy);
  auto beta = (C_2 + sqrt(reduced_energy)) / (C_3 + sqrt(reduced_energy));
  auto gamma = (C_4 + reduced_energy) / (C_5 + reduced_energy);
  auto B = impact_parameter / screening_length;
  auto A = 2 * alpha * reduced_energy * pow(B, beta);
  auto G = gamma / (sqrt(1 + A * A) - A);
  auto R_o = closest_approach / screening_length;
  auto result = A * (R_o - B) / (1 + G);
  VLOG(2) << "Magic formula correction parameter with reduced_energy "
          << reduced_energy << ", impact_parameter " << impact_parameter
          << ", screening_length " << screening_length << " result is "
          << result;
  return result;
}

angle_quantity
Ion::magic_formula_scattering_angle(length_quantity impact_parameter,
                                    length_quantity closest_approach,
                                    length_quantity radius_of_curvature,
                                    dimensionless_quantity correction_factor,
                                    length_quantity screening_distance) const {
  // TODO: Correct radius of curvature naming scheme
  const auto B = impact_parameter / screening_distance;
  const auto R_o = closest_approach / screening_distance;
  const auto R_c = radius_of_curvature / screening_distance;
  const auto result = acos((B + R_c + correction_factor) / (R_o + R_c));
  VLOG(2) << "Magic scattering angle with impact_parameter " << impact_parameter
          << " closest_approach " << closest_approach << " radius_of_curvature "
          << radius_of_curvature << " correction_factor " << correction_factor
          << " screening_distance " << screening_distance << " with result "
          << result;
  return result;
}

energy_quantity
Ion::nuclear_energy_loss(mass_quantity atom_mass, mass_quantity target_mass,
                         energy_quantity incident_energy,
                         angle_quantity scattering_angle) const {
  auto result = 4.0 * atom_mass * target_mass /
                pow<2>(atom_mass + target_mass) * incident_energy *
                pow<2>(sin(scattering_angle));
  VLOG(2) << "Nuclear scattering energy loss atom mass " << atom_mass
          << " target_mass " << target_mass << " incident_energy "
          << incident_energy << " scattering_angle " << scattering_angle
          << " has result " << result;
  return result;
}

angle_quantity
Ion::laboratory_scattering_angle(angle_quantity cm_scattering_angle,
                                 mass_quantity atom_mass,
                                 mass_quantity target_mass) const {
  auto result = atan(sin(cm_scattering_angle) /
                     (cos(cm_scattering_angle) + (atom_mass / target_mass)));

  VLOG(2) << "Laboratory scattering angle with cm_scattering_angle "
          << cm_scattering_angle << " atom_mass " << atom_mass
          << " target_mass " << target_mass << " with result " << result;
  return result;
}

// TODO: Switched to quaternions, likely to remove.{{{
// std::array<angle_quantity, 4>
// Ion::relative_to_absolute_angle(angle_quantity initial_altitude,
//                                 angle_quantity initial_azimuth,
//                                 angle_quantity incident_deflection,
//                                 angle_quantity target_deflection) const {
//   auto THETA1RELATIVE = incident_deflection;
//   auto THETA2RELATIVE = target_deflection;
//   auto ALPHA1RELATIVE =
//       uniform_distribution_(random_number_generator_) * 2.0 * pi;
//   auto ALPHA2RELATIVE = ALPHA1RELATIVE + pi;

//   auto THETAO = initial_altitude;
//   auto ALPHAO = initial_azimuth;

//   // Angle
//   auto X1 = sin(THETA1RELATIVE) * cos(ALPHA1RELATIVE);
//   auto Y1 = sin(THETA1RELATIVE) * sin(ALPHA1RELATIVE);
//   auto Z1 = cos(THETA1RELATIVE);

//   auto Y0 = Y1 * cos(THETAO) + Z1 * sin(THETAO);
//   auto Z0 = -Y1 * sin(THETAO) + Z1 * cos(THETAO);
//   auto X0 = X1;

//   auto Z = Z0;
//   auto X = X0 * sin(ALPHAO) + Y0 * cos(ALPHAO);
//   auto Y = -X0 * cos(ALPHAO) + Y0 * sin(ALPHAO);

//   angle_quantity THETA1;

//   if (Z > 0.0) {
//     THETA1 = atan(sqrt(X * X + Y * Y) / Z);
//   } else if (Z < 0.0) {
//     THETA1 = radian * pi + atan(sqrt(X * X + Y * Y) / Z);
//   } else {
//     THETA1 = radian * pi / 2.0;
//   }

//   angle_quantity ALPHA1;

//   if (sin(THETA1) != 0.0) {
//     if (X > 0.0) {
//       ALPHA1 = atan(Y / X);
//     } else if (X == 0.0) {
//       ALPHA1 = radian * (pi - (Y > 0.0 ? 1.0 : -1.0) * pi / 2.0);
//     } else {
//       ALPHA1 = radian * pi + atan(Y / X);
//     }
//   } else {
//     ALPHA1 = 0.0 * radian;
//   }

//   // Target
//   auto X3 = sin(THETA2RELATIVE) * cos(ALPHA2RELATIVE);
//   auto Y3 = sin(THETA2RELATIVE) * sin(ALPHA2RELATIVE);
//   auto Z3 = cos(THETA2RELATIVE);

//   auto Y2 = Y3 * cos(THETAO) + Z3 * sin(THETAO);
//   auto Z2 = -Y3 * sin(THETAO) + Z3 * cos(THETAO);
//   auto X2 = X3;

//   auto Z5 = Z2;
//   auto X5 = X2 * sin(ALPHAO) + Y2 * cos(ALPHAO);
//   auto Y5 = -X2 * cos(ALPHAO) + Y2 * sin(ALPHAO);

//   angle_quantity THETA2;
//   if (Z5 < 0.0) {
//     THETA2 = radian * pi + atan(sqrt(X5 * X5 + Y5 * Y5) / Z5);
//   } else if (Z5 > 0.0) {
//     THETA2 = atan(sqrt(X5 * X5 + Y5 * Y5) / Z5);
//   } else {
//     THETA2 = radian * pi / 2.0;
//   }

//   angle_quantity ALPHA2;

//   if (sin(THETA2) != 0.0) {
//     if (X5 < 0.0) {
//       ALPHA2 = radian * pi + atan(Y5 / X5);
//     } else if (X5 > 0.0) {
//       atan(Y5 / X5);
//     } else {
//       ALPHA2 = radian * (pi - 0.5 * ((Y5 > 0.0) ? 1 : -1) * pi);
//     }
//   } else {
//     ALPHA2 = radian * 0.0;
//   }
//   // newVelocity.zAngle = THETA1;
//   // newVelocity.xAngle = ALPHA1;
//   // targetVelocity.zAngle = THETA2;
//   // targetVelocity.xAngle = ALPHA2;
//   return {{THETA1, ALPHA1, THETA2, ALPHA2}};
// } // }}}

// Electron function implementations{{{
void Electron::fire() try {
  // TODO: Implement Electron simulation
  throw std::logic_error("Electron simulation is not implemented");
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Ion function implementations // {{{
void Ion::fire() try {
  VLOG(1) << "Firing ion with initial energy " << velocity_.energy_;
  auto ion_stopping_energy = simulation_.get_settings(
      &Simulation::Settings::
          ion_stopping_energy_); // TODO: Make this a user-configurable
                                 // element property. Each element may have a
                                 // different value for this in each material.
  VLOG(2) << "Ion stopping energy: " << ion_stopping_energy;

  // HACK: Temp check for bombardment counting
  size_t incident =
      (velocity_.energy_ ==
       simulation_.get_settings(&Simulation::Settings::incident_energy_))
          ? 0
          : 1;
 
          if(incident == 1) {
            VLOG(2) << "Incident particle is a knock-on particle.";
          } else {
            VLOG(2) << "Incident particle is the primary particle.";
          }

  // #ifdef NDEBUG [>{{{<]
  // #else
  //   // Grab an iterator to a view of the first (Debug code)
  //   auto comp_view = current_layer.get_relative_compositions().begin();
  //   for (size_t i = 0; i < current_layer.get_number_of_components(); i++) {
  //     auto direct_composition = current_layer.get_relative_composition(i);
  //     auto viewed_composition = *(comp_view++);
  //     VLOG(2) << "Index: " << i << " has directly measured composition of "
  //             << direct_composition << " view-based composition of "
  //             << "Relative (view): " << viewed_composition;
  //     assert(direct_composition == viewed_composition);
  //   }
  // #endif [>}}}<]
  // velocity_.energy_ =
  //     simulation_.get_settings(&Simulation::Settings::incident_energy_);
  // velocity_.quaternion_ = Quaternion(0.0, 0.0, 0.0, 1.0);

  // std::cout << "Incident velocity: (" << velocity_.quaternion_.r_ << ", "
  //           << velocity_.quaternion_.x_ << ", " << velocity_.quaternion_.y_
  //           << ", " << velocity_.quaternion_.z_ << ")" << std::endl;

  auto &energy = velocity_.energy_; // A reference is required so that energy
                                    // can changed in-place through the alias.
  VLOG(1) << "Initial energy: " << energy;
  auto max_depth = simulation_.get_volume().get_depth();
  VLOG(2) << "Max depth: " << max_depth;
  auto knock_on_threshold = simulation_.get_settings(
      &Simulation::Settings::
          ion_displacement_energy_); // TODO: Replace this with
                                     // element-specific displacement energy
                                     // (extracted from user-provided layer
                                     // information)
  VLOG(2) << "Knock-on threshold: " << knock_on_threshold;
  auto enable_knock_on =
      simulation_.get_settings(&Simulation::Settings::enable_damage_cascade_);
    VLOG(2) << "Knock-on enabled: " << enable_knock_on;

  while (velocity_.energy_ > ion_stopping_energy) {
    const auto current_layer_pair =
        simulation_.get_volume().get_layer(coordinate_.z_);
        VLOG(3) << "Current layer: " << current_layer_pair.second;
    const Layer &current_layer = *(
        current_layer_pair.first); // TODO: Check for next layer during history
    // Select atom type from current layer
    const size_t atom_index =
        current_layer.get_random_index(random_number_generator_);
    VLOG(3) << "Selected atom index: " << atom_index;
    const auto layer_properties = current_layer.get_property(atom_index);
    auto z_1 = properties_.charge_;
    VLOG(3) << "Ion charge: " << z_1;
    auto z_2 = layer_properties.charge_;
    VLOG(3) << "Target atom charge: " << z_2;
    auto m_1 = properties_.mass_;
    VLOG(3) << "Ion mass: " << m_1;
    auto m_2 = layer_properties.mass_;
    VLOG(3) << "Target atom mass: " << m_2;
    auto number_density = current_layer.get_number_density();
    VLOG(3) << "Layer number density: " << number_density;

    // Calculate electronic stopping energy
    auto e_energy_loss = electronic_stopping_energy(
        z_1, m_1, energy, number_density, layer_properties);
    VLOG(2) << "Electronic stopping energy loss: " << e_energy_loss;
    energy -= e_energy_loss;
    VLOG(2) << "Energy after electronic stopping: " << energy;
    // Record energy loss (electronic) at the actual particle depth.
    bombardment_.record_energy_loss((long double)e_energy_loss.value(),
                                    coordinate_.z_.value(),
                                    Bombardment::EnergyLossType::electronic);

    // Calculate nuclear scattering energy loss (energy, m_1, m_2, z_1, z_2,
    // number_density)
    auto speed = speed_from_energy(energy, m_1);
    VLOG(3) << "Speed from energy: " << speed;
    auto c_mass_energy = cm_energy(m_1, speed);
    VLOG(3) << "Center of mass energy: " << c_mass_energy;
    auto screen_length = screening_length(z_1, z_2);
    VLOG(3) << "Screening length: " << screen_length;
    auto rd_energy = reduced_energy(screen_length, c_mass_energy, z_1, z_2);
    VLOG(3) << "Reduced energy: " << rd_energy;
    auto path_length = free_flying_path_length(m_1, m_2, rd_energy,
                                               screen_length, number_density);
    VLOG(3) << "Free flying path length: " << path_length;
    auto i_param = impact_parameter(number_density, path_length, rd_energy);
    VLOG(3) << "Impact parameter: " << i_param;
    auto closest_appr = closest_approach(z_1, z_2, c_mass_energy, i_param);
    VLOG(3) << "Closest approach: " << closest_appr;
    auto rad_curv = radius_of_curvature(closest_appr, c_mass_energy, z_1, z_2);
    VLOG(3) << "Radius of curvature: " << rad_curv;
    auto magic_correction = magic_formula_correction_parameter(
        rd_energy, i_param, screen_length, closest_appr);
    VLOG(3) << "Magic formula correction parameter: " << magic_correction;
    auto cm_scatter_angle = magic_formula_scattering_angle(
        i_param, closest_appr, rad_curv, magic_correction, screen_length);
    VLOG(3) << "Center of mass scattering angle: " << cm_scatter_angle;
    auto nuclear_e_loss =
        nuclear_energy_loss(m_1, m_2, energy, cm_scatter_angle);
    VLOG(2) << "Nuclear energy loss: " << nuclear_e_loss;

    energy -= nuclear_e_loss;
    // Record nuclear and total energy loss at the actual particle depth.
    bombardment_.record_energy_loss((long double)nuclear_e_loss.value(),
                                    coordinate_.z_.value(),
                                    Bombardment::EnergyLossType::nuclear);
    bombardment_.record_energy_loss((long double)(e_energy_loss + nuclear_e_loss).value(),
                                    coordinate_.z_.value(),
                                    Bombardment::EnergyLossType::total);


    auto laboratory_scatter_angle =
        laboratory_scattering_angle(cm_scatter_angle, m_1, m_2);
    if (nuclear_e_loss > knock_on_threshold && enable_knock_on) {
      VLOG(2) << "Knock-on particle created with energy " << nuclear_e_loss;
      data_leaf_.add_child(ParticleDataLeaf());
      // Record displacement due to knock-on (top-level recoil)
      bombardment_.record_displacement(coordinate_.z_.value());
      auto knock_on_particle =
          Ion(z_2,
              static_cast<size_t>(std::round((m_2 / atomic_mass_unit).value())),
              simulation_, uniform_distribution_, random_number_generator_,
              data_leaf_.get_children().back(), bombardment_);
      auto scattering_angle_2 = asin(sin(laboratory_scatter_angle) *
                                     sqrt(energy * m_1 / nuclear_e_loss / m_2));
      auto velocity = Velocity(velocity_.quaternion_, nuclear_e_loss);
      knock_on_particle.set_velocity(velocity);
      knock_on_particle.set_coordinate(coordinate_);
      knock_on_particle.rotate_inplace(scattering_angle_2);
      knock_on_particle.fire();
    }
    rotate_inplace(laboratory_scatter_angle); // Changes velocity in-place
    auto &velocity = velocity_.quaternion_;

    VLOG(2) << "Original coordinate (" << coordinate_.x_ << ", "
            << coordinate_.y_ << ", " << coordinate_.z_ << ")";
    VLOG(2) << "Original velocity (" << velocity.x_ << ", " << velocity.y_
            << ", " << velocity.z_ << ")";

    data_leaf_.add_coordinate(coordinate_);
    coordinates_.emplace_back(coordinate_);
    coordinate_.x_ += path_length * velocity.x_;
    coordinate_.y_ += path_length * velocity.y_;
    coordinate_.z_ += path_length * velocity.z_;

    if (coordinate_.z_ < 0.0 * angstrom) {
      data_leaf_.sputter();
      bombardment_.record_sputtered();
      break;
    }

    if (coordinate_.z_.value() > max_depth.value()) {
      data_leaf_.escape();
      break;
    }

    VLOG(2) << "Final coordinate moved " << path_length << " away: ("
            << coordinate_.x_ << ", " << coordinate_.y_ << ", "
            << coordinate_.z_ << ")";

    // 2. Calculate recoil energy and velocity
    // 4. Determine if sputtering occurred
    // 5. Check if damage cascade occurs
    //    a. If so, if energy difference is larger than displacement
    //    energy, create a cascade particle and push it to the particles_
    //    vector
    // }
  }

  auto draw_function = bombardment_.get_simulation().get_draw_function();
  draw_function(coordinates_,
                incident); // HACK: Temp check for bombardment counting

  // // TODO: Remove
  // std::ofstream file_stream("output.csv");
  // file_stream << "scalar,x,y,z\n";
  // file_stream << string_stream.str();

} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

length_quantity
Ion::impact_parameter(number_density_quantity number_density,
                      length_quantity length,
                      dimensionless_quantity reduced_energy) const {
  length_quantity result;

  auto random_number = uniform_distribution_(random_number_generator_);
  // TODO: Fix this. Fix proper random numbers and calculation. High energy
  // gives ridiculous results. Research e >> 10 on 7-8 in SRIM book.
  if (reduced_energy > 10.0) {
    result = sqrt(-log(random_number) /
                  (pi * number_density * length * si::constants::codata::N_A));
    VLOG(2) << "Impact parameter at high energy for number density "
            << number_density << ", length " << length
            << ", and reduced energy " << reduced_energy << " is " << result;
  } else {
    result = sqrt(random_number /
                  (pi * pow<static_rational<2, 3>>(
                            number_density * si::constants::codata::N_A)));
    VLOG(2) << "Impact parameter at low energy for number density "
            << number_density << ", length " << length
            << ", and reduced energy " << reduced_energy << " is " << result;
  }
  return result;
}

length_quantity
Ion::interatomic_spacing(number_density_quantity number_density) const {
  return root<-3>(number_density * constants::N_A);
};

energy_quantity
Ion::electronic_stopping_energy(dimensionless_quantity charge,
                                mass_quantity mass, energy_quantity energy,
                                number_density_quantity number_density,
                                const Properties &properties) const try {
  VLOG(2) << "Calculating ion electronic stopping energy with charge " << charge
          << ", mass " << mass << ", and energy " << energy;
  auto c_1 =
      1.212 * root<2>(constants::m_u) * root<2>(electron_volt) / angstrom;
  VLOG(3) << "Value of c_1 " << c_1;
  const auto k_l = c_1 * pow<static_rational<7, 6>>(charge) *
                   properties.charge_ /
                   (pow<static_rational<3, 2>>(
                        pow<static_rational<2, 3>>(charge) +
                        pow<static_rational<2, 3>>(properties.charge_)) *
                    root<2>(mass));
  VLOG(3) << "Value of k_l " << k_l;
  const auto stopping_energy = k_l * root<2>(energy);
  VLOG(3) << "Value of stopping_energy " << stopping_energy;
  const auto c_2 = 1.59 * pow<3>(angstrom);
  VLOG(3) << "Value of c_2 " << c_2;
  const auto travel_length = interatomic_spacing(number_density);
  VLOG(3) << "Value of travel_length " << travel_length;
  const auto energy_loss =
      c_2 * travel_length * number_density * constants::N_A * stopping_energy;
  VLOG(3) << "Value of energy_loss " << energy_loss;
  return energy_loss;
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

// Quaternion function implementations{{{

Quaternion &Quaternion::operator+=(const Quaternion &b) {
  r_ += b.r_;
  x_ += b.x_;
  y_ += b.y_;
  z_ += b.z_;
  return *this;
}

Quaternion Quaternion::operator+(const Quaternion &b) const {
  return Quaternion(*this) += b;
}

Quaternion operator+(const Quaternion &a, const Quaternion &b) {
  return Quaternion(a) += b;
}

Quaternion &Quaternion::operator*=(const Quaternion &b) {
  auto r_1 = r_;
  auto r_2 = b.r_;
  auto x_1 = x_;
  auto y_1 = y_;
  auto z_1 = z_;
  auto x_2 = b.x_;
  auto y_2 = b.y_;
  auto z_2 = b.z_;

  r_ = r_1 * r_2 - x_1 * x_2 - y_1 * y_2 - z_1 * z_2;
  x_ = r_1 * x_2 + r_2 * x_1 + y_1 * z_2 - z_1 * y_2;
  y_ = r_1 * y_2 + r_2 * y_1 + z_1 * x_2 - x_1 * z_2;
  z_ = r_1 * z_2 + r_2 * z_1 + x_1 * y_2 - y_1 * x_2;
  return *this;
}

Quaternion Quaternion::operator*(const Quaternion &b) const {
  return Quaternion(*this) *= b;
}

// Quaternion operator*(const Quaternion &a, const Quaternion &b) {
//   return Quaternion(a) *= b;
// }

Quaternion &Quaternion::operator*=(double a) {
  r_ *= a;
  x_ *= a;
  y_ *= a;
  z_ *= a;
  return *this;
}

Quaternion Quaternion::operator*(const double b) const {
  return Quaternion(*this) *= b;
}

Quaternion operator*(double a, const Quaternion &b) {
  return Quaternion(b) *= a;
}

Quaternion &Quaternion::operator/=(double a) {
  r_ /= a;
  x_ /= a;
  y_ /= a;
  z_ /= a;
  return *this;
}

Quaternion Quaternion::operator/(const double b) const {
  return Quaternion(*this) /= b;
}

Quaternion operator/(const Quaternion &a, double b) {
  return Quaternion(a) /= b;
}

Quaternion &Quaternion::normalize_inplace() {
  *this /= (this->magnitude());
  return *this;
}

Quaternion &Quaternion::conjugate_inplace() {
  x_ = -x_;
  y_ = -y_;
  z_ = -z_;
  return *this;
}

Quaternion Quaternion::conjugate_copy() const {
  return Quaternion(r_, -x_, -y_, -z_);
}

double Quaternion::magnitude() const {
  return std::sqrt(r_ * r_ + x_ * x_ + y_ * y_ + z_ * z_);
} // }}}

boost::mt19937::result_type EntropySeeder::generate() {
  try {
    VLOG(1) << "Attempting to seed with hardware entropy source.";
    // Try hardware entropy first
    std::random_device hw_rng;
    // Test if hardware entropy is actually available
    if (hw_rng.entropy() > 0) {
      // Get multiple values to increase entropy
      std::array<std::random_device::result_type, 10> entropy_values;
      std::generate(entropy_values.begin(), entropy_values.end(),
                    std::ref(hw_rng));

      // Combine entropy values using XOR
      std::random_device::result_type seed = 0;
      for (const auto &value : entropy_values) {
        seed ^= value;
      }

      VLOG(1) << "Hardware entropy-generated seed value is " << seed;
      return static_cast<boost::mt19937::result_type>(seed);
    } else {
      VLOG(2) << "Hardware entropy generator not found";
    }
  } catch (const std::exception &e) {
    LOG(WARNING) << "Exception caught after attempting hardware entropy generation: "
            << e.what();
    // Hardware RNG might not be available, fall through to backup method
  }

  // Fallback to combining system clock with process/thread specific
  // information
  auto time_seed =
      std::chrono::high_resolution_clock::now().time_since_epoch().count();
  auto pid_seed = static_cast<boost::mt19937::result_type>(
      std::hash<std::thread::id>{}(std::this_thread::get_id()));
  LOG(WARNING) << "Falling back on time " << time_seed << " and pid " << pid_seed;
  auto seed = time_seed ^ pid_seed;
  VLOG(1) << "Calculated fallback seed as " << seed;

  return static_cast<boost::mt19937::result_type>(seed);
}

} // namespace common
} // namespace amcset