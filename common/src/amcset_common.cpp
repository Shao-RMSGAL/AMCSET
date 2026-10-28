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

#include <absl/log/log.h>

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
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <ranges>
#include <sstream>
#include <stdexcept>

#include "amcset_utilities.h" // }}}

namespace amcset {
namespace common {

// Simulation function implementations{{{
Simulation::Simulation(const Settings settings, Volume &&volume,
                       DrawFunction draw_function) try
    : settings_(settings), volume_(std::move(volume)),
      tallies_(settings.range_bin_count_, settings.range_.value(), settings.displacement_bin_count_,
               settings.electronic_energy_loss_bin_count_,
               settings.nuclear_energy_loss_bin_count_,
               settings.total_energy_loss_bin_count_),
      thread_pool_(settings.thread_count_), draw_function_(draw_function) {
  LOG(INFO) << "Creating simulation with settings:\n"
            << Simulation::print_settings();
                                             
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

void Simulation::push_bombardment_results(std::vector<CollisionPoint> &&vec) {
  std::lock_guard<std::mutex> lock(bombardment_queue_mutex_);
  bombardment_queue_.emplace_back(std::move(vec));
}

bool Simulation::try_pop_bombardment_results(std::vector<CollisionPoint> &out) {
  std::lock_guard<std::mutex> lock(bombardment_queue_mutex_);
  if (bombardment_queue_.empty()) return false;
  out = std::move(bombardment_queue_.front());
  bombardment_queue_.pop_front();
  return true;
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
  ss << "Screening model: "
      << (settings_.screening_model_ == Simulation::Settings::ScreeningModel::QuickBasic ? "QuickBasic" : "ThomasFermi")
      << std::endl;
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

  // Publish collected collision points for this bombardment into the
  // simulation queue so the frontend can consume them asynchronously.
  try {
    simulation_.push_bombardment_results(take_collision_points());
  } catch (...) {
    VLOG(1) << "Failed to publish bombardment collision points to queue";
  }

  // Post-bombardment processing.
  // This goes through the incident_particle_ data_leaf_ and
  // looks for knock-ons through the list of pointers
  // contained within the ParticleDataLeaf class. It is a tree-structure,
  // hence why it is called a "leaf".
  // TODO:Remove this comment


} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
} // }}}

} // namespace common
} // namespace amcset