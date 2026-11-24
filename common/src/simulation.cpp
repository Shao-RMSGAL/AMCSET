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

/*! \file simulation.cpp
 *  \brief Simulation::Settings validation, model selection, the
 *  multi-bombardment run loop, tally aggregation/output, and Bombardment.
 */

#include "amcset_common.h"
#include "diagnostics.h"
#include "element_data.h"

#include <absl/log/log.h>

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace common {

namespace {
// Map the Settings enums to the static model functions. The function
// pointers are stored in the Simulation so the transport loop never
// switches on the enums.
Simulation::ElectronicStoppingFunction select_electronic_stopping_function(
    Simulation::Settings::ElectronicStoppingModel model) {
  switch (model) {
  case Simulation::Settings::ElectronicStoppingModel::ZBL:
    return &Simulation::zbl_electronic_stopping;
  case Simulation::Settings::ElectronicStoppingModel::LSS:
    return &Simulation::lss_electronic_stopping;
  case Simulation::Settings::ElectronicStoppingModel::CustomPowerLaw:
    return &Simulation::custom_power_law_electronic_stopping;
  case Simulation::Settings::ElectronicStoppingModel::BetheBloch:
  case Simulation::Settings::ElectronicStoppingModel::Builtin:
  case Simulation::Settings::ElectronicStoppingModel::SrimTable:
    // Evaluated through Simulation::stopping_curve() (see Ion).
    return nullptr;
  }
  throw std::invalid_argument("Unsupported electronic stopping model.");
}

Simulation::NuclearScreeningFunction select_nuclear_screening_function(
    Simulation::Settings::ScreeningModel model) {
  switch (model) {
  case Simulation::Settings::ScreeningModel::ThomasFermi:
  case Simulation::Settings::ScreeningModel::QuickBasic:
  case Simulation::Settings::ScreeningModel::Custom:
    return &Simulation::thomas_fermi_screening;
  case Simulation::Settings::ScreeningModel::ZBL:
    return &Simulation::zbl_screening;
  case Simulation::Settings::ScreeningModel::Moliere:
    return &Simulation::moliere_screening;
  }
  throw std::invalid_argument("Unsupported nuclear screening model.");
}

Simulation::NuclearScreeningFunction
select_nuclear_screening_derivative_function(
    Simulation::Settings::ScreeningModel model) {
  switch (model) {
  case Simulation::Settings::ScreeningModel::ThomasFermi:
  case Simulation::Settings::ScreeningModel::QuickBasic:
  case Simulation::Settings::ScreeningModel::Custom:
    return &Simulation::thomas_fermi_screening_derivative;
  case Simulation::Settings::ScreeningModel::ZBL:
    return &Simulation::zbl_screening_derivative;
  case Simulation::Settings::ScreeningModel::Moliere:
    return &Simulation::moliere_screening_derivative;
  }
  throw std::invalid_argument("Unsupported nuclear screening model.");
}

Simulation::NuclearScreeningLengthFunction
select_nuclear_screening_length_function(
    Simulation::Settings::ScreeningModel model) {
  switch (model) {
  case Simulation::Settings::ScreeningModel::ThomasFermi:
    return &Simulation::thomas_fermi_screening_length;
  case Simulation::Settings::ScreeningModel::Custom:
    return &Simulation::custom_screening_length;
  case Simulation::Settings::ScreeningModel::QuickBasic:
    return &Simulation::firsov_screening_length;
  case Simulation::Settings::ScreeningModel::ZBL:
    return &Simulation::zbl_screening_length;
  case Simulation::Settings::ScreeningModel::Moliere:
    return &Simulation::moliere_screening_length;
  }
  throw std::invalid_argument("Unsupported nuclear screening model.");
}

const char *nuclear_model_name(Simulation::Settings::ScreeningModel model) {
  switch (model) {
  case Simulation::Settings::ScreeningModel::ThomasFermi:
    return "Thomas-Fermi";
  case Simulation::Settings::ScreeningModel::QuickBasic:
    return "Firsov";
  case Simulation::Settings::ScreeningModel::ZBL:
    return "ZBL";
  case Simulation::Settings::ScreeningModel::Moliere:
    return "Moliere";
  case Simulation::Settings::ScreeningModel::Custom:
    return "Custom Thomas-Fermi scale";
  }
  return "Unknown";
}

const char *electronic_model_name(
    Simulation::Settings::ElectronicStoppingModel model) {
  switch (model) {
  case Simulation::Settings::ElectronicStoppingModel::ZBL:
    return "ZBL";
  case Simulation::Settings::ElectronicStoppingModel::LSS:
    return "LSS";
  case Simulation::Settings::ElectronicStoppingModel::CustomPowerLaw:
    return "Custom power law";
  case Simulation::Settings::ElectronicStoppingModel::BetheBloch:
    return "Bethe-Bloch";
  case Simulation::Settings::ElectronicStoppingModel::Builtin:
    return "Built-in (Lindhard-Scharff / Bethe-Bloch)";
  case Simulation::Settings::ElectronicStoppingModel::SrimTable:
    return "SRIM-2013 tables";
  }
  return "Unknown";
}
} // namespace

// Settings overload (3): the full constructor. Overloads (1) and (2) below
// delegate here. This is the only place arguments are validated.
Simulation::Settings::Settings(
  energy_quantity electron_stopping_energy, size_t z_number,
  size_t mass_number, bool enable_damage_cascade,
  energy_quantity ion_stopping_energy,
  energy_quantity ion_displacement_energy,
  size_t divisor_angle_number, size_t flying_distance_number,
  length_quantity range, size_t bombardment_count, bool is_electron,
  energy_quantity incident_energy, size_t thread_count,
  angle_quantity incident_angle, size_t displacement_bin_count,
  size_t electronic_energy_loss_bin_count,
  size_t nuclear_energy_loss_bin_count,
  size_t total_energy_loss_bin_count,
  size_t range_bin_count,
  std::string displacement_histogram_output_file,
  std::string electronic_energy_loss_histogram_output_file,
  std::string nuclear_energy_loss_histogram_output_file,
  std::string total_energy_loss_histogram_output_file,
  bool save_per_bombardment_info,
  std::string per_bombardment_output_dir,
  std::string electron_screening_parameters_path,
  std::string mott_scattering_parameters_path,
  Simulation::Settings::ScreeningModel screening_model,
  std::string range_histogram_output_file,
  Simulation::Settings::ElectronicStoppingModel electronic_stopping_model,
  double custom_electronic_coefficient,
  double custom_electronic_exponent,
  double custom_recoil_electronic_coefficient,
  double custom_recoil_electronic_exponent,
  double custom_nuclear_screening_scale,
  std::string srim_stopping_table_path) try
  : electron_stopping_energy_(electron_stopping_energy), z_number_(z_number),
    mass_number_(mass_number),
    incident_particle_properties_(Particle::Properties(z_number, mass_number)),
    enable_damage_cascade_(enable_damage_cascade),
    ion_stopping_energy_(ion_stopping_energy),
    ion_displacement_energy_(ion_displacement_energy),
    divisor_angle_number_(divisor_angle_number),
    flying_distance_number_(flying_distance_number), range_(range),
    bombardment_count_(bombardment_count), is_electron_(is_electron),
    incident_energy_(incident_energy), thread_count_(thread_count),
    incident_angle_(incident_angle),
    displacement_bin_count_(displacement_bin_count),
    electronic_energy_loss_bin_count_(electronic_energy_loss_bin_count),
    nuclear_energy_loss_bin_count_(nuclear_energy_loss_bin_count),
    total_energy_loss_bin_count_(total_energy_loss_bin_count),
    range_bin_count_(range_bin_count),
    displacement_histogram_output_file_(std::move(displacement_histogram_output_file)),
    electronic_energy_loss_histogram_output_file_(
      std::move(electronic_energy_loss_histogram_output_file)),
    nuclear_energy_loss_histogram_output_file_(
      std::move(nuclear_energy_loss_histogram_output_file)),
    total_energy_loss_histogram_output_file_(
      std::move(total_energy_loss_histogram_output_file)),
    range_histogram_output_file_(
      range_histogram_output_file.empty()
        ? (std::filesystem::path(displacement_histogram_output_file_)
             .parent_path() /
           "range_histogram.csv")
            .string()
        : std::move(range_histogram_output_file)),
    save_per_bombardment_info_(save_per_bombardment_info),
    per_bombardment_output_dir_((per_bombardment_output_dir.empty()) ? std::string("output/bombardments/") : std::move(per_bombardment_output_dir)),
    electron_screening_parameters_path_(std::move(electron_screening_parameters_path)),
    mott_scattering_parameters_path_(std::move(mott_scattering_parameters_path)),
    screening_model_(screening_model),
    electronic_stopping_model_(electronic_stopping_model),
    custom_electronic_coefficient_(custom_electronic_coefficient),
    custom_electronic_exponent_(custom_electronic_exponent),
    custom_recoil_electronic_coefficient_(
        custom_recoil_electronic_coefficient),
    custom_recoil_electronic_exponent_(custom_recoil_electronic_exponent),
    custom_nuclear_screening_scale_(custom_nuclear_screening_scale),
    srim_stopping_table_path_(std::move(srim_stopping_table_path)) {
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
  if (electronic_stopping_model_ ==
          ElectronicStoppingModel::CustomPowerLaw &&
      (!std::isfinite(custom_electronic_coefficient_) ||
       custom_electronic_coefficient_ <= 0.0 ||
       !std::isfinite(custom_electronic_exponent_) ||
       custom_electronic_exponent_ < 0.0 ||
       !std::isfinite(custom_recoil_electronic_coefficient_) ||
       custom_recoil_electronic_coefficient_ <= 0.0 ||
       !std::isfinite(custom_recoil_electronic_exponent_) ||
       custom_recoil_electronic_exponent_ < 0.0)) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Custom electronic stopping coefficients must be finite and positive "
        "and their exponents must be finite and non-negative."));
  }
  if (electronic_stopping_model_ == ElectronicStoppingModel::SrimTable &&
      srim_stopping_table_path_.empty()) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "The SRIM-table electronic stopping model needs the path of the SRIM "
        "stopping tables."));
  }
  if (!std::isfinite(custom_nuclear_screening_scale_) ||
      custom_nuclear_screening_scale_ <= 0.0) {
    throw std::invalid_argument(EXCEPTION_MESSAGE(
        "Custom nuclear screening scale must be finite and positive."));
  }
  if (displacement_bin_count == 0 || electronic_energy_loss_bin_count == 0 ||
      nuclear_energy_loss_bin_count == 0 || total_energy_loss_bin_count == 0 ||
      range_bin_count == 0) {
    throw std::invalid_argument(
        EXCEPTION_MESSAGE("Histogram bin counts must be greater than 0"));
  }
} catch (...) {
  rethrow();
}

// Settings overload (1): optional per-bombardment output, default asset paths.
Simulation::Settings::Settings(
  energy_quantity electron_stopping_energy, size_t z_number,
  size_t mass_number, bool enable_damage_cascade,
  energy_quantity ion_stopping_energy,
  energy_quantity ion_displacement_energy,
  size_t divisor_angle_number, size_t flying_distance_number,
  length_quantity range, size_t bombardment_count, bool is_electron,
  energy_quantity incident_energy, size_t thread_count,
  angle_quantity incident_angle,
  size_t displacement_bin_count,
  size_t electronic_energy_loss_bin_count,
  size_t nuclear_energy_loss_bin_count,
  size_t total_energy_loss_bin_count,
  size_t range_bin_count,
  std::string displacement_histogram_output_file,
  std::string electronic_energy_loss_histogram_output_file,
  std::string nuclear_energy_loss_histogram_output_file,
  std::string total_energy_loss_histogram_output_file,
  bool save_per_bombardment_info,
  std::string per_bombardment_output_dir,
  Simulation::Settings::ScreeningModel screening_model,
  std::string range_histogram_output_file,
  Simulation::Settings::ElectronicStoppingModel electronic_stopping_model,
  double custom_electronic_coefficient,
  double custom_electronic_exponent,
  double custom_recoil_electronic_coefficient,
  double custom_recoil_electronic_exponent,
  double custom_nuclear_screening_scale,
  std::string srim_stopping_table_path)
  : Settings(electron_stopping_energy, z_number, mass_number,
         enable_damage_cascade, ion_stopping_energy,
         ion_displacement_energy,
         divisor_angle_number, flying_distance_number, range,
         bombardment_count, is_electron, incident_energy, thread_count,
         incident_angle, displacement_bin_count,
         electronic_energy_loss_bin_count,
         nuclear_energy_loss_bin_count,
         total_energy_loss_bin_count,
         range_bin_count,
         std::move(displacement_histogram_output_file),
         std::move(electronic_energy_loss_histogram_output_file),
         std::move(nuclear_energy_loss_histogram_output_file),
         std::move(total_energy_loss_histogram_output_file),
         save_per_bombardment_info, std::move(per_bombardment_output_dir),
         std::string("common/assets/electron_screening_parameters.csv"),
         std::string("common/assets/mott_scattering_parameters.csv"),
         screening_model, std::move(range_histogram_output_file),
         electronic_stopping_model, custom_electronic_coefficient,
         custom_electronic_exponent, custom_recoil_electronic_coefficient,
         custom_recoil_electronic_exponent, custom_nuclear_screening_scale,
         std::move(srim_stopping_table_path)) {}

// Settings overload (2): per-bombardment output disabled, default asset paths.
Simulation::Settings::Settings(
  energy_quantity electron_stopping_energy, size_t z_number,
  size_t mass_number, bool enable_damage_cascade,
  energy_quantity ion_stopping_energy,
  energy_quantity ion_displacement_energy,
  size_t divisor_angle_number, size_t flying_distance_number,
  length_quantity range, size_t bombardment_count, bool is_electron,
  energy_quantity incident_energy, size_t thread_count,
  angle_quantity incident_angle,
  size_t displacement_bin_count,
  size_t electronic_energy_loss_bin_count,
  size_t nuclear_energy_loss_bin_count,
  size_t total_energy_loss_bin_count,
  size_t range_bin_count,
  std::string displacement_histogram_output_file,
  std::string electronic_energy_loss_histogram_output_file,
  std::string nuclear_energy_loss_histogram_output_file,
  std::string total_energy_loss_histogram_output_file,
  Simulation::Settings::ScreeningModel screening_model,
  std::string range_histogram_output_file,
  Simulation::Settings::ElectronicStoppingModel electronic_stopping_model,
  double custom_electronic_coefficient,
  double custom_electronic_exponent,
  double custom_recoil_electronic_coefficient,
  double custom_recoil_electronic_exponent,
  double custom_nuclear_screening_scale,
  std::string srim_stopping_table_path)
  : Settings(electron_stopping_energy, z_number, mass_number,
         enable_damage_cascade, ion_stopping_energy,
         ion_displacement_energy,
         divisor_angle_number, flying_distance_number, range,
         bombardment_count, is_electron, incident_energy, thread_count,
         incident_angle, displacement_bin_count,
         electronic_energy_loss_bin_count,
         nuclear_energy_loss_bin_count,
         total_energy_loss_bin_count,
         range_bin_count,
         std::move(displacement_histogram_output_file),
         std::move(electronic_energy_loss_histogram_output_file),
         std::move(nuclear_energy_loss_histogram_output_file),
         std::move(total_energy_loss_histogram_output_file),
         false, std::string("output/bombardments/"),
         std::string("common/assets/electron_screening_parameters.csv"),
         std::string("common/assets/mott_scattering_parameters.csv"),
         screening_model, std::move(range_histogram_output_file),
         electronic_stopping_model, custom_electronic_coefficient,
         custom_electronic_exponent, custom_recoil_electronic_coefficient,
         custom_recoil_electronic_exponent, custom_nuclear_screening_scale,
         std::move(srim_stopping_table_path)) {}

Simulation::Simulation(const Settings settings, Volume &&volume) try
    : settings_(settings),
      electronic_stopping_function_(
          select_electronic_stopping_function(
              settings_.electronic_stopping_model_)),
      nuclear_stopping_function_(
          &Simulation::binary_collision_nuclear_stopping),
      nuclear_screening_function_(
          select_nuclear_screening_function(settings_.screening_model_)),
      nuclear_screening_derivative_function_(
          select_nuclear_screening_derivative_function(
              settings_.screening_model_)),
      nuclear_screening_length_function_(
          select_nuclear_screening_length_function(
              settings_.screening_model_)),
      volume_(std::move(volume)),
      tallies_(settings.range_bin_count_, settings.range_.value(), settings.displacement_bin_count_,
               settings.electronic_energy_loss_bin_count_,
               settings.nuclear_energy_loss_bin_count_,
               settings.total_energy_loss_bin_count_),
      thread_pool_(settings.thread_count_) {
  LOG(INFO) << "Creating simulation with settings:\n"
            << Simulation::print_settings();
  // Curve-based electronic stopping: one curve per (projectile, element).
  if (!settings_.is_electron_ && electronic_stopping_function_ == nullptr) {
    using Model = Settings::ElectronicStoppingModel;
    std::shared_ptr<const physics::SrimStoppingDatabase> db;
    physics::ElectronicStoppingModel model = physics::ElectronicStoppingModel::builtin;
    switch (settings_.electronic_stopping_model_) {
    case Model::BetheBloch: model = physics::ElectronicStoppingModel::bethe_bloch; break;
    case Model::SrimTable:
      model = physics::ElectronicStoppingModel::srim_tables;
      db = physics::SrimStoppingDatabase::load(settings_.srim_stopping_table_path_);
      break;
    default: break;
    }
    // Projectiles: the incident ion and every target isotope (recoils).
    std::vector<std::pair<unsigned, mass_quantity>> projectiles;
    projectiles.emplace_back(static_cast<unsigned>(settings_.z_number_),
                             settings_.incident_particle_properties_.mass_);
    std::vector<unsigned> elements;
    for (size_t li = 0; li < volume_.get_layer_count(); ++li) {
      const auto &layer = volume_.get_layer(li);
      for (size_t c = 0; c < layer.get_number_of_components(); ++c) {
        const auto &p = layer.get_property(c);
        const auto z = static_cast<unsigned>(p.charge_.value());
        projectiles.emplace_back(z, p.mass_);
        if (std::find(elements.begin(), elements.end(), z) == elements.end()) {
          elements.push_back(z);
        }
      }
    }
    for (const auto &[z1, m1] : projectiles) {
      for (unsigned z2 : elements) {
        const auto key = std::make_tuple(
            z1, std::llround((m1 / atomic_mass_unit).value() * 1e6), z2);
        if (stopping_curves_.count(key)) continue;
        // Representative atomic density of the element (only the
        // Bethe-Bloch density effect depends on it).
        const auto density = Element::by_z(z2).solid_atomic_density();
        stopping_curves_.emplace(
            key, physics::StoppingCurve(model, db.get(), z1, m1, {{z2, 1.0}}, false,
                                        1.0, density));
      }
    }
  }
  std::vector<Bombardment> bombardments;
  bombardments.reserve(settings.bombardment_count_);

  for (size_t i = 0; i < settings.bombardment_count_; i++) {
    VLOG(1) << "Creating bombardment " << i;
    bombardments.push_back(Bombardment(*this, i));
  }
  bombardments_ = std::move(bombardments);
  LOG(INFO) << settings.bombardment_count_ << " Bombardments Created";
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

const physics::StoppingCurve *Simulation::stopping_curve(unsigned z1, mass_quantity m1,
                                                       unsigned z2) const {
  const auto it = stopping_curves_.find(
      std::make_tuple(z1, std::llround((m1 / atomic_mass_unit).value() * 1e6), z2));
  return it == stopping_curves_.end() ? nullptr : &it->second;
}

int Simulation::run_simulation() try {
  LOG(INFO) << "Starting Simulation...";
  AMCSET_LOG(info) << "native engine run: " << bombardments_.size() << " bombardments, "
                   << settings_.thread_count_ << " threads";
  // An exception must not leave a pool thread (asio would let it reach
  // std::terminate and end the program): the first one is kept, the other
  // bombardments are cancelled, and it is rethrown after the join.
  std::exception_ptr worker_error;
  std::mutex worker_error_mutex;
  for (auto &bombardment : bombardments_) {
    boost::asio::post(thread_pool_, [this, &bombardment, &worker_error, &worker_error_mutex]() {
      if (cancellation_requested_.load(std::memory_order_acquire)) {
        return;
      }
      try {
        VLOG(1) << "Starting Bombardment...";
        bombardment.run_bombardment();
        increment_completed_bombardment_count();
        VLOG(1) << "...Bombardment Complete! Completed: "
                << get_completed_bombardment_count() << "/"
                << settings_.bombardment_count_;
      } catch (...) {
        std::lock_guard<std::mutex> lock(worker_error_mutex);
        if (!worker_error) worker_error = std::current_exception();
        cancellation_requested_.store(true, std::memory_order_release);
      }
    });
  }
  thread_pool_.join();
  if (worker_error) {
    AMCSET_LOG(error) << "native engine run failed after " << get_completed_bombardment_count()
                      << " bombardments";
    std::rethrow_exception(worker_error);
  }
  if (is_cancelled()) {
    LOG(INFO) << "...Simulation cancelled after "
              << get_completed_bombardment_count() << "/"
              << settings_.bombardment_count_ << " bombardments";
  } else {
    LOG(INFO) << "...Simulation complete";
  }
  try {
    // tallies.json goes next to the histogram files (the directory of the
    // displacement histogram; "output" if that path has no directory).
    std::filesystem::path output_root =
        std::filesystem::path(settings_.displacement_histogram_output_file_).parent_path();
    if (output_root.empty()) output_root = "output";
    std::filesystem::create_directories(output_root);
    // Remove the collisions.csv written by older versions so the GUI does not
    // display stale data.
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
    write_histogram_with_log("range", settings_.range_histogram_output_file_,
                 get_tallies().range_histogram, "count");
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
    LOG(ERROR) << "Failed to write simulation output: " << e.what();
  } catch (...) {
    LOG(ERROR) << "Failed to write simulation output";
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

void Simulation::clear_bombardment_results() {
  std::lock_guard<std::mutex> lock(bombardment_queue_mutex_);
  bombardment_queue_.clear();
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
  ss << "Range histogram bins: " << settings_.range_bin_count_ << std::endl;
  ss << "Nuclear stopping model: "
     << nuclear_model_name(settings_.screening_model_) << std::endl;
  ss << "Electronic stopping model: "
     << electronic_model_name(settings_.electronic_stopping_model_)
     << std::endl;
  if (settings_.electronic_stopping_model_ ==
      Settings::ElectronicStoppingModel::CustomPowerLaw) {
    ss << "Custom electronic stopping C: "
       << settings_.custom_electronic_coefficient_ << std::endl;
    ss << "Custom electronic stopping n: "
       << settings_.custom_electronic_exponent_ << std::endl;
    ss << "Custom recoil electronic stopping C: "
       << settings_.custom_recoil_electronic_coefficient_ << std::endl;
    ss << "Custom recoil electronic stopping n: "
       << settings_.custom_recoil_electronic_exponent_ << std::endl;
  }
  if (settings_.electronic_stopping_model_ ==
      Settings::ElectronicStoppingModel::SrimTable) {
    ss << "SRIM stopping tables: " << settings_.srim_stopping_table_path_ << std::endl;
  }
  if (settings_.screening_model_ == Settings::ScreeningModel::Custom) {
    ss << "Custom nuclear screening scale: "
       << settings_.custom_nuclear_screening_scale_ << std::endl;
  }
  return ss.str();
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

void Simulation::merge_tallies(const Tallies &tallies) {
  std::lock_guard<std::mutex> lk(tallies_mutex_);
  tallies_.merge(tallies);
}

Tallies Simulation::snapshot_tallies() const {
  std::lock_guard<std::mutex> lk(tallies_mutex_);
  return tallies_;
}

const Tallies &Simulation::get_tallies() const noexcept { return tallies_; }

void Simulation::write_histogram_csv(const std::string &output_file,
                                    const Histogram &histogram,
                                    const std::string &value_label) const {
  if (output_file.empty()) {
    throw std::invalid_argument("Histogram output file path cannot be empty.");
  }
  common::write_histogram_csv(std::filesystem::path(output_file), histogram, value_label);
}

void write_histogram_csv(const std::filesystem::path &output_path, const Histogram &histogram,
                         const std::string &value_label) {
  if (!output_path.parent_path().empty()) {
    std::filesystem::create_directories(output_path.parent_path());
  }

  std::ofstream output(output_path);
  if (!output.is_open()) {
    throw std::runtime_error("Failed to open histogram output: " + output_path.string());
  }

  // Histogram limits are stored in metres.
  const auto meters_to_angstrom = [](double meters) {
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
    output << meters_to_angstrom(left) << "," << meters_to_angstrom(right) << "," << histogram.bins[i] << "\n";
  }
  output.flush();
  LOG(INFO) << "Histogram written to " << output_path.string() << " with "
            << histogram.bins.size() << " bins";
}

void write_tallies_files(const std::filesystem::path &directory, const Tallies &tallies) {
  std::filesystem::create_directories(directory);
  write_histogram_csv(directory / "displacement_histogram.csv", tallies.displacement_histogram, "count");
  write_histogram_csv(directory / "range_histogram.csv", tallies.range_histogram, "count");
  write_histogram_csv(directory / "electronic_loss_histogram.csv",
                      tallies.electronic_energy_loss_histogram, "energy_loss_eV");
  write_histogram_csv(directory / "nuclear_loss_histogram.csv",
                      tallies.nuclear_energy_loss_histogram, "energy_loss_eV");
  write_histogram_csv(directory / "total_loss_histogram.csv",
                      tallies.total_energy_loss_histogram, "energy_loss_eV");
  std::ofstream json(directory / "tallies.json");
  if (!json) throw std::runtime_error("Failed to open " + (directory / "tallies.json").string());
  json << tallies.to_json();
}

// The depth histogram is a legacy combined histogram: it receives both the
// incident particle's final depth and every displacement depth.
void Bombardment::record_incident(double range_meters) {
  local_tallies_.incidents += 1;
  local_tallies_.range_stats.add((long double)range_meters);
  local_tallies_.range_histogram.add(range_meters);
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

Bombardment::Bombardment(Simulation &simulation, size_t id) try
    : simulation_(simulation),
      local_tallies_(simulation.get_settings(&Simulation::Settings::range_bin_count_),
              simulation.get_settings(&Simulation::Settings::range_).value(),
              simulation.get_settings(&Simulation::Settings::displacement_bin_count_),
              simulation.get_settings(&Simulation::Settings::electronic_energy_loss_bin_count_),
              simulation.get_settings(&Simulation::Settings::nuclear_energy_loss_bin_count_),
              simulation.get_settings(&Simulation::Settings::total_energy_loss_bin_count_)),
      random_number_generator_(EntropySeeder::generate()),
      uniform_distribution_(), id_(id), data_leaf_() {
  // Note: a former VLOG(1) here drew a random number, so enabling verbose
  // logging silently changed the random stream. It was removed.
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
  VLOG(1) << "...Bombardment complete.";

  try {
    // Range = depth of the incident particle's last recorded coordinate.
    // For ions that is the position before the final free flight; for
    // electrons the position after the last step. Sputtered and escaped
    // histories are included.
    auto &coords = incident_particle_->get_coordinates();
    if (!coords.empty()) {
      auto last_z = coords.back().z_.value();
      record_incident(last_z);
    }
  } catch (...) {
    VLOG(1) << "Failed to collect local tallies";
  }

  simulation_.merge_tallies(local_tallies_);

  try {
    // Before publishing collision points, optionally write them to
    // <per_bombardment_output_dir>/bombardment_<id>.csv. The column order is
    // relied upon by scripts/animate_cascade.py and scripts/julia.
    auto pts = take_collision_points();
    if (simulation_.get_settings(&Simulation::Settings::save_per_bombardment_info_)) {
      try {
        const auto outdir = simulation_.get_settings(&Simulation::Settings::per_bombardment_output_dir_);
        std::filesystem::create_directories(outdir);
        const auto filename = std::filesystem::path(outdir) / ("bombardment_" + std::to_string(id_) + ".csv");
        std::ofstream ofs(filename);
        if (ofs.is_open()) {
          ofs << "energy (eV),x position (A),y position (A),z position (A),x velocity component,y velocity component,z velocity component,z number,atomic mass (amu),generation,time (s)\n";
          for (const auto &p : pts) {
            const auto x = (p.coordinate.x_ / angstrom).value();
            const auto y = (p.coordinate.y_ / angstrom).value();
            const auto z = (p.coordinate.z_ / angstrom).value();
            ofs << p.energy_eV << "," << x << "," << y << "," << z << "," << p.vx << "," << p.vy << "," << p.vz << "," << p.z_number << "," << p.atomic_mass_number << "," << p.generation << "," << p.time << "\n";
          }
          ofs.flush();
        }
      } catch (...) {
        LOG(WARNING) << "Failed to write per-bombardment CSV for id " << id_;
      }
    }
    simulation_.push_bombardment_results(std::move(pts));
  } catch (...) {
    VLOG(1) << "Failed to publish bombardment collision points to queue";
  }
} catch (...) {
  rethrow(EXCEPTION_MESSAGE(""));
}

} // namespace common
} // namespace amcset
