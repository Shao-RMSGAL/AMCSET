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

#include <iostream>

#include <boost/interprocess/ipc/message_queue.hpp>
#include <boost/math/constants/constants.hpp>

#include <QApplication>
#include <QFrame>
#include <QLabel>
#include <QMainWindow>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "amcset_common.h"
#include "amcset_server.h"
#include "amcset_utilities.h"
#include "config.h"

DEFINE_bool(gui, true, "Enable or disable the graphical user interface (GUI)");

using namespace amcset::server;
using namespace amcset::common;

int main(int argc, char *argv[]) {
  try { // {{{

    // Initialize logger
    google::InitGoogleLogging(argv[0]);
    FLAGS_alsologtostderr = 1;

    std::string usage_message = "AMCSET. Call " + std::string(argv[0]);
    gflags::SetUsageMessage(usage_message);

    std::ostringstream version;
    version << amcset_VERSION_MAJOR << '.' << amcset_VERSION_MINOR << '.'
            << amcset_VERSION_PATCH;
    gflags::SetVersionString(version.str());
    // Process flags
    gflags::ParseCommandLineFlags(&argc, &argv, true);

    // Create AMCSET settings (TODO: Read this from a file and allow a user to
    // configure)
    energy_quantity electron_stopping_energy =
        40.0 * electron_volt;           // Energy at which electrons stop
    size_t z_number = 1;                // Hydrogen
    size_t mass_number = 1;             // Hydrogen
    bool enable_damage_cascade = false; // TEST: Update to test cascades
    energy_quantity ion_stopping_energy =
        40.0 * electron_volt; // Typical for iron
    energy_quantity ion_displacement_energy =
        40.0 * electron_volt; // Typical for iron
    bool log_single_displacement =
        false; // TEST: Test later, will consume a LOT of memory/disk space.
    size_t divisor_angle_number =
        1000; // Break angle range into 1000 segements.
              // TEST: Test with different numbers
    size_t flying_distance_number =
        1000; // Break flying distnaces into groups of 1000. TEST: Test with
              // different numbers
    length_quantity range =
        100000.0 * angstrom; // Target range. Should be unused. HACK: Remove
                             // once confirmed to be unused.
    size_t bombardment_count = 1; // TEST: Test with more later.
    bool is_electron =
        false; // Start with ions. TEST: Test with electrons later.
    energy_quantity incident_energy =
        100.0 * kilo_electron_volt; // TEST: Test with other energies
    size_t thread_count =
        1; // Try with single threading, PERF: Try with multithreading

    // PERF: Add an option to use GPU (Use Kokkos for cross-platform
    // compatability. Note that it is not available in conan. Either make a
    // conan recipe for it or use another installation method.

    auto settings = Simulation::Settings(
        electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
        ion_stopping_energy, ion_displacement_energy, log_single_displacement,
        divisor_angle_number, flying_distance_number, range, bombardment_count,
        is_electron, incident_energy, thread_count);

    // Create AMCSET volume (TODO: Read this from a file and allow a user to
    // configure)

    constexpr double natural_abundance_54 = 0.05845;
    constexpr double natural_abundance_56 = 0.91754;
    constexpr double natural_abundance_57 = 0.02119;
    constexpr double natural_abundance_58 = 0.00282;
    Layer::material_vector material{
        {natural_abundance_54, Particle::Properties(26, 54)},
        {natural_abundance_56, Particle::Properties(26, 56)},
        {natural_abundance_57, Particle::Properties(26, 57)},
        {natural_abundance_58, Particle::Properties(26, 58)},
    };
    auto depth = 10000.0 * angstrom;
    constexpr auto density = 7874.0 * kg_per_cubic_meter;

    std::vector<Layer> single_layer;
    single_layer.push_back(Layer(std::move(material), depth, density));
    auto volume = Volume(std::move(single_layer));

    // Create the AMCSET simulation object
    auto simulation = Simulation(settings, std::move(volume));
    int exit_code;

    if (FLAGS_gui) {
      LOG(INFO) << "Starting AMCSET GUI...";
      auto application = Application{argc, argv, &simulation};
      VLOG(1) << "Initializing window";
      auto main_window = MainWindow{&application};
      VLOG(1) << "Resizing window";
      main_window.resize(800, 600); // TODO: Make this a setting/ store it when
                                    // the user manually resizes.
      VLOG(1) << "Showing window";
      main_window.show();
      VLOG(1) << "Executing application";
      exit_code = application.exec();
      LOG(INFO) << "...Quitting AMCSET GUI. Exit code: " << exit_code;
      return exit_code;
    } else {
      LOG(INFO) << "Starting AMCSET in CLI mode...";
      exit_code = simulation.run_simulation();
      LOG(INFO) << "...Quitting AMCSET CLI.";
      return 0;
    }
  } catch (const std::exception &e) {
    LOG(FATAL) << (e.what());
  } // }}}
}
