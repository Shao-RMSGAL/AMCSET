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

using namespace amcset::server;
using namespace amcset::common;

int main(int argc, char *argv[]) {

  // boost::interprocess::message_queue::remove("server_to_client_message_queue");
  // boost::interprocess::message_queue::remove("client_to_server_message_queue");

  // boost::interprocess::message_queue server_to_client_message_queue(
  //     boost::interprocess::create_only, "server_to_client_message_queue",
  //     100, 256);
  // boost::interprocess::message_queue client_to_server_message_queue(
  //     boost::interprocess::create_only, "client_to_server_message_queue",
  //     100, 256);

  // size_t count = 0;

  // while (true) {
  //   char recieved_message[256];
  //   size_t recieved_size;
  //   unsigned int priority;
  //   client_to_server_message_queue.receive(
  //       &recieved_message, sizeof(recieved_message), recieved_size,
  //       priority);
  //   recieved_message[recieved_size] = '\0';

  //   if (std::string(recieved_message) == "exit") {
  //       std::cout << "Exiting..." << std::endl;
  //     return 0;
  //   }

  //   std::cout << "Recieved: " << recieved_message << std::endl;

  //   count++;

  //   std::string message = "Server call count: " + std::to_string(count);

  //   server_to_client_message_queue.send(message.c_str(), message.size(),
  //   0);
  // }

  // boost::interprocess::message_queue::remove("server_to_client_message_queue");
  // boost::interprocess::message_queue::remove("client_to_server_message_queue");

  // Initialize logger
  google::InitGoogleLogging(argv[0]);
  FLAGS_alsologtostderr = 1;

  // Process flags
  gflags::ParseCommandLineFlags(&argc, &argv, true);

  LOG(INFO) << "Starting AMCSET GUI...";

  // Create AMCSET settings (TODO: Read this from a file and allow a user to
  // configure)
  energy_quantity electron_stopping_energy =
      40.0 * electron_volt;           // Energy at which electrons stop
  size_t z_number = 1;                // Hydrogen
  size_t mass_number = 1;             // Hydrogen
  bool enable_damage_cascade = false; // TODO: Update to test cascades
  energy_quantity ion_stopping_energy =
      40.0 * electron_volt; // Typical for iron
  energy_quantity ion_displacement_energy =
      40.0 * electron_volt; // Typical for iron
  bool log_single_displacement =
      false; // TODO: Test later, will consume a LOT of memory/disk space.
  size_t divisor_angle_number = 1000; // Break angle range into 1000 segements.
                                      // TODO: Test with different numbers
  size_t flying_distance_number =
      1000; // Break flying distnaces into groups of 1000. TODO: Test with
            // different numbers
  length_quantity range =
      100000.0 * angstrom; // Target range. Should be unused. TODO: Remove once
                           // confirmed to be unused.
  size_t bombardment_count = 1; // TODO: Test with more later.
  bool is_electron = false; // Start with ions. TODO: Test with electrons later.
  energy_quantity incident_energy =
      100.0 * kilo_electron_volt; // TODO: Test with other energies
  size_t thread_count =
      1; // Try with single threading, TODO: Try with multithreading
  // TODO: Add an option to use GPU (Use Kokkos for cross-platform
  // compatability. Note that it is not available in conan. Either make a conan
  // recipe for it or use another installation method.

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

  bool no_gui = true;
  if (no_gui) {
    simulation.run_simulation();
    return 0;
  } else {
    auto application = Application{argc, argv, &simulation};
    auto window1 = Window1{&application};
    window1.resize(800, 600);
    window1.show();

    auto exit_code = application.exec();
    LOG(INFO) << "...Quitting AMCSET GUI. Exit code: " << exit_code;
    return exit_code;
  }
}
