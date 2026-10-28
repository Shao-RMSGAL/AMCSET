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

// EntropySeeder implementation separated from amcset_common.cpp
#include "amcset_common.h"

#include <array>
#include <algorithm>
#include <chrono>
#include <random>
#include <thread>
#include <absl/log/log.h>

namespace amcset {
namespace common {

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
