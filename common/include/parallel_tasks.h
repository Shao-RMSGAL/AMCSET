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

/*! \file parallel_tasks.h
 *  \brief Run a few independent tasks (e.g. output files) concurrently.
 */

#pragma once

#include <exception>
#include <functional>
#include <future>
#include <system_error>
#include <vector>

namespace amcset {
namespace common {

/*!
 * \brief Run every task, concurrently where threads are available, and wait
 * for all of them. The first exception (in task order) is rethrown after
 * every task has finished, so no task outlives the call. A task that cannot
 * get a thread runs on the calling thread.
 */
inline void run_tasks(const std::vector<std::function<void()>> &tasks) {
  std::vector<std::future<void>> running;
  std::vector<std::exception_ptr> errors(tasks.size());
  std::vector<std::size_t> index;
  running.reserve(tasks.size());
  for (std::size_t i = 0; i < tasks.size(); ++i) {
    if (i + 1 == tasks.size()) {
      // The last task runs here: the calling thread would only wait.
      try {
        tasks[i]();
      } catch (...) {
        errors[i] = std::current_exception();
      }
      break;
    }
    try {
      running.push_back(std::async(std::launch::async, tasks[i]));
      index.push_back(i);
    } catch (const std::system_error &) {
      try {
        tasks[i]();
      } catch (...) {
        errors[i] = std::current_exception();
      }
    }
  }
  for (std::size_t k = 0; k < running.size(); ++k) {
    try {
      running[k].get();
    } catch (...) {
      errors[index[k]] = std::current_exception();
    }
  }
  for (const auto &e : errors) {
    if (e) std::rethrow_exception(e);
  }
}

} // namespace common
} // namespace amcset
