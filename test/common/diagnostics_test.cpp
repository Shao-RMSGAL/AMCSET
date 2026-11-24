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

/*! \file diagnostics_test.cpp
 *  \brief The diagnostic log and the parallel task runner. (Crash reports
 *  are checked end to end by the cli.crash_report_* CTest tests.)
 */

#include <gtest/gtest.h>

#include <atomic>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "diagnostics.h"
#include "parallel_tasks.h"

using namespace amcset;

TEST(Diagnostics, LevelsParseAndFilter) {
  EXPECT_EQ(diag::parse_level("debug"), diag::Level::debug);
  EXPECT_EQ(diag::parse_level("warn"), diag::Level::warning);
  EXPECT_EQ(diag::parse_level("off"), diag::Level::off);
  EXPECT_THROW(diag::parse_level("loud"), std::invalid_argument);
  const auto before = diag::level();
  diag::set_level(diag::Level::warning);
  EXPECT_FALSE(diag::enabled(diag::Level::info));
  EXPECT_TRUE(diag::enabled(diag::Level::error));
  EXPECT_FALSE(diag::enabled(diag::Level::off));
  diag::set_level(before);
}

TEST(Diagnostics, LogFileGetsTheEnabledLines) {
  const auto dir = std::filesystem::temp_directory_path() / "amcset_test_diagnostics";
  std::filesystem::remove_all(dir);
  diag::Options o;
  o.program = "common-test";
  o.log_file = dir / "test.log";
  o.level = diag::Level::info;
  diag::initialize("common-test", o);
  EXPECT_EQ(diag::log_file(), o.log_file);
  EXPECT_EQ(std::filesystem::absolute(diag::crash_directory()), std::filesystem::absolute(dir));
  int evaluated = 0;
  AMCSET_LOG(info) << "visible line " << 42;
  AMCSET_LOG(debug) << "hidden line " << ++evaluated; // not evaluated when disabled
  EXPECT_EQ(evaluated, 0);
  // Back to memory-only logging (closes the file).
  diag::initialize("common-test", diag::Options{});
  std::ifstream in(o.log_file);
  std::stringstream text;
  text << in.rdbuf();
  EXPECT_NE(text.str().find("[info]"), std::string::npos);
  EXPECT_NE(text.str().find("visible line 42"), std::string::npos);
  EXPECT_EQ(text.str().find("hidden line"), std::string::npos);
}

TEST(ParallelTasks, RunsEveryTaskAndRethrowsTheFirstError) {
  std::atomic<int> ran{0};
  std::vector<std::function<void()>> tasks;
  for (int i = 0; i < 12; ++i) {
    tasks.push_back([&ran, i] {
      ++ran;
      if (i == 3) throw std::runtime_error("task 3");
      if (i == 7) throw std::logic_error("task 7");
    });
  }
  try {
    common::run_tasks(tasks);
    FAIL() << "expected an exception";
  } catch (const std::runtime_error &e) {
    EXPECT_STREQ(e.what(), "task 3"); // the first in task order
  }
  EXPECT_EQ(ran.load(), 12); // every task ran to the end
  ran = 0;
  common::run_tasks({});
  common::run_tasks({[&ran] { ++ran; }});
  EXPECT_EQ(ran.load(), 1);
}
