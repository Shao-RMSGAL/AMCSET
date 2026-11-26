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

#include "data_files.h"

#include <cstdlib>
#include <system_error>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#endif

namespace amcset::data {
namespace {

std::filesystem::path executable_directory() {
  std::error_code ec;
#if defined(_WIN32)
  std::wstring buf(32768, L'\0');
  const DWORD n = GetModuleFileNameW(nullptr, buf.data(), static_cast<DWORD>(buf.size()));
  if (n == 0 || n >= buf.size()) return {};
  buf.resize(n);
  return std::filesystem::path(buf).parent_path();
#elif defined(__APPLE__)
  uint32_t size = 0;
  _NSGetExecutablePath(nullptr, &size);
  std::string buf(size, '\0');
  if (_NSGetExecutablePath(buf.data(), &size) != 0) return {};
  return std::filesystem::weakly_canonical(std::filesystem::path(buf.c_str()), ec).parent_path();
#else
  const auto p = std::filesystem::read_symlink("/proc/self/exe", ec);
  return ec ? std::filesystem::path{} : p.parent_path();
#endif
}

} // namespace

std::vector<std::filesystem::path> search_directories() {
  std::vector<std::filesystem::path> dirs;
  if (const char *env = std::getenv("AMCSET_DATA_DIR"); env && *env) dirs.emplace_back(env);
  const auto exe = executable_directory();
  if (!exe.empty()) {
    dirs.push_back(exe / "data");
    dirs.push_back(exe.parent_path() / "data");
    dirs.push_back(exe.parent_path().parent_path() / "data");
    dirs.push_back(exe.parent_path() / "share" / "amcset");
  }
#ifdef AMCSET_SOURCE_DATA_DIR
  dirs.emplace_back(AMCSET_SOURCE_DATA_DIR);
#endif
  return dirs;
}

std::filesystem::path find(std::string_view relative) {
  std::error_code ec;
  for (const auto &d : search_directories()) {
    const auto p = d / std::filesystem::path(relative);
    if (std::filesystem::is_regular_file(p, ec)) return p;
  }
  return {};
}

std::string bundled_srim_stopping_table() { return find(srim_stopping_table).string(); }

std::filesystem::path find_electron_table(std::string_view file) {
  const auto relative = std::filesystem::path(electron_tables) / std::filesystem::path(file);
  if (auto p = find(relative.generic_string()); !p.empty()) return p;
#ifdef AMCSET_SOURCE_ASSETS_DIR
  std::error_code ec;
  const auto p = std::filesystem::path(AMCSET_SOURCE_ASSETS_DIR) / std::filesystem::path(file);
  if (std::filesystem::is_regular_file(p, ec)) return p;
#endif
  return {};
}

} // namespace amcset::data
