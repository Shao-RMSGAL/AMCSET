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

/*! \file dynlib.h
 *  \brief Loading a shared library and its functions at run time (the GPU
 *  frameworks come from the system's driver, not from AMCSET's build).
 */

#pragma once

#include <initializer_list>
#include <string>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace amcset {
namespace compute {

//! A shared library, closed when destroyed.
class DynamicLibrary {
public:
  DynamicLibrary() = default;
  //! The first of \p names that loads (empty if none).
  explicit DynamicLibrary(std::initializer_list<const char *> names) {
    for (const char *n : names) {
      if (open(n)) break;
    }
  }
  ~DynamicLibrary() { close(); }
  DynamicLibrary(const DynamicLibrary &) = delete;
  DynamicLibrary &operator=(const DynamicLibrary &) = delete;

  bool open(const std::string &name) {
    close();
#if defined(_WIN32)
    // A full path: its dependencies (hiprtc's builtins, comgr) load from its folder.
    const bool path = name.find_first_of("\\/") != std::string::npos;
    handle_ = reinterpret_cast<void *>(path ? LoadLibraryExA(name.c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH)
                                            : LoadLibraryA(name.c_str()));
#else
    handle_ = dlopen(name.c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
    if (handle_ != nullptr) name_ = name;
    return handle_ != nullptr;
  }
  void close() {
    if (handle_ == nullptr) return;
#if defined(_WIN32)
    FreeLibrary(reinterpret_cast<HMODULE>(handle_));
#else
    dlclose(handle_);
#endif
    handle_ = nullptr;
  }
  bool loaded() const noexcept { return handle_ != nullptr; }
  const std::string &name() const noexcept { return name_; }

  //! The function \p symbol, or null.
  template <class F> F get(const char *symbol) const {
    if (handle_ == nullptr) return nullptr;
#if defined(_WIN32)
    return reinterpret_cast<F>(reinterpret_cast<void *>(GetProcAddress(reinterpret_cast<HMODULE>(handle_), symbol)));
#else
    return reinterpret_cast<F>(dlsym(handle_, symbol));
#endif
  }

private:
  void *handle_ = nullptr;
  std::string name_;
};

} // namespace compute
} // namespace amcset
