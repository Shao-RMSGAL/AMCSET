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

/*! \file text.h
 *  \brief Fast text output of numbers: the printf conversions the output
 *  files use (%*.*E, %*.*f, %*zu), written into a caller's buffer without
 *  printf or allocation, and a buffered writer for long files.
 *
 * The output must not change: std::to_chars with a precision is specified as
 * printf with that precision in the "C" locale (both print the exactly
 * rounded decimal value), and the tests compare every helper with snprintf
 * (OutputText.*). Non-finite values go through snprintf itself, so even
 * their spelling ("inf", "-nan(ind)", ...) is printf's.
 */

#pragma once

#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

namespace amcset {
namespace output {
namespace text {

//! Room every helper below may use (fixed notation of 1e308 with 30
//! decimals is 341 characters).
inline constexpr std::size_t max_field = 400;

namespace detail {
//! Right-align [s, s + n) in a field of \p width (printf's default).
inline char *pad_left(char *p, const char *s, std::size_t n, int width, char fill = ' ') noexcept {
  for (int i = static_cast<int>(n); i < width; ++i) *p++ = fill;
  std::memcpy(p, s, n);
  return p + n;
}
} // namespace detail

/*!
 * \brief printf("%*.*E", width, precision, v) (or "%*.*e" with \p upper
 * false) at \p p; returns the end. \p p must have max_field characters.
 */
inline char *put_sci(char *p, double v, int width, int precision, bool upper = true) {
  char tmp[max_field];
  if (!std::isfinite(v)) {
    const int n = std::snprintf(tmp, sizeof tmp, upper ? "%.*E" : "%.*e", precision, v);
    return detail::pad_left(p, tmp, static_cast<std::size_t>(n > 0 ? n : 0), width);
  }
  const auto r = std::to_chars(tmp, tmp + sizeof tmp, v, std::chars_format::scientific, precision);
  if (upper) {
    for (char *q = tmp; q != r.ptr; ++q) {
      if (*q == 'e') *q = 'E';
    }
  }
  return detail::pad_left(p, tmp, static_cast<std::size_t>(r.ptr - tmp), width);
}

//! printf("%*.*f", width, precision, v) at \p p; returns the end.
inline char *put_fixed(char *p, double v, int width, int precision) {
  char tmp[max_field];
  if (!std::isfinite(v)) {
    const int n = std::snprintf(tmp, sizeof tmp, "%.*f", precision, v);
    return detail::pad_left(p, tmp, static_cast<std::size_t>(n > 0 ? n : 0), width);
  }
  const auto r = std::to_chars(tmp, tmp + sizeof tmp, v, std::chars_format::fixed, precision);
  return detail::pad_left(p, tmp, static_cast<std::size_t>(r.ptr - tmp), width);
}

//! printf("%*zu", width, v), or "%0*zu" with \p fill '0'; returns the end.
inline char *put_unsigned(char *p, unsigned long long v, int width, char fill = ' ') noexcept {
  char tmp[24];
  const auto r = std::to_chars(tmp, tmp + sizeof tmp, v);
  return detail::pad_left(p, tmp, static_cast<std::size_t>(r.ptr - tmp), width, fill);
}

/*!
 * \brief Collects text in a large buffer and writes it to a stream in
 * chunks (a stream insertion per number or line costs more than the
 * formatting). Typical use: `char *q = w.reserve(n); q = put_sci(q, ...);
 * w.commit(q);` with n the most the fields can take (max_field each).
 * Call finish() at the end; it throws std::runtime_error if the stream
 * failed. The destructor writes what is left but never throws.
 */
class BufferedWriter {
public:
  explicit BufferedWriter(std::ostream &o, std::size_t capacity = std::size_t{1} << 20)
      : o_(o), cap_(capacity < 64 * max_field ? 64 * max_field : capacity), buf_(new char[cap_]) {}
  BufferedWriter(const BufferedWriter &) = delete;
  BufferedWriter &operator=(const BufferedWriter &) = delete;
  ~BufferedWriter() {
    try {
      flush();
    } catch (...) {
      // finish() reports errors; a destructor must not throw.
    }
  }

  //! Room for \p n characters at the end of the buffer (written out first
  //! if needed); write at most n there, then commit() the end.
  char *reserve(std::size_t n) {
    if (n > cap_) throw std::length_error("BufferedWriter: more than its capacity at once");
    if (pos_ + n > cap_) flush();
    return buf_.get() + pos_;
  }
  //! The end of what was written into reserve()'s room.
  void commit(const char *end) noexcept { pos_ = static_cast<std::size_t>(end - buf_.get()); }
  void put(const char *s, std::size_t n) {
    if (n > cap_) {
      flush();
      o_.write(s, static_cast<std::streamsize>(n));
      return;
    }
    char *q = reserve(n);
    std::memcpy(q, s, n);
    commit(q + n);
  }
  void put(const std::string &s) { put(s.data(), s.size()); }
  void put(const char *s) { put(s, std::strlen(s)); }
  void put(char c) {
    char *q = reserve(1);
    *q = c;
    commit(q + 1);
  }
  void flush() {
    if (pos_ > 0) o_.write(buf_.get(), static_cast<std::streamsize>(pos_));
    pos_ = 0;
  }
  void finish() {
    flush();
    o_.flush();
    if (!o_) throw std::runtime_error("writing an output file failed (disk full?)");
  }

private:
  std::ostream &o_;
  std::size_t cap_;
  std::unique_ptr<char[]> buf_;
  std::size_t pos_ = 0;
};

} // namespace text
} // namespace output
} // namespace amcset
