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

/*! \file amcset_utilities.h
 *  \brief The utilities header file for the AMCSET common library
 *
 * This file contains useful constructs for AMCSET: it includes the
 * dimensional-analysis header amcset_units.h and defines the exception
 * helpers (EXCEPTION_MESSAGE, rethrow, print_exception).
 *
 * \note This header brings \c boost::units and
 * \c boost::math::double_constants into namespace amcset::common with
 * using-directives, so names such as \c pi, \c sqrt and \c pow are
 * visible unqualified throughout the library.
 *
 * Units, quantity types and physical constants live in amcset_units.h
 * (included here); see its file comment for the conventions.
 */

#pragma once

// Standard library
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// Dimensional analysis (quantity types, unit constants, physical constants)
#include "amcset_units.h"

#include <boost/math/constants/constants.hpp>

namespace amcset {
namespace common {

// The core library historically relies on unqualified Boost.Units functions
// (sqrt, pow, root, exp, ...) and Boost math constants (pi). New code should
// prefer the explicit names from amcset_units.h (pi_v, in_units, ...).
using namespace boost::units;
using namespace boost::math::double_constants;

//! Force (energy per length). Used for derivatives of interatomic potentials.
using force_quantity = boost::units::quantity<si::force>;

/*!
 * \brief A macro to provide debugging details for exception messages.
 *
 * This macro takes a message and appends the function name, file name, and line
 * number to the end of the message.
 */
#define EXCEPTION_MESSAGE(msg)                                                 \
  std::string(msg) + " [Function: " + __func__ + ", File: " + __FILE__ +       \
      ", Line: " + std::to_string(__LINE__) + "]"

/*!
 * \brief A nested exception handling throwing function.
 *
 * This will show the trace of functions which catch and throw exceptions as
 * they propagate through the call stack. It is useful for debugging, as it
 * shows the call stack of functions in a human-readable way.
 *
 * Code taken from Richard Hodges:
 * https://stackoverflow.com/questions/37227300/why-doesnt-c-use-stdnested-exception-to-allow-throwing-from-destructor/37227893#37227893
 *
 * Must be called from inside a catch block. The currently handled exception
 * is nested inside a new exception whose type mirrors it (invalid_argument,
 * out_of_range and logic_error are preserved; anything else becomes
 * runtime_error) and whose message is \p args joined with ", ".
 *
 * \param args Messages to include in the exception invocation.
 */
template <class... Args> [[noreturn]] void rethrow(Args &&...args) {
  std::ostringstream ss;
  using expand = int[];
  std::string sep;
  void(expand{0, ((ss << sep << args), sep = ", ", 0)...});
  try {
    std::rethrow_exception(std::current_exception());
  } catch (const std::invalid_argument &) {
    std::throw_with_nested(std::invalid_argument(ss.str()));
  } catch (const std::out_of_range &) {
    std::throw_with_nested(std::out_of_range(ss.str()));
  } catch (const std::logic_error &) {
    std::throw_with_nested(std::logic_error(ss.str()));
  } catch (...) {
    std::throw_with_nested(std::runtime_error(ss.str()));
  }
}

/*!
 * \brief A helper function that prints out nested exceptions.
 *
 * Code taken from Richard Hodges:
 * https://stackoverflow.com/questions/37227300/why-doesnt-c-use-stdnested-exception-to-allow-throwing-from-destructor/37227893#37227893
 *
 * \param e The exception to print out.
 *
 * \param depth The nesting depth of \p e (0 for the top-level exception).
 */
inline void print_exception(const std::exception &e, std::size_t depth = 0) {
  std::cerr << "exception: " << std::string(depth, ' ') << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &nested) {
    print_exception(nested, depth + 1);
  }
}

} // namespace common
} // namespace amcset
