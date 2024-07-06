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
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

/*! \file amcset_utilities.h
 *  \brief The utilites header file for the AMCSET common library
 *
 * This file contains useful constructs and includes for AMCSET,
 * such as Boost.units definitions, custom units and quantities, as well as
 * other useful constants.
 */

#pragma once

// Boost units
#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/units/io.hpp>
#include <boost/units/physical_dimensions/energy.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/systems/si/electric_potential.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/plane_angle.hpp>

// Boost math
#include <boost/math/constants/constants.hpp>
#include <iostream>
#include <string>

namespace amcset {
namespace common {

// Units definitions
using namespace boost::units;
namespace si = boost::units::si;
namespace metric = boost::units::metric;
namespace constants = boost::units::si::constants::codata;

using namespace boost::math::double_constants;

using length_quantity = quantity<si::length>;  //!< SI type for length (meter)
using energy_quantity = quantity<si::energy>;  //!< SI type for energy (Joule)
using angle_quantity =
    quantity<si::plane_angle>;             //!< SI type for plane angle (Radian)
using mass_quantity = quantity<si::mass>;  //!< SI type for mass (Kilogram)
using charge_quantity =
    quantity<si::electric_charge>;  //!< SI type for charge (Coulomb)
using dimensionless_quantity =
    quantity<si::dimensionless>;  //!< Type for dimensionless quantity

constexpr auto kilo_electron_volt =
    double(1000) * constants::e * si::volt;  //!< Quantity for kiloelectron volt
constexpr auto angstrom =
    metric::angstrom_base_unit::unit_type();  //!< Quantity for angstrom
constexpr auto radian =
    angle::radian_base_unit::unit_type();  //!< Quantity for radian
constexpr auto atomic_mass_unit =
    constants::m_u;  //!< Quantity for atomic mass unit
constexpr auto elementary_charge =
    constants::e;  //!< Quantity for the elementary charge

/*!
 * \brief A macro to provide debugging details for exception messages.
 *
 * This macro takes a message and appends the function name, file name, and line
 * number to the end of the message.
 */
#define EXCEPTION_MESSAGE(msg)                                           \
  std::string(msg) + " [Function: " + __func__ + ", File: " + __FILE__ + \
      ", Line: " + std::to_string(__LINE__) + "]"

/*!
 * \brief A helper function that converts a quantity to a string with the value
 * and units.
 *
 * \param quantity The quantity to convert to a string
 *
 * \return A string containing the value and unit.
 */
template <typename Unit, typename T>
std::string to_string_with_unit(const quantity<Unit, T>& quantity) {
  std::ostringstream oss;
  oss << quantity;
  return oss.str();
}

/*!
 * \brief A nested exception handling throwing function.
 *
 * This will show the trace of functions which catch and throw exceptions as
 * they propogate through the call stack. It is useful for debugging, as it
 * shows the call stack of functions in a human-readable way.
 *
 * Code taken from Richard Hodges:
 * https://stackoverflow.com/questions/37227300/why-doesnt-c-use-stdnested-exception-to-allow-throwing-from-destructor/37227893#37227893
 *
 * \param args Messages to include in the exception invocation.
 */
template <class... Args>
[[noreturn]] void rethrow(Args&&... args) {
  std::ostringstream ss;
  using expand = int[];
  std::string sep;
  void(expand{0, ((ss << sep << args), sep = ", ", 0)...});
  try {
    std::rethrow_exception(std::current_exception());
  } catch (const std::invalid_argument& e) {
    std::throw_with_nested(std::invalid_argument(ss.str()));
  } catch (const std::out_of_range& e) {
    std::throw_with_nested(std::out_of_range(ss.str()));
  } catch (const std::logic_error& e) {
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
 * \param depth The depth of exception messages. This starts out at for the
 * top-leve exception.
 */
inline void print_exception(const std::exception& e, std::size_t depth = 0) {
  std::cerr << "exception: " << std::string(depth, ' ') << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception& nested) {
    print_exception(nested, depth + 1);
  }
}

}  // namespace common
}  // namespace amcset
