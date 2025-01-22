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

// Standard library
#include <cmath>

// Boost units
#include <boost/mpl/list.hpp>
#include <boost/units/base_units/cgs/gram.hpp>
#include <boost/units/base_units/metric/angstrom.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
#include <boost/units/physical_dimensions/energy.hpp>
#include <boost/units/physical_dimensions/mass.hpp>
#include <boost/units/systems/cgs/length.hpp>
#include <boost/units/systems/cgs/mass.hpp>
#include <boost/units/systems/si/amount.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/electric_potential.hpp>
#include <boost/units/systems/si/length.hpp>
#include <boost/units/systems/si/mass.hpp>
#include <boost/units/systems/si/plane_angle.hpp>

// Boost math
#include <boost/math/constants/constants.hpp>
#include <boost/units/unit.hpp>
#include <iostream>
#include <string>

namespace amcset {
namespace common {

// Units definitions
using namespace boost::units;
namespace si = boost::units::si;
// namespace metric = boost::units::metric;
namespace constants = boost::units::si::constants::codata;

using namespace boost::math::double_constants;

using length_quantity = quantity<si::length>; //!< SI type for length (meter)
using energy_quantity = quantity<si::energy>; //!< SI type for energy (Joule)
using angle_quantity =
    quantity<si::plane_angle>;            //!< SI type for plane angle (Radian)
using mass_quantity = quantity<si::mass>; //!< SI type for mass (Kilogram)
using charge_quantity =
    quantity<si::electric_charge>; //!< SI type for charge (Coulomb)
using dimensionless_quantity =
    quantity<si::dimensionless>; //!< Type for dimensionless quantity
using volume_quantity = quantity<si::volume>; //!< Type for volume
using mass_density_dimension = make_dimension_list<
    boost::mpl::list2<dim<mass_base_dimension, static_rational<1>>,
                      dim<length_base_dimension, static_rational<-3>>>>::
    type; //!< Custom dimension for mass density
using statcoulomb_dimension = make_dimension_list<
    boost::mpl::list3<dim<length_base_dimension, static_rational<3, 2>>,
                      dim<mass_base_dimension, static_rational<1, 2>>,
                      dim<time_base_dimension, static_rational<-1>>>>::
    type; //!< Custom dimension statcoulombs

using mass_density_unit =
    unit<mass_density_dimension, si::system>; //!< Unit for mass density
using mass_density_quantity =
    quantity<mass_density_unit>; //!< Type for mass density (kg/m^3)
using number_density_unit =
    divide_typeof_helper<si::amount, si::volume>::type; //!< Custom unit
                                                        //!< for number density
using velocity_unit =
    divide_typeof_helper<si::length,
                         si::time>::type; //!< Custom unit for velocity
using voltage_unit =
    divide_typeof_helper<si::energy,
                         si::length>::type; //!< Custom unit for voltage
using number_density_quantity =
    quantity<number_density_unit>; //!< Unit for number density, (mol/m^3)
using number_density_quantity =
    quantity<number_density_unit>;                 //!< Number density
using velocity_quantity = quantity<velocity_unit>; //!< Velocity
using statcoulomb_quantity =
    quantity<statcoulomb_dimension>;             //!< Unit for statcoulomb,
                                                 //!< (m^(3/2)*kg^(1/2)/s)
using voltage_quantity = quantity<voltage_unit>; //!< Voltage quantity

constexpr auto electron_volt =
    constants::e * si::volt; //!< Quantity for electron volt
constexpr auto kilo_electron_volt =
    1000.0 * electron_volt; //!< Quantity for kiloelectron volt

constexpr auto angstrom = 1e-10 * si::meter; //!< Quantity for angstrom
constexpr auto radian =
    angle::radian_base_unit::unit_type(); //!< Quantity for radian
constexpr auto atomic_mass_unit =
    constants::m_u; //!< Quantity for atomic mass unit
constexpr auto elementary_charge =
    constants::e; //!< Quantity for the elementary charge
constexpr mass_density_quantity kg_per_cubic_meter =
    1.0 * si::kilogram /
    pow<3>(si::meter); //!< Unit quantity for mass density (kg/m^3)
constexpr number_density_quantity atoms_per_cubic_meter =
    1.0 * si::mole / pow<3>(si::meter); //!< Unit quantity for number density
                                        //!< (atoms/m^3)
constexpr velocity_quantity meters_per_second =
    1.0 * si::meter / si::second; //!< Unit quantity for velocity (m/s)

constexpr auto fine_structure_constant =
    si::constants::codata::e * si::constants::codata::e /
    (2.0 * si::constants::codata::epsilon_0 * si::constants::codata::h *
     si::constants::codata::c); //<! Fine structure constant

constexpr auto reduced_planck_constant =
    si::constants::codata::h / (2 * pi); //<! Reduced planck constant

constexpr length_quantity bohr_radius =
    reduced_planck_constant /
    (si::constants::codata::m_e * si::constants::codata::c *
     fine_structure_constant); //<! The bohr radius (m)

// sqrt is constexpr as of C++26, but clangd does not recognize this. TODO:
// Exchange const for constexpr when clang no longer complains about this. (Note
// that this will compile under gcc regardless.
const auto statcoulomb = pow<static_rational<3, 2>>(si::meter) / 1000.0 *
                         std::sqrt(1.0 / 1000.0) *
                         pow<static_rational<1, 2>>(si::kilogram) / si::second;
const auto e_statcoulomb = 4.80320425e-10 * statcoulomb;

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
 * \param args Messages to include in the exception invocation.
 */
template <class... Args> [[noreturn]] void rethrow(Args &&...args) {
  std::ostringstream ss;
  using expand = int[];
  std::string sep;
  void(expand{0, ((ss << sep << args), sep = ", ", 0)...});
  try {
    std::rethrow_exception(std::current_exception());
  } catch (const std::invalid_argument &e) {
    std::throw_with_nested(std::invalid_argument(ss.str()));
  } catch (const std::out_of_range &e) {
    std::throw_with_nested(std::out_of_range(ss.str()));
  } catch (const std::logic_error &e) {
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
