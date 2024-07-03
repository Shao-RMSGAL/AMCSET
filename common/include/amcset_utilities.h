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

}  // namespace common
}  // namespace amcset
