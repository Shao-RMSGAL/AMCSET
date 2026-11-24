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

/*! \file amcset_units.h
 *  \brief Dimensional analysis for AMCSET: Boost.Units quantity types, unit
 *  constants and physical constants.
 *
 * Conventions (read this before writing physics code):
 *
 *  - Every physical value is stored as a Boost.Units quantity in the SI
 *    system (`length_quantity` is metres, `energy_quantity` joules, ...).
 *    Mixing incompatible dimensions is a compile error.
 *  - Non-SI units (angstrom, eV, keV, amu, g/cm^3, ...) are provided as
 *    *unit quantities*: multiply to construct (`3.0 * angstrom`) and divide to
 *    extract a plain number (`(x / angstrom).value()` or `in_units(x,
 *    angstrom)`). Conversions therefore never appear as bare numeric factors.
 *  - Number densities are counts per volume (dimension length^-3). There is no
 *    "amount of substance" anywhere in the transport code.
 *  - Physical constants come from the CODATA values shipped with Boost.Units
 *    and are combined symbolically (e.g. `coulomb_energy_length` is
 *    e^2/(4 pi epsilon_0), the familiar 14.40 eV*Angstrom).
 *  - Dimensionless reduced quantities (reduced energy epsilon, reduced impact
 *    parameter b, ...) are plain `double` or `dimensionless_quantity`.
 *  - Empirical formulas whose published coefficients assume particular units
 *    (e.g. "T in eV") must express that assumption with these unit constants
 *    and say so in a comment.
 */

#pragma once

#include <boost/math/constants/constants.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
#include <boost/units/pow.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/physico-chemical_constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>

namespace amcset {
namespace common {

namespace si = boost::units::si;
namespace constants = boost::units::si::constants::codata;

// ---------------------------------------------------------------------------
// Units and quantity types
// ---------------------------------------------------------------------------

//! Unit of number density (per cubic metre).
using number_density_unit =
    boost::units::power_typeof_helper<si::length,
                                      boost::units::static_rational<-3>>::type;
//! Unit of an inverse length (per metre).
using inverse_length_unit =
    boost::units::power_typeof_helper<si::length,
                                      boost::units::static_rational<-1>>::type;
//! Unit of energy times length (e.g. e^2 / (4 pi epsilon_0)).
using energy_length_unit =
    boost::units::multiply_typeof_helper<si::energy, si::length>::type;
//! Unit of a stopping cross section (energy times area per atom).
using stopping_cross_section_unit =
    boost::units::multiply_typeof_helper<si::energy, si::area>::type;
//! Unit of a stopping power (energy per length).
using stopping_power_unit =
    boost::units::divide_typeof_helper<si::energy, si::length>::type;
//! Unit of an inverse energy (e.g. reduced energy per laboratory energy).
using inverse_energy_unit =
    boost::units::power_typeof_helper<si::energy,
                                      boost::units::static_rational<-1>>::type;
//! Unit of energy straggling per path length (energy^2 per length).
using energy_squared_per_length_unit = boost::units::divide_typeof_helper<
    boost::units::multiply_typeof_helper<si::energy, si::energy>::type,
    si::length>::type;
//! Unit of mass density.
using mass_density_unit = si::mass_density;
//! Unit of velocity.
using velocity_unit = si::velocity;

using length_quantity = boost::units::quantity<si::length>; //!< metres
using area_quantity = boost::units::quantity<si::area>;     //!< m^2
using volume_quantity = boost::units::quantity<si::volume>; //!< m^3
using energy_quantity = boost::units::quantity<si::energy>; //!< joules
using mass_quantity = boost::units::quantity<si::mass>;     //!< kilograms
using time_quantity = boost::units::quantity<si::time>;     //!< seconds
using angle_quantity = boost::units::quantity<si::plane_angle>; //!< radians
using temperature_quantity = boost::units::quantity<si::temperature>; //!< kelvin
using charge_quantity = boost::units::quantity<si::electric_charge>;
using dimensionless_quantity = boost::units::quantity<si::dimensionless>;
using velocity_quantity = boost::units::quantity<velocity_unit>; //!< m/s
using mass_density_quantity = boost::units::quantity<mass_density_unit>;
using number_density_quantity = boost::units::quantity<number_density_unit>;
using inverse_length_quantity = boost::units::quantity<inverse_length_unit>;
using inverse_energy_quantity = boost::units::quantity<inverse_energy_unit>;
using energy_length_quantity = boost::units::quantity<energy_length_unit>;
using stopping_cross_section_quantity =
    boost::units::quantity<stopping_cross_section_unit>;
using stopping_power_quantity = boost::units::quantity<stopping_power_unit>;
using energy_squared_per_length_quantity =
    boost::units::quantity<energy_squared_per_length_unit>;

// ---------------------------------------------------------------------------
// Unit quantities (multiply to construct, divide to extract)
// ---------------------------------------------------------------------------

inline constexpr double pi_v = boost::math::double_constants::pi;

inline constexpr length_quantity meter = 1.0 * si::meter;
inline constexpr length_quantity centimeter = 1e-2 * si::meter;
inline constexpr length_quantity micrometer = 1e-6 * si::meter;
inline constexpr length_quantity nanometer = 1e-9 * si::meter;
inline constexpr length_quantity angstrom = 1e-10 * si::meter;

inline constexpr energy_quantity joule = 1.0 * si::joule;
//! One electron volt, from the CODATA elementary charge.
inline constexpr energy_quantity electron_volt = constants::e * si::volt;
inline constexpr energy_quantity kilo_electron_volt = 1e3 * electron_volt;
inline constexpr energy_quantity mega_electron_volt = 1e6 * electron_volt;
inline constexpr energy_quantity giga_electron_volt = 1e9 * electron_volt;

inline constexpr mass_quantity kilogram = 1.0 * si::kilogram;
inline constexpr mass_quantity gram = 1e-3 * si::kilogram;
//! Unified atomic mass unit (CODATA).
inline constexpr mass_quantity atomic_mass_unit = constants::m_u;

inline constexpr angle_quantity radian = 1.0 * si::radian;
inline constexpr angle_quantity degree = (pi_v / 180.0) * si::radian;

inline constexpr temperature_quantity kelvin = 1.0 * si::kelvin;

inline constexpr volume_quantity cubic_centimeter =
    centimeter * centimeter * centimeter;
inline constexpr area_quantity square_centimeter = centimeter * centimeter;
inline constexpr area_quantity square_angstrom = angstrom * angstrom;

inline constexpr mass_density_quantity kg_per_cubic_meter =
    1.0 * si::kilogram / (si::meter * si::meter * si::meter);
inline constexpr mass_density_quantity gram_per_cubic_centimeter =
    gram / cubic_centimeter;

//! One atom per cubic metre.
inline constexpr number_density_quantity per_cubic_meter =
    1.0 / (si::meter * si::meter * si::meter);
//! Alias kept for readability where "atoms" is the natural word.
inline constexpr number_density_quantity atoms_per_cubic_meter =
    per_cubic_meter;
//! One atom per cubic centimetre.
inline constexpr number_density_quantity per_cubic_centimeter =
    1.0 / cubic_centimeter;
//! One atom per cubic angstrom.
inline constexpr number_density_quantity per_cubic_angstrom =
    1.0 / (angstrom * angstrom * angstrom);

inline constexpr velocity_quantity meters_per_second =
    1.0 * si::meter / si::second;

inline constexpr stopping_power_quantity electron_volt_per_angstrom =
    electron_volt / angstrom;
inline constexpr stopping_power_quantity kilo_electron_volt_per_micrometer =
    kilo_electron_volt / micrometer;

/*!
 * \brief The stopping-cross-section unit used by SRIM: eV / (1e15 atoms/cm^2).
 *
 * A stopping cross section S = -(1/N) dE/dx has dimension energy*area;
 * "per 1e15 atoms/cm^2" is the area 1e-15 cm^2.
 */
inline constexpr stopping_cross_section_quantity srim_stopping_unit =
    electron_volt * (1e-15 * square_centimeter);

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

inline constexpr charge_quantity elementary_charge = constants::e;
inline constexpr mass_quantity electron_mass = constants::m_e;
inline constexpr velocity_quantity speed_of_light = constants::c;

//! e^2 / (4 pi epsilon_0) = 14.3996 eV*Angstrom.
inline constexpr energy_length_quantity coulomb_energy_length =
    elementary_charge * elementary_charge /
    (4.0 * pi_v * constants::epsilon_0);

//! Fine-structure constant alpha = e^2 / (4 pi epsilon_0 hbar c).
inline constexpr dimensionless_quantity fine_structure_constant =
    coulomb_energy_length / (constants::hbar * constants::c);

//! Reduced Planck constant.
inline constexpr auto reduced_planck_constant = constants::hbar;

//! Boltzmann constant k_B.
inline constexpr auto boltzmann_constant = constants::k_B;

//! Bohr radius a_0 = hbar / (m_e c alpha) = 0.529177 Angstrom.
inline constexpr length_quantity bohr_radius =
    constants::hbar / (constants::m_e * constants::c * fine_structure_constant);

//! Bohr velocity v_0 = alpha c = e^2/(4 pi epsilon_0 hbar).
inline constexpr velocity_quantity bohr_velocity =
    fine_structure_constant * constants::c;

//! Electron rest energy m_e c^2 (0.511 MeV).
inline constexpr energy_quantity electron_rest_energy =
    constants::m_e * constants::c * constants::c;

//! Atomic mass unit rest energy m_u c^2 (931.494 MeV).
inline constexpr energy_quantity atomic_mass_unit_energy =
    constants::m_u * constants::c * constants::c;

//! Rydberg energy (13.606 eV) = e^2/(8 pi epsilon_0 a_0).
inline constexpr energy_quantity rydberg_energy =
    coulomb_energy_length / (2.0 * bohr_radius);

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/*!
 * \brief Express \p q as a plain number in the unit \p unit.
 *
 * `in_units(x, angstrom)` is `(x / angstrom).value()` with a compile-time
 * check that both have the same dimension.
 */
template <class Unit, class Y>
constexpr double in_units(const boost::units::quantity<Unit, Y> &q,
                          const boost::units::quantity<Unit, Y> &unit) {
  return static_cast<double>(q.value() / unit.value());
}

//! Plain value of a dimensionless quantity.
constexpr double value_of(const dimensionless_quantity &q) {
  return q.value();
}

} // namespace common
} // namespace amcset
