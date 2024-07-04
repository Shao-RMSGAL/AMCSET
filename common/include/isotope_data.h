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
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

/*! \file isotope_data.h
 * \brief Isotopic information for elements up to Z number 118.
 *
 * This information was taken from the National Institute of Science and
 * Technology (NIST) Atomic Weights and Isotopic Compositions with Relative
 * Atomic Masses database. At the time of writing, the data was aquired from
 * https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses.
 */

// Standard library
#include <array>
#include <stdexcept>
#include <string>

// AMCSET includes
#include "amcset_utilities.h"
#include "isotope_mass_data.h"

namespace amcset {
namespace common {
/*!
 * \class IsotopeData
 * \brief Class for allowing static access to isotopic data given a mass and z
 * number.
 *
 * This class enables static access to isotopic data with minimal overhead. By
 * calling the IsotopeData::getIsotopicMass() function, one can access the exact
 * mass of any known isotope (As tabulated by NIST)
 */
class IsotopeData {
 public:
  static constexpr size_t MAX_Z = 118;  //!< The maximum allowed atomic number

  static constexpr size_t MAX_A = 295;  //!< The maximum allowed mass number

  /*!
   * \brief A function that takes a Z number and mass number, and returns the
   * isotopic mass.
   *
   * This function will access tabulated data stored within IsotopeData and
   * statically return the value of an isotope's mass.
   *
   * \param atomic_number The atomic (Z) number of the isotope
   * \param mass_number The mass number of the isotope
   * \return The exact mass of the corresponding isotope in Relative Atomic Mass
   */
  static constexpr mass_quantity getIsotopeMass(size_t atomic_number,
                                                size_t mass_number) {
    if (atomic_number < 1 || atomic_number > MAX_Z) {
      throw std::out_of_range(
          "Invalid atomic number: " + std::to_string(atomic_number) +
          ". Allowed range is 1 to " + std::to_string(MAX_Z) + ".");
    }
    if (mass_number < 1 || mass_number > MAX_A) {
      throw std::out_of_range(
          "Invalid mass number: " + std::to_string(atomic_number) +
          ". Allowed range is 1 to " + std::to_string(MAX_Z) + ".");
    }

    mass_quantity mass = isotopic_masses[atomic_number - 1][mass_number - 1];

    if (mass == double(0.0) * atomic_mass_unit) {
      throw std::invalid_argument(
          "Isotope does not exist. Atomic number: " +
          std::to_string(atomic_number) +
          ". Mass number: " + std::to_string(mass_number) + ".");
    }
    return mass;
  };

 private:
  // Fancy way for instantiating a compile-time static data structure.
  static constexpr std::array<std::array<mass_quantity, MAX_A>, MAX_Z>
      isotopic_masses = []() {
        std::array<std::array<mass_quantity, MAX_A>, MAX_Z> masses = {};
        GENERATE_ISOTOPE_MASS_DATA(masses, atomic_mass_unit);
        return masses;
      }();
};
}  // namespace common
}  // namespace amcset
