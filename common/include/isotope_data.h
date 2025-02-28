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
 * Atomic Masses database. At the time of writing, the data was acquired from
 * https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses.
 */

#pragma once

// Boost units include
#include <boost/units/systems/si/codata/electron_constants.hpp>

// Standard library
#include <array>
#include <string>

// AMCSET includes
#include "amcset_utilities.h"
#include "isotope_database.h"

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
 *
 * This class also provides the natural isotopic abundances of various elements.
 */
class IsotopeData {
public:
  /*!
   * \struct Isotope
   * \brief Stores information about isotopes
   */
  struct Isotope {
    size_t atomic_number;
    std::string_view atomic_symbol;
    size_t mass_number;
    double relative_atomic_mass;
    double isotopic_composition;
    double standard_atomic_weight;
  };

  static constexpr size_t MAX_Z = 118; //!< The maximum allowed atomic number

  static constexpr size_t MAX_A = 295; //!< The maximum allowed mass number

  /*!
   * \brief A function that takes a Z number and mass number, and returns the
   * isotopic mass.
   *
   * This function will access tabulated data stored within IsotopeData and
   * statically return the value of an isotope's mass. If atomic_number and
   * mass_number are specifed as 0, the mass of an electron is returned.
   *
   * \param atomic_number The atomic (Z) number of the isotope. If this number
   * is zero and mass_number is 0, electron mass is returned.
   * \param mass_number The mass number of the isotope
   * \return The exact mass of the corresponding isotope in Relative Atomic Mass
   */
  inline static constexpr mass_quantity
  get_isotope_mass(size_t atomic_number, size_t mass_number) try {
    if (atomic_number == 0 && mass_number == 0) {
      return constants::m_e;
    }

    if (atomic_number == 0 || atomic_number > MAX_Z) {
      throw std::out_of_range(EXCEPTION_MESSAGE(
          "Invalid atomic number: " + std::to_string(atomic_number) +
          ". Allowed range is 1 to " + std::to_string(MAX_Z) +
          " for non-electrons."));
    }
    if (mass_number == 0 || mass_number > MAX_A) {
      throw std::out_of_range(EXCEPTION_MESSAGE(
          "Invalid mass number: " + std::to_string(atomic_number) +
          ". Allowed range is 1 to " + std::to_string(MAX_Z) +
          " for non-electrons."));
    }

    auto *isotope =
        find_isotope(atomic_number_to_element(atomic_number), mass_number);

    mass_quantity mass = isotope->relative_atomic_mass * atomic_mass_unit;
    // mass_quantity mass = isotopic_masses[atomic_number - 1][mass_number - 1];

    if (mass == double(0.0) * atomic_mass_unit) {
      throw std::invalid_argument(EXCEPTION_MESSAGE(
          "Isotope does not exist. Atomic number: " +
          std::to_string(atomic_number) +
          ". Mass number: " + std::to_string(mass_number) + "."));
    }
    return mass;
  } catch (...) {
    rethrow();
  }

  /*!
   * \brief A helper function for retrieving isotope data
   */
  inline static constexpr const Isotope *
  find_isotope(std::string_view atomic_symbol, size_t mass_number) {
    for (const auto &entry : database) {
      if (entry.atomic_symbol == atomic_symbol &&
          entry.mass_number == mass_number) {
        return &entry;
      }
    }
    throw std::out_of_range(
        EXCEPTION_MESSAGE("Invalid symbol or mass number "));
  }

  /*!
   * \brief Return a material_vector of a natural sample of an element.
   */
  static inline

#if __cplusplus >= 202302L
      constexpr
#endif
      std::vector<std::pair<size_t, double>>
      get_relative_abundance(size_t atomic_number) {
    size_t lowest_isotope = 0;
    auto composition_sum = 0.0;
    bool found_element = false; // Speed up evaluation.
    std::vector<std::pair<size_t, double>> material_vector = {};

    for (const auto &entry : database) {
      if (found_element && entry.atomic_number != atomic_number) {
        break;
      } // Assumes isotopic data is in increasing mass-number order.
      if (entry.atomic_number == atomic_number) {
        found_element = true;
        if (lowest_isotope == 0) {
          lowest_isotope = entry.mass_number;
        }
        if (entry.isotopic_composition != 0.0) {
          composition_sum += entry.isotopic_composition;
          material_vector.emplace_back(
              std::make_pair(entry.mass_number, entry.isotopic_composition));
        }
      }
    }

    // If there is no isotopic compositions, just return the lowest-mass
    // isotope.
    if (material_vector.empty()) {
      material_vector.emplace_back(std::make_pair(lowest_isotope, 1.0));
      composition_sum = 1.0;
    }

    // Ensure normalization
    if (composition_sum != 1.0) {
      for (auto &entry : material_vector) {
        entry.second /= composition_sum;
      }
    }
    return material_vector;

    // TODO: The user should be notified, and a suggestion made to manually
    // specify isotopes.
  }

  /*!
   *    \brief Converts atomic symbol to atomic number.
   */
  inline static constexpr size_t
  element_to_atomic_number(std::string_view symbol) {
    // Handle single character elements{{{
    if (symbol.length() == 1) {
      switch (symbol[0]) {
      case 'H':
        return 1; // Hydrogen
      case 'B':
        return 5; // Boron
      case 'C':
        return 6; // Carbon
      case 'N':
        return 7; // Nitrogen
      case 'O':
        return 8; // Oxygen
      case 'F':
        return 9; // Fluorine
      case 'P':
        return 15; // Phosphorus
      case 'S':
        return 16; // Sulfur
      case 'K':
        return 19; // Potassium
      case 'V':
        return 23; // Vanadium
      case 'Y':
        return 39; // Yttrium
      case 'I':
        return 53; // Iodine
      case 'W':
        return 74; // Tungsten
      case 'U':
        return 92; // Uranium
      default:
        return 0; // Unknown element
      }
    }
    // Handle two character elements
    else if (symbol.length() == 2) {
      char first = symbol[0];
      char second = symbol[1];

      if (first == 'H') {
        if (second == 'e')
          return 2; // Helium
        else if (second == 'g')
          return 80; // Mercury
        else if (second == 'o')
          return 67; // Holmium
        else if (second == 'f')
          return 72; // Hafnium
        else if (second == 's')
          return 108; // Hassium
      } else if (first == 'L') {
        if (second == 'i')
          return 3; // Lithium
        else if (second == 'a')
          return 57; // Lanthanum
        else if (second == 'u')
          return 71; // Lutetium
        else if (second == 'r')
          return 103; // Lawrencium
        else if (second == 'v')
          return 116; // Livermorium
      } else if (first == 'B') {
        if (second == 'e')
          return 4; // Beryllium
        else if (second == 'a')
          return 56; // Barium
        else if (second == 'i')
          return 83; // Bismuth
        else if (second == 'k')
          return 97; // Berkelium
        else if (second == 'r')
          return 35; // Bromine
        else if (second == 'h')
          return 107; // Bohrium
      } else if (first == 'N') {
        if (second == 'e')
          return 10; // Neon
        else if (second == 'a')
          return 11; // Sodium
        else if (second == 'i')
          return 28; // Nickel
        else if (second == 'b')
          return 41; // Niobium
        else if (second == 'd')
          return 60; // Neodymium
        else if (second == 'p')
          return 93; // Neptunium
        else if (second == 'o')
          return 102; // Nobelium
      } else if (first == 'M') {
        if (second == 'g')
          return 12; // Magnesium
        else if (second == 'n')
          return 25; // Manganese
        else if (second == 'o')
          return 42; // Molybdenum
        else if (second == 'd')
          return 101; // Mendelevium
        else if (second == 't')
          return 109; // Meitnerium
        else if (second == 'c')
          return 115; // Moscovium
      } else if (first == 'A') {
        if (second == 'l')
          return 13; // Aluminum
        else if (second == 'r')
          return 18; // Argon
        else if (second == 's')
          return 33; // Arsenic
        else if (second == 'g')
          return 47; // Silver
        else if (second == 'u')
          return 79; // Gold
        else if (second == 'c')
          return 89; // Actinium
        else if (second == 'm')
          return 95; // Americium
        else if (second == 't')
          return 85; // Astatine
      } else if (first == 'S') {
        if (second == 'i')
          return 14; // Silicon
        else if (second == 'c')
          return 21; // Scandium
        else if (second == 'e')
          return 34; // Selenium
        else if (second == 'r')
          return 38; // Strontium
        else if (second == 'b')
          return 51; // Antimony
        else if (second == 'n')
          return 50; // Tin
        else if (second == 'm')
          return 62; // Samarium
        else if (second == 'g')
          return 106; // Seaborgium
      } else if (first == 'C') {
        if (second == 'l')
          return 17; // Chlorine
        else if (second == 'a')
          return 20; // Calcium
        else if (second == 'r')
          return 24; // Chromium
        else if (second == 'o')
          return 27; // Cobalt
        else if (second == 'u')
          return 29; // Copper
        else if (second == 'd')
          return 48; // Cadmium
        else if (second == 's')
          return 55; // Cesium
        else if (second == 'e')
          return 58; // Cerium
        else if (second == 'f')
          return 98; // Californium
        else if (second == 'm')
          return 96; // Curium
        else if (second == 'n')
          return 112; // Copernicium
      } else if (first == 'T') {
        if (second == 'i')
          return 22; // Titanium
        else if (second == 'c')
          return 43; // Technetium
        else if (second == 'e')
          return 52; // Tellurium
        else if (second == 'a')
          return 73; // Tantalum
        else if (second == 'l')
          return 81; // Thallium
        else if (second == 'h')
          return 90; // Thorium
        else if (second == 'b')
          return 65; // Terbium
        else if (second == 'm')
          return 69; // Thulium
        else if (second == 's')
          return 117; // Tennessine
      } else if (first == 'F') {
        if (second == 'e')
          return 26; // Iron
        else if (second == 'r')
          return 87; // Francium
        else if (second == 'm')
          return 100; // Fermium
        else if (second == 'l')
          return 114; // Flerovium
      } else if (first == 'Z') {
        if (second == 'n')
          return 30; // Zinc
        else if (second == 'r')
          return 40; // Zirconium
      } else if (first == 'G') {
        if (second == 'a')
          return 31; // Gallium
        else if (second == 'e')
          return 32; // Germanium
        else if (second == 'd')
          return 64; // Gadolinium
      } else if (first == 'R') {
        if (second == 'b')
          return 37; // Rubidium
        else if (second == 'u')
          return 44; // Ruthenium
        else if (second == 'h')
          return 45; // Rhodium
        else if (second == 'e')
          return 75; // Rhenium
        else if (second == 'n')
          return 86; // Radon
        else if (second == 'a')
          return 88; // Radium
        else if (second == 'f')
          return 104; // Rutherfordium
        else if (second == 'g')
          return 111; // Roentgenium
      } else if (first == 'P') {
        if (second == 'd')
          return 46; // Palladium
        else if (second == 'r')
          return 59; // Praseodymium
        else if (second == 'm')
          return 61; // Promethium
        else if (second == 'o')
          return 84; // Polonium
        else if (second == 'a')
          return 91; // Protactinium
        else if (second == 'u')
          return 94; // Plutonium
        else if (second == 't')
          return 78; // Platinum
        else if (second == 'b')
          return 82; // Lead
      } else if (first == 'I') {
        if (second == 'n')
          return 49; // Indium
        else if (second == 'r')
          return 77; // Iridium
      } else if (first == 'X') {
        if (second == 'e')
          return 54; // Xenon
      } else if (first == 'E') {
        if (second == 'u')
          return 63; // Europium
        else if (second == 'r')
          return 68; // Erbium
        else if (second == 's')
          return 99; // Einsteinium
      } else if (first == 'D') {
        if (second == 'y')
          return 66; // Dysprosium
        else if (second == 'b')
          return 105; // Dubnium
        else if (second == 's')
          return 110; // Darmstadtium
      } else if (first == 'Y') {
        if (second == 'b')
          return 70; // Ytterbium
      } else if (first == 'O') {
        if (second == 's')
          return 76; // Osmium
        else if (second == 'g')
          return 118; // Oganesson
      } else if (first == 'K') {
        if (second == 'r')
          return 36; // Krypton
      } else if (first == 'J') {
        if (second == 'b')
          return 105; // Dubnium (alternative name Joliotium)
      } else if (first == 'U') {
        if (second == 'u')
          return 113; // Nihonium (alternative name Ununtrium)
        else if (second == 'u')
          return 114; // Flerovium (alternative name Ununquadium)
        else if (second == 'p')
          return 115; // Moscovium (alternative name Ununpentium)
        else if (second == 'h')
          return 116; // Livermorium (alternative name Ununhexium)
      }
    }

    // Handle three-letter symbols for newer elements or alternative names
    else if (symbol.length() == 3) {
      if (symbol == "Uut")
        return 113; // Nihonium (older name)
      if (symbol == "Uuq")
        return 114; // Flerovium (older name)
      if (symbol == "Uup")
        return 115; // Moscovium (older name)
      if (symbol == "Uuh")
        return 116; // Livermorium (older name)
      if (symbol == "Uus")
        return 117; // Tennessine (older name)
      if (symbol == "Uuo")
        return 118; // Oganesson (older name)
    }

    // std::ostringstream error_message;
    // error_message << "Invalid atomic symbol " << symbol;
    // throw std::logic_error(error_message.str());
    return 0;
  } // }}}

  /*!
   * \brief Converts atomic number to atomic symbol
   */
  inline static constexpr std::string_view
  atomic_number_to_element(int atomic_number) {
    switch (atomic_number) { // {{{
    case 1:
      return "H"; // Hydrogen
    case 2:
      return "He"; // Helium
    case 3:
      return "Li"; // Lithium
    case 4:
      return "Be"; // Beryllium
    case 5:
      return "B"; // Boron
    case 6:
      return "C"; // Carbon
    case 7:
      return "N"; // Nitrogen
    case 8:
      return "O"; // Oxygen
    case 9:
      return "F"; // Fluorine
    case 10:
      return "Ne"; // Neon
    case 11:
      return "Na"; // Sodium
    case 12:
      return "Mg"; // Magnesium
    case 13:
      return "Al"; // Aluminum
    case 14:
      return "Si"; // Silicon
    case 15:
      return "P"; // Phosphorus
    case 16:
      return "S"; // Sulfur
    case 17:
      return "Cl"; // Chlorine
    case 18:
      return "Ar"; // Argon
    case 19:
      return "K"; // Potassium
    case 20:
      return "Ca"; // Calcium
    case 21:
      return "Sc"; // Scandium
    case 22:
      return "Ti"; // Titanium
    case 23:
      return "V"; // Vanadium
    case 24:
      return "Cr"; // Chromium
    case 25:
      return "Mn"; // Manganese
    case 26:
      return "Fe"; // Iron
    case 27:
      return "Co"; // Cobalt
    case 28:
      return "Ni"; // Nickel
    case 29:
      return "Cu"; // Copper
    case 30:
      return "Zn"; // Zinc
    case 31:
      return "Ga"; // Gallium
    case 32:
      return "Ge"; // Germanium
    case 33:
      return "As"; // Arsenic
    case 34:
      return "Se"; // Selenium
    case 35:
      return "Br"; // Bromine
    case 36:
      return "Kr"; // Krypton
    case 37:
      return "Rb"; // Rubidium
    case 38:
      return "Sr"; // Strontium
    case 39:
      return "Y"; // Yttrium
    case 40:
      return "Zr"; // Zirconium
    case 41:
      return "Nb"; // Niobium
    case 42:
      return "Mo"; // Molybdenum
    case 43:
      return "Tc"; // Technetium
    case 44:
      return "Ru"; // Ruthenium
    case 45:
      return "Rh"; // Rhodium
    case 46:
      return "Pd"; // Palladium
    case 47:
      return "Ag"; // Silver
    case 48:
      return "Cd"; // Cadmium
    case 49:
      return "In"; // Indium
    case 50:
      return "Sn"; // Tin
    case 51:
      return "Sb"; // Antimony
    case 52:
      return "Te"; // Tellurium
    case 53:
      return "I"; // Iodine
    case 54:
      return "Xe"; // Xenon
    case 55:
      return "Cs"; // Cesium
    case 56:
      return "Ba"; // Barium
    case 57:
      return "La"; // Lanthanum
    case 58:
      return "Ce"; // Cerium
    case 59:
      return "Pr"; // Praseodymium
    case 60:
      return "Nd"; // Neodymium
    case 61:
      return "Pm"; // Promethium
    case 62:
      return "Sm"; // Samarium
    case 63:
      return "Eu"; // Europium
    case 64:
      return "Gd"; // Gadolinium
    case 65:
      return "Tb"; // Terbium
    case 66:
      return "Dy"; // Dysprosium
    case 67:
      return "Ho"; // Holmium
    case 68:
      return "Er"; // Erbium
    case 69:
      return "Tm"; // Thulium
    case 70:
      return "Yb"; // Ytterbium
    case 71:
      return "Lu"; // Lutetium
    case 72:
      return "Hf"; // Hafnium
    case 73:
      return "Ta"; // Tantalum
    case 74:
      return "W"; // Tungsten
    case 75:
      return "Re"; // Rhenium
    case 76:
      return "Os"; // Osmium
    case 77:
      return "Ir"; // Iridium
    case 78:
      return "Pt"; // Platinum
    case 79:
      return "Au"; // Gold
    case 80:
      return "Hg"; // Mercury
    case 81:
      return "Tl"; // Thallium
    case 82:
      return "Pb"; // Lead
    case 83:
      return "Bi"; // Bismuth
    case 84:
      return "Po"; // Polonium
    case 85:
      return "At"; // Astatine
    case 86:
      return "Rn"; // Radon
    case 87:
      return "Fr"; // Francium
    case 88:
      return "Ra"; // Radium
    case 89:
      return "Ac"; // Actinium
    case 90:
      return "Th"; // Thorium
    case 91:
      return "Pa"; // Protactinium
    case 92:
      return "U"; // Uranium
    case 93:
      return "Np"; // Neptunium
    case 94:
      return "Pu"; // Plutonium
    case 95:
      return "Am"; // Americium
    case 96:
      return "Cm"; // Curium
    case 97:
      return "Bk"; // Berkelium
    case 98:
      return "Cf"; // Californium
    case 99:
      return "Es"; // Einsteinium
    case 100:
      return "Fm"; // Fermium
    case 101:
      return "Md"; // Mendelevium
    case 102:
      return "No"; // Nobelium
    case 103:
      return "Lr"; // Lawrencium
    case 104:
      return "Rf"; // Rutherfordium
    case 105:
      return "Db"; // Dubnium
    case 106:
      return "Sg"; // Seaborgium
    case 107:
      return "Bh"; // Bohrium
    case 108:
      return "Hs"; // Hassium
    case 109:
      return "Mt"; // Meitnerium
    case 110:
      return "Ds"; // Darmstadtium
    case 111:
      return "Rg"; // Roentgenium
    case 112:
      return "Cn"; // Copernicium
    case 113:
      return "Nh"; // Nihonium
    case 114:
      return "Fl"; // Flerovium
    case 115:
      return "Mc"; // Moscovium
    case 116:
      return "Lv"; // Livermorium
    case 117:
      return "Ts"; // Tennessine
    case 118:
      return "Og"; // Oganesson
    default:
      // std::ostringstream error_message;
      // error_message << "Invalid atomic number " << atomic_number;
      // throw std::logic_error(error_message.str());
      return "";
    }
  } // }}}

#ifdef __clang__ // Temporarily disable clang to suppress annoying error
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Weverything"
#warning                                                                       \
    "Temporarily disable clang compilation to suppress annoying errors. Code will not behave properly."; // HACK: Fix this
#pragma clang diagnostic pop
  /*!
   * An array that stores isotope data.
   */
  inline static constexpr std::array<Isotope, 0> database = {}; // NOLINT
#else
  inline static constexpr std::array database = {
      GENERATE_ISOTOPE_DATA}; // NOLINT
#endif
};
} // namespace common
} // namespace amcset
