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

/*! \file element_data.cpp
 *  \brief The element table and its lookup functions.
 */

#include "element_data.h"

#include "element_database.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <utility>

namespace amcset {
namespace common {

namespace {

constexpr std::array<ElementRecord, Element::max_z> records = {
    AMCSET_ELEMENT_DATABASE};

template <std::size_t... I>
constexpr std::array<Element, sizeof...(I)>
make_elements(std::index_sequence<I...>) {
  return {Element(&records[I])...};
}

const std::array<Element, Element::max_z> elements =
    make_elements(std::make_index_sequence<Element::max_z>{});

bool iequals(std::string_view a, std::string_view b) {
  return a.size() == b.size() &&
         std::equal(a.begin(), a.end(), b.begin(), [](char x, char y) {
           return std::tolower(static_cast<unsigned char>(x)) ==
                  std::tolower(static_cast<unsigned char>(y));
         });
}

} // namespace

const std::array<Element, Element::max_z> &Element::all() { return elements; }

const Element &Element::by_z(unsigned z) {
  if (z < 1 || z > max_z) {
    throw std::out_of_range("No element data for Z = " + std::to_string(z) +
                            " (tabulated: 1.." + std::to_string(max_z) + ")");
  }
  return elements[z - 1];
}

std::optional<unsigned> Element::find_z(std::string_view symbol_or_name) {
  // Trim whitespace.
  while (!symbol_or_name.empty() &&
         std::isspace(static_cast<unsigned char>(symbol_or_name.front()))) {
    symbol_or_name.remove_prefix(1);
  }
  while (!symbol_or_name.empty() &&
         std::isspace(static_cast<unsigned char>(symbol_or_name.back()))) {
    symbol_or_name.remove_suffix(1);
  }
  for (const auto &r : records) {
    if (iequals(symbol_or_name, r.symbol) || iequals(symbol_or_name, r.name)) {
      return r.z;
    }
  }
  return std::nullopt;
}

const Element &Element::by_symbol(std::string_view symbol_or_name) {
  const auto z = find_z(symbol_or_name);
  if (!z) {
    throw std::invalid_argument("Unknown element: '" +
                                std::string(symbol_or_name) + "'");
  }
  return by_z(*z);
}

} // namespace common
} // namespace amcset
