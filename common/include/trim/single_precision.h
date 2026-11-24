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

/*! \file single_precision.h
 *  \brief TRIM.exe's single-precision depth tables.
 *
 * TRIM.exe adds 1 / (bin width) to a VB Single for every vacancy, knock-on,
 * replacement and stopped atom of its depth tables. Once a bin's sum is large,
 * each addition is rounded to the sum's precision, which shifts the tables by
 * up to about 1.5 % from the exact counts (e.g. +0.6 % for VACANCY.txt of
 * 2000 Fe ions at 2 MeV into Fe, while the totals in its header, kept in
 * double precision, are exact). These functions reproduce that rounding.
 */
#pragma once

#include <cstdint>

namespace amcset {
namespace trim {

/*!
 * \brief Single-precision sum of \p increments additions and \p decrements
 * subtractions of \p step, starting from zero, with the subtractions spread
 * evenly among the additions; each operation is rounded to float as in
 * TRIM.exe. Takes time proportional to the number of binary orders of
 * magnitude the sum passes through, not to the number of operations.
 */
double single_precision_sum(std::uint64_t increments, std::uint64_t decrements, float step);

/*!
 * \brief The value of a TRIM count table bin in counts: the count \p count
 * (of which \p removed were taken away again) accumulated as TRIM does, in
 * steps of 1 / \p bin_width_A in single precision, and converted back to a
 * count (so that the usual normalisations apply).
 */
double trim_table_count(double count, double removed, double bin_width_A);

} // namespace trim
} // namespace amcset
