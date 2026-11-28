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

/*! \file kokkos_device.h
 *  \brief The Kokkos backend's buffers (pointers into the execution space's
 *  memory) and its functions (kokkos_double.cpp, kokkos_mixed.cpp). Built
 *  only with AMCSET_WITH_KOKKOS.
 */

#pragma once

#include "trim_kernel.h"

#include <cstddef>

namespace amcset {
namespace compute {

struct KokkosBuffers {
  const KParams *params = nullptr;
  const KLayer *layers = nullptr;
  const KElement *elements = nullptr;
  const KPair *pairs = nullptr;
  const KCurve *curves = nullptr;
  const double *log_s = nullptr;
  const double *edges = nullptr;
  const KParticle *particles = nullptr;
  KMover *slots = nullptr;
  double *tallies = nullptr;
  double *lists = nullptr;
  KLane *lanes = nullptr;
  int *counter = nullptr;
  amk_i64 *slot_block = nullptr;
  amk_i64 ring = 1;
  amk_i64 tally_stride = 0;
};

void kokkos_lanes_double(const KokkosBuffers &b, std::size_t lanes, std::size_t first_block, std::size_t n_dev,
                         std::size_t limit, std::size_t steps, bool init);
void kokkos_lanes_mixed(const KokkosBuffers &b, std::size_t lanes, std::size_t first_block, std::size_t n_dev,
                        std::size_t limit, std::size_t steps, bool init);
std::size_t kokkos_scan(const KokkosBuffers &b, amk_i64 *out, std::size_t merged, std::size_t n_dev);
void kokkos_merge(const KokkosBuffers &b, double *total, std::size_t merged, std::size_t count);
void kokkos_put(const KokkosBuffers &b, const double *flat, std::size_t slot);
void kokkos_pack(const KokkosBuffers &b, const amk_i64 *offsets, double *positions, double *exits, std::size_t merged,
                 std::size_t count);

} // namespace compute
} // namespace amcset
