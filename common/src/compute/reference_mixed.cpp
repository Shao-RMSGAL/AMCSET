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

/*! \file reference_mixed.cpp
 *  \brief The kernel source compiled for the CPU in mixed precision
 *  (compute.precision = mixed on the reference backend): the same
 *  arithmetic as a GPU's mixed-precision kernel, up to its math library.
 */

#define AMK_MIXED
#include "trim_kernel.h"

#include "device.h"

namespace amcset {
namespace compute {

void reference_lane_mixed(const KParams *P, const KLayer *layers, const KElement *elements, const KPair *pairs,
                          const KCurve *curves, const double *log_s, const double *edges, const KParticle *particles,
                          KMover *slots, double *tallies, double *lists, KLane *lanes, int *counter,
                          amk_i64 *slot_block, amk_i64 lane, amk_i64 first_block, amk_i64 n_dev, amk_i64 ring,
                          amk_i64 limit, amk_i64 steps, amk_i64 init) {
  amk_run_lane(P, layers, elements, pairs, curves, log_s, edges, particles, slots, tallies, lists, lanes, counter,
               slot_block, lane, first_block, n_dev, ring, limit, steps, init);
}

} // namespace compute
} // namespace amcset
