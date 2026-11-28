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

// The Kokkos lanes, included by kokkos_double.cpp and kokkos_mixed.cpp with
// AMK_KOKKOS (and AMK_MIXED) defined and AMK_KOKKOS_LANES naming the function.
// The kernel source runs in Kokkos::DefaultExecutionSpace: a parallel_for
// over the lanes, each calling amk_run_lane as a GPU work item would.

#include "trim_kernel.h"

#include "kokkos_device.h"

namespace amcset {
namespace compute {

void AMK_KOKKOS_LANES(const KokkosBuffers &b, std::size_t lanes, std::size_t first_block, std::size_t n_dev,
                      std::size_t limit, std::size_t steps, bool init) {
  const KParams *P = b.params;
  const KLayer *layers = b.layers;
  const KElement *elements = b.elements;
  const KPair *pairs = b.pairs;
  const KCurve *curves = b.curves;
  const double *log_s = b.log_s;
  const double *edges = b.edges;
  const KParticle *particles = b.particles;
  KMover *slots = b.slots;
  double *tallies = b.tallies;
  double *lists = b.lists;
  KLane *lane_state = b.lanes;
  int *counter = b.counter;
  amk_i64 *slot_block = b.slot_block;
  const amk_i64 ring = b.ring;
  const auto fb = static_cast<amk_i64>(first_block);
  const auto nd = static_cast<amk_i64>(n_dev);
  const auto lim = static_cast<amk_i64>(limit);
  const auto st = static_cast<amk_i64>(steps);
  const amk_i64 in = init ? 1 : 0;
  Kokkos::parallel_for(
      "amcset_trim_lanes", Kokkos::RangePolicy<>(0, static_cast<amk_i64>(lanes)), KOKKOS_LAMBDA(const amk_i64 i) {
        amk_run_lane(P, layers, elements, pairs, curves, log_s, edges, particles, slots, tallies, lists, lane_state,
                     counter, slot_block, i, fb, nd, ring, lim, st, in);
      });
  Kokkos::fence();
}

#if defined(AMK_KOKKOS_AUX)
// The ring and the running total (the same in both precisions).
std::size_t kokkos_scan(const KokkosBuffers &b, amk_i64 *out, std::size_t merged, std::size_t n_dev) {
  const amk_i64 *slot_block = b.slot_block;
  const amk_i64 ring = b.ring;
  const auto m = static_cast<amk_i64>(merged);
  const auto nd = static_cast<amk_i64>(n_dev);
  Kokkos::parallel_for("amcset_trim_scan", Kokkos::RangePolicy<>(0, 1), KOKKOS_LAMBDA(const amk_i64) {
    amk_scan(slot_block, ring, m, nd, out);
  });
  Kokkos::fence();
  return 0;
}
void kokkos_merge(const KokkosBuffers &b, double *total, std::size_t merged, std::size_t count) {
  const KParams *P = b.params;
  const double *tallies = b.tallies;
  const amk_i64 ring = b.ring;
  const auto m = static_cast<amk_i64>(merged);
  const auto n = static_cast<amk_i64>(count);
  Kokkos::parallel_for("amcset_trim_merge", Kokkos::RangePolicy<>(0, b.tally_stride), KOKKOS_LAMBDA(const amk_i64 e) {
    amk_merge_element(P, tallies, total, ring, m, n, e);
  });
  Kokkos::fence();
}
void kokkos_put(const KokkosBuffers &b, const double *flat, std::size_t slot) {
  const KParams *P = b.params;
  double *tallies = b.tallies;
  const amk_i64 ring = b.ring;
  const auto s = static_cast<amk_i64>(slot);
  Kokkos::parallel_for("amcset_trim_put", Kokkos::RangePolicy<>(0, b.tally_stride), KOKKOS_LAMBDA(const amk_i64 e) {
    amk_put_element(P, tallies, flat, ring, s, e);
  });
  Kokkos::fence();
}
void kokkos_pack(const KokkosBuffers &b, const amk_i64 *offsets, double *positions, double *exits, std::size_t merged,
                 std::size_t count) {
  const KParams *P = b.params;
  const double *tallies = b.tallies;
  const double *lists = b.lists;
  const amk_i64 ring = b.ring;
  const auto m = static_cast<amk_i64>(merged);
  const auto n = static_cast<amk_i64>(count);
  Kokkos::parallel_for("amcset_trim_pack", Kokkos::RangePolicy<>(0, n), KOKKOS_LAMBDA(const amk_i64 k) {
    amk_pack_lists(P, tallies, lists, offsets, positions, exits, ring, m, n, k);
  });
  Kokkos::fence();
}
#endif

} // namespace compute
} // namespace amcset
