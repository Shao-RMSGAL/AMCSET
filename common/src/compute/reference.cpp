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

/*! \file reference.cpp
 *  \brief The "reference" backend: the kernel source compiled as C++ and
 *  run on the CPU's threads (one work item per block, as on a GPU).
 *
 * It reproduces the CPU transport bit for bit (same C library, no
 * contraction), which checks the kernel's port of the transport; a GPU
 * differs from it only by its math library.
 */

#include "device.h"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <thread>

namespace amcset {
namespace compute {

// reference_mixed.cpp: the same lane in mixed precision.
void reference_lane_mixed(const KParams *P, const KLayer *layers, const KElement *elements, const KPair *pairs,
                          const KCurve *curves, const double *log_s, const double *edges, const KParticle *particles,
                          KMover *slots, double *tallies, double *lists, KLane *lanes, int *counter,
                          amk_i64 *slot_block, amk_i64 lane, amk_i64 first_block, amk_i64 n_dev, amk_i64 ring,
                          amk_i64 limit, amk_i64 steps, amk_i64 init);

namespace {

void reference_lane_double(const KParams *P, const KLayer *layers, const KElement *elements, const KPair *pairs,
                           const KCurve *curves, const double *log_s, const double *edges, const KParticle *particles,
                           KMover *slots, double *tallies, double *lists, KLane *lanes, int *counter,
                           amk_i64 *slot_block, amk_i64 lane, amk_i64 first_block, amk_i64 n_dev, amk_i64 ring,
                           amk_i64 limit, amk_i64 steps, amk_i64 init) {
  amk_run_lane(P, layers, elements, pairs, curves, log_s, edges, particles, slots, tallies, lists, lanes, counter,
               slot_block, lane, first_block, n_dev, ring, limit, steps, init);
}

class ReferenceDevice final : public Device {
public:
  ReferenceDevice(unsigned threads, bool mixed)
      : threads_(std::max(1u, threads)), mixed_(mixed), lane_fn_(mixed ? reference_lane_mixed : reference_lane_double) {}

  std::string description() const override {
    return "reference: the kernel source on " + std::to_string(threads_) + " CPU threads" +
           (mixed_ ? " [mixed precision]" : "");
  }
  std::size_t preferred_lanes() const override { return threads_; }

  void load(const KernelInput &in, std::size_t ring, std::size_t max_lanes) override {
    in_ = &in;
    ring_ = ring;
    const KParams &P = in.params;
    slots_.assign(max_lanes * static_cast<std::size_t>(P.slot_capacity), KMover{});
    lanes_.assign(max_lanes, KLane{});
    tallies_.assign(ring * static_cast<std::size_t>(P.tally_stride), 0.0);
    lists_.assign(std::max<std::size_t>(1, ring * static_cast<std::size_t>(P.list_stride)), 0.0);
    slot_block_.assign(ring, 0);
    total_.assign(static_cast<std::size_t>(P.tally_stride), 0.0);
  }

  std::size_t launch(std::size_t first_block, std::size_t n_dev, std::size_t lanes, std::size_t limit,
                     std::size_t steps, bool init) override {
    // One thread per lane, as many steps as asked.
    const KernelInput &in = *in_;
    if (init) counter_[0] = counter_[1] = 0;
    const auto lane = [&](std::size_t i) {
      lane_fn_(&in.params, in.layers.data(), in.elements.data(), in.pairs.data(), in.curves.data(), in.log_s.data(),
               in.edges.data(), in.particles.data(), slots_.data(), tallies_.data(), lists_.data(), lanes_.data(),
               counter_, slot_block_.data(), static_cast<amk_i64>(i), static_cast<amk_i64>(first_block),
               static_cast<amk_i64>(n_dev), static_cast<amk_i64>(ring_), static_cast<amk_i64>(limit),
               static_cast<amk_i64>(steps), init ? 1 : 0);
    };
    std::vector<std::thread> pool;
    for (std::size_t t = 1; t < lanes; ++t) pool.emplace_back(lane, t);
    lane(0);
    for (auto &t : pool) t.join();
    return static_cast<std::size_t>(std::atomic_ref<int>(counter_[1]).load());
  }

  std::size_t scan(std::size_t merged, std::size_t n_dev) override {
    amk_i64 c = 0;
    amk_scan(slot_block_.data(), static_cast<amk_i64>(ring_), static_cast<amk_i64>(merged), static_cast<amk_i64>(n_dev), &c);
    return static_cast<std::size_t>(c);
  }
  void set_total(const double *flat) override { total_.assign(flat, flat + in_->params.tally_stride); }
  void merge(std::size_t merged, std::size_t count) override {
    for (amk_i64 e = 0; e < in_->params.tally_stride; ++e) {
      amk_merge_element(&in_->params, tallies_.data(), total_.data(), static_cast<amk_i64>(ring_),
                        static_cast<amk_i64>(merged), static_cast<amk_i64>(count), e);
    }
  }
  void read_total(double *flat) override { std::memcpy(flat, total_.data(), total_.size() * sizeof(double)); }
  void read_rows(const std::vector<std::size_t> &rows, double *out) override {
    for (std::size_t i = 0; i < rows.size(); ++i) {
      std::memcpy(out + i * ring_, tallies_.data() + rows[i] * ring_, ring_ * sizeof(double));
    }
  }
  void put_block(std::size_t slot, const double *flat) override {
    for (amk_i64 e = 0; e < in_->params.tally_stride; ++e) {
      amk_put_element(&in_->params, tallies_.data(), flat, static_cast<amk_i64>(ring_), static_cast<amk_i64>(slot), e);
    }
  }
  void pack_lists(std::size_t merged, std::size_t count, const std::vector<std::int64_t> &offsets,
                  std::size_t n_positions, std::size_t n_exits, double *positions, double *exits) override {
    std::vector<double> pp(std::max<std::size_t>(1, n_positions) * AMK_POSITION_DOUBLES);
    std::vector<double> pe(std::max<std::size_t>(1, n_exits));
    for (std::size_t k = 0; k < count; ++k) {
      amk_pack_lists(&in_->params, tallies_.data(), lists_.data(), reinterpret_cast<const amk_i64 *>(offsets.data()),
                     pp.data(), pe.data(), static_cast<amk_i64>(ring_), static_cast<amk_i64>(merged),
                     static_cast<amk_i64>(count), static_cast<amk_i64>(k));
    }
    std::memcpy(positions, pp.data(), n_positions * AMK_POSITION_DOUBLES * sizeof(double));
    std::memcpy(exits, pe.data(), n_exits * sizeof(double));
  }

private:
  unsigned threads_;
  bool mixed_;
  decltype(&reference_lane_double) lane_fn_;
  const KernelInput *in_ = nullptr;
  std::size_t ring_ = 1;
  std::vector<KMover> slots_;
  std::vector<KLane> lanes_;
  std::vector<double> tallies_;
  std::vector<double> lists_;
  std::vector<amk_i64> slot_block_;
  std::vector<double> total_;
  int counter_[2] = {0, 0};
};

} // namespace

std::unique_ptr<Device> open_reference(const Options &o) {
  return std::make_unique<ReferenceDevice>(std::thread::hardware_concurrency(), o.precision == Precision::mixed);
}

} // namespace compute
} // namespace amcset
