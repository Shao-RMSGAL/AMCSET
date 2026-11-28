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

/*! \file kokkos.cpp
 *  \brief The Kokkos backend (configure with -DAMCSET_WITH_KOKKOS=ON and a
 *  Kokkos installation on CMAKE_PREFIX_PATH).
 *
 * Unlike the other GPU backends, Kokkos is chosen when AMCSET is built: the
 * kernel source is compiled for Kokkos's default execution space (Serial,
 * Threads or OpenMP on the CPU; CUDA, HIP or SYCL on a GPU, with Kokkos's
 * compiler wrapper), and runs there as a parallel_for over the lanes. Without
 * AMCSET_WITH_KOKKOS this file only reports that the backend is not built in.
 */

#include "device.h"

#include <stdexcept>

#if defined(AMCSET_WITH_KOKKOS)
#include "kokkos_device.h"

#include <Kokkos_Core.hpp>

#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#endif

namespace amcset {
namespace compute {

#if defined(AMCSET_WITH_KOKKOS)

namespace {

/*!
 * Kokkos runs on one thread of its own: some execution spaces (Threads)
 * take work only from the thread that initialised them, and the engine
 * calls the device from its own worker. Kokkos is initialised there on
 * first use (all hardware threads) and finalised at exit.
 */
class Runtime {
public:
  Runtime() : thread_([this] { loop(); }) {
    call([] {
      if (!Kokkos::is_initialized()) {
        Kokkos::InitializationSettings s;
        s.set_disable_warnings(true);
        s.set_num_threads(static_cast<int>(std::max(1u, std::thread::hardware_concurrency())));
        Kokkos::initialize(s);
      }
    });
  }
  ~Runtime() {
    call([] {
      if (Kokkos::is_initialized()) Kokkos::finalize();
    });
    {
      std::lock_guard<std::mutex> lock(mutex_);
      stop_ = true;
    }
    cv_.notify_all();
    thread_.join();
  }
  //! Run \p f on the Kokkos thread and wait for it (its exception is rethrown here).
  void call(const std::function<void()> &f) {
    std::packaged_task<void()> task(f);
    auto done = task.get_future();
    {
      std::lock_guard<std::mutex> lock(mutex_);
      tasks_.push_back(std::move(task));
    }
    cv_.notify_all();
    done.get();
  }

private:
  void loop() {
    while (true) {
      std::packaged_task<void()> task;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return stop_ || !tasks_.empty(); });
        if (tasks_.empty()) return;
        task = std::move(tasks_.front());
        tasks_.pop_front();
      }
      task();
    }
  }
  std::mutex mutex_;
  std::condition_variable cv_;
  std::deque<std::packaged_task<void()>> tasks_;
  bool stop_ = false;
  std::thread thread_;
};
Runtime &runtime() {
  static Runtime r;
  return r;
}

using Space = Kokkos::DefaultExecutionSpace;
using Mem = Space::memory_space;

std::string space_name() { return std::string(Space::name()); }

int concurrency() {
  int n = 1;
  runtime().call([&] { n = Space().concurrency(); });
  return n;
}

class KokkosDevice final : public Device {
public:
  explicit KokkosDevice(bool mixed) : mixed_(mixed), concurrency_(concurrency()) {}
  ~KokkosDevice() override {
    // The views are freed on the Kokkos thread.
    runtime().call([this] {
      params_ = {};
      layers_ = {};
      elements_ = {};
      pairs_ = {};
      curves_ = {};
      log_s_ = {};
      edges_ = {};
      particles_ = {};
      slots_ = {};
      tallies_ = {};
      lists_ = {};
      lanes_ = {};
      counter_ = {};
      total_ = {};
      slot_block_ = {};
      scan_out_ = {};
      offsets_ = {};
      flat_ = {};
    });
  }

  std::string description() const override {
    return "Kokkos (" + space_name() + ", concurrency " + std::to_string(concurrency_) + ")" +
           (mixed_ ? " [mixed precision]" : "");
  }
  std::size_t preferred_lanes() const override { return static_cast<std::size_t>(std::max(1, concurrency_)); }

  void load(const KernelInput &in, std::size_t ring, std::size_t max_lanes) override {
    runtime().call([&] { load_here(in, ring, max_lanes); });
  }
  std::size_t launch(std::size_t first_block, std::size_t n_dev, std::size_t lanes, std::size_t limit,
                     std::size_t steps, bool init) override {
    std::size_t done = 0;
    runtime().call([&] {
      if (init) Kokkos::deep_copy(counter_, 0);
      (mixed_ ? kokkos_lanes_mixed : kokkos_lanes_double)(b_, lanes, first_block, n_dev, limit, steps, init);
      auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), counter_);
      done = static_cast<std::size_t>(host(1));
    });
    return done;
  }
  std::size_t scan(std::size_t merged, std::size_t n_dev) override {
    std::size_t c = 0;
    runtime().call([&] {
      kokkos_scan(b_, scan_out_.data(), merged, n_dev);
      auto host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), scan_out_);
      c = static_cast<std::size_t>(host(0));
    });
    return c;
  }
  void set_total(const double *flat) override {
    runtime().call([&] {
      Kokkos::View<const double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(flat, total_.extent(0));
      Kokkos::deep_copy(total_, h);
    });
  }
  void merge(std::size_t merged, std::size_t count) override {
    runtime().call([&] { kokkos_merge(b_, total_.data(), merged, count); });
  }
  void read_total(double *flat) override {
    runtime().call([&] {
      Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(flat, total_.extent(0));
      Kokkos::deep_copy(h, total_);
    });
  }
  void read_rows(const std::vector<std::size_t> &rows, double *out) override {
    runtime().call([&] {
      for (std::size_t i = 0; i < rows.size(); ++i) {
        Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(out + i * ring_, ring_);
        Kokkos::deep_copy(h, Kokkos::subview(tallies_, std::make_pair(rows[i] * ring_, (rows[i] + 1) * ring_)));
      }
    });
  }
  void put_block(std::size_t slot, const double *flat) override {
    runtime().call([&] {
      Kokkos::View<const double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(flat, flat_.extent(0));
      Kokkos::deep_copy(flat_, h);
      kokkos_put(b_, flat_.data(), slot);
    });
  }
  void pack_lists(std::size_t merged, std::size_t count, const std::vector<std::int64_t> &offsets,
                  std::size_t n_positions, std::size_t n_exits, double *positions, double *exits) override {
    runtime().call([&] {
      using HostD = Kokkos::View<double *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
      Kokkos::View<const amk_i64 *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> ho(
          reinterpret_cast<const amk_i64 *>(offsets.data()), offsets.size());
      Kokkos::deep_copy(Kokkos::subview(offsets_, std::make_pair(std::size_t{0}, offsets.size())), ho);
      Kokkos::View<double *, Mem> pp("positions", std::max<std::size_t>(1, n_positions) * AMK_POSITION_DOUBLES);
      Kokkos::View<double *, Mem> pe("exits", std::max<std::size_t>(1, n_exits));
      kokkos_pack(b_, offsets_.data(), pp.data(), pe.data(), merged, count);
      if (n_positions > 0) {
        HostD hp(positions, n_positions * AMK_POSITION_DOUBLES);
        Kokkos::deep_copy(hp, Kokkos::subview(pp, std::make_pair(std::size_t{0}, n_positions * AMK_POSITION_DOUBLES)));
      }
      if (n_exits > 0) {
        HostD he(exits, n_exits);
        Kokkos::deep_copy(he, Kokkos::subview(pe, std::make_pair(std::size_t{0}, n_exits)));
      }
    });
  }

private:
  void load_here(const KernelInput &in, std::size_t ring, std::size_t max_lanes) {
    in_ = &in;
    ring_ = ring;
    const KParams &P = in.params;
    params_ = upload("params", &in.params, 1);
    layers_ = upload("layers", in.layers.data(), in.layers.size());
    elements_ = upload("elements", in.elements.data(), in.elements.size());
    pairs_ = upload("pairs", in.pairs.data(), in.pairs.size());
    curves_ = upload("curves", in.curves.data(), in.curves.size());
    log_s_ = upload("log_s", in.log_s.data(), in.log_s.size());
    edges_ = upload("edges", in.edges.data(), in.edges.size());
    particles_ = upload("particles", in.particles.data(), in.particles.size());
    slots_ = Kokkos::View<KMover *, Mem>("slots", max_lanes * static_cast<std::size_t>(P.slot_capacity));
    tallies_ = Kokkos::View<double *, Mem>("tallies", ring * static_cast<std::size_t>(P.tally_stride));
    lists_ = Kokkos::View<double *, Mem>("lists", std::max<std::size_t>(1, ring * static_cast<std::size_t>(P.list_stride)));
    lanes_ = Kokkos::View<KLane *, Mem>("lanes", max_lanes);
    counter_ = Kokkos::View<int *, Mem>("counter", 2);
    slot_block_ = Kokkos::View<amk_i64 *, Mem>("slot_block", ring); // zero: no slot holds a finished block yet
    scan_out_ = Kokkos::View<amk_i64 *, Mem>("scan", 1);
    offsets_ = Kokkos::View<amk_i64 *, Mem>("offsets", 3 * ring + 2);
    total_ = Kokkos::View<double *, Mem>("total", static_cast<std::size_t>(P.tally_stride));
    flat_ = Kokkos::View<double *, Mem>("flat", static_cast<std::size_t>(P.tally_stride));
    b_.params = params_.data();
    b_.layers = layers_.data();
    b_.elements = elements_.data();
    b_.pairs = pairs_.data();
    b_.curves = curves_.data();
    b_.log_s = log_s_.data();
    b_.edges = edges_.data();
    b_.particles = particles_.data();
    b_.slots = slots_.data();
    b_.tallies = tallies_.data();
    b_.lists = lists_.data();
    b_.lanes = lanes_.data();
    b_.counter = counter_.data();
    b_.slot_block = slot_block_.data();
    b_.ring = static_cast<amk_i64>(ring);
    b_.tally_stride = P.tally_stride;
  }

  template <class T> Kokkos::View<T *, Mem> upload(const char *name, const T *data, std::size_t n) {
    Kokkos::View<T *, Mem> v(std::string(name), static_cast<std::size_t>(std::max<std::size_t>(n, 1)));
    if (n > 0) {
      Kokkos::View<const T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> h(data, n);
      Kokkos::deep_copy(Kokkos::subview(v, std::make_pair(std::size_t{0}, n)), h);
    }
    return v;
  }

  bool mixed_;
  int concurrency_ = 1;
  const KernelInput *in_ = nullptr;
  std::size_t ring_ = 1;
  KokkosBuffers b_;
  Kokkos::View<KParams *, Mem> params_;
  Kokkos::View<KLayer *, Mem> layers_;
  Kokkos::View<KElement *, Mem> elements_;
  Kokkos::View<KPair *, Mem> pairs_;
  Kokkos::View<KCurve *, Mem> curves_;
  Kokkos::View<double *, Mem> log_s_;
  Kokkos::View<double *, Mem> edges_;
  Kokkos::View<KParticle *, Mem> particles_;
  Kokkos::View<KMover *, Mem> slots_;
  Kokkos::View<double *, Mem> tallies_;
  Kokkos::View<double *, Mem> lists_;
  Kokkos::View<KLane *, Mem> lanes_;
  Kokkos::View<int *, Mem> counter_;
  Kokkos::View<amk_i64 *, Mem> slot_block_;
  Kokkos::View<amk_i64 *, Mem> scan_out_;
  Kokkos::View<amk_i64 *, Mem> offsets_;
  Kokkos::View<double *, Mem> total_;
  Kokkos::View<double *, Mem> flat_;
};

} // namespace

void list_kokkos(std::vector<DeviceInfo> &out, std::vector<std::string> *) {
  const int n = concurrency();
  DeviceInfo info;
  info.backend = Backend::kokkos;
  info.index = 0;
  info.name = "Kokkos " + space_name();
  info.details = "Kokkos " + std::to_string(KOKKOS_VERSION / 10000) + "." + std::to_string(KOKKOS_VERSION / 100 % 100) +
                 "." + std::to_string(KOKKOS_VERSION % 100) + ", default execution space " + space_name() +
                 ", concurrency " + std::to_string(n);
  out.push_back(info);
}

std::unique_ptr<Device> open_kokkos(const Options &o) {
  if (o.device != 0) throw std::runtime_error("Kokkos: one device (its default execution space), number 0");
  return std::make_unique<KokkosDevice>(o.precision == Precision::mixed);
}

#else

void list_kokkos(std::vector<DeviceInfo> &, std::vector<std::string> *notes) {
  if (notes) notes->push_back("Kokkos: not built in (configure with -DAMCSET_WITH_KOKKOS=ON)");
}

std::unique_ptr<Device> open_kokkos(const Options &) {
  throw std::runtime_error("Kokkos: AMCSET was built without Kokkos (configure with -DAMCSET_WITH_KOKKOS=ON)");
}

#endif

} // namespace compute
} // namespace amcset
