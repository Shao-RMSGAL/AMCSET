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

/*! \file device.h
 *  \brief The GPU backends' common interface (internal to amcset-common).
 *
 * TrimSimulation::run (trim_engine.cpp) builds a KernelInput from its
 * target model, opens a Device for TrimConfig::compute and lets it run
 * blocks of 16 ions; each block comes back as a flat tally array in the
 * layout of trim_kernel.h, which unpack_block turns into TrimTallies. The
 * blocks are merged in block order with those the CPU ran, exactly as the
 * CPU's own blocks are.
 */

#pragma once

#include "compute/compute.h"
#include "trim_kernel.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace amcset {
namespace trim {
struct TrimConfig;
struct TrimTallies;
struct IonPosition;
struct ExitRecord;
struct EnergyPositionRecord;
class TargetModel;
} // namespace trim

namespace compute {

//! Everything the kernel reads, in its layout.
struct KernelInput {
  KParams params{};
  std::vector<KLayer> layers;
  std::vector<KElement> elements;
  std::vector<KPair> pairs;
  std::vector<KCurve> curves;
  std::vector<double> log_s;
  std::vector<double> edges;
  std::vector<KParticle> particles; //!< at least one (a dummy for the beam)
  std::size_t blocks = 0;           //!< blocks of 16 ions in the run
};

/*!
 * \brief The kernel's input for a run of \p config on \p model (the
 * engine's), ions config.first_ion .. + \p ions, with \p tally_elements
 * target elements and the depth bins bin0 + k bin_w.
 */
KernelInput make_kernel_input(const trim::TrimConfig &config, const trim::TargetModel &model,
                              std::size_t ions, std::size_t tally_elements, double bin0, double bin_w);

//! Recoil slots a block needs (TrimConfig::compute.recoil_slots = 0): from
//! the CSDA path of the most energetic recoil, in monolayer steps.
std::size_t recoil_slots_needed(const trim::TrimConfig &config, const trim::TargetModel &model);

/*!
 * \brief Why \p config cannot run on a GPU backend (empty: it can). The
 * features the kernel leaves out are listed in trim_kernel.h.
 */
std::string unsupported_reason(const trim::TrimConfig &config);

//! Notes for the user about a run on a device (fast_math has no effect there).
std::vector<std::string> device_notes(const trim::TrimConfig &config);

/*!
 * \brief Whether block \p block (0-based in the run) must run on the CPU:
 * it holds ions whose trajectories are kept or that COLLISON.txt details.
 */
bool block_needs_cpu(const trim::TrimConfig &config, std::size_t block);

/*!
 * \brief The block's tallies from the kernel's flat arrays \p T (tally
 * layout) and \p lists (positions, exits) into \p out (sized for the target).
 */
void unpack_block(const KernelInput &in, const double *T, const double *lists, trim::TrimTallies &out);

//! Status of a block in its flat tallies (AMK_OK, AMK_SLOTS_FULL, ...).
int block_status(const KernelInput &in, const double *T);

//! \p t in the kernel's flat layout (the inverse of unpack_block, without
//! the lists): the running total starts from the host's tallies.
std::vector<double> pack_tallies(const KernelInput &in, const trim::TrimTallies &t);

//! Row index of scalar \p s (AMK_S_*) in the flat layout.
std::size_t scalar_row(const KernelInput &in, int s);

//! Packed position records (4 doubles each) and exit records (10 each,
//! the list kept to \p max_records as TrimTallies::add does) to the lists.
void append_positions(const double *packed, std::size_t n, std::vector<trim::IonPosition> &out);
void append_exits(const double *packed, std::size_t n, std::size_t max_records, std::vector<trim::ExitRecord> &out);
//! Packed EXYZ records (7 doubles each), kept to \p max_records.
void append_exyz(const double *packed, std::size_t n, std::size_t max_records,
                 std::vector<trim::EnergyPositionRecord> &out);

//! A compute device with the kernel loaded.
class Device {
public:
  virtual ~Device() = default;
  //! "OpenCL: gfx1100 (AMD Accelerated Parallel Processing)"
  virtual std::string description() const = 0;
  /*!
   * Upload the read-only data and size the buffers: \p ring block slots
   * (device block d keeps its tallies and lists in slot d % ring until it is
   * merged) and \p max_lanes lanes.
   */
  virtual void load(const KernelInput &in, std::size_t ring, std::size_t max_lanes) = 0;
  /*!
   * One launch: \p lanes lanes take up to \p steps steps each through the
   * device blocks 0 .. n_dev - 1 (the run's blocks first_block + d); a lane
   * takes block d only if d < \p limit (merged + ring: its slot is free).
   * \p init: the run's first launch. Synchronous; returns the lanes with no
   * block left (the run is done when every lane is and every block merged).
   */
  virtual std::size_t launch(std::size_t first_block, std::size_t n_dev, std::size_t lanes, std::size_t limit,
                             std::size_t steps, bool init) = 0;
  //! How many blocks from \p merged on have finished, in a row (the head of the ring).
  virtual std::size_t scan(std::size_t merged, std::size_t n_dev) = 0;

  // --- Results kept on the device ------------------------------------------
  //! Set the running total (tally_stride doubles, the flat layout).
  virtual void set_total(const double *flat) = 0;
  //! Add device blocks merged .. merged + count - 1 to the running total, in order.
  virtual void merge(std::size_t merged, std::size_t count) = 0;
  //! Copy the running total to the host.
  virtual void read_total(double *flat) = 0;
  //! Copy rows of the slots' tallies (row r: element r of every slot, ring doubles) to \p out.
  virtual void read_rows(const std::vector<std::size_t> &rows, double *out) = 0;
  //! Ring slot \p slot takes \p flat (tally_stride doubles; a block the CPU reran).
  virtual void put_block(std::size_t slot, const double *flat) = 0;
  /*!
   * The positions and exit records of blocks merged .. merged + count - 1
   * packed in block order: block merged + k's start at offsets[2k]
   * (positions) and offsets[2k + 1] (exits); blocks not AMK_OK are skipped.
   * With EXYZ records, offsets[2 count] is where they start in the exit
   * buffer (in doubles) and offsets[2 count + 1 + k] block k's first record.
   * \p positions gets n_positions * 4 doubles, \p exits \p n_exits doubles.
   */
  virtual void pack_lists(std::size_t merged, std::size_t count, const std::vector<std::int64_t> &offsets,
                          std::size_t n_positions, std::size_t n_exits, double *positions, double *exits) = 0;
  //! Lanes that fill the device (several wavefronts per SIMD).
  virtual std::size_t preferred_lanes() const = 0;
};

//! Bytes one ring slot takes on a device (tallies, lists).
std::size_t block_bytes(const KernelInput &in);
//! Bytes one lane takes (its recoil slots and state).
std::size_t lane_bytes(const KernelInput &in);

/*!
 * \brief Open device \p options.device of \p options.backend (not cpu).
 * \throws std::runtime_error with the reason (no runtime, no such device,
 * no double precision, the kernel does not compile)
 */
std::unique_ptr<Device> open_device(const Options &options);

// The backends (each in its own file; null when not built in).
std::unique_ptr<Device> open_reference(const Options &options);
std::unique_ptr<Device> open_opencl(const Options &options);
std::unique_ptr<Device> open_hip(const Options &options, bool cuda);
std::unique_ptr<Device> open_kokkos(const Options &options);
void list_opencl(std::vector<DeviceInfo> &out, std::vector<std::string> *notes);
void list_hip(std::vector<DeviceInfo> &out, std::vector<std::string> *notes, bool cuda);
void list_kokkos(std::vector<DeviceInfo> &out, std::vector<std::string> *notes);

//! The kernel source (trim_kernel.h as text, embedded at build time).
const std::string &kernel_source();

} // namespace compute
} // namespace amcset
