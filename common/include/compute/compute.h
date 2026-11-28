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

/*! \file compute.h
 *  \brief Where SRIM mode's histories run: the CPU or a GPU framework.
 *
 * The CPU (default) runs the trim2013 transport of trim_engine.cpp on
 * std::threads. The GPU backends run the same transport from one portable
 * kernel source (common/src/compute/trim_kernel.h) through OpenCL, HIP
 * (AMD), CUDA (NVIDIA) or Kokkos; "reference" runs that kernel source on the
 * CPU (it reproduces the CPU transport bit for bit and exists to check the
 * kernel). Every framework is loaded at run time from the system's driver
 * (OpenCL.dll / libOpenCL.so, amdhip64 / libamdhip64, nvcuda / libcuda), so
 * AMCSET runs on machines without them; Kokkos is chosen when AMCSET is
 * built (AMCSET_WITH_KOKKOS). docs/GPU.md.
 */

#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace compute {

//! The compute frameworks.
enum class Backend {
  cpu,       //!< The CPU transport on std::threads (default)
  opencl,    //!< OpenCL 1.2+ with double precision (AMD, NVIDIA, Intel, PoCL)
  hip,       //!< AMD HIP (ROCm / the AMD driver's HIP runtime)
  cuda,      //!< NVIDIA CUDA (driver API and NVRTC)
  kokkos,    //!< Kokkos, with the execution space AMCSET was built with
  reference  //!< The kernel source compiled for the CPU (checks the kernel)
};

//! "cpu", "opencl", "hip", "cuda", "kokkos", "reference".
std::string to_string(Backend b);
//! Inverse of to_string (case-insensitive; "gpu" = the first usable GPU
//! framework is not accepted here: see best_device()). \throws std::invalid_argument
Backend parse_backend(std::string_view s);
//! Every backend, in the order above.
std::vector<Backend> all_backends();

//! Arithmetic of the GPU kernel (compute.precision).
enum class Precision {
  double_precision, //!< "double": the CPU transport's arithmetic (default)
  /*!
   * "mixed": the magic formula and the stopping interpolation in single
   * precision, everything else in double. Consumer GPUs run single
   * precision 16 to 64 times faster; statistically the same results
   * (docs/GPU.md "Precision"). No effect on the CPU backend.
   */
  mixed
};
std::string to_string(Precision p);
Precision parse_precision(std::string_view s); //!< "double", "mixed"

//! Settings of a calculation (scenario key "compute"; SRIM mode).
struct Options {
  Backend backend = Backend::cpu;
  Precision precision = Precision::double_precision;
  //! Device number within the backend (amcset devices lists them).
  int device = 0;
  /*!
   * Upper limit of the device memory the run's buffers take (MiB): up to
   * 40 % for the lanes (their recoil stacks), the rest for the ring of block
   * slots. More memory, more lanes: the device is used better. 0: 2048 MiB.
   */
  std::size_t memory_mib = 0;
  /*!
   * Recoil slots per lane (0: from the path length of the most energetic
   * recoil, compute::recoil_slots_needed). TRIM's recoil stack gains a slot
   * at every collision of the atom on top; a block that needs more slots is
   * rerun on the CPU.
   */
  std::size_t recoil_slots = 0;
  //! Exit records (BACKSCAT/TRANSMIT/SPUTTER lists) per block of 16 ions
  //! (0: 256; none without output.exit_lists); a block with more is rerun
  //! on the CPU.
  std::size_t exit_records = 0;

  bool operator==(const Options &) const = default;
};

//! One device of one backend.
struct DeviceInfo {
  Backend backend = Backend::cpu;
  int index = 0;           //!< Options::device
  std::string name;        //!< "AMD Radeon RX 7900 XTX"
  std::string details;     //!< Driver, compute units, memory, architecture
  bool usable = true;      //!< Has what the kernel needs (double precision)
  std::string problem;     //!< Why not, if !usable
};

/*!
 * \brief The devices of every backend this machine has (probes the drivers;
 * a few hundred milliseconds). Backends without a runtime, or not built in,
 * are listed in \p notes with the reason.
 */
std::vector<DeviceInfo> list_devices(std::vector<std::string> *notes = nullptr);

//! Whether AMCSET was built with \p b (Kokkos only if AMCSET_WITH_KOKKOS).
bool built_with(Backend b);

} // namespace compute
} // namespace amcset
