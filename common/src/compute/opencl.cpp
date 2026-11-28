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

/*! \file opencl.cpp
 *  \brief The OpenCL backend: the kernel source built at run time by the
 *  platform's compiler (OpenCL 1.2 with cl_khr_fp64), loaded from the ICD
 *  loader of the system (OpenCL.dll, libOpenCL.so.1, the macOS framework).
 *
 * The OpenCL C API is declared here (the few calls and constants it needs,
 * from the Khronos specification): AMCSET needs no OpenCL SDK to build.
 */

#include "device.h"
#include "dynlib.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

namespace amcset {
namespace compute {

namespace {

#if defined(_WIN32)
#define AMK_CL_CALL __stdcall
#else
#define AMK_CL_CALL
#endif

using cl_int = std::int32_t;
using cl_uint = std::uint32_t;
using cl_ulong = std::uint64_t;
using cl_bitfield = cl_ulong;
using cl_platform_id = struct cl_platform_s *;
using cl_device_id = struct cl_device_s *;
using cl_context = struct cl_context_s *;
using cl_command_queue = struct cl_queue_s *;
using cl_mem = struct cl_mem_s *;
using cl_program = struct cl_program_s *;
using cl_kernel = struct cl_kernel_s *;
using cl_context_properties = std::intptr_t;

constexpr cl_bitfield CL_DEVICE_TYPE_GPU = 1 << 2;
constexpr cl_bitfield CL_DEVICE_TYPE_ACCELERATOR = 1 << 3;
constexpr cl_bitfield CL_DEVICE_TYPE_CPU = 1 << 1;
constexpr cl_uint CL_PLATFORM_NAME = 0x0902;
constexpr cl_uint CL_PLATFORM_VERSION = 0x0901;
constexpr cl_uint CL_DEVICE_TYPE = 0x1000;
constexpr cl_uint CL_DEVICE_MAX_COMPUTE_UNITS = 0x1002;
constexpr cl_uint CL_DEVICE_GLOBAL_MEM_SIZE = 0x101F;
constexpr cl_uint CL_DEVICE_MAX_MEM_ALLOC_SIZE = 0x1010;
constexpr cl_uint CL_DEVICE_NAME = 0x102B;
constexpr cl_uint CL_DEVICE_VERSION = 0x102F;
constexpr cl_uint CL_DEVICE_EXTENSIONS = 0x1030;
constexpr cl_uint CL_DRIVER_VERSION = 0x102D;
constexpr cl_uint CL_DEVICE_BOARD_NAME_AMD = 0x4038;
constexpr cl_uint CL_PROGRAM_BUILD_LOG = 0x1183;
constexpr cl_bitfield CL_MEM_READ_WRITE = 1 << 0;
constexpr cl_bitfield CL_MEM_READ_ONLY = 1 << 2;
constexpr cl_bitfield CL_MEM_COPY_HOST_PTR = 1 << 5;

struct Api {
  DynamicLibrary lib{
#if defined(_WIN32)
      "OpenCL.dll"
#elif defined(__APPLE__)
      "/System/Library/Frameworks/OpenCL.framework/OpenCL"
#else
      "libOpenCL.so.1", "libOpenCL.so"
#endif
  };
#define AMK_CL_FN(name, ret, ...) ret(AMK_CL_CALL *name)(__VA_ARGS__) = lib.get<ret(AMK_CL_CALL *)(__VA_ARGS__)>(#name);
  AMK_CL_FN(clGetPlatformIDs, cl_int, cl_uint, cl_platform_id *, cl_uint *)
  AMK_CL_FN(clGetPlatformInfo, cl_int, cl_platform_id, cl_uint, std::size_t, void *, std::size_t *)
  AMK_CL_FN(clGetDeviceIDs, cl_int, cl_platform_id, cl_bitfield, cl_uint, cl_device_id *, cl_uint *)
  AMK_CL_FN(clGetDeviceInfo, cl_int, cl_device_id, cl_uint, std::size_t, void *, std::size_t *)
  AMK_CL_FN(clCreateContext, cl_context, const cl_context_properties *, cl_uint, const cl_device_id *, void *, void *,
            cl_int *)
  AMK_CL_FN(clCreateCommandQueue, cl_command_queue, cl_context, cl_device_id, cl_bitfield, cl_int *)
  AMK_CL_FN(clCreateProgramWithSource, cl_program, cl_context, cl_uint, const char **, const std::size_t *, cl_int *)
  AMK_CL_FN(clBuildProgram, cl_int, cl_program, cl_uint, const cl_device_id *, const char *, void *, void *)
  AMK_CL_FN(clGetProgramBuildInfo, cl_int, cl_program, cl_device_id, cl_uint, std::size_t, void *, std::size_t *)
  AMK_CL_FN(clCreateKernel, cl_kernel, cl_program, const char *, cl_int *)
  AMK_CL_FN(clCreateBuffer, cl_mem, cl_context, cl_bitfield, std::size_t, void *, cl_int *)
  AMK_CL_FN(clSetKernelArg, cl_int, cl_kernel, cl_uint, std::size_t, const void *)
  AMK_CL_FN(clEnqueueNDRangeKernel, cl_int, cl_command_queue, cl_kernel, cl_uint, const std::size_t *,
            const std::size_t *, const std::size_t *, cl_uint, const void *, void *)
  AMK_CL_FN(clEnqueueReadBuffer, cl_int, cl_command_queue, cl_mem, cl_uint, std::size_t, std::size_t, void *, cl_uint,
            const void *, void *)
  AMK_CL_FN(clEnqueueWriteBuffer, cl_int, cl_command_queue, cl_mem, cl_uint, std::size_t, std::size_t, const void *,
            cl_uint, const void *, void *)
  AMK_CL_FN(clFinish, cl_int, cl_command_queue)
  AMK_CL_FN(clReleaseMemObject, cl_int, cl_mem)
  AMK_CL_FN(clReleaseKernel, cl_int, cl_kernel)
  AMK_CL_FN(clReleaseProgram, cl_int, cl_program)
  AMK_CL_FN(clReleaseCommandQueue, cl_int, cl_command_queue)
  AMK_CL_FN(clReleaseContext, cl_int, cl_context)
#undef AMK_CL_FN
  bool ok() const {
    return clGetPlatformIDs && clGetPlatformInfo && clGetDeviceIDs && clGetDeviceInfo && clCreateContext &&
           clCreateCommandQueue && clCreateProgramWithSource && clBuildProgram && clGetProgramBuildInfo &&
           clCreateKernel && clCreateBuffer && clSetKernelArg && clEnqueueNDRangeKernel && clEnqueueReadBuffer &&
           clEnqueueWriteBuffer &&
           clFinish && clReleaseMemObject && clReleaseKernel && clReleaseProgram && clReleaseCommandQueue &&
           clReleaseContext;
  }
};

Api &api() {
  static Api a;
  return a;
}

void check(cl_int err, const char *what) {
  if (err != 0) throw std::runtime_error(std::string("OpenCL: ") + what + " failed (error " + std::to_string(err) + ")");
}

std::string device_string(cl_device_id d, cl_uint what) {
  std::size_t n = 0;
  if (api().clGetDeviceInfo(d, what, 0, nullptr, &n) != 0 || n == 0) return {};
  std::string s(n, '\0');
  api().clGetDeviceInfo(d, what, n, s.data(), nullptr);
  while (!s.empty() && s.back() == '\0') s.pop_back();
  return s;
}
std::string platform_string(cl_platform_id p, cl_uint what) {
  std::size_t n = 0;
  if (api().clGetPlatformInfo(p, what, 0, nullptr, &n) != 0 || n == 0) return {};
  std::string s(n, '\0');
  api().clGetPlatformInfo(p, what, n, s.data(), nullptr);
  while (!s.empty() && s.back() == '\0') s.pop_back();
  return s;
}

struct ClDevice {
  cl_platform_id platform = nullptr;
  cl_device_id device = nullptr;
};

//! Every device of every platform (GPUs, accelerators and CPUs), in platform order.
std::vector<ClDevice> enumerate() {
  std::vector<ClDevice> out;
  Api &a = api();
  if (!a.lib.loaded() || !a.ok()) return out;
  cl_uint np = 0;
  if (a.clGetPlatformIDs(0, nullptr, &np) != 0 || np == 0) return out;
  std::vector<cl_platform_id> ps(np);
  a.clGetPlatformIDs(np, ps.data(), nullptr);
  for (auto p : ps) {
    cl_uint nd = 0;
    const cl_bitfield types = CL_DEVICE_TYPE_GPU | CL_DEVICE_TYPE_ACCELERATOR | CL_DEVICE_TYPE_CPU;
    if (a.clGetDeviceIDs(p, types, 0, nullptr, &nd) != 0 || nd == 0) continue;
    std::vector<cl_device_id> ds(nd);
    a.clGetDeviceIDs(p, types, nd, ds.data(), nullptr);
    for (auto d : ds) out.push_back({p, d});
  }
  return out;
}

bool has_fp64(cl_device_id d) { return device_string(d, CL_DEVICE_EXTENSIONS).find("cl_khr_fp64") != std::string::npos; }

std::string device_name(cl_device_id d) {
  // AMD reports the architecture (gfx1100) as the name and the product as the board name.
  std::string board = device_string(d, CL_DEVICE_BOARD_NAME_AMD);
  std::string name = device_string(d, CL_DEVICE_NAME);
  if (!board.empty() && board != name) return board + " (" + name + ")";
  return name;
}

class OpenClDevice final : public Device {
public:
  OpenClDevice(const ClDevice &d, std::string description, bool mixed)
      : d_(d), description_(std::move(description)) {
    Api &a = api();
    cl_int err = 0;
    ctx_ = a.clCreateContext(nullptr, 1, &d_.device, nullptr, nullptr, &err);
    check(err, "clCreateContext");
    queue_ = a.clCreateCommandQueue(ctx_, d_.device, 0, &err);
    check(err, "clCreateCommandQueue");
    const std::string &src = kernel_source();
    const char *text = src.c_str();
    const std::size_t len = src.size();
    program_ = a.clCreateProgramWithSource(ctx_, 1, &text, &len, &err);
    check(err, "clCreateProgramWithSource");
    std::string options = mixed ? "-cl-std=CL1.2 -DAMK_MIXED" : "-cl-std=CL1.2";
    // Developers' profiling switches (docs/GPU.md "Profile"), e.g. -DAMK_PROF_NO_TALLY.
    if (const char *extra = std::getenv("AMCSET_KERNEL_DEFINES"); extra != nullptr && *extra != '\0') {
      options += std::string(" ") + extra;
      description_ += std::string(" [profiling: ") + extra + "]";
    }
    err = a.clBuildProgram(program_, 1, &d_.device, options.c_str(), nullptr, nullptr);
    if (mixed) description_ += " [mixed precision]";
    if (err != 0) {
      std::size_t n = 0;
      a.clGetProgramBuildInfo(program_, d_.device, CL_PROGRAM_BUILD_LOG, 0, nullptr, &n);
      std::string log(n, '\0');
      a.clGetProgramBuildInfo(program_, d_.device, CL_PROGRAM_BUILD_LOG, n, log.data(), nullptr);
      if (log.size() > 4000) log.resize(4000);
      throw std::runtime_error("OpenCL: the kernel does not compile on " + description_ + " (error " +
                               std::to_string(err) + "):\n" + log);
    }
    kernel_ = a.clCreateKernel(program_, "amcset_trim_lanes", &err);
    check(err, "clCreateKernel");
    merge_kernel_ = a.clCreateKernel(program_, "amcset_trim_merge", &err);
    check(err, "clCreateKernel (merge)");
    put_kernel_ = a.clCreateKernel(program_, "amcset_trim_put", &err);
    check(err, "clCreateKernel (put)");
    pack_kernel_ = a.clCreateKernel(program_, "amcset_trim_pack", &err);
    check(err, "clCreateKernel (pack)");
    scan_kernel_ = a.clCreateKernel(program_, "amcset_trim_scan", &err);
    check(err, "clCreateKernel (scan)");
    cl_uint cu = 0;
    a.clGetDeviceInfo(d_.device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof cu, &cu, nullptr);
    compute_units_ = cu;
  }
  ~OpenClDevice() override {
    Api &a = api();
    release_buffers();
    if (kernel_) a.clReleaseKernel(kernel_);
    for (cl_kernel k : {merge_kernel_, put_kernel_, pack_kernel_, scan_kernel_}) {
      if (k) a.clReleaseKernel(k);
    }
    if (program_) a.clReleaseProgram(program_);
    if (queue_) a.clReleaseCommandQueue(queue_);
    if (ctx_) a.clReleaseContext(ctx_);
  }

  std::string description() const override { return description_; }

  //! Eight wavefronts of 64 per compute unit (AMD: two SIMDs, wave32 or 64).
  std::size_t preferred_lanes() const override { return std::max<std::size_t>(1, compute_units_) * 512; }

  void load(const KernelInput &in, std::size_t ring, std::size_t max_lanes) override {
    release_buffers();
    in_ = &in;
    ring_ = ring;
    const KParams &P = in.params;
    const auto ro = [&](const void *data, std::size_t bytes) {
      cl_int err = 0;
      cl_mem m = api().clCreateBuffer(ctx_, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, bytes, const_cast<void *>(data), &err);
      check(err, "clCreateBuffer");
      buffers_.push_back(m);
      return m;
    };
    const auto rw = [&](std::size_t bytes) {
      cl_int err = 0;
      cl_mem m = api().clCreateBuffer(ctx_, CL_MEM_READ_WRITE, std::max<std::size_t>(bytes, 8), nullptr, &err);
      check(err, "clCreateBuffer (device memory)");
      buffers_.push_back(m);
      return m;
    };
    cl_mem args[14] = {
        ro(&in.params, sizeof(KParams)),
        ro(in.layers.data(), in.layers.size() * sizeof(KLayer)),
        ro(in.elements.data(), in.elements.size() * sizeof(KElement)),
        ro(in.pairs.data(), in.pairs.size() * sizeof(KPair)),
        ro(in.curves.data(), in.curves.size() * sizeof(KCurve)),
        ro(in.log_s.data(), in.log_s.size() * sizeof(double)),
        ro(in.edges.data(), in.edges.size() * sizeof(double)),
        ro(in.particles.data(), in.particles.size() * sizeof(KParticle)),
        rw(max_lanes * static_cast<std::size_t>(P.slot_capacity) * sizeof(KMover)),
        rw(ring * static_cast<std::size_t>(P.tally_stride) * sizeof(double)),
        rw(ring * static_cast<std::size_t>(P.list_stride) * sizeof(double)),
        rw(max_lanes * sizeof(KLane)),
        rw(2 * sizeof(std::int32_t)),
        rw(ring * sizeof(std::int64_t)),
    };
    for (cl_uint i = 0; i < 14; ++i) check(api().clSetKernelArg(kernel_, i, sizeof(cl_mem), &args[i]), "clSetKernelArg");
    params_ = args[0];
    tallies_ = args[9];
    lists_ = args[10];
    counter_ = args[12];
    slot_block_ = args[13];
    total_ = rw(stride_bytes());
    flat_ = rw(stride_bytes());
    offsets_ = rw((3 * ring + 2) * sizeof(std::int64_t));
    scan_out_ = rw(sizeof(std::int64_t));
    // No slot holds a finished block yet.
    std::vector<std::int64_t> zeros(ring, 0);
    check(api().clEnqueueWriteBuffer(queue_, slot_block_, 1, 0, ring * sizeof(std::int64_t), zeros.data(), 0, nullptr,
                                     nullptr),
          "clEnqueueWriteBuffer");
  }

  std::size_t launch(std::size_t first_block, std::size_t n_dev, std::size_t lanes, std::size_t limit,
                     std::size_t steps, bool init) override {
    Api &a = api();
    if (init) {
      const std::int32_t zero[2] = {0, 0};
      check(a.clEnqueueWriteBuffer(queue_, counter_, 1, 0, sizeof zero, zero, 0, nullptr, nullptr),
            "clEnqueueWriteBuffer");
    }
    const std::int64_t v[7] = {static_cast<std::int64_t>(lanes), static_cast<std::int64_t>(first_block),
                               static_cast<std::int64_t>(n_dev),  static_cast<std::int64_t>(ring_),
                               static_cast<std::int64_t>(limit),  static_cast<std::int64_t>(steps),
                               init ? 1 : 0};
    for (cl_uint i = 0; i < 7; ++i) check(a.clSetKernelArg(kernel_, 14 + i, sizeof(std::int64_t), &v[i]), "clSetKernelArg");
    run_1d(kernel_, lanes);
    std::int32_t counts[2] = {0, 0};
    check(a.clEnqueueReadBuffer(queue_, counter_, 1, 0, sizeof counts, counts, 0, nullptr, nullptr),
          "the kernel (clEnqueueReadBuffer)");
    return static_cast<std::size_t>(counts[1]);
  }

  std::size_t scan(std::size_t merged, std::size_t n_dev) override {
    Api &a = api();
    set_args(scan_kernel_, {slot_block_});
    const std::int64_t v[3] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(merged),
                               static_cast<std::int64_t>(n_dev)};
    for (cl_uint i = 0; i < 3; ++i) check(a.clSetKernelArg(scan_kernel_, 1 + i, sizeof(std::int64_t), &v[i]), "clSetKernelArg");
    check(a.clSetKernelArg(scan_kernel_, 4, sizeof(cl_mem), &scan_out_), "clSetKernelArg");
    run_1d(scan_kernel_, 1);
    std::int64_t c = 0;
    check(a.clEnqueueReadBuffer(queue_, scan_out_, 1, 0, sizeof c, &c, 0, nullptr, nullptr), "clEnqueueReadBuffer");
    return static_cast<std::size_t>(c);
  }

  void set_total(const double *flat) override {
    check(api().clEnqueueWriteBuffer(queue_, total_, 1, 0, stride_bytes(), flat, 0, nullptr, nullptr), "clEnqueueWriteBuffer");
  }
  void merge(std::size_t merged, std::size_t count) override {
    const std::int64_t v[3] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(merged),
                               static_cast<std::int64_t>(count)};
    set_args(merge_kernel_, {params_, tallies_, total_});
    for (cl_uint i = 0; i < 3; ++i) check(api().clSetKernelArg(merge_kernel_, 3 + i, sizeof(std::int64_t), &v[i]), "clSetKernelArg");
    run_1d(merge_kernel_, static_cast<std::size_t>(in_->params.tally_stride));
  }
  void read_total(double *flat) override {
    check(api().clEnqueueReadBuffer(queue_, total_, 1, 0, stride_bytes(), flat, 0, nullptr, nullptr), "clEnqueueReadBuffer");
  }
  void read_rows(const std::vector<std::size_t> &rows, double *out) override {
    for (std::size_t i = 0; i < rows.size(); ++i) {
      check(api().clEnqueueReadBuffer(queue_, tallies_, 0, rows[i] * ring_ * sizeof(double), ring_ * sizeof(double),
                                      out + i * ring_, 0, nullptr, nullptr),
            "clEnqueueReadBuffer");
    }
    check(api().clFinish(queue_), "clFinish");
  }
  void put_block(std::size_t slot, const double *flat) override {
    check(api().clEnqueueWriteBuffer(queue_, flat_, 1, 0, stride_bytes(), flat, 0, nullptr, nullptr), "clEnqueueWriteBuffer");
    const std::int64_t v[2] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(slot)};
    set_args(put_kernel_, {params_, tallies_, flat_});
    for (cl_uint i = 0; i < 2; ++i) check(api().clSetKernelArg(put_kernel_, 3 + i, sizeof(std::int64_t), &v[i]), "clSetKernelArg");
    run_1d(put_kernel_, static_cast<std::size_t>(in_->params.tally_stride));
  }
  void pack_lists(std::size_t merged, std::size_t count, const std::vector<std::int64_t> &offsets,
                  std::size_t n_positions, std::size_t n_exits, double *positions, double *exits) override {
    Api &a = api();
    const std::size_t pos_bytes = std::max<std::size_t>(1, n_positions) * AMK_POSITION_DOUBLES * sizeof(double);
    const std::size_t exit_bytes = std::max<std::size_t>(1, n_exits) * sizeof(double);
    if (pos_bytes > packed_positions_bytes_) packed_positions_ = grow(packed_positions_, packed_positions_bytes_, pos_bytes);
    if (exit_bytes > packed_exits_bytes_) packed_exits_ = grow(packed_exits_, packed_exits_bytes_, exit_bytes);
    check(a.clEnqueueWriteBuffer(queue_, offsets_, 1, 0, offsets.size() * sizeof(std::int64_t), offsets.data(), 0, nullptr,
                                 nullptr),
          "clEnqueueWriteBuffer");
    set_args(pack_kernel_, {params_, tallies_, lists_, offsets_, packed_positions_, packed_exits_});
    const std::int64_t v[3] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(merged),
                               static_cast<std::int64_t>(count)};
    for (cl_uint i = 0; i < 3; ++i) check(a.clSetKernelArg(pack_kernel_, 6 + i, sizeof(std::int64_t), &v[i]), "clSetKernelArg");
    run_1d(pack_kernel_, count);
    if (n_positions > 0) {
      check(a.clEnqueueReadBuffer(queue_, packed_positions_, 1, 0, n_positions * AMK_POSITION_DOUBLES * sizeof(double),
                                  positions, 0, nullptr, nullptr),
            "clEnqueueReadBuffer");
    }
    if (n_exits > 0) {
      check(a.clEnqueueReadBuffer(queue_, packed_exits_, 1, 0, n_exits * sizeof(double), exits, 0,
                                  nullptr, nullptr),
            "clEnqueueReadBuffer");
    }
  }

private:
  void release_buffers() {
    for (cl_mem m : buffers_) api().clReleaseMemObject(m);
    buffers_.clear();
    for (cl_mem m : {packed_positions_, packed_exits_}) {
      if (m) api().clReleaseMemObject(m);
    }
    packed_positions_ = packed_exits_ = nullptr;
    packed_positions_bytes_ = packed_exits_bytes_ = 0;
  }
  std::size_t stride_bytes() const { return static_cast<std::size_t>(in_->params.tally_stride) * sizeof(double); }
  void set_args(cl_kernel k, std::initializer_list<cl_mem> mems) {
    cl_uint i = 0;
    for (cl_mem m : mems) check(api().clSetKernelArg(k, i++, sizeof(cl_mem), &m), "clSetKernelArg");
  }
  void run_1d(cl_kernel k, std::size_t n) {
    const std::size_t local = 64;
    const std::size_t global = std::max<std::size_t>(1, (n + local - 1) / local) * local;
    check(api().clEnqueueNDRangeKernel(queue_, k, 1, nullptr, &global, &local, 0, nullptr, nullptr),
          "clEnqueueNDRangeKernel");
  }
  cl_mem grow(cl_mem old, std::size_t &bytes, std::size_t need) {
    if (old) api().clReleaseMemObject(old);
    cl_int err = 0;
    bytes = std::max(need, bytes * 2);
    cl_mem m = api().clCreateBuffer(ctx_, CL_MEM_READ_WRITE, bytes, nullptr, &err);
    check(err, "clCreateBuffer (lists)");
    return m;
  }
  cl_kernel merge_kernel_ = nullptr, put_kernel_ = nullptr, pack_kernel_ = nullptr, scan_kernel_ = nullptr;
  cl_mem params_ = nullptr, total_ = nullptr, flat_ = nullptr, offsets_ = nullptr, slot_block_ = nullptr,
         scan_out_ = nullptr;
  std::size_t ring_ = 1;
  cl_mem packed_positions_ = nullptr, packed_exits_ = nullptr;
  std::size_t packed_positions_bytes_ = 0, packed_exits_bytes_ = 0;
  ClDevice d_;
  std::string description_;
  cl_context ctx_ = nullptr;
  cl_command_queue queue_ = nullptr;
  cl_program program_ = nullptr;
  cl_kernel kernel_ = nullptr;
  std::vector<cl_mem> buffers_;
  cl_mem tallies_ = nullptr;
  cl_mem lists_ = nullptr;
  cl_mem counter_ = nullptr;
  std::size_t compute_units_ = 0;
  const KernelInput *in_ = nullptr;
};

} // namespace

void list_opencl(std::vector<DeviceInfo> &out, std::vector<std::string> *notes) {
  Api &a = api();
  if (!a.lib.loaded()) {
    if (notes) notes->push_back("OpenCL: no OpenCL runtime (ICD loader) found");
    return;
  }
  if (!a.ok()) {
    if (notes) notes->push_back("OpenCL: the runtime " + a.lib.name() + " lacks OpenCL 1.2 functions");
    return;
  }
  const auto devices = enumerate();
  if (devices.empty() && notes) notes->push_back("OpenCL: no devices");
  int index = 0;
  for (const auto &d : devices) {
    DeviceInfo info;
    info.backend = Backend::opencl;
    info.index = index++;
    info.name = device_name(d.device);
    cl_uint cu = 0;
    cl_ulong mem = 0;
    a.clGetDeviceInfo(d.device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof cu, &cu, nullptr);
    a.clGetDeviceInfo(d.device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof mem, &mem, nullptr);
    info.details = platform_string(d.platform, CL_PLATFORM_NAME) + ", " + device_string(d.device, CL_DEVICE_VERSION) +
                   ", driver " + device_string(d.device, CL_DRIVER_VERSION) + ", " + std::to_string(cu) +
                   " compute units, " + std::to_string(mem >> 20) + " MiB";
    info.usable = has_fp64(d.device);
    if (!info.usable) info.problem = "no double precision (cl_khr_fp64)";
    out.push_back(info);
  }
}

std::unique_ptr<Device> open_opencl(const Options &o) {
  Api &a = api();
  if (!a.lib.loaded()) throw std::runtime_error("OpenCL: no OpenCL runtime (ICD loader) on this machine");
  if (!a.ok()) throw std::runtime_error("OpenCL: " + a.lib.name() + " lacks OpenCL 1.2 functions");
  const auto devices = enumerate();
  if (o.device < 0 || static_cast<std::size_t>(o.device) >= devices.size()) {
    throw std::runtime_error("OpenCL: no device " + std::to_string(o.device) + " (" + std::to_string(devices.size()) +
                             " found; amcset devices lists them)");
  }
  const ClDevice &d = devices[static_cast<std::size_t>(o.device)];
  const std::string desc = "OpenCL device " + std::to_string(o.device) + ": " + device_name(d.device) + " (" +
                           platform_string(d.platform, CL_PLATFORM_NAME) + ")";
  if (!has_fp64(d.device)) throw std::runtime_error(desc + " has no double precision (cl_khr_fp64)");
  return std::make_unique<OpenClDevice>(d, desc, o.precision == Precision::mixed);
}

} // namespace compute
} // namespace amcset
