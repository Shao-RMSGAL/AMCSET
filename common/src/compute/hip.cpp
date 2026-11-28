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

/*! \file hip.cpp
 *  \brief The HIP (AMD) and CUDA (NVIDIA) backends.
 *
 * Both load the vendor's driver API at run time and launch the kernel of
 * trim_kernel.h as a module:
 *
 * - HIP: the runtime of the AMD driver (amdhip64_7.dll, libamdhip64.so). The
 *   kernel comes from the code object compiled with hipcc when AMCSET was
 *   built (AMCSET_HIP_ARCHITECTURES; hip_code_object()), or, for another
 *   GPU, from hiprtc at run time (the HIP SDK's hiprtc*.dll, libhiprtc.so;
 *   AMCSET_HIPRTC names the library).
 * - CUDA: the driver API (nvcuda.dll, libcuda.so.1) and NVRTC from the CUDA
 *   toolkit (nvrtc64_*.dll, libnvrtc.so; AMCSET_NVRTC names it), which
 *   compiles the kernel to PTX for the device's compute capability.
 *
 * The functions and constants are declared here from the vendors' public
 * headers (hip_runtime_api.h of HIP 6/7, cuda.h, nvrtc.h): AMCSET builds
 * without either SDK. HIP is the CUDA driver API renamed, so one class
 * serves both. The CUDA backend has not been run on NVIDIA hardware yet.
 */

#include "device.h"
#include "dynlib.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace amcset {
namespace compute {

//! The HIP code objects compiled at build time, double and mixed precision
//! (empty if AMCSET was built without hipcc).
const std::string &hip_code_object();
const std::string &hip_code_object_mixed();

namespace {

#if defined(_WIN32)
#define AMK_GPU_CALL __stdcall
#else
#define AMK_GPU_CALL
#endif

using GpuResult = int; // hipError_t, CUresult, hiprtcResult, nvrtcResult: 0 = success
using Module = void *;
using Function = void *;
using DevicePtr = std::uint64_t; // hipDeviceptr_t / CUdeviceptr (64-bit)

//! The driver functions both vendors share (HIP's names in brackets).
struct DriverApi {
  bool cuda = false;
  DynamicLibrary lib;
  GpuResult(AMK_GPU_CALL *init)(unsigned) = nullptr;                          // cuInit [hipInit]
  GpuResult(AMK_GPU_CALL *device_count)(int *) = nullptr;                     // cuDeviceGetCount
  GpuResult(AMK_GPU_CALL *device_get)(int *, int) = nullptr;                  // cuDeviceGet
  GpuResult(AMK_GPU_CALL *device_name)(char *, int, int) = nullptr;           // cuDeviceGetName
  GpuResult(AMK_GPU_CALL *device_attribute)(int *, int, int) = nullptr;       // cuDeviceGetAttribute
  GpuResult(AMK_GPU_CALL *device_total_mem)(std::size_t *, int) = nullptr;    // cuDeviceTotalMem_v2
  GpuResult(AMK_GPU_CALL *primary_retain)(void **, int) = nullptr;            // cuDevicePrimaryCtxRetain
  GpuResult(AMK_GPU_CALL *primary_release)(int) = nullptr;                    // cuDevicePrimaryCtxRelease
  GpuResult(AMK_GPU_CALL *ctx_set_current)(void *) = nullptr;                 // cuCtxSetCurrent
  GpuResult(AMK_GPU_CALL *set_device)(int) = nullptr;                         // [hipSetDevice]
  GpuResult(AMK_GPU_CALL *module_load_data)(Module *, const void *) = nullptr;
  GpuResult(AMK_GPU_CALL *module_unload)(Module) = nullptr;
  GpuResult(AMK_GPU_CALL *module_get_function)(Function *, Module, const char *) = nullptr;
  GpuResult(AMK_GPU_CALL *launch)(Function, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned, unsigned,
                                  void *, void **, void **) = nullptr;
  GpuResult(AMK_GPU_CALL *mem_alloc)(DevicePtr *, std::size_t) = nullptr;
  GpuResult(AMK_GPU_CALL *mem_free)(DevicePtr) = nullptr;
  GpuResult(AMK_GPU_CALL *htod)(DevicePtr, const void *, std::size_t) = nullptr;
  GpuResult(AMK_GPU_CALL *dtoh)(void *, DevicePtr, std::size_t) = nullptr;
  GpuResult(AMK_GPU_CALL *synchronize)() = nullptr;                           // cuCtxSynchronize [hipDeviceSynchronize]
  GpuResult(AMK_GPU_CALL *driver_version)(int *) = nullptr;
  const char *(AMK_GPU_CALL *hip_error_string)(GpuResult) = nullptr;          // [hipGetErrorString]
  GpuResult(AMK_GPU_CALL *cu_error_string)(GpuResult, const char **) = nullptr; // cuGetErrorString
  // HIP only: the device properties (hipDeviceProp_tR0600, 1472 bytes in HIP
  // 6/7; gcnArchName, the "gfx1100:sramecc-:xnack-" target, at byte 1160).
  GpuResult(AMK_GPU_CALL *hip_properties)(void *, int) = nullptr;

  bool ok() const {
    return init && device_count && device_get && device_name && device_attribute && device_total_mem &&
           module_load_data && module_unload && module_get_function && launch && mem_alloc && mem_free && htod &&
           dtoh && synchronize && (cuda ? primary_retain && ctx_set_current : set_device != nullptr);
  }

  std::string error(GpuResult r) const {
    const char *s = nullptr;
    if (hip_error_string) s = hip_error_string(r);
    if (cu_error_string) cu_error_string(r, &s);
    return std::to_string(r) + (s ? std::string(" ") + s : std::string());
  }
};

template <class F> void bind(const DynamicLibrary &lib, F &f, const char *name) { f = lib.get<F>(name); }

std::unique_ptr<DriverApi> load_driver(bool cuda) {
  auto a = std::make_unique<DriverApi>();
  a->cuda = cuda;
  if (cuda) {
#if defined(_WIN32)
    a->lib.open("nvcuda.dll");
#else
    if (!a->lib.open("libcuda.so.1")) a->lib.open("libcuda.so");
#endif
    const auto &l = a->lib;
    bind(l, a->init, "cuInit");
    bind(l, a->device_count, "cuDeviceGetCount");
    bind(l, a->device_get, "cuDeviceGet");
    bind(l, a->device_name, "cuDeviceGetName");
    bind(l, a->device_attribute, "cuDeviceGetAttribute");
    bind(l, a->device_total_mem, "cuDeviceTotalMem_v2");
    bind(l, a->primary_retain, "cuDevicePrimaryCtxRetain");
    bind(l, a->primary_release, "cuDevicePrimaryCtxRelease_v2");
    bind(l, a->ctx_set_current, "cuCtxSetCurrent");
    bind(l, a->module_load_data, "cuModuleLoadData");
    bind(l, a->module_unload, "cuModuleUnload");
    bind(l, a->module_get_function, "cuModuleGetFunction");
    bind(l, a->launch, "cuLaunchKernel");
    bind(l, a->mem_alloc, "cuMemAlloc_v2");
    bind(l, a->mem_free, "cuMemFree_v2");
    bind(l, a->htod, "cuMemcpyHtoD_v2");
    bind(l, a->dtoh, "cuMemcpyDtoH_v2");
    bind(l, a->synchronize, "cuCtxSynchronize");
    bind(l, a->driver_version, "cuDriverGetVersion");
    bind(l, a->cu_error_string, "cuGetErrorString");
  } else {
#if defined(_WIN32)
    if (!a->lib.open("amdhip64_7.dll")) a->lib.open("amdhip64_6.dll");
#else
    if (!a->lib.open("libamdhip64.so.7") && !a->lib.open("libamdhip64.so.6")) a->lib.open("libamdhip64.so");
#endif
    const auto &l = a->lib;
    bind(l, a->init, "hipInit");
    bind(l, a->device_count, "hipGetDeviceCount");
    bind(l, a->device_get, "hipDeviceGet");
    bind(l, a->device_name, "hipDeviceGetName");
    bind(l, a->device_attribute, "hipDeviceGetAttribute");
    bind(l, a->device_total_mem, "hipDeviceTotalMem");
    bind(l, a->set_device, "hipSetDevice");
    bind(l, a->module_load_data, "hipModuleLoadData");
    bind(l, a->module_unload, "hipModuleUnload");
    bind(l, a->module_get_function, "hipModuleGetFunction");
    bind(l, a->launch, "hipModuleLaunchKernel");
    bind(l, a->mem_alloc, "hipMalloc");
    bind(l, a->mem_free, "hipFree");
    bind(l, a->htod, "hipMemcpyHtoD");
    bind(l, a->dtoh, "hipMemcpyDtoH");
    bind(l, a->synchronize, "hipDeviceSynchronize");
    bind(l, a->driver_version, "hipDriverGetVersion");
    bind(l, a->hip_error_string, "hipGetErrorString");
    bind(l, a->hip_properties, "hipGetDevicePropertiesR0600");
  }
  return a;
}

//! Device attributes: CUDA's (cuda.h) and HIP's (hip_runtime_api.h) numbers.
int attr_multiprocessors(bool cuda) { return cuda ? 16 : 63; }
int attr_cc_major(bool cuda) { return cuda ? 75 : 23; }
int attr_cc_minor(bool cuda) { return cuda ? 76 : 61; }

//! One shared driver per vendor (initialised once).
DriverApi *driver(bool cuda, std::string *problem) {
  static std::once_flag once[2];
  static std::unique_ptr<DriverApi> d[2];
  static std::string why[2];
  const int i = cuda ? 1 : 0;
  std::call_once(once[i], [&] {
    d[i] = load_driver(cuda);
    if (!d[i]->lib.loaded()) {
      why[i] = cuda ? "no NVIDIA driver (nvcuda / libcuda)" : "no AMD HIP runtime (amdhip64 / libamdhip64)";
    } else if (!d[i]->ok()) {
      why[i] = d[i]->lib.name() + " lacks functions the backend needs";
    } else if (const GpuResult r = d[i]->init(0); r != 0) {
      why[i] = d[i]->lib.name() + ": initialisation failed (" + d[i]->error(r) + ")";
    }
  });
  if (problem) *problem = why[i];
  return why[i].empty() ? d[i].get() : nullptr;
}

std::string hip_arch(DriverApi &a, int device) {
  if (!a.hip_properties) return {};
  std::vector<char> buf(8192, '\0');
  if (a.hip_properties(buf.data(), device) != 0) return {};
  std::string arch(buf.data() + 1160, ::strnlen(buf.data() + 1160, 256));
  return arch;
}

// --- Runtime compilers --------------------------------------------------------

//! hiprtc and NVRTC share their API (hiprtcX = nvrtcX).
struct RtcApi {
  DynamicLibrary lib;
  GpuResult(AMK_GPU_CALL *create)(void **, const char *, const char *, int, const char *const *,
                                  const char *const *) = nullptr;
  GpuResult(AMK_GPU_CALL *compile)(void *, int, const char *const *) = nullptr;
  GpuResult(AMK_GPU_CALL *log_size)(void *, std::size_t *) = nullptr;
  GpuResult(AMK_GPU_CALL *log)(void *, char *) = nullptr;
  GpuResult(AMK_GPU_CALL *code_size)(void *, std::size_t *) = nullptr;
  GpuResult(AMK_GPU_CALL *code)(void *, char *) = nullptr;
  GpuResult(AMK_GPU_CALL *destroy)(void **) = nullptr;
  bool ok() const { return create && compile && log_size && log && code_size && code && destroy; }
};

//! Folders to search for a runtime compiler: $env, the SDK's, next to AMCSET.
std::vector<std::filesystem::path> rtc_folders(bool cuda) {
  std::vector<std::filesystem::path> dirs;
  const auto env = [&](const char *v, const char *sub) {
    if (const char *p = std::getenv(v); p != nullptr && *p != '\0') dirs.push_back(std::filesystem::path(p) / sub);
  };
  if (cuda) {
    env("CUDA_PATH", "bin");
    env("CUDA_PATH", "bin/x64");
    env("CUDA_HOME", "lib64");
#if !defined(_WIN32)
    dirs.emplace_back("/usr/local/cuda/lib64");
#endif
  } else {
    env("HIP_PATH", "bin");
    env("ROCM_PATH", "bin");
    env("HIP_PATH", "lib");
    env("ROCM_PATH", "lib");
#if defined(_WIN32)
    std::error_code ec;
    for (const auto &e : std::filesystem::directory_iterator("C:/Program Files/AMD/ROCm", ec)) dirs.push_back(e.path() / "bin");
#else
    dirs.emplace_back("/opt/rocm/lib");
#endif
  }
  return dirs;
}

std::unique_ptr<RtcApi> load_rtc(bool cuda, std::string &tried) {
  auto r = std::make_unique<RtcApi>();
  const char *override_var = cuda ? "AMCSET_NVRTC" : "AMCSET_HIPRTC";
  const std::string prefix = cuda ? "nvrtc" : "hiprtc";
  if (const char *p = std::getenv(override_var); p != nullptr && *p != '\0') {
    r->lib.open(p);
    tried += std::string(p) + " ";
  }
  if (!r->lib.loaded()) {
#if defined(_WIN32)
    // Versioned names (hiprtc0605.dll, hiprtc0717.dll, nvrtc64_120_0.dll):
    // the newest in the folders above, then the plain search path.
    for (const auto &dir : rtc_folders(cuda)) {
      std::error_code ec;
      std::vector<std::filesystem::path> found;
      for (const auto &e : std::filesystem::directory_iterator(dir, ec)) {
        const std::string n = e.path().filename().string();
        if (n.rfind(prefix, 0) == 0 && e.path().extension() == ".dll" && n.find("builtins") == std::string::npos) {
          found.push_back(e.path());
        }
      }
      std::sort(found.rbegin(), found.rend());
      for (const auto &f : found) {
        tried += f.string() + " ";
        if (r->lib.open(f.string())) break;
      }
      if (r->lib.loaded()) break;
    }
#else
    for (const auto &dir : rtc_folders(cuda)) {
      const auto f = dir / ("lib" + prefix + ".so");
      tried += f.string() + " ";
      if (r->lib.open(f.string())) break;
    }
    if (!r->lib.loaded()) r->lib.open("lib" + prefix + ".so");
#endif
  }
  const std::string p = cuda ? "nvrtc" : "hiprtc";
  bind(r->lib, r->create, (p + "CreateProgram").c_str());
  bind(r->lib, r->compile, (p + "CompileProgram").c_str());
  bind(r->lib, r->log_size, (p + "GetProgramLogSize").c_str());
  bind(r->lib, r->log, (p + "GetProgramLog").c_str());
  bind(r->lib, r->code_size, (cuda ? "nvrtcGetPTXSize" : "hiprtcGetCodeSize"));
  bind(r->lib, r->code, (cuda ? "nvrtcGetPTX" : "hiprtcGetCode"));
  bind(r->lib, r->destroy, (p + "DestroyProgram").c_str());
  return r;
}

//! The kernel compiled at run time for \p arch (gfx1100, compute_86).
std::string compile_at_run_time(bool cuda, const std::string &arch, bool mixed) {
  std::string tried;
  const auto rtc = load_rtc(cuda, tried);
  const char *name = cuda ? "NVRTC" : "hiprtc";
  if (!rtc->lib.loaded() || !rtc->ok()) {
    throw std::runtime_error(std::string(name) + " (the runtime compiler) was not found" +
                             (tried.empty() ? std::string() : " (tried " + tried + ")") + "; install the " +
                             (cuda ? "CUDA toolkit or set AMCSET_NVRTC" : "HIP SDK or set AMCSET_HIPRTC") +
                             " to the library");
  }
  void *prog = nullptr;
  const std::string &src = kernel_source();
  if (rtc->create(&prog, src.c_str(), "trim_kernel.cu", 0, nullptr, nullptr) != 0) {
    throw std::runtime_error(std::string(name) + ": cannot create the program");
  }
  std::vector<std::string> opts;
  if (cuda) {
    opts = {"--gpu-architecture=" + arch, "--fmad=false", "-default-device"};
  } else {
    opts = {"-O3", "-ffp-contract=off"};
    if (!arch.empty()) opts.push_back("--offload-arch=" + arch);
  }
  if (mixed) opts.push_back("-DAMK_MIXED");
  std::vector<const char *> o;
  for (const auto &s : opts) o.push_back(s.c_str());
  const GpuResult r = rtc->compile(prog, static_cast<int>(o.size()), o.data());
  std::size_t n = 0;
  rtc->log_size(prog, &n);
  std::string log(n, '\0');
  if (n > 0) rtc->log(prog, log.data());
  if (r != 0) {
    rtc->destroy(&prog);
    if (log.size() > 4000) log.resize(4000);
    throw std::runtime_error(std::string(name) + ": the kernel does not compile (" + std::to_string(r) + "):\n" + log);
  }
  rtc->code_size(prog, &n);
  std::string code(n, '\0');
  rtc->code(prog, code.data());
  rtc->destroy(&prog);
  return code;
}

class GpuDevice final : public Device {
public:
  GpuDevice(DriverApi &a, int index, std::string description, bool mixed)
      : a_(a), index_(index), description_(std::move(description)) {
    if (a_.cuda) {
      int dev = 0;
      check(a_.device_get(&dev, index_), "cuDeviceGet");
      check(a_.primary_retain(&ctx_, dev), "cuDevicePrimaryCtxRetain");
      check(a_.ctx_set_current(ctx_), "cuCtxSetCurrent");
      int major = 0, minor = 0;
      a_.device_attribute(&major, attr_cc_major(true), dev);
      a_.device_attribute(&minor, attr_cc_minor(true), dev);
      code_ = compile_at_run_time(true, "compute_" + std::to_string(major) + std::to_string(minor), mixed);
      origin_ = "NVRTC";
    } else {
      check(a_.set_device(index_), "hipSetDevice");
      // The code object built with AMCSET first (if it has this GPU), then hiprtc.
      const std::string &prebuilt = mixed ? hip_code_object_mixed() : hip_code_object();
      if (!prebuilt.empty() && a_.module_load_data(&module_, prebuilt.data()) == 0) {
        origin_ = "built-in code object";
      } else {
        module_ = nullptr;
        std::string arch = hip_arch(a_, index_);
        if (const auto colon = arch.find(':'); colon != std::string::npos) arch.resize(colon);
        code_ = compile_at_run_time(false, arch, mixed);
        origin_ = "hiprtc";
      }
    }
    if (mixed) origin_ += ", mixed precision";
    if (module_ == nullptr) check(a_.module_load_data(&module_, code_.data()), "loading the kernel module");
    check(a_.module_get_function(&fn_, module_, "amcset_trim_lanes"), "finding the kernel");
    check(a_.module_get_function(&merge_fn_, module_, "amcset_trim_merge"), "finding the merge kernel");
    check(a_.module_get_function(&put_fn_, module_, "amcset_trim_put"), "finding the put kernel");
    check(a_.module_get_function(&pack_fn_, module_, "amcset_trim_pack"), "finding the pack kernel");
    check(a_.module_get_function(&scan_fn_, module_, "amcset_trim_scan"), "finding the scan kernel");
    description_ += " [" + origin_ + "]";
    int dev = index_;
    if (a_.cuda) a_.device_get(&dev, index_);
    int mp = 0;
    a_.device_attribute(&mp, attr_multiprocessors(a_.cuda), dev);
    compute_units_ = static_cast<std::size_t>(std::max(mp, 1));
  }

  //! Eight wavefronts of 64 per compute unit (AMD) or 16 warps per multiprocessor (NVIDIA).
  std::size_t preferred_lanes() const override { return compute_units_ * 512; }
  ~GpuDevice() override {
    activate();
    free_buffers();
    if (module_) a_.module_unload(module_);
    if (a_.cuda && a_.primary_release) {
      int dev = 0;
      if (a_.device_get(&dev, index_) == 0) a_.primary_release(dev);
    }
  }

  std::string description() const override { return description_; }

  void load(const KernelInput &in, std::size_t ring, std::size_t max_lanes) override {
    activate();
    free_buffers();
    in_ = &in;
    ring_ = ring;
    const KParams &P = in.params;
    const auto upload = [&](const void *data, std::size_t bytes) {
      DevicePtr p = 0;
      check(a_.mem_alloc(&p, std::max<std::size_t>(bytes, 8)), "allocating device memory");
      buffers_.push_back(p);
      if (bytes > 0) check(a_.htod(p, data, bytes), "copying to the device");
      return p;
    };
    const auto alloc = [&](std::size_t bytes) {
      DevicePtr p = 0;
      check(a_.mem_alloc(&p, std::max<std::size_t>(bytes, 8)), "allocating device memory");
      buffers_.push_back(p);
      return p;
    };
    args_[0] = upload(&in.params, sizeof(KParams));
    args_[1] = upload(in.layers.data(), in.layers.size() * sizeof(KLayer));
    args_[2] = upload(in.elements.data(), in.elements.size() * sizeof(KElement));
    args_[3] = upload(in.pairs.data(), in.pairs.size() * sizeof(KPair));
    args_[4] = upload(in.curves.data(), in.curves.size() * sizeof(KCurve));
    args_[5] = upload(in.log_s.data(), in.log_s.size() * sizeof(double));
    args_[6] = upload(in.edges.data(), in.edges.size() * sizeof(double));
    args_[7] = upload(in.particles.data(), in.particles.size() * sizeof(KParticle));
    args_[8] = alloc(max_lanes * static_cast<std::size_t>(P.slot_capacity) * sizeof(KMover));
    args_[9] = alloc(ring * static_cast<std::size_t>(P.tally_stride) * sizeof(double));
    args_[10] = alloc(ring * static_cast<std::size_t>(P.list_stride) * sizeof(double));
    args_[11] = alloc(max_lanes * sizeof(KLane));
    args_[12] = alloc(2 * sizeof(std::int32_t));
    std::vector<std::int64_t> zeros(ring, 0); // no slot holds a finished block yet
    args_[13] = upload(zeros.data(), ring * sizeof(std::int64_t));
    total_ = alloc(stride_bytes());
    flat_ = alloc(stride_bytes());
    offsets_ = alloc((3 * ring + 2) * sizeof(std::int64_t));
    scan_out_ = alloc(sizeof(std::int64_t));
    packed_positions_ = packed_exits_ = 0;
    packed_positions_bytes_ = packed_exits_bytes_ = 0;
  }

  std::size_t launch(std::size_t first_block, std::size_t n_dev, std::size_t lanes, std::size_t limit,
                     std::size_t steps, bool init) override {
    activate();
    if (init) {
      const std::int32_t zero[2] = {0, 0};
      check(a_.htod(args_[12], zero, sizeof zero), "copying to the device");
    }
    std::int64_t v[7] = {static_cast<std::int64_t>(lanes), static_cast<std::int64_t>(first_block),
                         static_cast<std::int64_t>(n_dev),  static_cast<std::int64_t>(ring_),
                         static_cast<std::int64_t>(limit),  static_cast<std::int64_t>(steps),
                         init ? 1 : 0};
    void *params[21];
    for (int i = 0; i < 14; ++i) params[i] = &args_[i];
    for (int i = 0; i < 7; ++i) params[14 + i] = &v[i];
    run_1d(fn_, lanes, params);
    std::int32_t counts[2] = {0, 0};
    check(a_.dtoh(counts, args_[12], sizeof counts), "copying from the device");
    return static_cast<std::size_t>(counts[1]);
  }

  std::size_t scan(std::size_t merged, std::size_t n_dev) override {
    activate();
    std::int64_t v[3] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(merged),
                         static_cast<std::int64_t>(n_dev)};
    void *params[5] = {&args_[13], &v[0], &v[1], &v[2], &scan_out_};
    run_1d(scan_fn_, 1, params);
    std::int64_t c = 0;
    check(a_.dtoh(&c, scan_out_, sizeof c), "copying from the device");
    return static_cast<std::size_t>(c);
  }
  void set_total(const double *flat) override {
    activate();
    check(a_.htod(total_, flat, stride_bytes()), "copying to the device");
  }
  void merge(std::size_t merged, std::size_t count) override {
    activate();
    std::int64_t v[3] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(merged),
                         static_cast<std::int64_t>(count)};
    void *params[6] = {&args_[0], &args_[9], &total_, &v[0], &v[1], &v[2]};
    run_1d(merge_fn_, static_cast<std::size_t>(in_->params.tally_stride), params);
  }
  void read_total(double *flat) override {
    activate();
    check(a_.dtoh(flat, total_, stride_bytes()), "copying from the device");
  }
  void read_rows(const std::vector<std::size_t> &rows, double *out) override {
    activate();
    for (std::size_t i = 0; i < rows.size(); ++i) {
      check(a_.dtoh(out + i * ring_, args_[9] + rows[i] * ring_ * sizeof(double), ring_ * sizeof(double)),
            "copying from the device");
    }
  }
  void put_block(std::size_t slot, const double *flat) override {
    activate();
    check(a_.htod(flat_, flat, stride_bytes()), "copying to the device");
    std::int64_t v[2] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(slot)};
    void *params[5] = {&args_[0], &args_[9], &flat_, &v[0], &v[1]};
    run_1d(put_fn_, static_cast<std::size_t>(in_->params.tally_stride), params);
  }
  void pack_lists(std::size_t merged, std::size_t count, const std::vector<std::int64_t> &offsets,
                  std::size_t n_positions, std::size_t n_exits, double *positions, double *exits) override {
    activate();
    const std::size_t pos_bytes = std::max<std::size_t>(1, n_positions) * AMK_POSITION_DOUBLES * sizeof(double);
    const std::size_t exit_bytes = std::max<std::size_t>(1, n_exits) * sizeof(double);
    if (pos_bytes > packed_positions_bytes_) grow(packed_positions_, packed_positions_bytes_, pos_bytes);
    if (exit_bytes > packed_exits_bytes_) grow(packed_exits_, packed_exits_bytes_, exit_bytes);
    check(a_.htod(offsets_, offsets.data(), offsets.size() * sizeof(std::int64_t)), "copying to the device");
    std::int64_t v[3] = {static_cast<std::int64_t>(ring_), static_cast<std::int64_t>(merged),
                         static_cast<std::int64_t>(count)};
    void *params[9] = {&args_[0], &args_[9], &args_[10], &offsets_, &packed_positions_, &packed_exits_,
                       &v[0], &v[1], &v[2]};
    run_1d(pack_fn_, count, params);
    if (n_positions > 0) {
      check(a_.dtoh(positions, packed_positions_, n_positions * AMK_POSITION_DOUBLES * sizeof(double)),
            "copying from the device");
    }
    if (n_exits > 0) {
      check(a_.dtoh(exits, packed_exits_, n_exits * sizeof(double)), "copying from the device");
    }
  }

private:
  void check(GpuResult r, const char *what) const {
    if (r != 0) throw std::runtime_error((a_.cuda ? "CUDA: " : "HIP: ") + std::string(what) + " failed (" + a_.error(r) + ")");
  }
  //! The calling thread's context/device (the engine's device thread is not the one that opened it).
  void activate() const {
    if (a_.cuda) a_.ctx_set_current(ctx_);
    else a_.set_device(index_);
  }
  void free_buffers() {
    for (DevicePtr p : buffers_) a_.mem_free(p);
    buffers_.clear();
  }
  std::size_t stride_bytes() const { return static_cast<std::size_t>(in_->params.tally_stride) * sizeof(double); }
  void run_1d(Function f, std::size_t n, void **params) {
    const unsigned threads = 64;
    const auto grid = static_cast<unsigned>(std::max<std::size_t>(1, (n + threads - 1) / threads));
    check(a_.launch(f, grid, 1, 1, threads, 1, 1, 0, nullptr, params, nullptr), "launching a kernel");
    check(a_.synchronize(), "a kernel");
  }
  //! A list buffer of at least \p need bytes (kept in buffers_, freed with them).
  void grow(DevicePtr &p, std::size_t &bytes, std::size_t need) {
    bytes = std::max(need, bytes * 2);
    DevicePtr q = 0;
    check(a_.mem_alloc(&q, bytes), "allocating device memory");
    buffers_.push_back(q);
    p = q;
  }
  Function merge_fn_ = nullptr, put_fn_ = nullptr, pack_fn_ = nullptr, scan_fn_ = nullptr;
  DevicePtr total_ = 0, flat_ = 0, offsets_ = 0, packed_positions_ = 0, packed_exits_ = 0, scan_out_ = 0;
  std::size_t ring_ = 1;
  std::size_t packed_positions_bytes_ = 0, packed_exits_bytes_ = 0;
  DriverApi &a_;
  int index_ = 0;
  std::string description_;
  std::string origin_;
  void *ctx_ = nullptr;
  Module module_ = nullptr;
  Function fn_ = nullptr;
  std::string code_;
  std::vector<DevicePtr> buffers_;
  DevicePtr args_[14] = {};
  std::size_t compute_units_ = 1;
  const KernelInput *in_ = nullptr;
};

std::string device_label(DriverApi &a, int i) {
  char name[256] = {};
  a.device_name(name, sizeof name, i);
  return name;
}

} // namespace

void list_hip(std::vector<DeviceInfo> &out, std::vector<std::string> *notes, bool cuda) {
  std::string problem;
  DriverApi *a = driver(cuda, &problem);
  const char *what = cuda ? "CUDA" : "HIP";
  if (a == nullptr) {
    if (notes) notes->push_back(std::string(what) + ": " + problem);
    return;
  }
  int n = 0;
  if (a->device_count(&n) != 0) n = 0;
  if (n == 0 && notes) notes->push_back(std::string(what) + ": no devices");
  for (int i = 0; i < n; ++i) {
    DeviceInfo info;
    info.backend = cuda ? Backend::cuda : Backend::hip;
    info.index = i;
    info.name = device_label(*a, i);
    int dev = i;
    if (cuda) a->device_get(&dev, i);
    int mp = 0;
    std::size_t mem = 0;
    a->device_attribute(&mp, attr_multiprocessors(cuda), dev);
    a->device_total_mem(&mem, dev);
    int ver = 0;
    if (a->driver_version) a->driver_version(&ver);
    info.details = a->lib.name() + ", driver " + std::to_string(ver) + ", " + std::to_string(mp) +
                   (cuda ? " multiprocessors, " : " compute units, ") + std::to_string(mem >> 20) + " MiB";
    if (cuda) {
      int major = 0, minor = 0;
      a->device_attribute(&major, attr_cc_major(true), dev);
      a->device_attribute(&minor, attr_cc_minor(true), dev);
      info.details += ", compute capability " + std::to_string(major) + "." + std::to_string(minor);
    } else {
      const std::string arch = hip_arch(*a, i);
      if (!arch.empty()) info.details += ", " + arch;
    }
    out.push_back(info);
  }
}

std::unique_ptr<Device> open_hip(const Options &o, bool cuda) {
  std::string problem;
  DriverApi *a = driver(cuda, &problem);
  const std::string what = cuda ? "CUDA" : "HIP";
  if (a == nullptr) throw std::runtime_error(what + ": " + problem);
  int n = 0;
  if (a->device_count(&n) != 0) n = 0;
  if (o.device < 0 || o.device >= n) {
    throw std::runtime_error(what + ": no device " + std::to_string(o.device) + " (" + std::to_string(n) +
                             " found; amcset devices lists them)");
  }
  const std::string desc = what + " device " + std::to_string(o.device) + ": " + device_label(*a, o.device);
  return std::make_unique<GpuDevice>(*a, o.device, desc, o.precision == Precision::mixed);
}

} // namespace compute
} // namespace amcset
