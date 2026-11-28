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

/*! \file trim_kernel.h
 *  \brief The trim2013 transport as one portable kernel source (GPU backends).
 *
 * This file is compiled four ways from the same text: as OpenCL C (the
 * OpenCL backend builds it at run time), as HIP or CUDA (hiprtc / NVRTC at
 * run time, or hipcc ahead of time), and as C++ (the reference backend that
 * runs it on the CPU, and the host side, which shares its data layout). It
 * is therefore written in the common subset of C99 and C++: no references,
 * templates, classes, standard headers or recursion, only plain structs of
 * doubles and integers whose layout is the same on every compiler.
 *
 * It reproduces Transporter::run_ion_2013 and the functions it calls
 * (common/src/trim/trim_engine.cpp) for the beam and TRIM.DAT ion sources,
 * expression by expression and in the same order of random draws. Lengths
 * are metres and energies joules, as the Boost.Units quantities of the CPU
 * transport hold them, and every conversion is the same division or
 * multiplication by the unit's value (KParams::ang, ::ev), so that the C++
 * build gives bit-identical tallies (test ComputeKernel.*). On a GPU the
 * vendor's exp/log/pow/sin/cos differ from the C library's in the last bit,
 * so histories part after a while: the results are statistically the same
 * and, for a given device and driver, deterministic (every block of 16 ions
 * is one work item that writes its own tallies; there are no atomics).
 *
 * What it leaves out runs on the CPU instead (compute::device_supports):
 * trajectories, MD events, COLLISON.txt details, EXYZ.txt, TRIM.DAT
 * recoils, the refined and crystal transports and fast_math's tables. A
 * block whose recoil stack or exit list outgrows its fixed capacity stops
 * with a status and is rerun on the CPU.
 *
 * Keep it in step with trim_engine.cpp: a change to the trim2013 transport
 * makes ComputeKernel.MatchesTheCpuTransport fail until it is made here too.
 */

#ifndef AMCSET_TRIM_KERNEL_H
#define AMCSET_TRIM_KERNEL_H

// ---------------------------------------------------------------------------
// Language layer
// ---------------------------------------------------------------------------
#if defined(__OPENCL_VERSION__)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL FP_CONTRACT OFF
#define AMK_FN static inline
#define AMK_G __global
#define AMK_KERNEL __kernel
#define AMK_GLOBAL_ID() ((long)get_global_id(0))
typedef long amk_i64;
typedef ulong amk_u64;
typedef int amk_i32;
#define AMK_SQRT sqrt
#define AMK_EXP exp
#define AMK_LOG log
#define AMK_LOG10 log10
#define AMK_POW pow
#define AMK_COS cos
#define AMK_SIN sin
#define AMK_FABS fabs
#define AMK_FLOOR floor
#define AMK_COPYSIGN copysign
#define AMK_RINT rint
#define AMK_ISFINITE(x) isfinite(x)
#elif defined(__HIPCC__) || defined(__CUDACC__) || defined(__HIPCC_RTC__) || defined(__CUDACC_RTC__)
#define AMK_FN static __device__ inline
#define AMK_G
#define AMK_KERNEL extern "C" __global__
#define AMK_GLOBAL_ID() ((long long)blockIdx.x * (long long)blockDim.x + (long long)threadIdx.x)
typedef long long amk_i64;
typedef unsigned long long amk_u64;
typedef int amk_i32;
#define AMK_SQRT sqrt
#define AMK_EXP exp
#define AMK_LOG log
#define AMK_LOG10 log10
#define AMK_POW pow
#define AMK_COS cos
#define AMK_SIN sin
#define AMK_FABS fabs
#define AMK_FLOOR floor
#define AMK_COPYSIGN copysign
#define AMK_RINT rint
#define AMK_ISFINITE(x) isfinite(x)
#elif defined(AMK_KOKKOS)
// C++ through Kokkos (kokkos_impl.inc): host or device execution spaces.
#include <Kokkos_Core.hpp>
#include <cstdint>
#define AMK_FN static KOKKOS_INLINE_FUNCTION // static: the double and mixed TUs differ
#define AMK_G
typedef std::int64_t amk_i64;
typedef std::uint64_t amk_u64;
typedef std::int32_t amk_i32;
#define AMK_SQRT Kokkos::sqrt
#define AMK_EXP Kokkos::exp
#define AMK_LOG Kokkos::log
#define AMK_LOG10 Kokkos::log10
#define AMK_POW Kokkos::pow
#define AMK_COS Kokkos::cos
#define AMK_SIN Kokkos::sin
#define AMK_FABS Kokkos::fabs
#define AMK_FLOOR Kokkos::floor
#define AMK_COPYSIGN Kokkos::copysign
#define AMK_RINT Kokkos::nearbyint
#define AMK_ISFINITE(x) Kokkos::isfinite(x)
#else
// C++: the reference backend and the host side.
#include <cmath>
#include <cstdint>
#define AMK_FN static inline
#define AMK_G
typedef std::int64_t amk_i64;
typedef std::uint64_t amk_u64;
typedef std::int32_t amk_i32;
#define AMK_SQRT std::sqrt
#define AMK_EXP std::exp
#define AMK_LOG std::log
#define AMK_LOG10 std::log10
#define AMK_POW std::pow
#define AMK_COS std::cos
#define AMK_SIN std::sin
#define AMK_FABS std::fabs
#define AMK_FLOOR std::floor
#define AMK_COPYSIGN std::copysign
#define AMK_RINT std::nearbyint
#define AMK_ISFINITE(x) std::isfinite(x)
#endif

// ---------------------------------------------------------------------------
// Data shared by host and device (8-byte members only: one layout everywhere)
// ---------------------------------------------------------------------------

//! TrimBugs flags the trim2013 ion source paths read (KParams::bugs).
#define AMK_BUG_STALE_PRIMARY 1
#define AMK_BUG_RECOIL_ENERGY_CREATED 2
#define AMK_BUG_STRUCK_ATOM_DIRECTION 4
#define AMK_BUG_REPLACEMENT_INTERSTITIAL 8
#define AMK_BUG_SUBTHRESHOLD_SPUTTER 16
#define AMK_BUG_ION_PHONONS_TWICE 32
#define AMK_BUG_CASCADE_REPLAY 64
#define AMK_BUG_RECOIL_EVENT_TWICE 128
#define AMK_BUG_EXYZ_SKIPS_FIRST 256

//! KParams::mode (DamageCalculation).
#define AMK_MODE_KP 1
#define AMK_MODE_FULL 2
#define AMK_MODE_MONOLAYER 3

//! Status of a block (tally scalar AMK_S_STATUS).
#define AMK_OK 0
#define AMK_SLOTS_FULL 1   //!< recoil stack beyond KParams::slot_capacity
#define AMK_EXITS_FULL 2   //!< exit list beyond KParams::exit_capacity
#define AMK_NOT_A_NUMBER 3 //!< an energy or depth became NaN (the CPU throws)

/*
 * Unit conversions of the SI quantities (in_units). Double precision
 * divides by the unit's value as in_units does (bit for bit the CPU);
 * mixed precision multiplies by its reciprocal (a double division takes
 * tens of instructions at 1/32 rate on consumer GPUs).
 */
#if defined(AMK_MIXED)
#define AMK_SQRT_FAST(x) ((double)AMK_SQRT((float)(x)))
#define AMK_POW_FAST(x, y) ((double)AMK_POW((float)(x), (float)(y)))
#define AMK_A(P, x) ((x) * (P)->inv_ang)
#define AMK_EV(P, e) ((e) * (P)->inv_ev)
#else
#define AMK_SQRT_FAST(x) AMK_SQRT(x)
#define AMK_POW_FAST(x, y) AMK_POW(x, y)
#define AMK_A(P, x) ((x) / (P)->ang)
#define AMK_EV(P, e) ((e) / (P)->ev)
#endif

typedef struct {
  // Unit values of the SI quantities (amcset_units.h) and constants.
  double ang;        //!< angstrom in metres
  double ev;         //!< electron_volt in joules
  double ev_per_ang; //!< electron_volt_per_angstrom in J/m
  double pi;         //!< pi_v
  double stop_floor; //!< 0.001 * electron_volt (SELOSS's lowest energy)
  // Depth window and bins (Angstrom).
  double bin0, inv_w, w0, w1;
  double thickness;   //!< metres
  double thickness_A; //!< in_units(thickness, angstrom)
  // The ion.
  double ion_energy; //!< joules
  double ion_cx, ion_cy;
  double ion_stop, recoil_stop; //!< joules
  double first_monolayer;       //!< layer 0's monolayer, metres
  double inv_ang, inv_ev; //!< 1 / ang, 1 / ev (mixed precision)
  double exyz_interval_eV; //!< EXYZ.txt's energy step (0: off)
  double spare_d[3];
  amk_i64 nbins, n_layers, n_species, n_elements_total, n_tally;
  amk_i64 ion_z;
  amk_i64 mode, beam, bugs, straggling, integral;
  amk_i64 seed, max_recoils;
  amk_i64 record_positions, record_exits, max_exit_records;
  amk_i64 slot_capacity, exit_capacity;
  amk_i64 ion0, n_ions, particle_repeat;
  amk_i64 tally_stride, list_stride;
  amk_i64 recoil_events; //!< TRIM.DAT recoils (code 6): every particle is a cascade's primary
  amk_i64 exyz_capacity; //!< EXYZ records a block holds
  amk_i64 spare_i[4];
} KParams;

typedef struct {
  double top, bottom, density, monolayer; // SI
  double monolayer_A, bottom_A, density_A3;
  double pmax0, screening, fi, epsdg, kp_eps_per_eV, kp_k, k_low, k_bohr;
  amk_i64 first_element, n_elements;
} KLayer;

typedef struct {
  double fraction, ed, eb, es;   // joules
  double ed_eV, eb_eV;
  double arc_b, arc_c, rpa_c, arc_scale, rpa_bc, rpa_scale;
  amk_i64 z, species, tally_index, has_efficiency;
} KElement;

//! One collision pair: species s on element e of the target is pairs[s * n_elements_total + e].
typedef struct {
  double screening_length, rep, mass_ratio, gamma;
  double c[4];
  double d[4];
} KPair;

//! Stopping curve of species s in layer l: curves[s * n_layers + l].
typedef struct {
  double log_e0, inv_step, mc2_eV, spare;
  amk_i64 offset, n;
} KCurve;

//! A TRIM.DAT ion (SI).
typedef struct {
  double x, y, z, cx, cy, cz, energy;
  amk_i64 layer;         //!< layer at the starting depth
  amk_i64 element_layer; //!< recoils: the layer whose element entry describes the atom
  amk_i64 element;       //!< recoils: its element index there
} KParticle;

// ---------------------------------------------------------------------------
// Tally layout of one block (doubles); the host reads it back in this order.
// ---------------------------------------------------------------------------
#define AMK_H_ION_STOPPED 0
#define AMK_H_KNOCK_ONS 1
#define AMK_H_VAC_BY_ION 2
#define AMK_H_VAC_BY_RECOILS 3
#define AMK_H_REPLACEMENTS 4
#define AMK_H_DISPLACEMENTS 5
#define AMK_H_INTERSTITIALS 6
#define AMK_H_ION_BY_ION 7
#define AMK_H_ION_BY_RECOILS 8
#define AMK_H_PHON_BY_ION 9
#define AMK_H_PHON_BY_RECOILS 10
#define AMK_H_E2RECOIL 11
#define AMK_H_PKA_TD 12
#define AMK_H_NRT 13
#define AMK_H_KP 14
#define AMK_H_ARC 15
#define AMK_H_RPA 16
#define AMK_SINGLE_HISTOGRAMS 17
// Per element: kind * n_tally + element.
#define AMK_E_RECOIL_STOPPED 0
#define AMK_E_VACANCIES 1
#define AMK_E_VAC_REMOVED 2
#define AMK_E_ABSORBED 3
#define AMK_ELEMENT_HISTOGRAMS 4
// Scalars.
#define AMK_S_IONS 0
#define AMK_S_PKA_TD 1
#define AMK_S_NRT 2
#define AMK_S_KP 3
#define AMK_S_ARC 4
#define AMK_S_RPA 5
#define AMK_S_PKA_NO_EFF 6
#define AMK_S_DISPLACEMENTS 7
#define AMK_S_VACANCIES 8
#define AMK_S_REPLACEMENTS 9
#define AMK_S_INTERSTITIALS 10
#define AMK_S_IONIZ_ION 11
#define AMK_S_IONIZ_RECOILS 12
#define AMK_S_PHON_ION 13
#define AMK_S_PHON_RECOILS 14
#define AMK_S_E2RECOIL 15
#define AMK_S_BINDING 16
#define AMK_S_BACKSCATTERED 17
#define AMK_S_BACKSCATTERED_E 18
#define AMK_S_TRANSMITTED 19
#define AMK_S_TRANSMITTED_E 20
#define AMK_S_RECOILS_FOLLOWED 21
#define AMK_S_RECOIL_LIMIT 22
#define AMK_S_POSITIONS 23
#define AMK_S_EXITS 24
#define AMK_S_STATUS 25
#define AMK_S_STATUS_ION 26
#define AMK_S_MAX_K 27 // deepest recoil slot the block used
#define AMK_S_STEPS 28 // steps of the step machine the block took (profiling)
#define AMK_S_EXYZ 29  // EXYZ records of the block
#define AMK_SCALARS 32
// Moments (n, s1..s4 each): range, lateral_y, lateral_abs_y, radial.
#define AMK_M_RANGE 0
#define AMK_M_LAT_Y 5
#define AMK_M_LAT_ABS_Y 10
#define AMK_M_RADIAL 15
#define AMK_MOMENTS 20
// Per element arrays: kind * n_tally + element.
#define AMK_A_SPUTTERED 0
#define AMK_A_SPUTTERED_E 1
#define AMK_A_TRANSMITTED_RECOILS 2
#define AMK_ELEMENT_ARRAYS 3
#define AMK_SPECTRUM_BINS 100   // TrimTallies::spectrum_bins_per_decade * spectrum_decades
#define AMK_SURFACE_BINS 200    // TrimTallies::surface_energy_bins
#define AMK_POSITION_DOUBLES 4  // ion (1-based), x, y, z
#define AMK_EXIT_DOUBLES 10     // kind, ion, z, energy, x, y, z, cx, cy, cz
#define AMK_EXYZ_DOUBLES 7      // ion, energy (keV), x, y, z, Se (eV/A), last T (eV)
#define AMK_BLOCK_IONS 16       // TrimSimulation::run's block

AMK_FN amk_i64 amk_off_single(amk_i64 nbins, amk_i64 h) { return h * nbins; }
AMK_FN amk_i64 amk_off_element(amk_i64 nbins, amk_i64 ne, amk_i64 kind, amk_i64 e) {
  return (AMK_SINGLE_HISTOGRAMS + kind * ne + e) * nbins;
}
AMK_FN amk_i64 amk_off_lateral(amk_i64 nbins, amk_i64 ne) {
  return (AMK_SINGLE_HISTOGRAMS + AMK_ELEMENT_HISTOGRAMS * ne) * nbins;
}
AMK_FN amk_i64 amk_off_scalars(amk_i64 nbins, amk_i64 ne) { return amk_off_lateral(nbins, ne) + 5 * nbins; }
AMK_FN amk_i64 amk_off_moments(amk_i64 nbins, amk_i64 ne) { return amk_off_scalars(nbins, ne) + AMK_SCALARS; }
AMK_FN amk_i64 amk_off_arrays(amk_i64 nbins, amk_i64 ne) { return amk_off_moments(nbins, ne) + AMK_MOMENTS; }
AMK_FN amk_i64 amk_off_spectra(amk_i64 nbins, amk_i64 ne) {
  return amk_off_arrays(nbins, ne) + AMK_ELEMENT_ARRAYS * ne;
}
AMK_FN amk_i64 amk_off_surface(amk_i64 nbins, amk_i64 ne) { return amk_off_spectra(nbins, ne) + 2 * AMK_SPECTRUM_BINS; }
AMK_FN amk_i64 amk_tally_size(amk_i64 nbins, amk_i64 ne) { return amk_off_surface(nbins, ne) + AMK_SURFACE_BINS * ne; }
AMK_FN amk_i64 amk_list_size(amk_i64 exit_capacity, amk_i64 exyz_capacity) {
  return AMK_BLOCK_IONS * AMK_POSITION_DOUBLES + exit_capacity * AMK_EXIT_DOUBLES + exyz_capacity * AMK_EXYZ_DOUBLES;
}

// ---------------------------------------------------------------------------
// Device-side state
// ---------------------------------------------------------------------------

//! Transporter's Mover, the members the trim2013 tallies need.
typedef struct {
  double x, y, zz, cx, cy, cz, energy;
  amk_i32 species, z, tally_element, layer_hint, element_hint, origin_layer;
  amk_i32 primary; //!< a TRIM.DAT recoil event's own atom (the CPU's track 1), carried by the moving copy
  amk_i32 pad;
} KMover;

//! LecuyerRandom (trim/random.h).
typedef struct {
  amk_i64 s1, s2;
} KRng;

//! What a lane is doing (KLane::phase).
#define AMK_PHASE_NEXT_ION 0      //!< start the block's next ion, or the next block
#define AMK_PHASE_ION 1           //!< one flight and collision of the ion
#define AMK_PHASE_CASCADE 2       //!< one step of the recoil stack (RECOILLOOP / NEXTRECOIL)
#define AMK_PHASE_AFTER_CASCADE 3 //!< the ion's loop after its cascade
#define AMK_PHASE_DONE 4          //!< no block left

/*!
 * A lane's state between kernel launches (one per work item, in device
 * memory). A lane runs blocks one after the other, each from its first to
 * its last ion, in steps; what the CPU's Transporter keeps between the ions
 * of a block, and the ion's loop state across its cascades, are here.
 */
typedef struct {
  KMover ion;
  KRng rng;   //!< the ion's stream
  KRng crng;  //!< the cascade's stream (a copy of the ion's, or its own)
  amk_i64 phase, block, ion_in_block, ion_index;
  amk_i64 k, pop, cascade_on_ion_rng;
  amk_i64 status, flight_counter, primary_displaced, sub_threshold;
  amk_i64 cascades_this_ion, recoils_this_ion, max_k;
  amk_i64 spare[2]; // [0]: steps; [1]: TRIM.DAT recoil event flags
  double primary_stop_x;
  double last_T_eV;
  amk_i64 exyz_next;
  amk_i64 spare2;
} KLane;

//! The data of one work item.
typedef struct {
  AMK_G const KParams *P;
  AMK_G const KLayer *layers;
  AMK_G const KElement *elements;
  AMK_G const KPair *pairs;
  AMK_G const KCurve *curves;
  AMK_G const double *log_s;
  AMK_G const double *edges;
  AMK_G const KParticle *particles;
  AMK_G KMover *slots;
  AMK_G double *T;     // this block's tallies
  AMK_G double *lists; // this block's positions and exits
  amk_i64 status;
  // Transporter members that persist between the ions of a block.
  amk_i64 flight_counter;
  amk_i64 primary_displaced;
  amk_i64 sub_threshold;
  amk_i64 cascades_this_ion;
  amk_i64 recoils_this_ion;
  amk_i64 ion_index;
  // The step machine (KLane).
  KMover mv;
  KRng rng, crng;
  amk_i64 phase, block, ion_in_block, k, pop, cascade_on_ion_rng, max_k, steps;
  amk_i64 ns; // ring slots: the tallies' stride (element e of this block at T[e * ns])
  amk_i64 ring, first_block;
  amk_i64 track_primary, has_primary_stop; // TRIM.DAT recoil events
  double primary_stop_x;
  amk_i64 exyz_next;  // TRIM's G360: the energy step of the last EXYZ record
  double last_T_eV;  // TRIM's G5DC: the energy given to the last struck atom
} KCtx;

// ---------------------------------------------------------------------------
// std:: semantics of the CPU transport
// ---------------------------------------------------------------------------
AMK_FN double amk_max(double a, double b) { return (a < b) ? b : a; }  // std::max
AMK_FN double amk_min(double a, double b) { return (b < a) ? b : a; }  // std::min
AMK_FN double amk_clamp(double v, double lo, double hi) { return (v < lo) ? lo : (hi < v) ? hi : v; }
AMK_FN amk_i64 amk_clamp_i(amk_i64 v, amk_i64 lo, amk_i64 hi) { return (v < lo) ? lo : (hi < v) ? hi : v; }

// ---------------------------------------------------------------------------
// Random numbers (trim/random.h)
// ---------------------------------------------------------------------------
AMK_FN amk_i64 amk_mul_mod(amk_i64 a, amk_i64 s, amk_i64 m) {
  const amk_i64 x = a * s;
  const amk_i64 c = (((amk_i64)1) << 31) - m;
  amk_i64 r = (x & 0x7fffffff) + c * (x >> 31);
  if (r >= m) r -= m;
  return r;
}
AMK_FN double amk_rand(KRng *g) {
  g->s1 = amk_mul_mod(40014, g->s1, 2147483563);
  g->s2 = amk_mul_mod(40692, g->s2, 2147483399);
  amk_i64 z = g->s1 - g->s2;
  if (z < 1) z += 2147483563 - 1;
  return (double)z * (1.0 / (double)2147483563);
}
AMK_FN amk_u64 amk_splitmix64(amk_u64 x) {
  x += 0x9E3779B97F4A7C15ULL;
  x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31);
}
AMK_FN KRng amk_make_stream(amk_u64 seed, amk_u64 ion, amk_u64 stream) {
  const amk_u64 base = amk_splitmix64((seed * 0x100000001B3ULL) ^ amk_splitmix64(ion * 4 + stream));
  const amk_u64 a = amk_splitmix64(base);
  const amk_u64 b = amk_splitmix64(base ^ 0xD1B54A32D192ED03ULL);
  KRng g;
  g.s1 = (amk_i64)(a % (amk_u64)(2147483563 - 1)) + 1;
  g.s2 = (amk_i64)(b % (amk_u64)(2147483399 - 1)) + 1;
  return g;
}

// ---------------------------------------------------------------------------
// Tallies (Transporter::bin_of, add, deposit, spectrum_bin)
// ---------------------------------------------------------------------------
/*!
 * Tally element e of the lane's block. The blocks of a batch are
 * interleaved (element e of block j at e * n_blocks + j), so that the lanes
 * of a wavefront, which run neighbouring blocks, update neighbouring
 * addresses: their writes to the same quantity coalesce.
 */
AMK_FN AMK_G double *amk_t(KCtx *c, amk_i64 e) { return c->T + e * c->ns; }

AMK_FN amk_i64 amk_bin_of(AMK_G const KParams *P, double x) {
  const double f = (x - P->bin0) * P->inv_w;
  if (!(f >= 0.0) || f >= (double)P->nbins) return -1;
  return (amk_i64)f;
}
AMK_FN void amk_add(KCtx *c, amk_i64 off, double x, double v) {
#if defined(AMK_PROF_NO_TALLY) // profiling only: no depth-histogram writes
  return;
#endif
  const amk_i64 b = amk_bin_of(c->P, x);
  if (b >= 0) *amk_t(c, off + b) += v;
}
AMK_FN void amk_hadd(KCtx *c, amk_i64 h, double x, double v) { amk_add(c, amk_off_single(c->P->nbins, h), x, v); }
AMK_FN void amk_eadd(KCtx *c, amk_i64 kind, amk_i64 e, double x, double v) {
  amk_add(c, amk_off_element(c->P->nbins, c->P->n_tally, kind, e), x, v);
}
AMK_FN AMK_G double *amk_scalar(KCtx *c, amk_i64 s) {
  return amk_t(c, amk_off_scalars(c->P->nbins, c->P->n_tally) + s);
}
AMK_FN AMK_G double *amk_array(KCtx *c, amk_i64 kind, amk_i64 e) {
  return amk_t(c, amk_off_arrays(c->P->nbins, c->P->n_tally) + kind * c->P->n_tally + e);
}
AMK_FN void amk_deposit(KCtx *c, amk_i64 h, double x0, double x1, double amount) {
  AMK_G const KParams *P = c->P;
#if defined(AMK_PROF_NO_TALLY)
  return;
#endif
  if (amount == 0.0) return;
  if (x1 < x0) {
    const double t = x0;
    x0 = x1;
    x1 = t;
  }
  const double len = x1 - x0;
  if (len <= 1e-12) {
    amk_add(c, h, x0, amount);
    return;
  }
  const double lo = amk_max(x0, P->bin0);
  const double hi = amk_min(x1, c->edges[P->nbins]);
  if (!(hi > lo)) return;
  amk_i64 b0 = (amk_i64)((lo - P->bin0) * P->inv_w);
  amk_i64 b1 = (amk_i64)((hi - P->bin0) * P->inv_w);
  b0 = amk_clamp_i(b0, 0, P->nbins - 1);
  b1 = amk_clamp_i(b1, 0, P->nbins - 1);
  for (amk_i64 b = b0; b <= b1; ++b) {
    const double e0 = amk_max(lo, c->edges[b]);
    const double e1 = amk_min(hi, c->edges[b + 1]);
    if (e1 > e0) *amk_t(c, h + b) += amount * (e1 - e0) / len;
  }
}
AMK_FN amk_i64 amk_spectrum_bin(double er) {
#if defined(AMK_MIXED)
  const double l = (double)AMK_LOG10((float)er);
#else
  const double l = AMK_LOG10(er);
#endif
  const double f = AMK_FLOOR(l * 10);
  return f >= 0.0 && f < 1e6 ? (amk_i64)f : -1;
}
AMK_FN int amk_in_window(AMK_G const KParams *P, double x) { return x >= P->w0 && x <= P->w1; }

// ---------------------------------------------------------------------------
// Physics (physics/scattering.cpp, physics/stopping.cpp)
// ---------------------------------------------------------------------------
AMK_FN void amk_phi_dphi(AMK_G const KPair *f, double x, double *phi_x, double *dphi_x) {
  const double e0 = AMK_EXP(-f->d[0] * x);
  const double e1 = AMK_EXP(-f->d[1] * x);
  const double e2 = AMK_EXP(-f->d[2] * x);
  const double e3 = AMK_EXP(-f->d[3] * x);
  *phi_x = f->c[0] * e0 + f->c[1] * e1 + f->c[2] * e2 + f->c[3] * e3;
  *dphi_x = -(f->c[0] * f->d[0] * e0 + f->c[1] * f->d[1] * e1 + f->c[2] * f->d[2] * e2 + f->c[3] * f->d[3] * e3);
}
AMK_FN double amk_phi(AMK_G const KPair *f, double x) {
  return f->c[0] * AMK_EXP(-f->d[0] * x) + f->c[1] * AMK_EXP(-f->d[1] * x) + f->c[2] * AMK_EXP(-f->d[2] * x) +
         f->c[3] * AMK_EXP(-f->d[3] * x);
}
AMK_FN double amk_closest_approach(AMK_G const KPair *f, double eps, double b) {
  double r = b;
  if (eps * b > 0.0) {
    double rr = -2.7 * AMK_LOG(eps * b);
    if (rr >= b) {
      rr = -2.7 * AMK_LOG(eps * rr);
      if (rr >= b) r = rr;
    }
  }
  if (!(r > 0.0)) r = amk_max(b, 1e-6);
  double lo = 0.0;
  double hi = amk_max(amk_max(2.0 / eps, 2.0 * b), 1e-6);
  for (int iter = 0; iter < 100; ++iter) {
    double ph = 0.0;
    double dph = 0.0;
    amk_phi_dphi(f, r, &ph, &dph);
    const double F = r - ph / eps - b * b / r;
    const double dF = 1.0 - dph / eps + b * b / (r * r);
    if (F < 0.0) lo = amk_max(lo, r);
    else hi = amk_min(hi, r);
    double next = r - F / dF;
    if (!(next > lo) || !(next < hi) || !AMK_ISFINITE(next)) next = 0.5 * (lo + hi);
    if (AMK_FABS(next - r) <= 1e-9 * next) return next;
    r = next;
  }
  return r;
}
AMK_FN double amk_magic_sin2(AMK_G const KPair *f, double eps, double b) {
  const double C1 = 0.99229;
  const double C2 = 0.011615;
  const double C3 = 0.007122;
  const double C4 = 9.3066;
  const double C5 = 14.813;
  if (b <= 0.0) return 1.0;
  const double r0 = amk_closest_approach(f, eps, b);
  double ph0 = 0.0;
  double dph0 = 0.0;
  amk_phi_dphi(f, r0, &ph0, &dph0);
  const double v = ph0 / r0;
  const double v1 = dph0 / r0 - v / r0;
  const double rho = -2.0 * (eps - v) / v1;
  if (!(AMK_FABS(rho) <= 1.7976931348623157e308)) return 0.0;
  const double sqe = AMK_SQRT(eps);
  const double alpha = 1.0 + C1 / sqe;
  const double beta = (C2 + sqe) / (C3 + sqe);
  const double gam = (C4 + eps) / (C5 + eps);
  const double A = 2.0 * alpha * eps * AMK_POW_FAST(b, beta);
  const double G = 1.0 / (gam * (AMK_SQRT(1.0 + A * A) - A));
  const double delta = A * (r0 - b) / (1.0 + G);
  const double cos_half = (b + rho + delta) / (r0 + rho);
  const double s2 = 1.0 - cos_half * cos_half;
  return amk_clamp(s2, 0.0, 1.0);
}
AMK_FN double amk_integral_sin2(AMK_G const KPair *f, double pi, double eps, double b) {
  if (b <= 0.0) return 1.0;
  const double r0 = amk_closest_approach(f, eps, b);
  const int n = 100; // integral_sin2_half_theta's default nodes (even, >= 4)
  double sum = 0.0;
  for (int j = 1; j <= n / 2; ++j) {
    const double u = AMK_COS((2.0 * j - 1.0) * pi / (2.0 * n));
    const double r = r0 / u;
    const double g = 1.0 - amk_phi(f, r) / (eps * r) - b * b / (r * r);
    sum += AMK_SQRT((1.0 - u * u) / amk_max(g, 1e-300));
  }
  const double theta = pi * (1.0 - 2.0 * b / (n * r0) * sum);
  const double s = AMK_SIN(0.5 * theta);
  return amk_clamp(s * s, 0.0, 1.0);
}
#if defined(AMK_MIXED)
/*
 * Mixed precision (compute.precision = mixed): the magic formula in single
 * precision, the rest of the transport in double. Consumer GPUs run single
 * precision 16 to 64 times faster than double. The closest approach is
 * converged in single precision and finished with one Newton step in
 * double, and 1 - cos^2(theta/2) is formed in double: in single precision
 * its cancellation biased small-angle collisions (docs/GPU.md "Precision").
 */
AMK_FN void amk_phi_dphi_f(AMK_G const KPair *f, float x, float *phi_x, float *dphi_x) {
  const float c0 = (float)f->c[0], c1 = (float)f->c[1], c2 = (float)f->c[2], c3 = (float)f->c[3];
  const float d0 = (float)f->d[0], d1 = (float)f->d[1], d2 = (float)f->d[2], d3 = (float)f->d[3];
  const float e0 = AMK_EXP(-d0 * x);
  const float e1 = AMK_EXP(-d1 * x);
  const float e2 = AMK_EXP(-d2 * x);
  const float e3 = AMK_EXP(-d3 * x);
  *phi_x = c0 * e0 + c1 * e1 + c2 * e2 + c3 * e3;
  *dphi_x = -(c0 * d0 * e0 + c1 * d1 * e1 + c2 * d2 * e2 + c3 * d3 * e3);
}
AMK_FN double amk_magic_sin2_f(AMK_G const KPair *f, double eps_d, double b_d) {
  if (b_d <= 0.0) return 1.0;
  const float eps = (float)eps_d;
  const float b = (float)b_d;
  float r = b;
  if (eps * b > 0.0f) {
    float rr = -2.7f * AMK_LOG(eps * b);
    if (rr >= b) {
      rr = -2.7f * AMK_LOG(eps * rr);
      if (rr >= b) r = rr;
    }
  }
  if (!(r > 0.0f)) r = b > 1e-6f ? b : 1e-6f;
  float lo = 0.0f;
  float hi = 2.0f / eps > 2.0f * b ? 2.0f / eps : 2.0f * b;
  if (hi < 1e-6f) hi = 1e-6f;
  for (int iter = 0; iter < 60; ++iter) {
    float ph = 0.0f, dph = 0.0f;
    amk_phi_dphi_f(f, r, &ph, &dph);
    const float F = r - ph / eps - b * b / r;
    const float dF = 1.0f - dph / eps + b * b / (r * r);
    if (F < 0.0f) lo = lo > r ? lo : r;
    else hi = hi < r ? hi : r;
    float next = r - F / dF;
    if (!(next > lo) || !(next < hi) || !AMK_ISFINITE(next)) next = 0.5f * (lo + hi);
    const float dr = next - r;
    r = next;
    if (AMK_FABS(dr) <= 4e-7f * next) break;
  }
  // sin^2(theta/2) = 1 - cos^2(theta/2) = (1 - cos)(1 + cos) with
  // 1 - cos = (r0 - b - delta)/(r0 + rho): for distant collisions r0 - b is
  // a small difference and needs r0 in double. One Newton step in double from
  // the single-precision root (good to 4e-7) makes r0 good to about 1e-13.
  // rho enters only through r0 + rho, so the potential of that step, at the
  // single-precision root, serves for it (4e-7 relative). A and G (fit
  // constants) enter through delta only and stay single precision.
  const double r1 = (double)r;
  double ph0 = 0.0;
  double dph0 = 0.0;
  amk_phi_dphi(f, r1, &ph0, &dph0);
  double r0 = r1;
  {
    const double F = r1 - ph0 / eps_d - b_d * b_d / r1;
    const double dF = 1.0 - dph0 / eps_d + b_d * b_d / (r1 * r1);
    const double next = r1 - F / dF;
    if (next > 0.0 && AMK_ISFINITE(next) && AMK_FABS(next - r1) <= 1e-3 * r1) r0 = next;
  }
  const double v = ph0 / r1;
  const double v1 = dph0 / r1 - v / r1;
  const double rho = -2.0 * (eps_d - v) / v1;
  if (!(AMK_FABS(rho) <= 1.7976931348623157e308)) return 0.0;
  const float sqe = AMK_SQRT(eps);
  const float alpha = 1.0f + 0.99229f / sqe;
  const float beta = (0.011615f + sqe) / (0.007122f + sqe);
  const float gam = (9.3066f + eps) / (14.813f + eps);
  const float A = 2.0f * alpha * eps * AMK_POW(b, beta);
  const float G = 1.0f / (gam * (AMK_SQRT(1.0f + A * A) - A));
  const double delta = (double)A * (r0 - b_d) / (1.0 + (double)G);
  const double cos_half = (b_d + rho + delta) / (r0 + rho);
  return amk_clamp(1.0 - cos_half * cos_half, 0.0, 1.0);
}
#endif

AMK_FN double amk_sin2_half_theta(AMK_G const KParams *P, AMK_G const KPair *f, double eps, double b) {
  if (eps > 10.0) {
    const double x = 2.0 * eps * b;
    return 1.0 / (1.0 + (1.0 + b * (1.0 + b)) * x * x);
  }
  if (P->integral) return amk_integral_sin2(f, P->pi, eps, b);
#if defined(AMK_PROF_CHEAP_SCATTER) // profiling only: screened Rutherford for every collision
  {
    const double x = 2.0 * eps * b;
    return 1.0 / (1.0 + (1.0 + b * (1.0 + b)) * x * x);
  }
#endif
#if defined(AMK_MIXED)
  return amk_magic_sin2_f(f, eps, b);
#else
  return amk_magic_sin2(f, eps, b);
#endif
}

typedef struct {
  double cos_projectile, sin_projectile, cos_recoil, sin_recoil, transferred, s2;
} KKin;

//! collision_kinematics(pair, energy, trim_sin2_half_theta(...)).
AMK_FN KKin amk_scatter(AMK_G const KParams *P, AMK_G const KPair *pair, double energy, double eps, double b) {
  KKin k;
  const double s2 = amk_clamp(amk_sin2_half_theta(P, pair, eps, b), 0.0, 1.0);
  const double c2 = 1.0 - s2;
  k.s2 = s2;
#if defined(AMK_MIXED)
  {
    // Single-precision kinematics (the energy transfer stays double).
    const float s2f = (float)s2, c2f = (float)c2;
    const float cos_theta_cm = (float)(1.0 - 2.0 * s2);
    const float sin_theta = 2.0f * AMK_SQRT(s2f * c2f);
    k.transferred = pair->gamma * s2 * energy;
    const float mu = (float)pair->mass_ratio;
    const float den = AMK_SQRT(1.0f + 2.0f * mu * cos_theta_cm + mu * mu);
    if (den < 1e-6f) {
      k.cos_projectile = 0.0;
      k.sin_projectile = 1.0;
    } else {
      const float inv = 1.0f / den;
      k.cos_projectile = (double)((cos_theta_cm + mu) * inv);
      k.sin_projectile = (double)(sin_theta * inv);
    }
    k.cos_recoil = (double)AMK_SQRT(s2f);
    k.sin_recoil = (double)AMK_SQRT(c2f);
    return k;
  }
#endif
  const double cos_theta_cm = 1.0 - 2.0 * s2;
  const double sin_theta = 2.0 * AMK_SQRT(s2 * c2);
  k.transferred = pair->gamma * s2 * energy;
  const double mu = pair->mass_ratio;
  const double den = AMK_SQRT(1.0 + 2.0 * mu * cos_theta_cm + mu * mu);
  if (den < 1e-12) {
    k.cos_projectile = 0.0;
    k.sin_projectile = 1.0;
  } else {
    k.cos_projectile = (cos_theta_cm + mu) / den;
    k.sin_projectile = sin_theta / den;
  }
  k.cos_recoil = AMK_SQRT(s2);
  k.sin_recoil = AMK_SQRT(c2);
  return k;
}

//! StoppingCurve::stopping_power (SI, J/m).
AMK_FN double amk_stopping(KCtx *c, amk_i64 species, amk_i64 layer, double energy) {
  AMK_G const KCurve *cu = c->curves + species * c->P->n_layers + layer;
  if (cu->n == 0 || energy <= 0.0 * c->P->ev) return 0.0 * c->P->ev_per_ang;
  AMK_G const double *ls = c->log_s + cu->offset;
#if defined(AMK_PROF_CHEAP_STOPPING) // profiling only: no logarithm, exponential or table
  return ls[cu->n / 2] * 0.0 + 30.0 * c->P->ev_per_ang * AMK_SQRT_FAST(energy * 1e16);
#endif
#if defined(AMK_MIXED)
  {
    // Single precision: ln E and exp of the interpolated ln S (1e-7 relative).
    const float xf = ((float)AMK_LOG((float)(AMK_EV(c->P, energy))) - (float)cu->log_e0) * (float)cu->inv_step;
    float vf = 0.0f;
    if (xf <= 0.0f) {
      vf = (float)ls[0] + 0.5f * xf / (float)cu->inv_step;
    } else {
      const float lastf = (float)(cu->n - 1);
      if (xf >= lastf) {
        vf = (float)ls[cu->n - 1] + (xf - lastf) * (float)(ls[cu->n - 1] - ls[cu->n - 2]);
      } else {
        const amk_i64 i = (amk_i64)xf;
        const float t = xf - (float)i;
        vf = (float)ls[i] + t * (float)(ls[i + 1] - ls[i]);
      }
    }
    return (double)AMK_EXP(vf) * c->P->ev_per_ang;
  }
#endif
  const double x = (AMK_LOG(AMK_EV(c->P, energy)) - cu->log_e0) * cu->inv_step;
  double value = 0.0;
  if (x <= 0.0) {
    value = ls[0] + 0.5 * x / cu->inv_step;
  } else {
    const double last = (double)(cu->n - 1);
    if (x >= last) {
      value = ls[cu->n - 1] + (x - last) * (ls[cu->n - 1] - ls[cu->n - 2]);
    } else {
      const amk_i64 i = (amk_i64)x;
      const double t = x - (double)i;
      value = ls[i] + t * (ls[i + 1] - ls[i]);
    }
  }
  return AMK_EXP(value) * c->P->ev_per_ang;
}

// Direction cosines (rotate, rotate_about_x_2013).
#if defined(AMK_MIXED)
// Mixed precision: the rotations in single precision (one reciprocal square
// root instead of a double square root and five double divisions).
AMK_FN void amk_rotate(double *cx, double *cy, double *cz, double cos_t, double sin_t, double cos_phi,
                       double sin_phi) {
  const float x = (float)*cx, y = (float)*cy, z = (float)*cz;
  const float ct = (float)cos_t, st = (float)sin_t, cp = (float)cos_phi, sp = (float)sin_phi;
  float px = 0.0f, py = cp, pz = sp;
  // sqrt(cy^2 + cz^2), not sqrt(1 - cx^2): near the depth axis cx rounds to
  // steps of 6e-8 in single precision, which would quantize s to 3e-4.
  const float s2 = y * y + z * z;
  if (s2 > 1e-20f) {
    const float s = AMK_SQRT(s2);
    const float inv_s = 1.0f / s;
    px = -s * cp;
    py = (x * y * cp - z * sp) * inv_s;
    pz = (x * z * cp + y * sp) * inv_s;
  }
  const float nx = x * ct + px * st;
  const float ny = y * ct + py * st;
  const float nz = z * ct + pz * st;
  const float inv = 1.0f / AMK_SQRT(nx * nx + ny * ny + nz * nz);
  *cx = (double)(nx * inv);
  *cy = (double)(ny * inv);
  *cz = (double)(nz * inv);
}
AMK_FN void amk_rotate_about_x(double *cx, double *cy, double *cz, double cos_t, double sin_t, double cos_phi,
                               double sin_phi) {
  const float nx = (float)*cx * (float)cos_t;
  const float ny = -(float)sin_t * (float)cos_phi;
  const float nz = -(float)sin_t * (float)sin_phi;
  const float n2 = nx * nx + ny * ny + nz * nz;
  if (n2 == 0.0f) return;
  const float inv = 1.0f / AMK_SQRT(n2);
  *cx = (double)(nx * inv);
  *cy = (double)(ny * inv);
  *cz = (double)(nz * inv);
}
AMK_FN void amk_azimuth(KCtx *c, KRng *rng, double *cos_phi, double *sin_phi) {
  const float phi = 6.2831853f * (float)amk_rand(rng);
  *cos_phi = (double)AMK_COS(phi);
  *sin_phi = (double)AMK_SIN(phi);
  (void)c;
}
#else
AMK_FN void amk_rotate(double *cx, double *cy, double *cz, double cos_t, double sin_t, double cos_phi,
                       double sin_phi) {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  const double s = AMK_SQRT(amk_max(0.0, 1.0 - *cx * *cx));
  if (s < 1e-10) {
    px = 0.0;
    py = cos_phi;
    pz = sin_phi;
  } else {
    px = -s * cos_phi;
    py = (*cx * *cy * cos_phi - *cz * sin_phi) / s;
    pz = (*cx * *cz * cos_phi + *cy * sin_phi) / s;
  }
  const double nx = *cx * cos_t + px * sin_t;
  const double ny = *cy * cos_t + py * sin_t;
  const double nz = *cz * cos_t + pz * sin_t;
  const double norm = AMK_SQRT(nx * nx + ny * ny + nz * nz);
  *cx = nx / norm;
  *cy = ny / norm;
  *cz = nz / norm;
}
AMK_FN void amk_rotate_about_x(double *cx, double *cy, double *cz, double cos_t, double sin_t, double cos_phi,
                               double sin_phi) {
  const double nx = *cx * cos_t;
  const double ny = -sin_t * cos_phi;
  const double nz = -sin_t * sin_phi;
  const double norm = AMK_SQRT(nx * nx + ny * ny + nz * nz);
  if (norm == 0.0) return;
  *cx = nx / norm;
  *cy = ny / norm;
  *cz = nz / norm;
}
AMK_FN void amk_azimuth(KCtx *c, KRng *rng, double *cos_phi, double *sin_phi) {
  const double u = amk_rand(rng);
  const double phi = 2.0 * c->P->pi * u;
  *cos_phi = AMK_COS(phi);
  *sin_phi = AMK_SIN(phi);
}
#endif

// ---------------------------------------------------------------------------
// The trim2013 transport (Transporter::*_2013)
// ---------------------------------------------------------------------------
AMK_FN amk_i64 amk_layer_index(KCtx *c, double x) {
  for (amk_i64 i = 0; i < c->P->n_layers; ++i) {
    if (x < c->layers[i].bottom) return i;
  }
  return c->P->n_layers - 1;
}

AMK_FN amk_i64 amk_pick_element(AMK_G const KLayer *layer, AMK_G const KElement *els, KRng *rng) {
  double r = amk_rand(rng);
  if (layer->n_elements == 1) return 0;
  for (amk_i64 j = 0; j < layer->n_elements; ++j) {
    r -= els[layer->first_element + j].fraction;
    if (r < 0.0) return j;
  }
  return layer->n_elements - 1;
}

//! electronic_loss_2013 (\p known_se >= 0: the stopping power at the mover's energy).
AMK_FN double amk_electronic_loss(KCtx *c, amk_i64 species, double energy, amk_i64 li, double L, KRng *rng,
                                  int has_se, double known_se) {
  AMK_G const KParams *P = c->P;
  const double e = amk_max(energy, P->stop_floor);
  const double mean = (has_se && e == energy ? known_se : amk_stopping(c, species, li, e)) * L;
  double de = mean;
  const double r = amk_rand(rng);
  if (P->straggling) {
    AMK_G const KLayer *t = c->layers + li;
    const double e_eV = AMK_EV(P, e);
    double omega = t->k_low * e_eV;
    if (omega > t->k_bohr) {
      const double mc2_eV = c->curves[species * P->n_layers + li].mc2_eV;
      const double g = 1.0 + e_eV / mc2_eV;
      const double beta2 = 1.0 - 1.0 / (g * g);
      const double inv = 1.0 / beta2;
      omega = t->k_bohr * (inv - 0.5) / (inv - 1.0);
    }
#if defined(AMK_MIXED)
    double strag = (r - 0.5) * omega * (double)AMK_SQRT((float)AMK_A(P, L));
#else
    double strag = (r - 0.5) * omega * AMK_SQRT(AMK_A(P, L));
#endif
    const double m_eV = AMK_EV(P, mean);
    if (AMK_FABS(strag) > m_eV) strag = AMK_COPYSIGN(m_eV, strag);
    de = mean + strag * P->ev;
  }
  if (de > e) de = e;
  return de;
}

AMK_FN double amk_edge(KCtx *c, const KMover *mv, amk_i64 li, double L, KRng *rng) {
  AMK_G const KParams *P = c->P;
  AMK_G const KLayer *ly = c->layers + li;
  const double x = AMK_A(P, mv->x);
  const double l0 = ly->monolayer_A;
  const double binw = (P->w1 - P->w0) / (double)P->nbins;
  double D = 0.0;
  if (mv->cx >= 0.0) {
    D = P->w1 - x;
    if (D < 0.0) D = ly->bottom_A - x;
  } else {
    D = x - P->w0;
  }
  const int quick = P->mode == AMK_MODE_KP;
  double l = AMK_A(P, L);
  while (AMK_FABS(l * mv->cx) > D) {
    l *= amk_rand(rng);
    if (quick) {
      if (AMK_FABS(l * mv->cx) < binw) {
        l = AMK_FABS(binw / mv->cx);
        break;
      }
    } else if (AMK_FABS(l * mv->cx) < l0) {
      l = l0;
      break;
    }
  }
  if (mv->cx > 0.0 && l * mv->cx > D) l = (D + ly->density_A3) / mv->cx;
  c->flight_counter = 2;
  return l * P->ang;
}

AMK_FN double amk_interface(KCtx *c, const KMover *mv, amk_i64 li, double L, KRng *rng) {
  AMK_G const KParams *P = c->P;
  const double x = AMK_A(P, mv->x);
  const double l0 = c->layers[li].monolayer_A;
  const double binw = (P->w1 - P->w0) / (double)P->nbins;
  double D = 0.0;
  if (mv->cx >= 0.0) {
    D = c->layers[li].bottom_A - x;
  } else {
    D = li == 0 ? x : x - c->layers[li - 1].bottom_A;
  }
  if (D < l0) D = l0;
  double l = AMK_A(P, L);
  while (AMK_FABS(l * mv->cx) > D) {
    l *= amk_rand(rng);
    const double g = AMK_FABS(l * mv->cx);
    if (binw > 3.0 * l0) {
      if (g < binw && x > 10.0 * l0) {
        l = amk_min(binw, D / amk_max(AMK_FABS(mv->cx), 0.01));
        break;
      }
    } else if (g < l0) {
      l = l0;
      break;
    }
  }
  c->flight_counter = 2;
  return l * P->ang;
}

AMK_FN double amk_flight_path(KCtx *c, const KMover *mv, amk_i64 li, KRng *rng) {
  AMK_G const KParams *P = c->P;
  if (mv->x < P->first_monolayer && mv->energy == P->ion_energy) {
    c->flight_counter = 2;
    return P->first_monolayer;
  }
  AMK_G const KLayer *layer = c->layers + li;
  if (P->mode == AMK_MODE_MONOLAYER) return layer->monolayer;
  const double eps = mv->energy * layer->fi;
  const double eeg = AMK_SQRT_FAST(eps * layer->epsdg);
  const double pmax = layer->screening / (eeg + AMK_SQRT_FAST(eeg) + 0.125 * AMK_POW_FAST(eeg, 0.1));
  double L = 1.0 / (P->pi * pmax * pmax * layer->density);
  const double xn = mv->x + L * mv->cx;
  if (amk_layer_index(c, xn) != li) L = amk_interface(c, mv, li, L, rng);
  const double xna = AMK_A(P, xn);
  if (xna < P->w0 || xna > P->w1) L = amk_edge(c, mv, li, L, rng);
  if (mv->x < pmax && mv->cx < 0.0) L = amk_edge(c, mv, li, L, rng);
  if (P->mode == AMK_MODE_FULL && P->beam) return amk_max(layer->monolayer, amk_rand(rng) * L);
  if (c->flight_counter > 0) {
    L = amk_max(layer->monolayer, amk_rand(rng) * L);
    --c->flight_counter;
  }
  return L;
}

AMK_FN void amk_book_pka(KCtx *c, AMK_G const KElement *el, double x, double td_eV) {
  const double ed = el->ed_eV;
  if (!(td_eV > 0.0) || !(ed > 0.0)) return;
  const double knee = 2.0 * ed / 0.8;
  const double nrt = td_eV < ed ? 0.0 : td_eV < knee ? 1.0 : 0.8 * td_eV / (2.0 * ed);
  const double kp = td_eV < ed ? 0.0 : td_eV < 2.0 * ed ? 1.0 : td_eV / (2.0 * ed);
  amk_hadd(c, AMK_H_PKA_TD, x, td_eV);
  amk_hadd(c, AMK_H_NRT, x, nrt);
  amk_hadd(c, AMK_H_KP, x, kp);
  *amk_scalar(c, AMK_S_PKA_TD) += td_eV;
  *amk_scalar(c, AMK_S_NRT) += nrt;
  *amk_scalar(c, AMK_S_KP) += kp;
  if (!el->has_efficiency) {
    *amk_scalar(c, AMK_S_PKA_NO_EFF) += 1.0;
    return;
  }
  double arc = nrt;
  double rpa = nrt;
  if (td_eV >= knee) {
    const double xi_arc = el->arc_scale * AMK_POW_FAST(td_eV, el->arc_b) + el->arc_c;
    const double tc = AMK_POW_FAST(td_eV, el->rpa_c);
    const double xi_rpa = el->rpa_scale * tc / (el->rpa_bc + tc);
    arc = nrt * xi_arc;
    rpa = nrt * xi_rpa;
  }
  amk_hadd(c, AMK_H_ARC, x, arc);
  amk_hadd(c, AMK_H_RPA, x, rpa);
  *amk_scalar(c, AMK_S_ARC) += arc;
  *amk_scalar(c, AMK_S_RPA) += rpa;
}

AMK_FN double amk_layer_damage_energy(AMK_G const KLayer *t, double T_eV) {
  const double eps = t->kp_eps_per_eV * T_eV;
  const double g = eps + 0.4 * AMK_POW_FAST(eps, 0.75) + 3.4 * AMK_POW_FAST(eps, 1.0 / 6.0);
  return T_eV / (1.0 + t->kp_k * g);
}

AMK_FN void amk_spectra(KCtx *c, double er, int pka) {
  const amk_i64 bin = amk_spectrum_bin(er);
  if (bin >= 0 && bin < AMK_SPECTRUM_BINS) {
    const amk_i64 s = amk_off_spectra(c->P->nbins, c->P->n_tally);
    if (pka) *amk_t(c, s + bin) += 1.0;
    *amk_t(c, s + AMK_SPECTRUM_BINS + bin) += 1.0;
  }
}

AMK_FN void amk_displace(KCtx *c, AMK_G const KElement *el, double x, int by_ion, double T_eV) {
  amk_hadd(c, AMK_H_DISPLACEMENTS, x, 1.0);
  *amk_scalar(c, AMK_S_DISPLACEMENTS) += 1.0;
  amk_eadd(c, AMK_E_VACANCIES, el->tally_index, x, 1.0);
  *amk_scalar(c, AMK_S_VACANCIES) += 1.0;
  if (by_ion) amk_hadd(c, AMK_H_KNOCK_ONS, x, 1.0);
  const double er = amk_max(T_eV - el->eb_eV, 1.0);
  amk_spectra(c, er, by_ion);
}

AMK_FN void amk_kp(KCtx *c, AMK_G const KElement *el, AMK_G const KLayer *t, double T_eV, double x) {
  const double eps = t->kp_eps_per_eV * T_eV;
  const double g = eps + 0.4 * AMK_POW_FAST(eps, 0.75) + 3.4 * AMK_POW_FAST(eps, 1.0 / 6.0);
  const double ev = T_eV / (1.0 + t->kp_k * g);
  const double ed = el->ed_eV;
  const double eb = el->eb_eV;
  amk_book_pka(c, el, x, ev);
  amk_hadd(c, AMK_H_ION_BY_RECOILS, x, T_eV - ev);
  *amk_scalar(c, AMK_S_IONIZ_RECOILS) += T_eV - ev;
  if (ev < ed) {
    amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, ev);
    *amk_scalar(c, AMK_S_PHON_RECOILS) += ev;
    return;
  }
  const double nu = ev > 2.5 * ed ? 0.8 * ev / (2.0 * ed) : 1.0;
  amk_hadd(c, AMK_H_DISPLACEMENTS, x, 1.0);
  *amk_scalar(c, AMK_S_DISPLACEMENTS) += 1.0;
  amk_hadd(c, AMK_H_KNOCK_ONS, x, 1.0);
  amk_hadd(c, AMK_H_VAC_BY_ION, x, 1.0);
  amk_hadd(c, AMK_H_VAC_BY_RECOILS, x, nu - 1.0);
  amk_eadd(c, AMK_E_VACANCIES, el->tally_index, x, nu);
  *amk_scalar(c, AMK_S_VACANCIES) += nu;
  {
    const double er = amk_max(T_eV - eb, 1.0);
    amk_spectra(c, er, 1);
  }
  amk_hadd(c, AMK_H_E2RECOIL, x, -eb);
  *amk_scalar(c, AMK_S_E2RECOIL) -= eb;
  amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, ev - nu * eb);
  *amk_scalar(c, AMK_S_PHON_RECOILS) += ev - nu * eb;
  *amk_scalar(c, AMK_S_BINDING) += nu * eb;
}

AMK_FN void amk_recoil_stop(KCtx *c, amk_i64 k) {
  AMK_G const KParams *P = c->P;
  AMK_G const KMover *r = c->slots + k;
  const double x = AMK_A(P, r->x);
  const double e = AMK_EV(P, r->energy);
  if (c->track_primary && r->primary) {
    c->primary_stop_x = x;
    c->has_primary_stop = 1;
  }
  if (amk_in_window(P, x)) {
    amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, e);
    *amk_scalar(c, AMK_S_PHON_RECOILS) += e;
  }
  if (!c->primary_displaced) {
    c->primary_displaced = 1;
    return;
  }
  if (c->sub_threshold) {
    c->sub_threshold = 0;
    return;
  }
  int replaced = 0;
  if (k > 1) {
    AMK_G const KMover *below = c->slots + (k - 1);
    AMK_G const KElement *bl =
        c->elements + c->layers[below->origin_layer].first_element + below->element_hint;
    if (below->energy > bl->ed && below->species == r->species) {
      replaced = 1;
      amk_eadd(c, AMK_E_VACANCIES, bl->tally_index, x, -1.0);
      amk_eadd(c, AMK_E_VAC_REMOVED, bl->tally_index, x, 1.0);
      *amk_scalar(c, AMK_S_VACANCIES) -= 1.0;
      amk_hadd(c, AMK_H_REPLACEMENTS, x, 1.0);
      *amk_scalar(c, AMK_S_REPLACEMENTS) += 1.0;
      const double eb = bl->eb_eV;
      amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, eb);
      *amk_scalar(c, AMK_S_PHON_RECOILS) += eb;
    }
  }
  amk_eadd(c, AMK_E_RECOIL_STOPPED, r->tally_element, x, 1.0);
  if (replaced && !(P->bugs & AMK_BUG_REPLACEMENT_INTERSTITIAL)) return;
  amk_hadd(c, AMK_H_INTERSTITIALS, x, 1.0);
  *amk_scalar(c, AMK_S_INTERSTITIALS) += 1.0;
}

AMK_FN void amk_surface_arrival(KCtx *c, AMK_G const KMover *r) {
  const double en = AMK_EV(c->P, r->energy) * r->cx * r->cx;
  double b = AMK_RINT(en / 0.1);
  if (!(b >= 1.0)) b = 1.0;
  if (b > (double)AMK_SURFACE_BINS) b = (double)AMK_SURFACE_BINS;
  const amk_i64 h = amk_off_surface(c->P->nbins, c->P->n_tally) + (amk_i64)r->tally_element * AMK_SURFACE_BINS;
  *amk_t(c, h + (amk_i64)b - 1) += 1.0;
}

//! Append an exit record; false if the block's list is full (status set).
AMK_FN int amk_exit(KCtx *c, double kind, double z, double energy, double x, double y, double zz, double cx,
                    double cy, double cz) {
  AMK_G const KParams *P = c->P;
  AMK_G double *n = amk_scalar(c, AMK_S_EXITS);
  if (!P->record_exits || *n >= (double)P->max_exit_records) return 1;
  if (*n >= (double)P->exit_capacity) {
    c->status = AMK_EXITS_FULL;
    return 0;
  }
  AMK_G double *e = c->lists + AMK_BLOCK_IONS * AMK_POSITION_DOUBLES + (amk_i64)(*n) * AMK_EXIT_DOUBLES;
  e[0] = kind;
  e[1] = (double)(c->ion_index + 1);
  e[2] = z;
  e[3] = energy;
  e[4] = x;
  e[5] = y;
  e[6] = zz;
  e[7] = cx;
  e[8] = cy;
  e[9] = cz;
  *n += 1.0;
  return 1;
}

// ---------------------------------------------------------------------------
// The step machine: the trim2013 loops of run_ion_2013 and cascade_2013 cut
// into steps (one flight of the ion, one step of the recoil stack), so that
// the lanes of a wavefront do the same kind of work at a time and a launch
// can stop after any step. The order of every operation within a block is
// the CPU transport's.
// ---------------------------------------------------------------------------

#if defined(__OPENCL_VERSION__)
#define AMK_ATOMIC_INC(p) atomic_inc(p)
#define AMK_ATOMIC_CAS(p, cmp, val) atomic_cmpxchg((p), (cmp), (val))
#define AMK_ATOMIC_LOAD(p) (*(p))
#define AMK_COUNTER volatile int
#elif defined(__HIPCC__) || defined(__CUDACC__) || defined(__HIPCC_RTC__) || defined(__CUDACC_RTC__)
#define AMK_ATOMIC_INC(p) atomicAdd((p), 1)
#define AMK_ATOMIC_CAS(p, cmp, val) atomicCAS((p), (cmp), (val))
#define AMK_ATOMIC_LOAD(p) atomicAdd((p), 0)
#define AMK_COUNTER int
#elif defined(AMK_KOKKOS)
#define AMK_ATOMIC_INC(p) Kokkos::atomic_fetch_add((p), 1)
#define AMK_ATOMIC_CAS(p, cmp, val) Kokkos::atomic_compare_exchange((p), (cmp), (val))
#define AMK_ATOMIC_LOAD(p) Kokkos::atomic_load(p)
#define AMK_COUNTER int
#else
#include <atomic>
#define AMK_ATOMIC_INC(p) std::atomic_ref<int>(*(p)).fetch_add(1)
static inline int amk_cas_host(int *p, int cmp, int val) {
  int expected = cmp;
  std::atomic_ref<int>(*p).compare_exchange_strong(expected, val);
  return expected; // the old value (cmp on success)
}
#define AMK_ATOMIC_CAS(p, cmp, val) amk_cas_host((p), (cmp), (val))
#define AMK_ATOMIC_LOAD(p) std::atomic_ref<int>(*(p)).load()
#define AMK_COUNTER int
#endif

//! One iteration of cascade_2013's loop: NEXTRECOIL (pop) or RECOILLOOP.
AMK_FN void amk_cascade_step(KCtx *c) {
  AMK_G const KParams *P = c->P;
  KRng *rng = c->cascade_on_ion_rng ? &c->rng : &c->crng;
  const double stop = P->recoil_stop;
  const double th = P->thickness_A;
  const double mono1 = P->first_monolayer;
  amk_i64 k = c->k;
  if (c->pop) {
    // NEXTRECOIL.
    if (--k == 0) {
      c->k = 0;
      c->phase = AMK_PHASE_AFTER_CASCADE;
      return;
    }
    c->k = k;
    c->sub_threshold = 0;
    AMK_G KMover *r = c->slots + k;
    AMK_G const KElement *el = c->elements + c->layers[r->origin_layer].first_element + r->element_hint;
    const double x = AMK_A(P, r->x);
    const double rn = amk_rand(rng);
    if (r->energy < el->ed && r->x > 9.0 * mono1 * rn) {
      if (r->energy != 0.0 * P->ev && amk_in_window(P, x)) {
        const double e = AMK_EV(P, r->energy);
        amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, e);
        *amk_scalar(c, AMK_S_PHON_RECOILS) += e;
      }
      return;
    }
    if (r->energy >= el->ed) {
      if (amk_in_window(P, x)) {
        amk_displace(c, el, x, 0, AMK_EV(P, r->energy));
      } else {
        *amk_scalar(c, AMK_S_DISPLACEMENTS) += 1.0;
        *amk_scalar(c, AMK_S_VACANCIES) += 1.0;
      }
      *amk_scalar(c, AMK_S_RECOILS_FOLLOWED) += 1.0;
      if (r->x > 2.0 * mono1) {
        r->energy -= el->eb;
        *amk_scalar(c, AMK_S_BINDING) += el->eb_eV;
      }
    } else {
      c->sub_threshold = 1;
    }
    if (r->energy > stop) {
      c->pop = 0;
    } else {
      amk_recoil_stop(c, k);
      c->pop = 1;
    }
    return;
  }

  // RECOILLOOP.
  AMK_G KMover *r = c->slots + k;
  if (!(r->energy >= 0.0 * P->ev)) {
    c->status = AMK_NOT_A_NUMBER;
    return;
  }
  const double l0 = c->layers[r->layer_hint].monolayer;
  r->x += l0 * r->cx;
  r->y += l0 * r->cy;
  r->zz += l0 * r->cz;
  double x = AMK_A(P, r->x);
  if (x > th) {
    *amk_array(c, AMK_A_TRANSMITTED_RECOILS, r->tally_element) += 1.0;
    c->pop = 1;
    return;
  }
  if (x < 0.0) {
    AMK_G const KElement *el = c->elements + c->layers[r->origin_layer].first_element + r->element_hint;
    const double normal = r->energy * r->cx * r->cx;
    amk_surface_arrival(c, r);
    if (normal < el->es) {
      r->cx = -r->cx;
      r->x = 1e-4 * P->ang;
      x = 1e-4;
    } else {
      const double out = amk_max(r->energy - el->es, 0.0 * P->ev);
      const double eo = AMK_EV(P, out);
      *amk_array(c, AMK_A_SPUTTERED, r->tally_element) += 1.0;
      *amk_array(c, AMK_A_SPUTTERED_E, r->tally_element) += eo;
      if (c->sub_threshold && (P->bugs & AMK_BUG_SUBTHRESHOLD_SPUTTER)) {
        *amk_scalar(c, AMK_S_DISPLACEMENTS) += 1.0;
        *amk_scalar(c, AMK_S_VACANCIES) += 1.0;
      }
      if (!amk_exit(c, 83.0, (double)r->z, AMK_EV(P, r->energy), 0.0, AMK_A(P, r->y), AMK_A(P, r->zz),
                    amk_clamp(r->cx, -0.9999999, 0.9999999), amk_clamp(r->cy, -0.9999999, 0.9999999),
                    amk_clamp(r->cz, -0.9999999, 0.9999999))) {
        return;
      }
      r->energy = 0.0 * P->ev;
      c->pop = 1;
      return;
    }
  }
  const amk_i64 li = amk_layer_index(c, r->x);
  r->layer_hint = (amk_i32)li;
  double de = amk_electronic_loss(c, r->species, r->energy, li, c->layers[li].monolayer, rng, 0, 0.0);
  if (de > r->energy) de = r->energy;
  {
    const double ev = AMK_EV(P, de);
    if (amk_in_window(P, x)) {
      amk_hadd(c, AMK_H_ION_BY_RECOILS, x, ev);
      *amk_scalar(c, AMK_S_IONIZ_RECOILS) += ev;
    }
  }
  if (r->energy < stop) {
    amk_recoil_stop(c, k);
    c->pop = 1;
    return;
  }
  AMK_G const KLayer *tl = c->layers + li;
  const amk_i64 j = amk_pick_element(tl, c->elements, rng);
  AMK_G const KElement *el = c->elements + tl->first_element + j;
  double p = AMK_SQRT_FAST(amk_rand(rng)) * tl->pmax0;
  if (r->z < 21 && r->x < 2.0 * mono1) {
    if (r->z < 14) {
      p *= amk_rand(rng);
    } else {
      const double p1 = AMK_SQRT_FAST(amk_rand(rng)) * tl->pmax0;
      const double p2 = p1 * AMK_SQRT_FAST(amk_rand(rng));
      p = p2 + (p1 - p2) * ((double)r->z - 13.0) / 7.0;
    }
  }
  AMK_G const KPair *pair = c->pairs + (amk_i64)r->species * P->n_elements_total + tl->first_element + j;
  const double eps = r->energy * pair->rep;
  const KKin kin = amk_scatter(P, pair, r->energy, eps, p / pair->screening_length);
  double T = amk_min(kin.transferred, r->energy);
  if (!(P->bugs & AMK_BUG_RECOIL_ENERGY_CREATED) && T > r->energy - de) {
    T = amk_max(r->energy - de, 0.0 * P->ev);
  }
  double cphi = 1.0;
  double sphi = 0.0;
  amk_azimuth(c, rng, &cphi, &sphi);

  if (k + 1 >= P->slot_capacity) {
    c->status = AMK_SLOTS_FULL;
    return;
  }
  KMover struck = *r;
  KMover moved = *r;
  moved.energy = amk_max(r->energy - T - de, 0.0 * P->ev);
  amk_rotate(&moved.cx, &moved.cy, &moved.cz, kin.cos_projectile, kin.sin_projectile, -cphi, -sphi);
  struck.species = (amk_i32)el->species;
  struck.z = (amk_i32)el->z;
  struck.tally_element = (amk_i32)el->tally_index;
  struck.layer_hint = (amk_i32)li;
  struck.element_hint = (amk_i32)j;
  struck.origin_layer = (amk_i32)li;
  struck.energy = T;
  struck.primary = 0;
  if (P->bugs & AMK_BUG_STRUCK_ATOM_DIRECTION) {
    amk_rotate_about_x(&struck.cx, &struck.cy, &struck.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
  } else {
    amk_rotate(&struck.cx, &struck.cy, &struck.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
  }
  moved.layer_hint = (amk_i32)li;
  c->slots[k] = struck;
  c->slots[k + 1] = moved;
  ++k;
  c->k = k;
  if (k > c->max_k) c->max_k = k;
  if (P->max_recoils != 0 && ++c->recoils_this_ion > P->max_recoils) {
    *amk_scalar(c, AMK_S_RECOIL_LIMIT) += 1.0;
  }
  if (c->slots[k].energy >= stop) {
    c->pop = 0;
  } else {
    amk_recoil_stop(c, k);
    c->pop = 1;
  }
}

//! start_cascade_2013: the primary takes slot 1 at the collision site; the
//! lane then steps through the cascade.
AMK_FN void amk_start_cascade(KCtx *c, AMK_G const KElement *el, amk_i64 li, amk_i64 j, const KMover *site,
                              double T, int restore_rng) {
  AMK_G const KParams *P = c->P;
  const double x1 = AMK_A(P, site->x);
  const double T_eV = AMK_EV(P, T);
  AMK_G KMover *r = c->slots + 1;
  r->species = (amk_i32)el->species;
  r->z = (amk_i32)el->z;
  r->tally_element = (amk_i32)el->tally_index;
  r->layer_hint = (amk_i32)li;
  r->element_hint = (amk_i32)j;
  r->origin_layer = (amk_i32)li;
  r->x = site->x;
  r->y = site->y;
  r->zz = site->zz;
  r->cx = site->cx;
  r->cy = site->cy;
  r->cz = site->cz;
  r->primary = (amk_i32)c->track_primary;
  c->sub_threshold = 0;
  if (T > el->ed) {
    c->primary_displaced = 1;
    if (amk_in_window(P, x1)) {
      amk_displace(c, el, x1, 1, T_eV);
    } else {
      *amk_scalar(c, AMK_S_DISPLACEMENTS) += 1.0;
      *amk_scalar(c, AMK_S_VACANCIES) += 1.0;
    }
    *amk_scalar(c, AMK_S_RECOILS_FOLLOWED) += 1.0;
    r->energy = amk_max(T - el->eb, 0.0 * P->ev);
    *amk_scalar(c, AMK_S_BINDING) += el->eb_eV;
  } else {
    c->primary_displaced = 0; // the slot keeps its energy from the previous cascade
  }
  if (restore_rng && !(P->bugs & AMK_BUG_CASCADE_REPLAY)) {
    c->crng = amk_make_stream((amk_u64)P->seed, (amk_u64)c->ion_index, (amk_u64)(2 + c->cascades_this_ion++));
    c->cascade_on_ion_rng = 0;
  } else if (restore_rng) {
    c->crng = c->rng; // TRIM saves the generator and restores it after the cascade
    c->cascade_on_ion_rng = 0;
  } else {
    c->cascade_on_ion_rng = 1; // TRIM.DAT: the ion continues after the cascade's numbers
  }
  c->k = 1;
  c->pop = 0;
  c->phase = AMK_PHASE_CASCADE;
}

AMK_FN void amk_moments(KCtx *c, amk_i64 m, double x) {
  const amk_i64 s = amk_off_moments(c->P->nbins, c->P->n_tally) + m;
  const double x2 = x * x;
  *amk_t(c, s) += 1.0;
  *amk_t(c, s + 1) += x;
  *amk_t(c, s + 2) += x2;
  *amk_t(c, s + 3) += x2 * x;
  *amk_t(c, s + 4) += x2 * x2;
}

//! stop_mover for the ion.
AMK_FN void amk_stop_ion(KCtx *c, const KMover *mv) {
  AMK_G const KParams *P = c->P;
  const double x = AMK_A(P, mv->x);
  const double e = AMK_EV(P, mv->energy);
  const double y = AMK_A(P, mv->y);
  const double z = AMK_A(P, mv->zz);
  const double r = AMK_SQRT(y * y + z * z);
  amk_hadd(c, AMK_H_ION_STOPPED, x, 1.0);
  amk_moments(c, AMK_M_RANGE, x);
  amk_moments(c, AMK_M_LAT_Y, y);
  amk_moments(c, AMK_M_LAT_ABS_Y, AMK_FABS(y));
  amk_moments(c, AMK_M_RADIAL, r);
  const amk_i64 b = amk_bin_of(P, x);
  if (b >= 0) {
    const amk_i64 lb = amk_off_lateral(P->nbins, P->n_tally) + 5 * b;
    *amk_t(c, lb) += 1.0;
    *amk_t(c, lb + 1) += AMK_FABS(y);
    *amk_t(c, lb + 2) += y * y;
    *amk_t(c, lb + 3) += r;
    *amk_t(c, lb + 4) += r * r;
  }
  amk_hadd(c, AMK_H_PHON_BY_ION, x, e);
  *amk_scalar(c, AMK_S_PHON_ION) += e;
  if (P->record_positions) {
    AMK_G double *n = amk_scalar(c, AMK_S_POSITIONS);
    AMK_G double *p = c->lists + (amk_i64)(*n) * AMK_POSITION_DOUBLES;
    p[0] = (double)(c->ion_index + 1);
    p[1] = x;
    p[2] = y;
    p[3] = z;
    *n += 1.0;
  }
}

//! record_exyz: the ion's state (EXYZ.txt), up to the block's capacity.
AMK_FN void amk_exyz(KCtx *c, double se_eV_per_A) {
  AMK_G const KParams *P = c->P;
  AMK_G double *n = amk_scalar(c, AMK_S_EXYZ);
  if (*n >= (double)P->max_exit_records) return;
  if (*n >= (double)P->exyz_capacity) {
    c->status = AMK_EXITS_FULL;
    return;
  }
  AMK_G double *r = c->lists + AMK_BLOCK_IONS * AMK_POSITION_DOUBLES + P->exit_capacity * AMK_EXIT_DOUBLES +
                    (amk_i64)(*n) * AMK_EXYZ_DOUBLES;
  const double e = AMK_EV(P, c->mv.energy);
  r[0] = (double)(c->ion_index + 1);
  r[1] = e < AMK_EV(P, P->ion_stop) ? 0.0 : e / 1000.0;
  r[2] = AMK_A(P, c->mv.x);
  r[3] = AMK_A(P, c->mv.y);
  r[4] = AMK_A(P, c->mv.zz);
  r[5] = se_eV_per_A;
  r[6] = c->last_T_eV;
  *n += 1.0;
}

//! The end of run_ion_2013.
AMK_FN void amk_ion_end(KCtx *c) {
  *amk_scalar(c, AMK_S_IONS) += 1.0;
  c->phase = AMK_PHASE_NEXT_ION;
}

//! The end of run_recoil_event_2013.
AMK_FN void amk_recoil_event_end(KCtx *c) {
  c->track_primary = 0;
  if (c->has_primary_stop && (c->P->bugs & AMK_BUG_RECOIL_EVENT_TWICE)) {
    amk_hadd(c, AMK_H_ION_STOPPED, c->primary_stop_x, 1.0);
  }
  c->has_primary_stop = 0;
  amk_ion_end(c);
}

/*!
 * run_recoil_event_2013: TRIM.DAT codes 6 (0x449330). An atom of the target
 * is set in motion at the given point and direction, without an incident
 * ion; TRIM books the event as an ion that stops where the atom starts.
 */
AMK_FN void amk_recoil_event_begin(KCtx *c, amk_i64 ion_index) {
  AMK_G const KParams *P = c->P;
  AMK_G const KParticle *sp = c->particles + ion_index / P->particle_repeat;
  c->ion_index = ion_index;
  c->recoils_this_ion = 0;
  c->cascades_this_ion = 0;
  c->rng = amk_make_stream((amk_u64)P->seed, (amk_u64)ion_index, 0);
  AMK_G const KLayer *layer = c->layers + sp->element_layer;
  AMK_G const KElement *el = c->elements + layer->first_element + sp->element;
  KMover site;
  site.x = sp->x;
  site.y = sp->y;
  site.zz = sp->z;
  site.cx = sp->cx;
  site.cy = sp->cy;
  site.cz = sp->cz;
  site.energy = 0.0 * P->ev;
  site.species = 0;
  site.z = 0;
  site.tally_element = 0;
  site.layer_hint = 0;
  site.element_hint = 0;
  site.origin_layer = 0;
  site.primary = 0;
  site.pad = 0;
  const double x = AMK_A(P, sp->x);
  const double T = sp->energy;
  const double T_eV = AMK_EV(P, T);
  amk_stop_ion(c, &site);
  c->has_primary_stop = 0;
  c->track_primary = 1;
  if (T >= P->recoil_stop && P->mode != AMK_MODE_KP) {
    amk_book_pka(c, el, x, amk_layer_damage_energy(c->layers + sp->layer, T_eV));
  }
  if (T < P->recoil_stop) {
    amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, T_eV);
    *amk_scalar(c, AMK_S_PHON_RECOILS) += T_eV;
  } else if (P->mode == AMK_MODE_KP) {
    amk_kp(c, el, c->layers + sp->layer, T_eV, x);
  } else if (T <= el->ed && !(P->bugs & AMK_BUG_STALE_PRIMARY)) {
    amk_hadd(c, AMK_H_PHON_BY_RECOILS, x, T_eV);
    *amk_scalar(c, AMK_S_PHON_RECOILS) += T_eV;
  } else {
    amk_start_cascade(c, el, sp->element_layer, sp->element, &site, T, 0);
    return; // the lane runs the cascade, then amk_recoil_event_end
  }
  amk_recoil_event_end(c);
}

//! The start of run_ion_2013 (beam and TRIM.DAT ions).
AMK_FN void amk_ion_begin(KCtx *c, amk_i64 ion_index) {
  if (c->P->recoil_events) {
    amk_recoil_event_begin(c, ion_index);
    return;
  }
  AMK_G const KParams *P = c->P;
  KMover *mv = &c->mv;
  mv->species = 0;
  mv->z = (amk_i32)P->ion_z;
  mv->tally_element = 0;
  mv->layer_hint = 0;
  mv->element_hint = 0;
  mv->origin_layer = 0;
  mv->primary = 0;
  mv->pad = 0;
  c->ion_index = ion_index;
  c->recoils_this_ion = 0;
  c->cascades_this_ion = 0;
  c->rng = amk_make_stream((amk_u64)P->seed, (amk_u64)ion_index, 0);
  KRng *rng = &c->rng;
  const double l0 = P->first_monolayer;
  if (P->beam) {
    mv->x = 0.0 * P->ang;
    mv->cx = P->ion_cx;
    mv->cy = P->ion_cy;
    mv->cz = 0.0;
    mv->energy = P->ion_energy;
    c->flight_counter = 0;
    mv->y = l0 * mv->cy * amk_rand(rng);
    mv->zz = l0 * mv->cz * amk_rand(rng);
  } else {
    AMK_G const KParticle *sp = c->particles + ion_index / P->particle_repeat;
    mv->x = sp->x;
    mv->cx = sp->cx;
    mv->cy = sp->cy;
    mv->cz = sp->cz;
    const double norm = AMK_SQRT(mv->cx * mv->cx + mv->cy * mv->cy + mv->cz * mv->cz);
    if (norm > 0.0) {
      mv->cx /= norm;
      mv->cy /= norm;
      mv->cz /= norm;
    }
    mv->energy = sp->energy;
    if (sp->x == 0.0 * P->ang) c->flight_counter = 0;
    mv->y = sp->y == 0.0 * P->ang ? l0 * mv->cy * amk_rand(rng) : sp->y;
    mv->zz = sp->z == 0.0 * P->ang ? l0 * mv->cz * amk_rand(rng) : sp->z;
  }
  c->last_T_eV = 0.0;
  if (P->exyz_interval_eV > 0.0) {
    c->exyz_next = (amk_i64)AMK_FLOOR(AMK_EV(P, mv->energy) / P->exyz_interval_eV) + 1;
  }
  c->phase = AMK_PHASE_ION;
}

//! The ion's loop after a collision (and its cascade): it stops or flies on.
AMK_FN void amk_ion_after_collision(KCtx *c) {
  if (c->mv.energy <= c->P->ion_stop) {
    amk_stop_ion(c, &c->mv);
    if (c->P->exyz_interval_eV > 0.0) amk_exyz(c, 0.0);
    amk_ion_end(c);
  } else {
    c->phase = AMK_PHASE_ION;
  }
}

//! One iteration of run_ion_2013's loop: a flight and a collision.
AMK_FN void amk_ion_step(KCtx *c) {
  AMK_G const KParams *P = c->P;
  KMover *mvp = &c->mv;
  KRng *rng = &c->rng;
  const double stop = P->ion_stop;
  const int full = P->mode != AMK_MODE_KP;
  const int restore_rng = (int)P->beam;
  if (!(mvp->energy >= 0.0 * P->ev) || !AMK_ISFINITE(AMK_A(P, mvp->x))) {
    c->status = AMK_NOT_A_NUMBER;
    return;
  }
  amk_i64 li = amk_layer_index(c, mvp->x);
  const double se = amk_stopping(c, 0, li, mvp->energy);
  if (P->exyz_interval_eV > 0.0) {
    // 0x43460b; a record at the TRIM.IN energy sets the next one a step lower (0x434745).
    const double e = AMK_EV(P, mvp->energy);
    const amk_i64 step = (amk_i64)AMK_FLOOR(e / P->exyz_interval_eV);
    const int at_start_energy = mvp->energy == P->ion_energy;
    if (at_start_energy || step < c->exyz_next) {
      amk_exyz(c, se / P->ev_per_ang);
      c->exyz_next = at_start_energy && (P->bugs & AMK_BUG_EXYZ_SKIPS_FIRST) ? step - 1 : step;
    }
  }
  const double L = amk_flight_path(c, mvp, li, rng);
  const double de = amk_electronic_loss(c, 0, mvp->energy, li, L, rng, 1, se);
  mvp->energy -= de;
  if (mvp->energy > P->ion_energy) mvp->energy = P->ion_energy;
  if (mvp->energy < stop) mvp->energy = stop;
  const double x0 = AMK_A(P, mvp->x);
  mvp->x += L * mvp->cx;
  mvp->y += L * mvp->cy;
  mvp->zz += L * mvp->cz;
  const double x1 = AMK_A(P, mvp->x);
  const double de_eV = AMK_EV(P, de);
  const int crosses = !((x0 < P->w0 && x1 < P->w0) || (x0 > P->w1 && x1 > P->w1));
  const amk_i64 h_ion = amk_off_single(P->nbins, AMK_H_ION_BY_ION);
  if (mvp->x < 0.0 * P->ang) {
    if (crosses) {
      amk_deposit(c, h_ion, amk_max(x0, 0.0), 0.0, de_eV * (x0 / amk_max(x0 - x1, 1e-12)));
      *amk_scalar(c, AMK_S_IONIZ_ION) += de_eV * (x0 / amk_max(x0 - x1, 1e-12));
    }
    // leave_front (the ion).
    if (mvp->cx < 0.0) {
      const double s = AMK_A(P, mvp->x) / mvp->cx;
      mvp->y -= s * mvp->cy * P->ang;
      mvp->zz -= s * mvp->cz * P->ang;
    }
    mvp->x = 0.0 * P->ang;
    const double e = AMK_EV(P, mvp->energy);
    *amk_scalar(c, AMK_S_BACKSCATTERED) += 1.0;
    *amk_scalar(c, AMK_S_BACKSCATTERED_E) += e;
    amk_exit(c, 66.0, (double)mvp->z, e, 0.0, AMK_A(P, mvp->y), AMK_A(P, mvp->zz), mvp->cx, mvp->cy, mvp->cz);
    if (P->exyz_interval_eV > 0.0) amk_exyz(c, 0.0);
    amk_ion_end(c);
    return;
  }
  if (mvp->x > P->thickness) {
    const double th = P->thickness_A;
    if (crosses) {
      const double frac = (th - x0) / amk_max(x1 - x0, 1e-12);
      amk_deposit(c, h_ion, x0, th, de_eV * frac);
      *amk_scalar(c, AMK_S_IONIZ_ION) += de_eV * frac;
    }
    // leave_back (the ion).
    const double e = AMK_EV(P, mvp->energy);
    *amk_scalar(c, AMK_S_TRANSMITTED) += 1.0;
    *amk_scalar(c, AMK_S_TRANSMITTED_E) += e;
    amk_exit(c, 84.0, (double)mvp->z, e, P->thickness_A, AMK_A(P, mvp->y), AMK_A(P, mvp->zz), mvp->cx, mvp->cy,
             mvp->cz);
    if (P->exyz_interval_eV > 0.0) amk_exyz(c, 0.0);
    amk_ion_end(c);
    return;
  }
  li = amk_layer_index(c, mvp->x);
  AMK_G const KLayer *layer = c->layers + li;
  const amk_i64 j = amk_pick_element(layer, c->elements, rng);
  AMK_G const KElement *el = c->elements + layer->first_element + j;
  const double p = AMK_SQRT_FAST(amk_rand(rng) / (P->pi * L * layer->density));
  AMK_G const KPair *pair = c->pairs + layer->first_element + j; // species 0
  const double eps = mvp->energy * pair->rep;
  const KKin kin = amk_scatter(P, pair, mvp->energy, eps, p / pair->screening_length);
  const double T = amk_min(kin.transferred, mvp->energy);
  double cphi = 1.0;
  double sphi = 0.0;
  amk_azimuth(c, rng, &cphi, &sphi);
  amk_rotate(&mvp->cx, &mvp->cy, &mvp->cz, kin.cos_projectile, kin.sin_projectile, -cphi, -sphi);
  mvp->energy -= T;
  const double T_eV = AMK_EV(P, T);
  c->last_T_eV = T_eV;

  if (crosses) {
    const double a = amk_clamp(x0, P->w0, P->w1);
    const double b = amk_clamp(x1, P->w0, P->w1);
    amk_deposit(c, h_ion, a, b, de_eV);
    *amk_scalar(c, AMK_S_IONIZ_ION) += de_eV;
    if (T >= el->ed) {
      amk_deposit(c, amk_off_single(P->nbins, AMK_H_E2RECOIL), a, b, T_eV);
      amk_deposit(c, amk_off_element(P->nbins, P->n_tally, AMK_E_ABSORBED, el->tally_index), a, b, T_eV);
      *amk_scalar(c, AMK_S_E2RECOIL) += T_eV;
    } else {
      amk_deposit(c, amk_off_single(P->nbins, AMK_H_PHON_BY_ION), a, b, T_eV);
      *amk_scalar(c, AMK_S_PHON_ION) += T_eV;
    }
  }

  if (amk_in_window(P, x1)) {
    if (T < P->recoil_stop) {
      if (P->bugs & AMK_BUG_ION_PHONONS_TWICE) {
        amk_hadd(c, AMK_H_PHON_BY_RECOILS, x1, T_eV);
        *amk_scalar(c, AMK_S_PHON_RECOILS) += T_eV;
      }
    } else if (!full) {
      amk_kp(c, el, layer, T_eV, x1);
    } else {
      KMover site = *mvp;
      amk_rotate(&site.cx, &site.cy, &site.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
      amk_book_pka(c, el, x1, amk_layer_damage_energy(layer, T_eV));
      if (T > el->ed || (P->bugs & AMK_BUG_STALE_PRIMARY)) {
        amk_start_cascade(c, el, li, j, &site, T, restore_rng);
        return; // the lane runs the cascade, then amk_ion_after_collision
      }
    }
  }
  amk_ion_after_collision(c);
}

//! Ions in block \p block (absolute block number).
AMK_FN amk_i64 amk_block_ions(AMK_G const KParams *P, amk_i64 block) {
  amk_i64 count = P->n_ions - block * AMK_BLOCK_IONS;
  return count > AMK_BLOCK_IONS ? AMK_BLOCK_IONS : count;
}

//! The end of a block: its status (and how deep its recoil stack went).
AMK_FN void amk_block_end(KCtx *c) {
  if (c->status != AMK_OK) *amk_scalar(c, AMK_S_STATUS_ION) = (double)(c->ion_index + 1);
  *amk_scalar(c, AMK_S_STATUS) = (double)c->status;
  *amk_scalar(c, AMK_S_MAX_K) = (double)c->max_k;
  *amk_scalar(c, AMK_S_STEPS) = (double)c->steps;
}

/*!
 * Device block \p d (the run's block first_block + d) in its ring slot
 * d % ring: a new Transporter, as TrimSimulation::run_block makes for every
 * block (zero tallies, default recoil slots, its counters).
 */
AMK_FN void amk_block_begin(KCtx *c, AMK_G double *tallies, AMK_G double *lists, amk_i64 d) {
  AMK_G const KParams *P = c->P;
  const amk_i64 slot = d % c->ring;
  c->block = d;
  c->T = tallies + slot; // the slot's column (amk_t)
  c->lists = lists + slot * P->list_stride;
  for (amk_i64 i = 0; i < P->tally_stride; ++i) *amk_t(c, i) = 0.0;
  for (amk_i64 s = 0; s < 2 && s < P->slot_capacity; ++s) {
    KMover m;
    m.x = 0.0;
    m.y = 0.0;
    m.zz = 0.0;
    m.cx = 1.0;
    m.cy = 0.0;
    m.cz = 0.0;
    m.energy = 0.0;
    m.species = 0;
    m.z = 0;
    m.tally_element = 0;
    m.layer_hint = 0;
    m.element_hint = 0;
    m.origin_layer = 0;
    m.primary = 0;
    m.pad = 0;
    c->slots[s] = m;
  }
  c->status = AMK_OK;
  c->flight_counter = 0;
  c->primary_displaced = 0;
  c->sub_threshold = 0;
  c->ion_in_block = 0;
  c->max_k = 0;
  c->steps = 0;
}

/*!
 * One step of a lane; returns 1 when the lane must wait for the next
 * launch (the ring of block slots is full). A lane that has finished an
 * ion starts the next ion of its block; at the end of a block it marks the
 * slot finished (slot_block[slot] = d + 1) and takes the next device block
 * from counter[0], if its slot is free (d < limit: the host merged the
 * block that held it). Which lane runs a block changes nothing.
 */
AMK_FN int amk_step(KCtx *c, AMK_G AMK_COUNTER *counter, AMK_G amk_i64 *slot_block, AMK_G double *tallies,
                    AMK_G double *lists, amk_i64 n_dev, amk_i64 limit) {
  if (c->status != AMK_OK) c->phase = AMK_PHASE_NEXT_ION; // the block is given back to the CPU
  switch (c->phase) {
  case AMK_PHASE_ION:
    ++c->steps;
    amk_ion_step(c);
    return 0;
  case AMK_PHASE_CASCADE:
    ++c->steps;
    amk_cascade_step(c);
    return 0;
  case AMK_PHASE_AFTER_CASCADE:
    ++c->steps;
    if (c->P->recoil_events) amk_recoil_event_end(c);
    else amk_ion_after_collision(c);
    return 0;
  case AMK_PHASE_NEXT_ION: {
    if (c->block >= 0) {
      if (c->status == AMK_OK && c->ion_in_block < amk_block_ions(c->P, c->first_block + c->block)) {
        ++c->steps;
        amk_ion_begin(c, c->P->ion0 + (c->first_block + c->block) * AMK_BLOCK_IONS + c->ion_in_block);
        ++c->ion_in_block;
        return 0;
      }
      amk_block_end(c);
      slot_block[c->block % c->ring] = c->block + 1;
      c->block = -1;
    }
    int n = AMK_ATOMIC_LOAD(counter);
    while (1) {
      if ((amk_i64)n >= n_dev) {
        c->phase = AMK_PHASE_DONE;
        return 0;
      }
      if ((amk_i64)n >= limit) return 1;
      const int old = AMK_ATOMIC_CAS(counter, n, n + 1);
      if (old == n) break;
      n = old;
    }
    amk_block_begin(c, tallies, lists, (amk_i64)n);
    return 0; // the block's first ion starts with the next step
  }
  default: return 1;
  }
}

AMK_FN void amk_lane_load(KCtx *c, AMK_G const KLane *s) {
  c->mv = s->ion;
  c->rng = s->rng;
  c->crng = s->crng;
  c->phase = s->phase;
  c->block = s->block;
  c->ion_in_block = s->ion_in_block;
  c->ion_index = s->ion_index;
  c->k = s->k;
  c->pop = s->pop;
  c->cascade_on_ion_rng = s->cascade_on_ion_rng;
  c->status = s->status;
  c->flight_counter = s->flight_counter;
  c->primary_displaced = s->primary_displaced;
  c->sub_threshold = s->sub_threshold;
  c->cascades_this_ion = s->cascades_this_ion;
  c->recoils_this_ion = s->recoils_this_ion;
  c->max_k = s->max_k;
  c->steps = s->spare[0];
  c->track_primary = s->spare[1] & 1;
  c->has_primary_stop = (s->spare[1] >> 1) & 1;
  c->primary_stop_x = s->primary_stop_x;
  c->last_T_eV = s->last_T_eV;
  c->exyz_next = s->exyz_next;
}
AMK_FN void amk_lane_store(const KCtx *c, AMK_G KLane *s) {
  s->ion = c->mv;
  s->rng = c->rng;
  s->crng = c->crng;
  s->phase = c->phase;
  s->block = c->block;
  s->ion_in_block = c->ion_in_block;
  s->ion_index = c->ion_index;
  s->k = c->k;
  s->pop = c->pop;
  s->cascade_on_ion_rng = c->cascade_on_ion_rng;
  s->status = c->status;
  s->flight_counter = c->flight_counter;
  s->primary_displaced = c->primary_displaced;
  s->sub_threshold = c->sub_threshold;
  s->cascades_this_ion = c->cascades_this_ion;
  s->recoils_this_ion = c->recoils_this_ion;
  s->max_k = c->max_k;
  s->spare[0] = c->steps;
  s->spare[1] = c->track_primary | (c->has_primary_stop << 1);
  s->primary_stop_x = c->primary_stop_x;
  s->last_T_eV = c->last_T_eV;
  s->exyz_next = c->exyz_next;
}

/*!
 * Lane \p lane: up to \p steps steps (a launch), from its state in
 * \p lanes (\p init: the run's first launch; every lane starts without a
 * block). The device runs blocks first_block .. + n_dev; device block d
 * keeps its tallies in ring slot d % ring (element e at e * ring + slot) until
 * the host merges it; \p limit (merged + ring) is the first block without a
 * free slot. A lane that runs out of blocks counts itself in counter[1].
 */
AMK_FN void amk_run_lane(AMK_G const KParams *P, AMK_G const KLayer *layers, AMK_G const KElement *elements,
                         AMK_G const KPair *pairs, AMK_G const KCurve *curves, AMK_G const double *log_s,
                         AMK_G const double *edges, AMK_G const KParticle *particles, AMK_G KMover *slots,
                         AMK_G double *tallies, AMK_G double *lists, AMK_G KLane *lanes, AMK_G AMK_COUNTER *counter,
                         AMK_G amk_i64 *slot_block, amk_i64 lane, amk_i64 first_block, amk_i64 n_dev, amk_i64 ring,
                         amk_i64 limit, amk_i64 steps, amk_i64 init) {
  KCtx c;
  c.P = P;
  c.layers = layers;
  c.elements = elements;
  c.pairs = pairs;
  c.curves = curves;
  c.log_s = log_s;
  c.edges = edges;
  c.particles = particles;
  c.slots = slots + lane * P->slot_capacity;
  c.T = tallies;
  c.lists = lists;
  c.ns = ring;
  c.ring = ring;
  c.first_block = first_block;
  if (init) {
    c.phase = AMK_PHASE_NEXT_ION;
    c.block = -1;
    c.status = AMK_OK;
    c.ion_in_block = 0;
    c.ion_index = 0;
    c.k = 0;
    c.pop = 0;
    c.cascade_on_ion_rng = 0;
    c.flight_counter = 0;
    c.primary_displaced = 0;
    c.sub_threshold = 0;
    c.cascades_this_ion = 0;
    c.recoils_this_ion = 0;
    c.max_k = 0;
    c.steps = 0;
    c.track_primary = 0;
    c.has_primary_stop = 0;
    c.primary_stop_x = 0.0;
    c.last_T_eV = 0.0;
    c.exyz_next = 0;
    c.rng.s1 = 1;
    c.rng.s2 = 1;
    c.crng = c.rng;
  } else {
    amk_lane_load(&c, lanes + lane);
    if (c.phase == AMK_PHASE_DONE) return;
  }
  if (c.block >= 0) {
    const amk_i64 slot = c.block % ring;
    c.T = tallies + slot;
    c.lists = lists + slot * P->list_stride;
  }
  for (amk_i64 s = 0; s < steps && c.phase != AMK_PHASE_DONE; ++s) {
    if (amk_step(&c, counter, slot_block, tallies, lists, n_dev, limit)) break;
  }
  if (c.phase == AMK_PHASE_DONE) AMK_ATOMIC_INC(counter + 1);
  amk_lane_store(&c, lanes + lane);
}

// ---------------------------------------------------------------------------
// Keeping the results on the device: the running total, block by block in
// block order (TrimTallies::add element by element: the CPU's sums), and the
// lists packed to their length.
// ---------------------------------------------------------------------------

//! The finished blocks at the head of the ring: merged, merged + 1, ... (out[0]).
AMK_FN void amk_scan(AMK_G const amk_i64 *slot_block, amk_i64 ring, amk_i64 merged, amk_i64 n_dev,
                     AMK_G amk_i64 *out) {
  amk_i64 c = 0;
  while (c < ring && merged + c < n_dev && slot_block[(merged + c) % ring] == merged + c + 1) ++c;
  out[0] = c;
}

//! Element e of the running total: add device blocks merged .. merged + count - 1, in order.
AMK_FN void amk_merge_element(AMK_G const KParams *P, AMK_G const double *tallies, AMK_G double *total,
                              amk_i64 ring, amk_i64 merged, amk_i64 count, amk_i64 e) {
  if (e >= P->tally_stride) return;
  AMK_G const double *row = tallies + e * ring;
  double acc = total[e];
  for (amk_i64 k = 0; k < count; ++k) acc += row[(merged + k) % ring];
  total[e] = acc;
}

//! Ring slot \p slot takes the flat tallies \p flat (a block the CPU reran).
AMK_FN void amk_put_element(AMK_G const KParams *P, AMK_G double *tallies, AMK_G const double *flat, amk_i64 ring,
                            amk_i64 slot, amk_i64 e) {
  if (e >= P->tally_stride) return;
  tallies[e * ring + slot] = flat[e];
}

/*!
 * The positions and exit records of device block merged + k to their place
 * in the packed lists (offsets[2k], offsets[2k + 1]: the records before it,
 * in block order); blocks the CPU reran (status not AMK_OK) are skipped.
 */
AMK_FN void amk_pack_lists(AMK_G const KParams *P, AMK_G const double *tallies, AMK_G const double *lists,
                           AMK_G const amk_i64 *offsets, AMK_G double *packed_positions, AMK_G double *packed_exits,
                           amk_i64 ring, amk_i64 merged, amk_i64 count, amk_i64 k) {
  if (k >= count) return;
  const amk_i64 slot = (merged + k) % ring;
  const amk_i64 s = amk_off_scalars(P->nbins, P->n_tally);
  if (tallies[(s + AMK_S_STATUS) * ring + slot] != (double)AMK_OK) return;
  AMK_G const double *src = lists + slot * P->list_stride;
  const amk_i64 np = (amk_i64)tallies[(s + AMK_S_POSITIONS) * ring + slot];
  const amk_i64 ne = (amk_i64)tallies[(s + AMK_S_EXITS) * ring + slot];
  AMK_G double *pp = packed_positions + offsets[2 * k] * AMK_POSITION_DOUBLES;
  for (amk_i64 i = 0; i < np * AMK_POSITION_DOUBLES; ++i) pp[i] = src[i];
  AMK_G double *pe = packed_exits + offsets[2 * k + 1] * AMK_EXIT_DOUBLES;
  AMK_G const double *se = src + AMK_BLOCK_IONS * AMK_POSITION_DOUBLES;
  for (amk_i64 i = 0; i < ne * AMK_EXIT_DOUBLES; ++i) pe[i] = se[i];
  // EXYZ records follow the exits in the packed exit buffer, after all of the exits (offset from the host).
  if (P->exyz_capacity > 0) {
    const amk_i64 nx = (amk_i64)tallies[(s + AMK_S_EXYZ) * ring + slot];
    AMK_G double *px = packed_exits + offsets[2 * count] + offsets[2 * count + 1 + k] * AMK_EXYZ_DOUBLES;
    AMK_G const double *sx = se + P->exit_capacity * AMK_EXIT_DOUBLES;
    for (amk_i64 i = 0; i < nx * AMK_EXYZ_DOUBLES; ++i) px[i] = sx[i];
  }
}

#if defined(__OPENCL_VERSION__) || defined(__HIPCC__) || defined(__CUDACC__) || defined(__HIPCC_RTC__) || \
    defined(__CUDACC_RTC__)
AMK_KERNEL void amcset_trim_scan(AMK_G const amk_i64 *slot_block, amk_i64 ring, amk_i64 merged, amk_i64 n_dev,
                                 AMK_G amk_i64 *out) {
  if (AMK_GLOBAL_ID() == 0) amk_scan(slot_block, ring, merged, n_dev, out);
}
AMK_KERNEL void amcset_trim_merge(AMK_G const KParams *P, AMK_G const double *tallies, AMK_G double *total,
                                  amk_i64 ring, amk_i64 merged, amk_i64 count) {
  amk_merge_element(P, tallies, total, ring, merged, count, AMK_GLOBAL_ID());
}
AMK_KERNEL void amcset_trim_put(AMK_G const KParams *P, AMK_G double *tallies, AMK_G const double *flat, amk_i64 ring,
                                amk_i64 slot) {
  amk_put_element(P, tallies, flat, ring, slot, AMK_GLOBAL_ID());
}
AMK_KERNEL void amcset_trim_pack(AMK_G const KParams *P, AMK_G const double *tallies, AMK_G const double *lists,
                                 AMK_G const amk_i64 *offsets, AMK_G double *packed_positions,
                                 AMK_G double *packed_exits, amk_i64 ring, amk_i64 merged, amk_i64 count) {
  amk_pack_lists(P, tallies, lists, offsets, packed_positions, packed_exits, ring, merged, count, AMK_GLOBAL_ID());
}
//! Work item i is lane i (of n_lanes).
AMK_KERNEL void amcset_trim_lanes(AMK_G const KParams *P, AMK_G const KLayer *layers,
                                  AMK_G const KElement *elements, AMK_G const KPair *pairs,
                                  AMK_G const KCurve *curves, AMK_G const double *log_s, AMK_G const double *edges,
                                  AMK_G const KParticle *particles, AMK_G KMover *slots, AMK_G double *tallies,
                                  AMK_G double *lists, AMK_G KLane *lanes, AMK_G AMK_COUNTER *counter,
                                  AMK_G amk_i64 *slot_block, amk_i64 n_lanes, amk_i64 first_block, amk_i64 n_dev,
                                  amk_i64 ring, amk_i64 limit, amk_i64 steps, amk_i64 init) {
  const amk_i64 i = AMK_GLOBAL_ID();
  if (i >= n_lanes) return;
  amk_run_lane(P, layers, elements, pairs, curves, log_s, edges, particles, slots, tallies, lists, lanes, counter,
               slot_block, i, first_block, n_dev, ring, limit, steps, init);
}
#endif

#endif // AMCSET_TRIM_KERNEL_H
