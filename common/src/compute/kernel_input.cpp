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

/*! \file kernel_input.cpp
 *  \brief The kernel's arrays from SRIM mode's target model, the block
 *  tallies back, and the backend registry.
 *
 * Every number the kernel reads is computed here with the expression the
 * CPU transport evaluates (Transporter's constructor and the *_2013
 * functions), so that the kernel compiled as C++ reproduces it bit for bit.
 */

#include "device.h"

#include "../trim/target_model.h"
#include "trim/trim.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

static_assert(sizeof(KParams) % 8 == 0 && sizeof(KLayer) == 17 * 8 && sizeof(KElement) == 16 * 8 &&
                  sizeof(KPair) == 12 * 8 && sizeof(KCurve) == 6 * 8 && sizeof(KParticle) == 10 * 8 &&
                  sizeof(KMover) == 88 && sizeof(KLane) % 8 == 0,
              "kernel structs must have the same layout on every compiler");

namespace amcset {
namespace compute {

using namespace common;
using trim::TrimConfig;

std::string to_string(Backend b) {
  switch (b) {
  case Backend::cpu: return "cpu";
  case Backend::opencl: return "opencl";
  case Backend::hip: return "hip";
  case Backend::cuda: return "cuda";
  case Backend::kokkos: return "kokkos";
  case Backend::reference: return "reference";
  }
  return "cpu";
}

Backend parse_backend(std::string_view s) {
  std::string t(s);
  std::transform(t.begin(), t.end(), t.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
  for (Backend b : all_backends()) {
    if (t == to_string(b)) return b;
  }
  throw std::invalid_argument("unknown compute backend '" + std::string(s) +
                              "' (cpu, opencl, hip, cuda, kokkos or reference)");
}

std::string to_string(Precision p) { return p == Precision::mixed ? "mixed" : "double"; }

Precision parse_precision(std::string_view s) {
  if (s == "double") return Precision::double_precision;
  if (s == "mixed") return Precision::mixed;
  throw std::invalid_argument("unknown precision '" + std::string(s) + "' (double or mixed)");
}

std::vector<Backend> all_backends() {
  return {Backend::cpu, Backend::opencl, Backend::hip, Backend::cuda, Backend::kokkos, Backend::reference};
}

bool built_with(Backend b) {
#if defined(AMCSET_WITH_KOKKOS)
  constexpr bool kokkos = true;
#else
  constexpr bool kokkos = false;
#endif
  return b != Backend::kokkos || kokkos;
}

std::vector<DeviceInfo> list_devices(std::vector<std::string> *notes) {
  std::vector<DeviceInfo> out;
  DeviceInfo cpu;
  cpu.backend = Backend::cpu;
  cpu.name = "CPU";
  cpu.details = "the CPU transport on all cores";
  out.push_back(cpu);
  list_opencl(out, notes);
  list_hip(out, notes, false);
  list_hip(out, notes, true);
  list_kokkos(out, notes);
  DeviceInfo ref;
  ref.backend = Backend::reference;
  ref.name = "CPU (kernel source)";
  ref.details = "the GPU kernel compiled for the CPU; reproduces the CPU transport (for checking)";
  out.push_back(ref);
  return out;
}

std::unique_ptr<Device> open_device(const Options &options) {
  std::unique_ptr<Device> d;
  switch (options.backend) {
  case Backend::cpu: throw std::logic_error("open_device: the CPU backend has no device");
  case Backend::reference: d = open_reference(options); break;
  case Backend::opencl: d = open_opencl(options); break;
  case Backend::hip: d = open_hip(options, false); break;
  case Backend::cuda: d = open_hip(options, true); break;
  case Backend::kokkos: d = open_kokkos(options); break;
  }
  if (!d) throw std::runtime_error("compute backend " + to_string(options.backend) + " is not available");
  return d;
}

// ---------------------------------------------------------------------------
// What runs where
// ---------------------------------------------------------------------------

namespace {
//! EXYZ records of a block: a record at the start, one per energy step and one at the end, per ion.
std::size_t exyz_records_per_block(const TrimConfig &c) {
  if (!(c.exyz_interval > 0.0 * electron_volt)) return 0;
  const double steps = std::floor(in_units(c.ion.energy, electron_volt) / in_units(c.exyz_interval, electron_volt));
  return static_cast<std::size_t>(std::min(steps + 3.0, 1e9)) * AMK_BLOCK_IONS;
}
constexpr std::size_t max_exyz_per_block = 4096;
//! Below this many ions a GPU run is usually slower than the CPU's (a note).
constexpr std::size_t small_run_ions = 100000;
} // namespace

/*
 * What the GPU kernel leaves to the CPU, and why (docs/GPU.md "What runs
 * where"): each was weighed against what a GPU would gain.
 */
std::string unsupported_reason(const TrimConfig &c) {
  if (c.transport == TrimConfig::Transport::crystal) {
    return "the crystal transport runs on the CPU only (its lattice search, thermal displacements, emptied sites "
           "and grains are not in the GPU kernel, which has the trim2013 transport)";
  }
  if (c.transport == TrimConfig::Transport::refined) {
    return "the refined transport runs on the CPU only (AMCSET's earlier transport, kept for comparison; the GPU "
           "kernel has the trim2013 transport)";
  }
  if (c.md_events.enabled()) {
    return "MD events (md_events) are recorded on the CPU only (they serve to hand a modest number of cascades "
           "to molecular dynamics)";
  }
  if (exyz_records_per_block(c) > max_exyz_per_block) {
    return "EXYZ.txt with " + std::to_string(exyz_records_per_block(c) / AMK_BLOCK_IONS) +
           " records per ion (exyz_interval) is recorded on the CPU (at most " +
           std::to_string(max_exyz_per_block / AMK_BLOCK_IONS) + " on a GPU; writing that many records, not the "
           "transport, would take the time)";
  }
  if (c.trajectory_ions == TrimConfig::all_ions) {
    return "trajectories of every ion are recorded on the CPU only (writing the points, not the transport, takes "
           "the time)";
  }
  return {};
}

std::vector<std::string> device_notes(const TrimConfig &c) {
  std::vector<std::string> notes;
  if (c.fast_math && c.compute.precision == Precision::double_precision) {
    notes.push_back("physics.fast_math selects compute.precision mixed on a GPU (its tables are the CPU's)");
  }
  // The 16 ions of a block run in order, so a run has ion_count / 16
  // independent work items; a GPU needs some ten thousand of them, and each
  // of its lanes is far slower than a CPU core (docs/GPU.md "Which to
  // choose": 2,000 ions of 1 MeV Fe took 3 s on the CPU and 40 to 64 s on an
  // RX 7900 XTX).
  const bool gpu = c.compute.backend == Backend::opencl || c.compute.backend == Backend::hip ||
                   c.compute.backend == Backend::cuda;
  if (gpu && c.ion_count < small_run_ions) {
    notes.push_back(std::to_string(c.ion_count) + " ions are " +
                    std::to_string((c.ion_count + AMK_BLOCK_IONS - 1) / AMK_BLOCK_IONS) +
                    " blocks of 16, too few to fill a GPU: below about " + std::to_string(small_run_ions) +
                    " ions the CPU backend is usually faster");
  }
  return notes;
}

bool block_needs_cpu(const TrimConfig &c, std::size_t block) {
  const std::size_t first = block * AMK_BLOCK_IONS; // relative to first_ion
  if (first < c.trajectory_ions) return true;
  if (c.first_ion + first < c.collision_details_ions) return true;
  return false;
}

// ---------------------------------------------------------------------------
// Kernel input
// ---------------------------------------------------------------------------

std::size_t recoil_slots_needed(const TrimConfig &cfg, const trim::TargetModel &m) {
  // TRIM's recoil stack grows by one slot with every collision of the atom
  // on top (cascade_2013): its depth is about the number of monolayer steps
  // along a recoil's path. The longest path is that of the target atom with
  // the most energy any collision can give (gamma E0 from the ion), as its
  // CSDA path length; nested recoils start higher up the stack, so twice
  // that, plus a margin. A block that still needs more is rerun on the CPU.
  if (cfg.damage == trim::DamageCalculation::kinchin_pease) return 4;
  double mono_A = 1e300;
  for (const auto &l : m.layers) mono_A = std::min(mono_A, in_units(l.monolayer, angstrom));
  double t_max_eV = 0.0;
  std::vector<std::pair<unsigned, mass_quantity>> atoms;
  for (const auto &l : m.layers) {
    for (const auto &e : l.elements) {
      const double m1 = in_units(m.species.front().mass, atomic_mass_unit);
      const double m2 = in_units(e.mass, atomic_mass_unit);
      t_max_eV = std::max(t_max_eV, 4.0 * m1 * m2 / ((m1 + m2) * (m1 + m2)) * in_units(cfg.ion.energy, electron_volt));
      atoms.emplace_back(e.z, e.mass);
    }
  }
  double path_A = 0.0;
  for (const auto &[z, mass] : atoms) {
    TrimConfig c = cfg;
    c.ion.z = z;
    c.ion.mass = mass;
    c.ion.energy = t_max_eV * electron_volt;
    c.source = TrimConfig::Source::beam;
    try {
      path_A = std::max(path_A, in_units(trim::estimate_penetration(c), angstrom));
    } catch (const std::exception &) {
      path_A = std::max(path_A, in_units(m.thickness, angstrom));
    }
  }
  const double slots = 2.0 * path_A / std::max(mono_A, 0.1) + 64.0;
  return static_cast<std::size_t>(std::clamp(slots, 64.0, 65536.0));
}

KernelInput make_kernel_input(const TrimConfig &cfg, const trim::TargetModel &m, std::size_t ions,
                              std::size_t tally_elements, double bin0, double bin_w) {
  KernelInput in;
  KParams &P = in.params;
  P.ang = angstrom.value();
  P.ev = electron_volt.value();
  P.ev_per_ang = electron_volt_per_angstrom.value();
  P.pi = pi_v;
  P.stop_floor = (0.001 * electron_volt).value();
  P.inv_ang = 1.0 / P.ang;
  P.inv_ev = 1.0 / P.ev;
  // Transporter's constructor.
  const double inv_w = 1.0 / bin_w;
  P.bin0 = bin0;
  P.inv_w = inv_w;
  P.w0 = in_units(cfg.window_start, angstrom);
  P.w1 = cfg.window_end > 0.0 * angstrom ? in_units(cfg.window_end, angstrom) : in_units(m.thickness, angstrom);
  P.thickness = m.thickness.value();
  P.thickness_A = in_units(m.thickness, angstrom);
  P.ion_energy = cfg.ion.energy.value();
  P.ion_cx = std::cos(in_units(cfg.ion.angle, radian));
  P.ion_cy = std::sin(in_units(cfg.ion.angle, radian));
  P.ion_stop = (cfg.ion_stop_energy > 0.0 * electron_volt ? cfg.ion_stop_energy : cfg.final_energy).value();
  P.recoil_stop = (cfg.recoil_stop_energy > 0.0 * electron_volt ? cfg.recoil_stop_energy : cfg.final_energy).value();
  P.first_monolayer = m.layers.front().monolayer.value();
  P.nbins = static_cast<amk_i64>(cfg.depth_bins);
  P.n_layers = static_cast<amk_i64>(m.layers.size());
  P.n_species = static_cast<amk_i64>(m.species.size());
  P.n_tally = static_cast<amk_i64>(tally_elements);
  P.ion_z = cfg.ion.z;
  P.mode = static_cast<amk_i64>(cfg.damage);
  P.beam = cfg.source == TrimConfig::Source::beam ? 1 : 0;
  const trim::TrimBugs &b = cfg.bugs;
  P.bugs = (b.stale_primary_energy ? AMK_BUG_STALE_PRIMARY : 0) | (b.recoil_energy_created ? AMK_BUG_RECOIL_ENERGY_CREATED : 0) |
           (b.struck_atom_direction ? AMK_BUG_STRUCK_ATOM_DIRECTION : 0) |
           (b.replacement_interstitial ? AMK_BUG_REPLACEMENT_INTERSTITIAL : 0) |
           (b.subthreshold_sputter_counted ? AMK_BUG_SUBTHRESHOLD_SPUTTER : 0) |
           (b.ion_phonons_counted_twice ? AMK_BUG_ION_PHONONS_TWICE : 0) |
           (b.cascade_random_replay ? AMK_BUG_CASCADE_REPLAY : 0) |
           (b.recoil_event_counted_twice ? AMK_BUG_RECOIL_EVENT_TWICE : 0) |
           (b.exyz_skips_first_step ? AMK_BUG_EXYZ_SKIPS_FIRST : 0);
  P.recoil_events = cfg.source == TrimConfig::Source::particle_recoils && cfg.damage != trim::DamageCalculation::monolayer ? 1 : 0;
  P.straggling = cfg.electronic_straggling ? 1 : 0;
  P.integral = cfg.scattering == physics::ScatteringMethod::integral ? 1 : 0;
  P.seed = cfg.seed;
  P.max_recoils = static_cast<amk_i64>(cfg.max_recoils_per_ion);
  P.record_positions = cfg.record_ion_positions ? 1 : 0;
  P.record_exits = cfg.record_exiting_atoms ? 1 : 0;
  P.max_exit_records = static_cast<amk_i64>(std::min<std::size_t>(cfg.max_exit_records, std::size_t{1} << 62));
  P.slot_capacity = static_cast<amk_i64>(cfg.compute.recoil_slots != 0 ? cfg.compute.recoil_slots
                                                                       : recoil_slots_needed(cfg, m));
  // Per block (deterministic: a shared pool would make which blocks overflow
  // depend on the lanes' timing); none when the lists are not kept.
  P.exit_capacity = !cfg.record_exiting_atoms ? 0
                    : static_cast<amk_i64>(cfg.compute.exit_records != 0 ? cfg.compute.exit_records : 256);
  if (P.slot_capacity < 4) P.slot_capacity = 4;
  P.ion0 = static_cast<amk_i64>(cfg.first_ion);
  P.n_ions = static_cast<amk_i64>(ions);
  P.particle_repeat = static_cast<amk_i64>(std::max<std::size_t>(cfg.particle_repeat, 1));
  P.exyz_interval_eV = in_units(cfg.exyz_interval, electron_volt);
  P.exyz_capacity = static_cast<amk_i64>(exyz_records_per_block(cfg));
  P.tally_stride = amk_tally_size(P.nbins, P.n_tally);
  P.list_stride = amk_list_size(P.exit_capacity, P.exyz_capacity);

  in.edges.resize(cfg.depth_bins + 1);
  for (std::size_t k = 0; k <= cfg.depth_bins; ++k) in.edges[k] = bin0 + static_cast<double>(k) / inv_w;

  for (const auto &l : m.layers) {
    KLayer k{};
    k.top = l.top.value();
    k.bottom = l.bottom.value();
    k.density = l.density.value();
    k.monolayer = l.monolayer.value();
    k.monolayer_A = in_units(l.monolayer, angstrom);
    k.bottom_A = in_units(l.bottom, angstrom);
    k.density_A3 = in_units(l.density, per_cubic_angstrom);
    k.pmax0 = l.t13.pmax0.value();
    k.screening = l.t13.screening.value();
    k.fi = l.t13.fi.value();
    k.epsdg = l.t13.epsdg;
    k.kp_eps_per_eV = l.t13.kp_eps_per_eV;
    k.kp_k = l.t13.kp_k;
    k.k_low = l.t13.k_low;
    k.k_bohr = l.t13.k_bohr;
    k.first_element = static_cast<amk_i64>(in.elements.size());
    k.n_elements = static_cast<amk_i64>(l.elements.size());
    in.layers.push_back(k);
    for (const auto &e : l.elements) {
      KElement ke{};
      ke.fraction = e.fraction;
      ke.ed = e.ed.value();
      ke.eb = e.eb.value();
      ke.es = e.es.value();
      ke.ed_eV = in_units(e.ed, electron_volt);
      ke.eb_eV = in_units(e.eb, electron_volt);
      ke.arc_b = e.arc_b;
      ke.arc_c = e.arc_c;
      ke.rpa_c = e.rpa_c;
      ke.arc_scale = e.arc_scale;
      ke.rpa_bc = e.rpa_bc;
      ke.rpa_scale = e.rpa_scale;
      ke.z = e.z;
      ke.species = static_cast<amk_i64>(e.species);
      ke.tally_index = static_cast<amk_i64>(e.tally_index);
      ke.has_efficiency = e.has_efficiency ? 1 : 0;
      in.elements.push_back(ke);
    }
  }
  P.n_elements_total = static_cast<amk_i64>(in.elements.size());

  for (std::size_t s = 0; s < m.species.size(); ++s) {
    for (std::size_t li = 0; li < m.layers.size(); ++li) {
      const trim::SpeciesLayer &d = m.data[s][li];
      for (std::size_t j = 0; j < m.layers[li].elements.size(); ++j) {
        KPair kp{};
        if (j < d.pairs.size()) {
          const auto &p = d.pairs[j];
          kp.screening_length = p.screening_length.value();
          kp.rep = p.reduced_energy_per_energy.value();
          kp.mass_ratio = p.mass_ratio;
          kp.gamma = p.gamma;
          for (int i = 0; i < 4; ++i) {
            kp.c[i] = p.screening.c[static_cast<std::size_t>(i)];
            kp.d[i] = p.screening.d[static_cast<std::size_t>(i)];
          }
        }
        in.pairs.push_back(kp);
      }
    }
  }
  for (std::size_t s = 0; s < m.species.size(); ++s) {
    const double mc2_eV = in_units(m.species[s].mass, atomic_mass_unit) * 931494.0e3;
    for (std::size_t li = 0; li < m.layers.size(); ++li) {
      const auto &curve = m.data[s][li].stopping;
      KCurve k{};
      k.mc2_eV = mc2_eV;
      k.log_e0 = curve.log_e0();
      k.inv_step = curve.inv_step();
      k.offset = static_cast<amk_i64>(in.log_s.size());
      k.n = static_cast<amk_i64>(curve.log_table().size());
      in.log_s.insert(in.log_s.end(), curve.log_table().begin(), curve.log_table().end());
      in.curves.push_back(k);
    }
  }
  if (in.log_s.empty()) in.log_s.push_back(0.0);
  for (const auto &rp : m.sources) {
    const auto &sp = rp.particle;
    KParticle k{};
    k.x = sp.x.value();
    k.y = sp.y.value();
    k.z = sp.z_position.value();
    k.cx = sp.cos_x;
    k.cy = sp.cos_y;
    k.cz = sp.cos_z;
    k.energy = sp.energy.value();
    k.layer = static_cast<amk_i64>(rp.layer);
    k.element_layer = static_cast<amk_i64>(rp.element_layer);
    k.element = static_cast<amk_i64>(rp.element);
    in.particles.push_back(k);
  }
  if (in.particles.empty()) in.particles.push_back(KParticle{});
  in.blocks = (ions + AMK_BLOCK_IONS - 1) / AMK_BLOCK_IONS;
  return in;
}

std::size_t block_bytes(const KernelInput &in) {
  const KParams &P = in.params;
  return static_cast<std::size_t>(P.tally_stride + P.list_stride) * sizeof(double);
}

std::size_t lane_bytes(const KernelInput &in) {
  return static_cast<std::size_t>(in.params.slot_capacity) * sizeof(KMover) + sizeof(KLane);
}


// ---------------------------------------------------------------------------
// Block tallies back
// ---------------------------------------------------------------------------

int block_status(const KernelInput &in, const double *T) {
  return static_cast<int>(T[amk_off_scalars(in.params.nbins, in.params.n_tally) + AMK_S_STATUS]);
}

namespace {

/*!
 * Every double of TrimTallies that the flat layout holds, with its index:
 * f(double &field, std::size_t index). unpack_block and pack_tallies both
 * use it, so the two directions cannot disagree.
 */
template <class F> void for_each_field(const KernelInput &in, trim::TrimTallies &t, F &&f) {
  const amk_i64 nb = in.params.nbins;
  const amk_i64 ne = in.params.n_tally;
  const auto nbins = static_cast<std::size_t>(nb);
  const auto hist = [&](trim::DepthHistogram &h, amk_i64 off) {
    for (std::size_t b = 0; b < nbins; ++b) f(h.values[b], static_cast<std::size_t>(off) + b);
  };
  const std::pair<trim::DepthHistogram *, int> singles[] = {
      {&t.ion_stopped, AMK_H_ION_STOPPED},       {&t.knock_ons, AMK_H_KNOCK_ONS},
      {&t.vacancies_by_ion, AMK_H_VAC_BY_ION},   {&t.vacancies_by_recoils, AMK_H_VAC_BY_RECOILS},
      {&t.replacements, AMK_H_REPLACEMENTS},     {&t.displacements, AMK_H_DISPLACEMENTS},
      {&t.interstitials, AMK_H_INTERSTITIALS},   {&t.ionization_by_ion, AMK_H_ION_BY_ION},
      {&t.ionization_by_recoils, AMK_H_ION_BY_RECOILS}, {&t.phonons_by_ion, AMK_H_PHON_BY_ION},
      {&t.phonons_by_recoils, AMK_H_PHON_BY_RECOILS}, {&t.energy_to_recoils, AMK_H_E2RECOIL},
      {&t.pka_damage_energy, AMK_H_PKA_TD},      {&t.nrt_displacements, AMK_H_NRT},
      {&t.kp_displacements, AMK_H_KP},           {&t.arc_displacements, AMK_H_ARC},
      {&t.rpa_replacements, AMK_H_RPA}};
  for (const auto &[h, k] : singles) hist(*h, amk_off_single(nb, k));
  for (amk_i64 e = 0; e < ne; ++e) {
    const auto ei = static_cast<std::size_t>(e);
    hist(t.recoil_stopped[ei], amk_off_element(nb, ne, AMK_E_RECOIL_STOPPED, e));
    hist(t.vacancies[ei], amk_off_element(nb, ne, AMK_E_VACANCIES, e));
    hist(t.vacancies_removed[ei], amk_off_element(nb, ne, AMK_E_VAC_REMOVED, e));
    hist(t.energy_absorbed[ei], amk_off_element(nb, ne, AMK_E_ABSORBED, e));
  }
  const auto lat = static_cast<std::size_t>(amk_off_lateral(nb, ne));
  for (std::size_t b = 0; b < nbins; ++b) {
    for (std::size_t k = 0; k < 5; ++k) f(t.lateral_bins[b][k], lat + 5 * b + k);
  }
  const auto S = static_cast<std::size_t>(amk_off_scalars(nb, ne));
  const std::pair<double *, int> scalars[] = {
      {&t.total_pka_damage_energy, AMK_S_PKA_TD}, {&t.total_nrt, AMK_S_NRT}, {&t.total_kp, AMK_S_KP},
      {&t.total_arc, AMK_S_ARC}, {&t.total_rpa, AMK_S_RPA}, {&t.pka_without_efficiency, AMK_S_PKA_NO_EFF},
      {&t.total_displacements, AMK_S_DISPLACEMENTS}, {&t.total_vacancies, AMK_S_VACANCIES},
      {&t.total_replacements, AMK_S_REPLACEMENTS}, {&t.total_interstitials, AMK_S_INTERSTITIALS},
      {&t.total_ionization_ion, AMK_S_IONIZ_ION}, {&t.total_ionization_recoils, AMK_S_IONIZ_RECOILS},
      {&t.total_phonons_ion, AMK_S_PHON_ION}, {&t.total_phonons_recoils, AMK_S_PHON_RECOILS},
      {&t.total_energy_to_recoils, AMK_S_E2RECOIL}, {&t.total_binding_energy_loss, AMK_S_BINDING},
      {&t.backscattered_energy, AMK_S_BACKSCATTERED_E}, {&t.transmitted_energy, AMK_S_TRANSMITTED_E},
      {&t.recoils_followed, AMK_S_RECOILS_FOLLOWED}, {&t.recoil_limit_hits, AMK_S_RECOIL_LIMIT}};
  for (const auto &[v, k] : scalars) f(*v, S + static_cast<std::size_t>(k));
  const auto M = static_cast<std::size_t>(amk_off_moments(nb, ne));
  const std::pair<trim::Moments *, int> moments[] = {{&t.range, AMK_M_RANGE},
                                                      {&t.lateral_y, AMK_M_LAT_Y},
                                                      {&t.lateral_abs_y, AMK_M_LAT_ABS_Y},
                                                      {&t.radial, AMK_M_RADIAL}};
  for (const auto &[m, k] : moments) {
    const std::size_t o = M + static_cast<std::size_t>(k);
    f(m->n, o);
    f(m->s1, o + 1);
    f(m->s2, o + 2);
    f(m->s3, o + 3);
    f(m->s4, o + 4);
  }
  const auto A = static_cast<std::size_t>(amk_off_arrays(nb, ne));
  for (std::size_t e = 0; e < static_cast<std::size_t>(ne); ++e) {
    f(t.sputtered[e], A + AMK_A_SPUTTERED * static_cast<std::size_t>(ne) + e);
    f(t.sputtered_energy[e], A + AMK_A_SPUTTERED_E * static_cast<std::size_t>(ne) + e);
    f(t.transmitted_recoils[e], A + AMK_A_TRANSMITTED_RECOILS * static_cast<std::size_t>(ne) + e);
  }
  const auto sp = static_cast<std::size_t>(amk_off_spectra(nb, ne));
  for (std::size_t i = 0; i < AMK_SPECTRUM_BINS && i < t.pka_spectrum.size(); ++i) {
    f(t.pka_spectrum[i], sp + i);
    f(t.recoil_spectrum[i], sp + AMK_SPECTRUM_BINS + i);
  }
  const auto su = static_cast<std::size_t>(amk_off_surface(nb, ne));
  for (std::size_t e = 0; e < static_cast<std::size_t>(ne); ++e) {
    for (std::size_t k = 0; k < AMK_SURFACE_BINS; ++k) f(t.surface_energy[e][k], su + e * AMK_SURFACE_BINS + k);
  }
}

} // namespace

std::size_t scalar_row(const KernelInput &in, int s) {
  return static_cast<std::size_t>(amk_off_scalars(in.params.nbins, in.params.n_tally) + s);
}

std::vector<double> pack_tallies(const KernelInput &in, const trim::TrimTallies &t) {
  std::vector<double> flat(static_cast<std::size_t>(in.params.tally_stride), 0.0);
  trim::TrimTallies copy = t.without_lists();
  if (copy.recoil_stopped.size() != static_cast<std::size_t>(in.params.n_tally)) {
    copy.resize(static_cast<std::size_t>(in.params.nbins), static_cast<std::size_t>(in.params.n_tally));
  }
  for_each_field(in, copy, [&](double &v, std::size_t i) { flat[i] = v; });
  flat[scalar_row(in, AMK_S_IONS)] = static_cast<double>(copy.ions);
  flat[scalar_row(in, AMK_S_BACKSCATTERED)] = static_cast<double>(copy.backscattered);
  flat[scalar_row(in, AMK_S_TRANSMITTED)] = static_cast<double>(copy.transmitted);
  return flat;
}

void unpack_block(const KernelInput &in, const double *T, const double *lists, trim::TrimTallies &t) {
  t.resize(static_cast<std::size_t>(in.params.nbins), static_cast<std::size_t>(in.params.n_tally));
  for_each_field(in, t, [&](double &v, std::size_t i) { v = T[i]; });
  const double *S = T + amk_off_scalars(in.params.nbins, in.params.n_tally);
  t.ions = static_cast<std::size_t>(S[AMK_S_IONS]);
  t.backscattered = static_cast<std::size_t>(S[AMK_S_BACKSCATTERED]);
  t.transmitted = static_cast<std::size_t>(S[AMK_S_TRANSMITTED]);
  if (lists == nullptr) return;
  const auto n_pos = static_cast<std::size_t>(S[AMK_S_POSITIONS]);
  for (std::size_t i = 0; i < n_pos; ++i) {
    const double *p = lists + i * AMK_POSITION_DOUBLES;
    t.positions.push_back({static_cast<std::size_t>(p[0]), p[1], p[2], p[3]});
  }
  const auto n_exit = static_cast<std::size_t>(S[AMK_S_EXITS]);
  for (std::size_t i = 0; i < n_exit; ++i) {
    const double *e = lists + AMK_BLOCK_IONS * AMK_POSITION_DOUBLES + i * AMK_EXIT_DOUBLES;
    trim::ExitRecord r;
    r.kind = static_cast<char>(e[0]);
    r.ion = static_cast<std::size_t>(e[1]);
    r.z = static_cast<unsigned>(e[2]);
    r.energy_eV = e[3];
    r.x = e[4];
    r.y = e[5];
    r.z_pos = e[6];
    r.cx = e[7];
    r.cy = e[8];
    r.cz = e[9];
    t.exits.push_back(r);
  }
}

void append_positions(const double *packed, std::size_t n, std::vector<trim::IonPosition> &out) {
  for (std::size_t i = 0; i < n; ++i) {
    const double *p = packed + i * AMK_POSITION_DOUBLES;
    out.push_back({static_cast<std::size_t>(p[0]), p[1], p[2], p[3]});
  }
}

void append_exyz(const double *packed, std::size_t n, std::size_t max_records,
                 std::vector<trim::EnergyPositionRecord> &out) {
  for (std::size_t i = 0; i < n && out.size() < max_records; ++i) {
    const double *r = packed + i * AMK_EXYZ_DOUBLES;
    out.push_back({static_cast<std::size_t>(r[0]), r[1], r[2], r[3], r[4], r[5], r[6]});
  }
}

void append_exits(const double *packed, std::size_t n, std::size_t max_records, std::vector<trim::ExitRecord> &out) {
  for (std::size_t i = 0; i < n && out.size() < max_records; ++i) {
    const double *e = packed + i * AMK_EXIT_DOUBLES;
    trim::ExitRecord r;
    r.kind = static_cast<char>(e[0]);
    r.ion = static_cast<std::size_t>(e[1]);
    r.z = static_cast<unsigned>(e[2]);
    r.energy_eV = e[3];
    r.x = e[4];
    r.y = e[5];
    r.z_pos = e[6];
    r.cx = e[7];
    r.cy = e[8];
    r.cz = e[9];
    out.push_back(r);
  }
}

} // namespace compute
} // namespace amcset
