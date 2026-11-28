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

/*!
 * \file compute_test.cpp
 *
 * \brief Tests of the GPU backends (compute/): the kernel source compiled as
 * C++ ("reference") reproduces the CPU transport's tallies exactly; blocks
 * that outgrow the device's buffers or keep trajectories run on the CPU;
 * the GPU frameworks present on the machine agree with the CPU statistically.
 */

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <string>

#include "compute/compute.h"
#include "run_config.h"
#include "trim/trim.h"

using namespace amcset;
using namespace amcset::common;
using namespace amcset::trim;

namespace {

TrimConfig fe_on_fe(double energy_keV, DamageCalculation damage, std::size_t ions) {
  TrimConfig c;
  c.ion = IonBeam::defaults(26, energy_keV * kilo_electron_volt);
  TargetLayer layer;
  layer.name = "Fe";
  layer.density = 7.8658 * gram_per_cubic_centimeter;
  layer.width = (20.0 * energy_keV + 500.0) * angstrom;
  auto fe = TargetElement::defaults(26);
  fe.displacement_energy = 40.0 * electron_volt;
  fe.lattice_binding_energy = 0.0 * electron_volt;
  layer.elements.push_back(fe);
  c.layers.push_back(layer);
  c.damage = damage;
  c.ion_count = ions;
  c.threads = 2;
  return c;
}

//! Ar on SiO2 on Si: two layers, a compound, a light and a heavy recoil.
TrimConfig ar_on_oxide(double energy_keV, DamageCalculation damage, std::size_t ions) {
  TrimConfig c;
  c.ion = IonBeam::defaults(18, energy_keV * kilo_electron_volt);
  c.ion.angle = 7.0 * amcset::common::degree;
  TargetLayer ox;
  ox.name = "SiO2";
  ox.density = 2.2 * gram_per_cubic_centimeter;
  ox.width = 150.0 * angstrom;
  ox.elements = {TargetElement::defaults(14, 1.0), TargetElement::defaults(8, 2.0)};
  TargetLayer sub;
  sub.name = "Si";
  sub.density = 2.3212 * gram_per_cubic_centimeter;
  sub.width = 2000.0 * angstrom;
  sub.elements = {TargetElement::defaults(14)};
  c.layers = {ox, sub};
  c.damage = damage;
  c.ion_count = ions;
  c.threads = 2;
  return c;
}

void expect_same_hist(const DepthHistogram &a, const DepthHistogram &b, const char *what) {
  ASSERT_EQ(a.values.size(), b.values.size()) << what;
  for (std::size_t i = 0; i < a.values.size(); ++i) EXPECT_EQ(a.values[i], b.values[i]) << what << " bin " << i;
}

//! Every tally bit for bit.
void expect_same(const TrimTallies &a, const TrimTallies &b) {
  EXPECT_EQ(a.ions, b.ions);
  const std::pair<const DepthHistogram TrimTallies::*, const char *> hs[] = {
      {&TrimTallies::ion_stopped, "ion_stopped"},
      {&TrimTallies::knock_ons, "knock_ons"},
      {&TrimTallies::vacancies_by_ion, "vacancies_by_ion"},
      {&TrimTallies::vacancies_by_recoils, "vacancies_by_recoils"},
      {&TrimTallies::replacements, "replacements"},
      {&TrimTallies::displacements, "displacements"},
      {&TrimTallies::interstitials, "interstitials"},
      {&TrimTallies::ionization_by_ion, "ionization_by_ion"},
      {&TrimTallies::ionization_by_recoils, "ionization_by_recoils"},
      {&TrimTallies::phonons_by_ion, "phonons_by_ion"},
      {&TrimTallies::phonons_by_recoils, "phonons_by_recoils"},
      {&TrimTallies::energy_to_recoils, "energy_to_recoils"},
      {&TrimTallies::pka_damage_energy, "pka_damage_energy"},
      {&TrimTallies::nrt_displacements, "nrt"},
      {&TrimTallies::kp_displacements, "kp"},
      {&TrimTallies::arc_displacements, "arc"},
      {&TrimTallies::rpa_replacements, "rpa"}};
  for (const auto &[h, name] : hs) expect_same_hist(a.*h, b.*h, name);
  ASSERT_EQ(a.vacancies.size(), b.vacancies.size());
  for (std::size_t e = 0; e < a.vacancies.size(); ++e) {
    expect_same_hist(a.vacancies[e], b.vacancies[e], "vacancies");
    expect_same_hist(a.vacancies_removed[e], b.vacancies_removed[e], "vacancies_removed");
    expect_same_hist(a.recoil_stopped[e], b.recoil_stopped[e], "recoil_stopped");
    expect_same_hist(a.energy_absorbed[e], b.energy_absorbed[e], "energy_absorbed");
    EXPECT_EQ(a.sputtered[e], b.sputtered[e]);
    EXPECT_EQ(a.sputtered_energy[e], b.sputtered_energy[e]);
    EXPECT_EQ(a.transmitted_recoils[e], b.transmitted_recoils[e]);
    EXPECT_EQ(a.surface_energy[e], b.surface_energy[e]);
  }
  EXPECT_EQ(a.lateral_bins, b.lateral_bins);
  for (const auto m : {&TrimTallies::range, &TrimTallies::lateral_y, &TrimTallies::lateral_abs_y, &TrimTallies::radial}) {
    EXPECT_EQ((a.*m).n, (b.*m).n);
    EXPECT_EQ((a.*m).s1, (b.*m).s1);
    EXPECT_EQ((a.*m).s2, (b.*m).s2);
    EXPECT_EQ((a.*m).s3, (b.*m).s3);
    EXPECT_EQ((a.*m).s4, (b.*m).s4);
  }
  const std::pair<double TrimTallies::*, const char *> ss[] = {
      {&TrimTallies::total_pka_damage_energy, "total_pka_damage_energy"},
      {&TrimTallies::total_nrt, "total_nrt"},
      {&TrimTallies::total_kp, "total_kp"},
      {&TrimTallies::total_arc, "total_arc"},
      {&TrimTallies::total_rpa, "total_rpa"},
      {&TrimTallies::pka_without_efficiency, "pka_without_efficiency"},
      {&TrimTallies::total_displacements, "total_displacements"},
      {&TrimTallies::total_vacancies, "total_vacancies"},
      {&TrimTallies::total_replacements, "total_replacements"},
      {&TrimTallies::total_interstitials, "total_interstitials"},
      {&TrimTallies::total_ionization_ion, "total_ionization_ion"},
      {&TrimTallies::total_ionization_recoils, "total_ionization_recoils"},
      {&TrimTallies::total_phonons_ion, "total_phonons_ion"},
      {&TrimTallies::total_phonons_recoils, "total_phonons_recoils"},
      {&TrimTallies::total_energy_to_recoils, "total_energy_to_recoils"},
      {&TrimTallies::total_binding_energy_loss, "total_binding_energy_loss"},
      {&TrimTallies::backscattered_energy, "backscattered_energy"},
      {&TrimTallies::transmitted_energy, "transmitted_energy"},
      {&TrimTallies::recoils_followed, "recoils_followed"},
      {&TrimTallies::recoil_limit_hits, "recoil_limit_hits"}};
  for (const auto &[s, name] : ss) EXPECT_EQ(a.*s, b.*s) << name;
  EXPECT_EQ(a.backscattered, b.backscattered);
  EXPECT_EQ(a.transmitted, b.transmitted);
  EXPECT_EQ(a.pka_spectrum, b.pka_spectrum);
  EXPECT_EQ(a.recoil_spectrum, b.recoil_spectrum);
  ASSERT_EQ(a.positions.size(), b.positions.size());
  for (std::size_t i = 0; i < a.positions.size(); ++i) {
    EXPECT_EQ(a.positions[i].ion, b.positions[i].ion);
    EXPECT_EQ(a.positions[i].x, b.positions[i].x);
    EXPECT_EQ(a.positions[i].y, b.positions[i].y);
    EXPECT_EQ(a.positions[i].z, b.positions[i].z);
  }
  ASSERT_EQ(a.exits.size(), b.exits.size());
  for (std::size_t i = 0; i < a.exits.size(); ++i) {
    const auto &x = a.exits[i];
    const auto &y = b.exits[i];
    EXPECT_EQ(x.kind, y.kind);
    EXPECT_EQ(x.ion, y.ion);
    EXPECT_EQ(x.z, y.z);
    EXPECT_EQ(x.energy_eV, y.energy_eV);
    EXPECT_EQ(x.x, y.x);
    EXPECT_EQ(x.y, y.y);
    EXPECT_EQ(x.z_pos, y.z_pos);
    EXPECT_EQ(x.cx, y.cx);
    EXPECT_EQ(x.cy, y.cy);
    EXPECT_EQ(x.cz, y.cz);
  }
}

TrimTallies run(TrimConfig c, compute::Backend b, TrimSimulation::ComputeReport *report = nullptr) {
  c.compute.backend = b;
  TrimSimulation sim(c);
  sim.run();
  if (report) *report = sim.compute_report();
  return sim.tallies();
}

void expect_reference_matches(const TrimConfig &c) {
  TrimSimulation::ComputeReport rep;
  const TrimTallies cpu = run(c, compute::Backend::cpu);
  const TrimTallies ref = run(c, compute::Backend::reference, &rep);
  EXPECT_FALSE(rep.device.empty());
  EXPECT_GT(rep.device_blocks, 0u);
  expect_same(cpu, ref);
}

} // namespace

TEST(ComputeBackend, NamesRoundTrip) {
  for (auto b : compute::all_backends()) EXPECT_EQ(compute::parse_backend(compute::to_string(b)), b);
  EXPECT_EQ(compute::parse_backend("OpenCL"), compute::Backend::opencl);
  EXPECT_THROW(compute::parse_backend("vulkan"), std::invalid_argument);
}

TEST(ComputeBackend, ScenarioKey) {
  const auto rc = config::parse_config(R"({"ion": {"element": "Fe", "energy": "10 keV"},
      "layers": [{"thickness": "100 nm", "elements": [{"element": "Fe"}]}],
      "compute": {"backend": "opencl", "device": 1, "memory_mib": 256}})");
  EXPECT_EQ(rc.trim.compute.backend, compute::Backend::opencl);
  EXPECT_EQ(rc.trim.compute.device, 1);
  EXPECT_EQ(rc.trim.compute.memory_mib, 256u);
  const auto back = config::parse_config(config::to_json(rc));
  EXPECT_EQ(back.trim.compute, rc.trim.compute);
  const auto short_form = config::parse_config(R"({"ion": {"element": "Fe", "energy": "10 keV"},
      "layers": [{"thickness": "100 nm", "elements": [{"element": "Fe"}]}], "compute": "hip"})");
  EXPECT_EQ(short_form.trim.compute.backend, compute::Backend::hip);
  EXPECT_THROW(config::parse_config(R"({"ion": {"element": "Fe", "energy": "10 keV"},
      "layers": [{"thickness": "100 nm", "elements": [{"element": "Fe"}]}], "compute": {"backend": "metal"}})"),
               std::invalid_argument);
}

TEST(ComputeBackend, ListsTheCpu) {
  std::vector<std::string> notes;
  const auto devices = compute::list_devices(&notes);
  ASSERT_FALSE(devices.empty());
  EXPECT_EQ(devices.front().backend, compute::Backend::cpu);
}

// The kernel source compiled as C++ reproduces the CPU transport bit for bit.
TEST(ComputeKernel, MatchesTheCpuTransportFullCascades) {
  expect_reference_matches(fe_on_fe(20.0, DamageCalculation::full_cascade, 70));
}

TEST(ComputeKernel, MatchesTheCpuTransportQuickDamage) {
  expect_reference_matches(fe_on_fe(50.0, DamageCalculation::kinchin_pease, 100));
}

TEST(ComputeKernel, MatchesTheCpuTransportCompoundLayers) {
  expect_reference_matches(ar_on_oxide(30.0, DamageCalculation::full_cascade, 48));
  expect_reference_matches(ar_on_oxide(30.0, DamageCalculation::kinchin_pease, 64));
}

TEST(ComputeKernel, MatchesTheCpuTransportMonolayerSputtering) {
  auto c = ar_on_oxide(2.0, DamageCalculation::monolayer, 48);
  c.ion.angle = 60.0 * amcset::common::degree;
  expect_reference_matches(c);
}

TEST(ComputeKernel, MatchesTheCpuTransportWithEveryBugFixed) {
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 48);
  c.bugs.set_all(false);
  c.max_recoils_per_ion = 50;
  expect_reference_matches(c);
}

TEST(ComputeKernel, MatchesTheCpuTransportInADepthWindowAndThinTarget) {
  auto c = fe_on_fe(100.0, DamageCalculation::full_cascade, 40);
  c.layers[0].width = 300.0 * angstrom; // most ions go through
  c.window_start = 50.0 * angstrom;
  c.window_end = 250.0 * angstrom;
  c.depth_bins = 37;
  c.first_ion = 1000;
  expect_reference_matches(c);
}

TEST(ComputeKernel, MatchesTheCpuTransportWithTheScatteringIntegral) {
  auto c = fe_on_fe(5.0, DamageCalculation::full_cascade, 16);
  c.scattering = physics::ScatteringMethod::integral;
  expect_reference_matches(c);
}

TEST(ComputeKernel, MatchesTheCpuTransportForTrimDatIons) {
  auto c = fe_on_fe(20.0, DamageCalculation::full_cascade, 0);
  c.source = TrimConfig::Source::particle_ions;
  for (int i = 0; i < 6; ++i) {
    SourceParticle p;
    p.z = 26;
    p.energy = (5.0 + 2.0 * i) * kilo_electron_volt;
    p.x = (i % 2 == 0 ? 0.0 : 30.0 * i) * angstrom;
    p.y = (i % 3 == 0 ? 0.0 : 5.0) * angstrom;
    p.cos_x = 1.0;
    p.cos_y = 0.1 * i;
    p.cos_z = 0.0;
    c.particles.push_back(p);
  }
  c.particle_repeat = 5;
  expect_reference_matches(c);
}

// EXYZ.txt: the ion's state at its start, at every energy step and at its end.
TEST(ComputeKernel, MatchesTheCpuTransportWithExyz) {
  for (const bool bug : {true, false}) {
    auto c = fe_on_fe(50.0, DamageCalculation::full_cascade, 40);
    c.exyz_interval = 5.0 * kilo_electron_volt;
    c.bugs.exyz_skips_first_step = bug;
    expect_reference_matches(c);
    const TrimTallies t = run(c, compute::Backend::reference);
    EXPECT_GT(t.exyz.size(), 40u * 3);
  }
  auto huge = fe_on_fe(50.0, DamageCalculation::kinchin_pease, 16);
  huge.exyz_interval = 1.0 * electron_volt; // 50 000 records per ion: the CPU
  huge.compute.backend = compute::Backend::reference;
  TrimSimulation sim(huge);
  bool warned = false;
  for (const auto &w : sim.warnings()) warned = warned || w.find("EXYZ") != std::string::npos;
  EXPECT_TRUE(warned);
}

// TRIM.DAT recoils (cascade code 6): one atom of the target set in motion
// per event, booked as an ion that stops where it starts.
TEST(ComputeKernel, MatchesTheCpuTransportForTrimDatRecoils) {
  for (const auto damage : {DamageCalculation::full_cascade, DamageCalculation::kinchin_pease}) {
    for (const bool bugs : {true, false}) {
      auto c = fe_on_fe(20.0, damage, 0);
      c.source = TrimConfig::Source::particle_recoils;
      c.bugs.set_all(bugs);
      for (int i = 0; i < 5; ++i) {
        SourceParticle p;
        p.z = 26;
        p.energy = (i == 4 ? 30.0 : 0.5 + 3.0 * i) * kilo_electron_volt; // includes one below Ed? no: 0.5 keV
        p.x = (50.0 + 40.0 * i) * angstrom;
        p.y = 3.0 * i * angstrom;
        p.cos_x = 0.6;
        p.cos_y = 0.8;
        p.cos_z = 0.0;
        c.particles.push_back(p);
      }
      SourceParticle slow;
      slow.z = 26;
      slow.energy = 20.0 * electron_volt; // below Ed (40 eV): the stale-energy branch
      slow.x = 100.0 * angstrom;
      slow.cos_x = 1.0;
      c.particles.push_back(slow);
      c.particle_repeat = 7;
      expect_reference_matches(c);
    }
  }
}

// A block whose recoil stack outgrows the device's slots runs on the CPU
// (the results stay the CPU's), as do the blocks that keep trajectories.
TEST(ComputeKernel, BlocksTheDeviceCannotHoldRunOnTheCpu) {
  auto c = fe_on_fe(20.0, DamageCalculation::full_cascade, 64);
  c.compute.recoil_slots = 6;
  c.trajectory_ions = 20; // blocks 0 and 1 keep trajectories
  TrimSimulation::ComputeReport rep;
  const TrimTallies cpu = run(c, compute::Backend::cpu);
  const TrimTallies ref = run(c, compute::Backend::reference, &rep);
  expect_same(cpu, ref);
  EXPECT_GT(rep.reruns, 0u);
  EXPECT_EQ(rep.cpu_blocks, 2u + rep.reruns); // the trajectory blocks and the reruns
  EXPECT_EQ(rep.device_blocks + rep.reruns, 2u);
}

TEST(ComputeKernel, UnsupportedRunsGoToTheCpuWithAWarning) {
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 16);
  c.transport = TrimConfig::Transport::refined;
  c.compute.backend = compute::Backend::reference;
  TrimSimulation sim(c);
  bool warned = false;
  for (const auto &w : sim.warnings()) warned = warned || w.find("runs on the CPU") != std::string::npos;
  EXPECT_TRUE(warned);
  sim.run();
  EXPECT_TRUE(sim.compute_report().device.empty());
}

// The report both front ends show (CLI summary, GUI log): none on the CPU, four lines on a device.
TEST(ComputeKernel, ReportLines) {
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 64);
  TrimSimulation cpu(c);
  cpu.run();
  EXPECT_TRUE(cpu.compute_report().lines().empty());
  EXPECT_EQ(cpu.summary_json().find("\"compute\""), std::string::npos);
  c.compute.backend = compute::Backend::reference;
  TrimSimulation dev(c);
  dev.run();
  const auto lines = dev.compute_report().lines();
  ASSERT_EQ(lines.size(), 4u);
  EXPECT_EQ(lines[0].rfind("Ran on ", 0), 0u);
  EXPECT_NE(lines[3].find("steps per block"), std::string::npos);
  EXPECT_NE(dev.summary_json().find("\"compute\": {\"backend\": \"reference\""), std::string::npos);
}

// A GPU run of few blocks is slower than the CPU's: a note, the run works.
TEST(ComputeKernel, SmallGpuRunsGetANote) {
  const auto has_note = [](const TrimConfig &c) {
    TrimSimulation sim(c); // the warnings are made here; no device is opened
    for (const auto &n : sim.warnings())
      if (n.find("too few to fill a GPU") != std::string::npos) return true;
    return false;
  };
  auto c = fe_on_fe(10.0, DamageCalculation::full_cascade, 2000);
  c.compute.backend = compute::Backend::opencl;
  EXPECT_TRUE(has_note(c));
  c.ion_count = 1000000;
  EXPECT_FALSE(has_note(c));
  c.ion_count = 2000;
  c.compute.backend = compute::Backend::reference; // runs on the CPU's cores
  EXPECT_FALSE(has_note(c));
}

// Kokkos on a CPU execution space (Serial, Threads, OpenMP) computes with
// the C library's functions: bit for bit the CPU transport.
TEST(ComputeKernel, KokkosOnTheHostMatchesTheCpuTransport) {
  if (!compute::built_with(compute::Backend::kokkos)) GTEST_SKIP() << "AMCSET was built without Kokkos";
  std::string space;
  for (const auto &d : compute::list_devices()) {
    if (d.backend == compute::Backend::kokkos) space = d.name;
  }
  if (space.find("Serial") == std::string::npos && space.find("Threads") == std::string::npos &&
      space.find("OpenMP") == std::string::npos) {
    GTEST_SKIP() << "Kokkos runs on a GPU here (" << space << ")";
  }
  auto c = fe_on_fe(20.0, DamageCalculation::full_cascade, 70);
  const TrimTallies cpu = run(c, compute::Backend::cpu);
  const TrimTallies kk = run(c, compute::Backend::kokkos);
  expect_same(cpu, kk);
}

// Mixed precision (single-precision scattering angle and stopping) gives
// the same physics within the statistics: the reference backend's mixed
// kernel against the CPU, 1600 ions in a compound.
TEST(ComputeKernel, MixedPrecisionAgreesWithinStatistics) {
  auto c = ar_on_oxide(30.0, DamageCalculation::full_cascade, 1600);
  c.compute.precision = compute::Precision::mixed;
  const TrimTallies cpu = run(c, compute::Backend::cpu);
  const TrimTallies mix = run(c, compute::Backend::reference);
  const double n = static_cast<double>(cpu.ions);
  ASSERT_EQ(mix.ions, cpu.ions);
  const double sd = cpu.range.stddev();
  EXPECT_NEAR(mix.range.mean(), cpu.range.mean(), 3.0 * sd * std::sqrt(2.0 / n));
  EXPECT_NEAR(mix.total_vacancies / n, cpu.total_vacancies / n, 0.02 * cpu.total_vacancies / n);
  EXPECT_NEAR(mix.total_replacements / n, cpu.total_replacements / n, 0.05 * cpu.total_replacements / n);
  EXPECT_NEAR(mix.total_ionization_ion / n, cpu.total_ionization_ion / n, 0.01 * cpu.total_ionization_ion / n);
  // Not number for number: the single-precision angle changes some histories.
  EXPECT_NE(mix.total_vacancies, cpu.total_vacancies);
}

// A GPU framework present on this machine (AMCSET_TEST_GPU=opencl, hip, cuda
// or kokkos, optional device AMCSET_TEST_GPU_DEVICE): deterministic, and the
// same physics as the CPU within the statistics of 2000 ions.
TEST(ComputeGpu, AgreesWithTheCpu) {
  const char *which = std::getenv("AMCSET_TEST_GPU");
  if (which == nullptr || *which == '\0') GTEST_SKIP() << "set AMCSET_TEST_GPU=opencl|hip|cuda|kokkos";
  auto c = fe_on_fe(50.0, DamageCalculation::full_cascade, 2000);
  c.threads = 0;
  if (const char *d = std::getenv("AMCSET_TEST_GPU_DEVICE")) c.compute.device = std::atoi(d);
  c.compute.memory_mib = 128;
  const auto backend = compute::parse_backend(which);
  TrimSimulation::ComputeReport rep;
  const TrimTallies gpu = run(c, backend, &rep);
  const TrimTallies again = run(c, backend);
  const TrimTallies cpu = run(c, compute::Backend::cpu);
  std::cout << "[ device ] " << rep.device << ": " << rep.device_blocks << " blocks, " << rep.reruns
            << " reran on the CPU, " << rep.device_seconds << " s\n";
  expect_same(gpu, again); // reproducible
  ASSERT_EQ(gpu.ions, cpu.ions);
  // Mean range and its straggling, vacancies per ion: within 4 sigma.
  const double n = static_cast<double>(cpu.ions);
  const double sd = cpu.range.stddev();
  EXPECT_NEAR(gpu.range.mean(), cpu.range.mean(), 4.0 * sd * std::sqrt(2.0 / n));
  EXPECT_NEAR(gpu.range.stddev(), sd, 4.0 * sd / std::sqrt(n));
  const double vac_cpu = cpu.total_vacancies / n;
  const double vac_gpu = gpu.total_vacancies / n;
  EXPECT_NEAR(vac_gpu, vac_cpu, 0.05 * vac_cpu);
  EXPECT_NEAR(gpu.total_ionization_ion / n, cpu.total_ionization_ion / n, 0.03 * cpu.total_ionization_ion / n);
}
