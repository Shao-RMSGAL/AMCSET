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

/*! \file trim_tables.cpp
 *  \brief SRIM-mode results as output tables (TRIM's units).
 */

#include "trim/trim_tables.h"

#include "config.h"
#include "element_data.h"

#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <sstream>

namespace amcset {
namespace trim {

using namespace common;
using output::Column;
using output::Dimension;
using output::Scalar;
using output::Table;

namespace {

std::string key_of(const std::string &name, std::size_t index,
                   const std::vector<std::string> &all) {
  // Column-key-safe element name; repeated names get their index appended.
  std::string k;
  for (char c : name) {
    if (std::isalnum(static_cast<unsigned char>(c))) k += c;
    else k += '_';
  }
  if (k.empty()) k = "element";
  std::size_t same = 0;
  for (const auto &n : all) same += n == name ? 1 : 0;
  if (same > 1) k += "_" + std::to_string(index + 1);
  return k;
}

std::string iso_time_now() {
  const auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::tm tm{};
#if defined(_WIN32)
  gmtime_s(&tm, &t);
#else
  gmtime_r(&t, &tm);
#endif
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", &tm);
  return buf;
}

std::string fmt_g(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%g", v);
  return buf;
}

} // namespace

output::Dataset make_dataset(const TrimSimulation &sim) { return make_dataset(sim, sim.tallies()); }

output::Dataset make_dataset(const TrimSimulation &sim, const TrimTallies &t) {
  const TrimConfig &c = sim.config();
  const auto &names = sim.element_names();
  const auto &zs = sim.element_z();
  const std::size_t ne = names.size();
  const std::size_t nb = c.depth_bins;
  const double ions = std::max<double>(1.0, static_cast<double>(t.ions));
  const double w = sim.bin_width();
  const std::string ion_sym(Element::by_z(c.ion.z).symbol());
  const bool kp = c.damage == DamageCalculation::kinchin_pease;
  const bool recoil_source = c.source == TrimConfig::Source::particle_recoils;

  std::vector<std::string> keys;
  for (std::size_t i = 0; i < ne; ++i) keys.push_back(key_of(names[i], i, names));

  output::Dataset d;
  d.title = c.title;
  {
    std::ostringstream target;
    for (std::size_t i = 0; i < c.layers.size(); ++i) {
      const auto &l = c.layers[i];
      if (i) target << " + ";
      target << l.name << " (" << fmt_g(in_units(l.width, angstrom)) << " A, "
             << fmt_g(in_units(l.density, gram_per_cubic_centimeter)) << " g/cm3"
             << (l.gas ? ", gas" : "") << ")";
    }
    d.metadata = {
        {"program", std::string("AMCSET ") + AMCSET_VERSION + " (SRIM mode)"},
        {"created", iso_time_now()},
        {"ion", ion_sym + " (Z=" + std::to_string(c.ion.z) + ", " +
                    fmt_g((c.ion.mass / atomic_mass_unit).value()) + " amu)"},
        {"energy", fmt_g(in_units(c.ion.energy, kilo_electron_volt)) + " keV"},
        {"angle", fmt_g(in_units(c.ion.angle, degree)) + " deg"},
        {"target", target.str()},
        {"damage", to_string(c.damage)},
        {"cascade_code", std::to_string(c.trim_code())},
        {"source", c.source == TrimConfig::Source::beam
                       ? std::string("beam")
                       : (recoil_source ? "TRIM.DAT recoils" : "TRIM.DAT ions") +
                             (c.particles_title.empty() ? std::string() : ": " + c.particles_title)},
        {"transport", to_string(c.transport)},
        {"stopping", physics::to_string(c.stopping_model)},
        {"seed", std::to_string(c.seed)},
        {"ions", std::to_string(t.ions)},
        {"depth_bins", std::to_string(nb) + " bins of " + fmt_g(w) + " A from " +
                           fmt_g(sim.bin_start(0)) + " A"},
    };
  }

  // --- Scalars ---------------------------------------------------------------
  const double e0 = in_units(c.ion.energy, electron_volt);
  auto scalar = [&](std::string key, std::string label, Dimension dim, std::string unit, double v,
                    std::string desc = {}) {
    d.scalars.push_back(Scalar{std::move(key), std::move(label), dim, std::move(unit), v, std::move(desc)});
  };
  scalar("ions", "Incident particles", Dimension::dimensionless, "1", static_cast<double>(t.ions));
  scalar("range_mean", "Ion average range (projected)", Dimension::length, "A", t.range.mean());
  scalar("range_straggling", "Range straggling", Dimension::length, "A", t.range.stddev());
  scalar("range_skewness", "Range skewness", Dimension::dimensionless, "1", t.range.skewness());
  scalar("range_kurtosis", "Range kurtosis", Dimension::dimensionless, "1", t.range.kurtosis());
  scalar("lateral_range", "Lateral projected range", Dimension::length, "A", t.lateral_abs_y.mean());
  scalar("lateral_straggling", "Lateral straggling", Dimension::length, "A", t.lateral_y.rms());
  scalar("radial_range", "Radial range", Dimension::length, "A", t.radial.mean());
  scalar("radial_straggling", "Radial straggling", Dimension::length, "A", t.radial.stddev());
  scalar("stopped_ions", "Ions stopped in the target", Dimension::dimensionless, "1", t.range.n);
  scalar("backscattered_fraction", "Backscattered ions per ion", Dimension::dimensionless, "1",
         static_cast<double>(t.backscattered) / ions);
  scalar("backscattered_energy_fraction", "Energy backscattered (fraction of incident)",
         Dimension::dimensionless, "1", t.backscattered_energy / ions / e0);
  scalar("transmitted_fraction", "Transmitted ions per ion", Dimension::dimensionless, "1",
         static_cast<double>(t.transmitted) / ions);
  scalar("transmitted_energy_fraction", "Energy transmitted (fraction of incident)",
         Dimension::dimensionless, "1", t.transmitted_energy / ions / e0);
  scalar("vacancies_per_ion", "Vacancies per ion", Dimension::dimensionless, "1", t.total_vacancies / ions);
  // The quick calculation follows no recoils: every vacancy is a displaced
  // atom, there are no replacement collisions, and its displacement count
  // is the ion's own knock-ons.
  const bool quick = c.damage == DamageCalculation::kinchin_pease;
  scalar("displacements_per_ion", "Displacements per ion", Dimension::dimensionless, "1",
         (quick ? t.total_vacancies : t.total_displacements) / ions);
  if (quick) {
    scalar("knock_ons_per_ion", "Knock-ons (atoms displaced by the ion) per ion", Dimension::dimensionless,
           "1", t.total_displacements / ions);
  } else {
    scalar("replacements_per_ion", "Replacement collisions per ion", Dimension::dimensionless, "1",
           t.total_replacements / ions);
    scalar("interstitials_per_ion", "Interstitials per ion", Dimension::dimensionless, "1",
           t.total_interstitials / ions);
  }
  double sputtered = 0.0;
  for (double s : t.sputtered) sputtered += s;
  scalar("sputtering_yield", "Sputtered atoms per ion", Dimension::dimensionless, "1", sputtered / ions);
  scalar("ionization_ion_percent", "Energy loss: ionization by ions", Dimension::dimensionless, "%",
         t.total_ionization_ion / ions / e0, "fraction of the incident energy");
  scalar("ionization_recoils_percent", "Energy loss: ionization by recoils", Dimension::dimensionless,
         "%", t.total_ionization_recoils / ions / e0, "fraction of the incident energy");
  scalar("phonons_ion_percent", "Energy loss: phonons by ions", Dimension::dimensionless, "%",
         t.total_phonons_ion / ions / e0, "fraction of the incident energy");
  scalar("phonons_recoils_percent", "Energy loss: phonons by recoils", Dimension::dimensionless, "%",
         t.total_phonons_recoils / ions / e0, "fraction of the incident energy");
  scalar("vacancies_energy_percent", "Energy loss: vacancies (binding energy)", Dimension::dimensionless,
         "%", t.total_binding_energy_loss / ions / e0, "fraction of the incident energy");
  // The '%' unit stores fractions scaled to percent: convert the fractions.
  for (auto &s : d.scalars) {
    if (s.unit == "%") s.value *= 100.0;
  }
  scalar("elapsed_time", "Calculation time", Dimension::time, "s", sim.elapsed_seconds());

  // --- Depth tables ---------------------------------------------------------
  std::vector<double> depth(nb);
  for (std::size_t i = 0; i < nb; ++i) depth[i] = sim.bin_start(i) + w;
  const auto depth_col = [&]() {
    return Column::real("depth", "Depth", Dimension::length, "A", depth,
                        "end of the depth bin (TRIM's convention)");
  };
  const auto scaled = [&](const DepthHistogram &h, double norm) {
    std::vector<double> v(nb, 0.0);
    for (std::size_t i = 0; i < nb && i < h.values.size(); ++i) v[i] = h.values[i] * norm;
    return v;
  };
  const double per_cm = 1.0 / (ions * in_units(w * angstrom, centimeter)); // (atoms/cm3)/(atoms/cm2)
  const double per_A = 1.0 / (ions * w);

  {
    Table tb;
    tb.key = "range";
    tb.title = "Ion and recoil atom distributions";
    tb.description = "Final positions of the ions and of the target atoms that came to rest "
                     "(interstitials) per unit depth and per incident ion; multiplied by the "
                     "fluence (atoms/cm2) this gives the concentration (atoms/cm3).";
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("ions", ion_sym + " ions", Dimension::inverse_length, "1/cm",
                                      scaled(t.ion_stopped, per_cm)));
    for (std::size_t e = 0; e < ne; ++e) {
      tb.columns.push_back(Column::real("recoils_" + keys[e], names[e] + " recoils",
                                        Dimension::inverse_length, "1/cm",
                                        scaled(t.recoil_stopped[e], per_cm)));
    }
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "lateral";
    tb.title = "Lateral spread of the ions";
    tb.description = "For the ions stopped in each depth bin: mean |y| (lateral projected range) "
                     "and its straggling sqrt(<y^2>), mean radial distance from the beam axis and "
                     "its standard deviation.";
    std::vector<double> lr(nb), ls(nb), rr(nb), rs(nb);
    for (std::size_t i = 0; i < nb && i < t.lateral_bins.size(); ++i) {
      const auto &b = t.lateral_bins[i];
      if (b[0] > 0) {
        lr[i] = b[1] / b[0];
        ls[i] = std::sqrt(b[2] / b[0]);
        rr[i] = b[3] / b[0];
        rs[i] = std::sqrt(std::max(0.0, b[4] / b[0] - rr[i] * rr[i]));
      }
    }
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("lateral_range", "Lateral projected range", Dimension::length, "A", lr));
    tb.columns.push_back(Column::real("lateral_straggling", "Projected straggling", Dimension::length, "A", ls));
    tb.columns.push_back(Column::real("radial_range", "Radial range", Dimension::length, "A", rr));
    tb.columns.push_back(Column::real("radial_straggling", "Radial straggling", Dimension::length, "A", rs));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "vacancies";
    tb.title = "Vacancy production";
    tb.columns.push_back(depth_col());
    std::vector<double> total(nb, 0.0);
    if (kp) {
      tb.description = "Vacancies per unit depth and per ion, from the quick (Kinchin-Pease) "
                       "estimate: those made by the ion (primary knock-ons) and by its recoils.";
      auto a = scaled(t.vacancies_by_ion, per_A);
      auto b = scaled(t.vacancies_by_recoils, per_A);
      for (std::size_t i = 0; i < nb; ++i) total[i] = a[i] + b[i];
      tb.columns.push_back(Column::real("by_ions", "Vacancies by ions", Dimension::inverse_length, "1/A", a));
      tb.columns.push_back(Column::real("by_recoils", "Vacancies by recoils", Dimension::inverse_length, "1/A", b));
    } else {
      tb.description = "Primary knock-ons (target atoms displaced by the ion) and vacancies of "
                       "each target element per unit depth and per ion.";
      tb.columns.push_back(Column::real("knock_ons", "Knock-ons", Dimension::inverse_length, "1/A",
                                        scaled(t.knock_ons, per_A)));
      for (std::size_t e = 0; e < ne; ++e) {
        auto v = scaled(t.vacancies[e], per_A);
        for (std::size_t i = 0; i < nb; ++i) total[i] += v[i];
        tb.columns.push_back(Column::real("vacancies_" + keys[e], names[e] + " vacancies",
                                          Dimension::inverse_length, "1/A", std::move(v)));
      }
    }
    tb.columns.push_back(Column::real("total", "Total vacancies", Dimension::inverse_length, "1/A", total));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "replacements";
    tb.title = "Replacement collisions";
    tb.description = "Replacement collisions per unit depth and per ion.";
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("replacements", "Replacement collisions", Dimension::inverse_length,
                                      "1/A", scaled(t.replacements, per_A)));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "damage";
    tb.title = "Collision events";
    tb.description = "Displacements (= vacancies + replacement collisions), vacancies, "
                     "replacement collisions and interstitials per unit depth and per ion.";
    std::vector<double> vac(nb, 0.0);
    if (kp) {
      auto a = scaled(t.vacancies_by_ion, per_A);
      auto b = scaled(t.vacancies_by_recoils, per_A);
      for (std::size_t i = 0; i < nb; ++i) vac[i] = a[i] + b[i];
    } else {
      for (std::size_t e = 0; e < ne; ++e) {
        auto v = scaled(t.vacancies[e], per_A);
        for (std::size_t i = 0; i < nb; ++i) vac[i] += v[i];
      }
    }
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("displacements", "Displacements", Dimension::inverse_length, "1/A",
                                      scaled(t.displacements, per_A)));
    tb.columns.push_back(Column::real("vacancies", "Vacancies", Dimension::inverse_length, "1/A", vac));
    tb.columns.push_back(Column::real("replacements", "Replacement collisions", Dimension::inverse_length,
                                      "1/A", scaled(t.replacements, per_A)));
    tb.columns.push_back(Column::real("interstitials", "Interstitials", Dimension::inverse_length, "1/A",
                                      scaled(t.interstitials, per_A)));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "ionization";
    tb.title = "Ionization (energy loss to target electrons)";
    tb.description = "Energy given to electrons per unit depth and per ion, by the ions and by the recoils.";
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("by_ions", "Ionization by ions", Dimension::energy_per_length, "eV/A",
                                      scaled(t.ionization_by_ion, per_A)));
    tb.columns.push_back(Column::real("by_recoils", "Ionization by recoils", Dimension::energy_per_length,
                                      "eV/A", scaled(t.ionization_by_recoils, per_A)));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "phonons";
    tb.title = "Phonons (energy loss to target atoms below the displacement threshold)";
    tb.description = "Energy left as lattice vibrations per unit depth and per ion.";
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("by_ions", "Phonons by ions", Dimension::energy_per_length, "eV/A",
                                      scaled(t.phonons_by_ion, per_A)));
    tb.columns.push_back(Column::real("by_recoils", "Phonons by recoils", Dimension::energy_per_length,
                                      "eV/A", scaled(t.phonons_by_recoils, per_A)));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "energy_to_recoils";
    tb.title = "Energy transferred to recoils";
    tb.description = "Energy the ions give to target atoms that are displaced, and the energy "
                     "absorbed by each target element, per unit depth and per ion.";
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("from_ions", "Energy from ions", Dimension::energy_per_length, "eV/A",
                                      scaled(t.energy_to_recoils, per_A)));
    for (std::size_t e = 0; e < ne; ++e) {
      tb.columns.push_back(Column::real("absorbed_" + keys[e], "Absorbed by " + names[e],
                                        Dimension::energy_per_length, "eV/A",
                                        scaled(t.energy_absorbed[e], per_A)));
    }
    d.tables.push_back(std::move(tb));
  }

  // --- Per-element results --------------------------------------------------
  {
    Table tb;
    tb.key = "sputtering";
    tb.title = "Sputtering and transmitted recoils by target element";
    std::vector<std::string> el;
    std::vector<double> z, yield, mean_e, trans;
    for (std::size_t e = 0; e < ne; ++e) {
      el.push_back(names[e]);
      z.push_back(zs[e]);
      yield.push_back(t.sputtered[e] / ions);
      mean_e.push_back(t.sputtered[e] > 0 ? t.sputtered_energy[e] / t.sputtered[e] : 0.0);
      trans.push_back(t.transmitted_recoils[e] / ions);
    }
    tb.columns.push_back(Column::strings("element", "Element", el));
    tb.columns.push_back(Column::integer("z", "Z", z));
    tb.columns.push_back(Column::real("yield", "Sputtering yield", Dimension::dimensionless, "1", yield,
                                      "sputtered atoms per incident ion"));
    tb.columns.push_back(Column::real("mean_energy", "Mean energy of sputtered atoms", Dimension::energy, "eV",
                                      mean_e, "after the surface barrier"));
    tb.columns.push_back(Column::real("transmitted_recoils", "Transmitted recoils", Dimension::dimensionless,
                                      "1", trans, "recoil atoms leaving through the back surface, per ion"));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "surface_energy";
    tb.title = "Atoms reaching the surface, by normal energy";
    tb.description = "Target atoms reaching the front surface from inside, per incident ion, by "
                     "the energy of their motion normal to the surface (E cos^2 of the angle to "
                     "the normal) in 0.1 eV bins (the last bin collects everything above). Atoms "
                     "above the surface binding energy are sputtered. 'yield_above' is the "
                     "number per ion at or above the bin, i.e. the sputtering yield the target "
                     "would have with that surface binding energy.";
    const std::size_t m = TrimTallies::surface_energy_bins;
    std::vector<double> e(m);
    for (std::size_t k = 0; k < m; ++k) e[k] = static_cast<double>(k + 1) * TrimTallies::surface_energy_bin_eV;
    tb.columns.push_back(Column::real("energy", "Normal energy", Dimension::energy, "eV", e));
    for (std::size_t el = 0; el < ne; ++el) {
      std::vector<double> diff(m, 0.0), integral(m, 0.0);
      if (el < t.surface_energy.size()) {
        for (std::size_t k = 0; k < m && k < t.surface_energy[el].size(); ++k) diff[k] = t.surface_energy[el][k] / ions;
      }
      double acc = 0.0;
      for (std::size_t k = m; k-- > 0;) {
        acc += diff[k];
        integral[k] = acc;
      }
      tb.columns.push_back(Column::real("arrivals_" + keys[el], names[el] + " atoms per bin",
                                        Dimension::dimensionless, "1", std::move(diff)));
      tb.columns.push_back(Column::real("yield_above_" + keys[el], names[el] + " at or above",
                                        Dimension::dimensionless, "1", std::move(integral)));
    }
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "recoil_spectrum";
    tb.title = "Energy spectrum of displaced atoms";
    tb.description = "Displaced target atoms per ion by recoil energy (transferred energy minus the "
                     "lattice binding energy), in logarithmic bins (10 per decade): primary "
                     "knock-ons made by the ion, and all displaced atoms (full cascades only).";
    const int ppd = TrimTallies::spectrum_bins_per_decade;
    std::vector<double> lo, hi, pka, all;
    for (std::size_t i = 0; i < t.pka_spectrum.size(); ++i) {
      lo.push_back(std::pow(10.0, static_cast<double>(i) / ppd));
      hi.push_back(std::pow(10.0, static_cast<double>(i + 1) / ppd));
      pka.push_back(t.pka_spectrum[i] / ions);
      all.push_back(i < t.recoil_spectrum.size() ? t.recoil_spectrum[i] / ions : 0.0);
    }
    tb.columns.push_back(Column::real("energy_low", "Energy from", Dimension::energy, "eV", lo));
    tb.columns.push_back(Column::real("energy_high", "Energy to", Dimension::energy, "eV", hi));
    tb.columns.push_back(Column::real("primary_knock_ons", "Primary knock-ons", Dimension::dimensionless, "1", pka));
    tb.columns.push_back(Column::real("all_recoils", "All displaced atoms", Dimension::dimensionless, "1", all));
    d.tables.push_back(std::move(tb));
  }

  if (c.transport == TrimConfig::Transport::trim2013 && !t.nrt_displacements.values.empty()) {
    // Displacements per atom by depth, by several damage models (docs/OUTPUT.md "dpa").
    Table tb;
    tb.key = "dpa";
    tb.title = "Displacements per atom (dpa)";
    char fluence[64];
    std::snprintf(fluence, sizeof fluence, "%g", c.dpa_fluence_per_cm2);
    tb.description =
        std::string("dpa at a fluence of ") + fluence +
        " ions/cm2 (output.dpa_fluence; dpa scale linearly with it), from the displacements per ion "
        "and unit depth divided by the layer's atomic density. trim_vacancies: TRIM's vacancies (" +
        (kp ? "quick calculation: the NRT rule with TRIM's damage energy" : "full cascades") +
        "). The others apply a damage function to every primary knock-on atom, at the depth where it "
        "starts, with its damage energy T_d (TRIM's Robinson partition): nrt (Norgett-Robinson-Torrens, "
        "0.8 T_d/(2 Ed)), kinchin_pease (T_d/(2 Ed)), arc_dpa (NRT times the athermal-recombination "
        "efficiency, Nordlund et al. 2018) and rpa (replaced atoms per atom, the same paper). "
        "damage_energy: the primaries' T_d per ion and depth.";
    std::vector<double> td, trimv, nrt, kpo, arc, rpa;
    // Atomic density (per cubic Angstrom) of the layer at each bin centre.
    const auto &layers = sim.effective_layers();
    std::vector<double> edge{0.0};
    std::vector<double> density;
    for (const auto &l : layers) {
      edge.push_back(edge.back() + in_units(l.width, angstrom));
      double total = 0.0;
      mass_quantity mean = 0.0 * kilogram;
      for (const auto &e : l.elements) total += e.stoichiometry;
      for (const auto &e : l.elements) mean += (e.stoichiometry / total) * e.mass;
      density.push_back(in_units(l.density / mean, per_cubic_angstrom));
    }
    const double fluence_per_A2 = c.dpa_fluence_per_cm2 * 1e-16; // 1 cm2 = 1e16 A2
    for (std::size_t i = 0; i < nb; ++i) {
      const double centre = sim.bin_start(i) + 0.5 * w;
      std::size_t li = 0;
      while (li + 1 < density.size() && centre >= edge[li + 1]) ++li;
      const double scale = fluence_per_A2 / (ions * w * density[li]); // dpa per event in the bin
      double vac = 0.0;
      for (const auto &h : t.vacancies) vac += i < h.values.size() ? h.values[i] : 0.0;
      const auto v = [i](const DepthHistogram &h) { return i < h.values.size() ? h.values[i] : 0.0; };
      trimv.push_back(vac * scale);
      nrt.push_back(v(t.nrt_displacements) * scale);
      kpo.push_back(v(t.kp_displacements) * scale);
      arc.push_back(v(t.arc_displacements) * scale);
      rpa.push_back(v(t.rpa_replacements) * scale);
      td.push_back(v(t.pka_damage_energy) / (ions * w));
    }
    tb.columns.push_back(depth_col());
    tb.columns.push_back(Column::real("trim_vacancies", kp ? "TRIM vacancies (quick, = NRT)" : "TRIM vacancies (full cascades)",
                                      Dimension::dimensionless, "1", trimv));
    tb.columns.push_back(Column::real("nrt", "NRT-dpa", Dimension::dimensionless, "1", nrt));
    tb.columns.push_back(Column::real("kinchin_pease", "Kinchin-Pease", Dimension::dimensionless, "1", kpo));
    // arc-dpa and rpa need the constants of every primary's element.
    if (t.pka_without_efficiency == 0.0) {
      tb.columns.push_back(Column::real("arc_dpa", "arc-dpa", Dimension::dimensionless, "1", arc));
      tb.columns.push_back(Column::real("rpa", "Replacements per atom (rpa)", Dimension::dimensionless, "1", rpa));
    } else {
      tb.description += " arc_dpa and rpa are left out: some primary knock-on atoms are of elements without "
                        "arc-dpa/rpa constants (layers.N.elements.M.damage_efficiency).";
    }
    tb.columns.push_back(Column::real("damage_energy", "Damage energy of the primaries", Dimension::energy_per_length,
                                      "eV/A", td));
    d.tables.push_back(std::move(tb));
  }
  if (!t.md.empty()) {
    // Hybrid BCA-MD events (TrimConfig::md_events), one row per E_c.
    Table tb;
    tb.key = "md_events_summary";
    tb.title = "Recoil events for MD";
    tb.description = std::string("Per ion, for each transition energy E_c: events (recoils set in motion below "
                                 "E_c, and moving atoms whose energy fell below E_c), the energy they carry, "
                                 "and recoils set in motion below E_min (left to the BCA). Mode: ") +
                     (c.md_events.handoff ? "handoff (the BCA did not follow the events)"
                                          : "record (the transport and tallies are unchanged)") +
                     ". The events themselves are in the MD event file.";
    std::vector<double> ec, ev, born, cross, energy, share, below, below_e, handed;
    const double ion_eV = in_units(c.ion.energy, electron_volt);
    for (const auto &m : t.md) {
      ec.push_back(m.transition_energy_eV);
      ev.push_back(m.events / ions);
      born.push_back(m.birth_events / ions);
      cross.push_back(m.crossing_events / ions);
      energy.push_back(m.energy_eV / ions);
      share.push_back(100.0 * m.energy_eV / ions / ion_eV);
      below.push_back(m.below_min / ions);
      below_e.push_back(m.below_min_energy_eV / ions);
      handed.push_back(m.ions_handed);
    }
    tb.columns.push_back(Column::real("E_c", "Transition energy E_c", Dimension::energy, "eV", ec));
    tb.columns.push_back(Column::real("events", "Events per ion", Dimension::dimensionless, "1", ev));
    tb.columns.push_back(Column::real("birth_events", "Born below E_c, per ion", Dimension::dimensionless, "1", born));
    tb.columns.push_back(Column::real("crossing_events", "Fell below E_c, per ion", Dimension::dimensionless, "1", cross));
    tb.columns.push_back(Column::real("energy", "Energy to MD per ion", Dimension::energy, "eV", energy));
    tb.columns.push_back(Column::real("energy_percent", "Energy to MD", Dimension::dimensionless, "%", share));
    tb.columns.push_back(Column::real("below_min", "Recoils below E_min per ion", Dimension::dimensionless, "1", below));
    tb.columns.push_back(Column::real("below_min_energy", "Their energy per ion", Dimension::energy, "eV", below_e));
    tb.columns.push_back(Column::real("ions_handed", "Ions that became events", Dimension::dimensionless, "1", handed));
    d.tables.push_back(std::move(tb));
  }

  // --- Per-particle lists ---------------------------------------------------
  {
    Table tb;
    tb.key = "exits";
    tb.title = "Atoms leaving the target";
    tb.description = "Backscattered (B) and transmitted (T) ions and sputtered (S) atoms: ion "
                     "number, atom, energy, last position and direction cosines (TRIM's "
                     "BACKSCAT/TRANSMIT/SPUTTER/TRIMOUT files; the same columns as TRIM.DAT). " +
                     std::string(c.transport == TrimConfig::Transport::trim2013
                                     ? "As in TRIM, sputtered atoms are listed with their energy "
                                       "and direction inside the surface (Es not subtracted)."
                                     : "Sputtered atoms are listed after the surface barrier.");
    std::vector<std::string> kind, el;
    std::vector<double> ion, z, e, x, y, zz, cx, cy, cz;
    for (const auto &r : t.exits) {
      kind.emplace_back(1, r.kind);
      ion.push_back(static_cast<double>(r.ion));
      z.push_back(r.z);
      el.emplace_back(Element::by_z(r.z).symbol());
      e.push_back(r.energy_eV);
      x.push_back(r.x);
      y.push_back(r.y);
      zz.push_back(r.z_pos);
      cx.push_back(r.cx);
      cy.push_back(r.cy);
      cz.push_back(r.cz);
    }
    tb.columns.push_back(Column::strings("kind", "Kind", kind, "B backscattered, T transmitted, S sputtered"));
    tb.columns.push_back(Column::integer("ion", "Ion number", ion));
    tb.columns.push_back(Column::integer("z", "Z", z));
    tb.columns.push_back(Column::strings("element", "Element", el));
    tb.columns.push_back(Column::real("energy", "Energy", Dimension::energy, "eV", e));
    tb.columns.push_back(Column::real("x", "Depth X", Dimension::length, "A", x));
    tb.columns.push_back(Column::real("y", "Lateral Y", Dimension::length, "A", y));
    tb.columns.push_back(Column::real("z_position", "Lateral Z", Dimension::length, "A", zz));
    tb.columns.push_back(Column::real("cos_x", "cos(X)", Dimension::dimensionless, "1", cx));
    tb.columns.push_back(Column::real("cos_y", "cos(Y)", Dimension::dimensionless, "1", cy));
    tb.columns.push_back(Column::real("cos_z", "cos(Z)", Dimension::dimensionless, "1", cz));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "ion_positions";
    tb.title = "Final positions of the stopped ions";
    tb.description = "Depth X and lateral Y, Z of every ion that stopped in the target (RANGE_3D.txt).";
    std::vector<double> ion, x, y, z;
    for (const auto &p : t.positions) {
      ion.push_back(static_cast<double>(p.ion));
      x.push_back(p.x);
      y.push_back(p.y);
      z.push_back(p.z);
    }
    tb.columns.push_back(Column::integer("ion", "Ion number", ion));
    tb.columns.push_back(Column::real("x", "Depth X", Dimension::length, "A", x));
    tb.columns.push_back(Column::real("y", "Lateral Y", Dimension::length, "A", y));
    tb.columns.push_back(Column::real("z", "Lateral Z", Dimension::length, "A", z));
    d.tables.push_back(std::move(tb));
  }
  if (c.exyz_interval > 0.0 * electron_volt || !t.exyz.empty()) {
    Table tb;
    tb.key = "exyz";
    tb.title = "Ion energy versus position";
    tb.description = "The ion's position, electronic stopping and the energy given to the last "
                     "struck atom each time its energy has dropped by the EXYZ interval, at its "
                     "start and where it stops or leaves the target (EXYZ.txt).";
    std::vector<double> ion, e, x, y, z, se, last;
    for (const auto &r : t.exyz) {
      ion.push_back(static_cast<double>(r.ion));
      e.push_back(r.energy_keV);
      x.push_back(r.x);
      y.push_back(r.y);
      z.push_back(r.z);
      se.push_back(r.se_eV_per_A);
      last.push_back(r.last_recoil_eV);
    }
    tb.columns.push_back(Column::integer("ion", "Ion number", ion));
    tb.columns.push_back(Column::real("energy", "Energy", Dimension::energy, "keV", e));
    tb.columns.push_back(Column::real("x", "Depth X", Dimension::length, "A", x));
    tb.columns.push_back(Column::real("y", "Lateral Y", Dimension::length, "A", y));
    tb.columns.push_back(Column::real("z", "Lateral Z", Dimension::length, "A", z));
    tb.columns.push_back(Column::real("electronic_stopping", "Electronic stopping",
                                      Dimension::energy_per_length, "eV/A", se));
    tb.columns.push_back(Column::real("last_recoil_energy", "Energy to last recoil", Dimension::energy,
                                      "eV", last));
    d.tables.push_back(std::move(tb));
  }
  if (c.collision_details_ions > 0 || !t.cascades.empty()) {
    Table tb;
    tb.key = "collisions";
    tb.title = "Ion collisions that displaced target atoms";
    tb.description = "Every collision of an ion that displaced a target atom, with the summary of "
                     "the cascade it started (COLLISON.txt).";
    std::vector<double> ion, e, x, y, z, se, zh, t_eV, disp, vac, repl, inter;
    std::vector<std::string> atom;
    for (const auto &r : t.cascades) {
      ion.push_back(static_cast<double>(r.ion));
      e.push_back(r.ion_energy_keV);
      x.push_back(r.x);
      y.push_back(r.y);
      z.push_back(r.z);
      se.push_back(r.se_eV_per_A);
      zh.push_back(r.z_hit);
      atom.emplace_back(Element::by_z(r.z_hit).symbol());
      t_eV.push_back(r.recoil_energy_eV);
      disp.push_back(r.displacements);
      vac.push_back(r.vacancies);
      repl.push_back(r.replacements);
      inter.push_back(r.interstitials);
    }
    tb.columns.push_back(Column::integer("ion", "Ion number", ion));
    tb.columns.push_back(Column::real("ion_energy", "Ion energy", Dimension::energy, "keV", e));
    tb.columns.push_back(Column::real("x", "Depth X", Dimension::length, "A", x));
    tb.columns.push_back(Column::real("y", "Lateral Y", Dimension::length, "A", y));
    tb.columns.push_back(Column::real("z", "Lateral Z", Dimension::length, "A", z));
    tb.columns.push_back(Column::real("electronic_stopping", "Electronic stopping",
                                      Dimension::energy_per_length, "eV/A", se));
    tb.columns.push_back(Column::strings("atom", "Atom hit", atom));
    tb.columns.push_back(Column::integer("atom_z", "Z of atom hit", zh));
    tb.columns.push_back(Column::real("recoil_energy", "Recoil energy", Dimension::energy, "eV", t_eV));
    tb.columns.push_back(Column::real("displacements", "Displacements", Dimension::dimensionless, "1", disp));
    tb.columns.push_back(Column::real("vacancies", "Vacancies", Dimension::dimensionless, "1", vac));
    tb.columns.push_back(Column::real("replacements", "Replacements", Dimension::dimensionless, "1", repl));
    tb.columns.push_back(Column::real("interstitials", "Interstitials", Dimension::dimensionless, "1", inter));
    d.tables.push_back(std::move(tb));
  }
  if (!t.ion_summaries.empty()) {
    Table tb;
    tb.key = "ion_damage";
    tb.title = "Damage per ion";
    tb.description = "Displacements, replacements, vacancies, interstitials and sputtered and "
                     "transmitted target atoms of each ion with collision details (COLLISON.txt).";
    std::vector<double> ion, disp, repl, vac, inter, sput, trans;
    for (const auto &r : t.ion_summaries) {
      ion.push_back(static_cast<double>(r.ion));
      disp.push_back(r.displacements);
      repl.push_back(r.replacements);
      vac.push_back(r.vacancies);
      inter.push_back(r.interstitials);
      sput.push_back(r.sputtered);
      trans.push_back(r.transmitted);
    }
    tb.columns.push_back(Column::integer("ion", "Ion number", ion));
    tb.columns.push_back(Column::real("displacements", "Displacements", Dimension::dimensionless, "1", disp));
    tb.columns.push_back(Column::real("replacements", "Replacements", Dimension::dimensionless, "1", repl));
    tb.columns.push_back(Column::real("vacancies", "Vacancies", Dimension::dimensionless, "1", vac));
    tb.columns.push_back(Column::real("interstitials", "Interstitials", Dimension::dimensionless, "1", inter));
    tb.columns.push_back(Column::real("sputtered", "Sputtered atoms", Dimension::dimensionless, "1", sput));
    tb.columns.push_back(Column::real("transmitted", "Transmitted atoms", Dimension::dimensionless, "1", trans));
    d.tables.push_back(std::move(tb));
  }
  bool any_atoms = false;
  for (const auto &r : t.cascades) any_atoms = any_atoms || !r.atoms.empty();
  if (any_atoms) {
    Table tb;
    tb.key = "cascade_atoms";
    tb.title = "Displaced and replacing atoms of the cascades";
    tb.description = "Every displacement (vacancy = 1, with the atom's starting energy) and "
                     "replacement (replacement = 1) of the cascades listed in the collisions "
                     "table, in the order they happen (COLLISON.txt with recoils).";
    std::vector<double> ion, cascade, num, zz, e, x, y, z, vac, repl;
    std::vector<std::string> atom;
    for (std::size_t k = 0; k < t.cascades.size(); ++k) {
      for (const auto &a : t.cascades[k].atoms) {
        ion.push_back(static_cast<double>(t.cascades[k].ion));
        cascade.push_back(static_cast<double>(k + 1));
        num.push_back(a.number);
        zz.push_back(a.z);
        atom.emplace_back(Element::by_z(a.z).symbol());
        e.push_back(a.energy_eV);
        x.push_back(a.x);
        y.push_back(a.y);
        z.push_back(a.z_pos);
        vac.push_back(a.vacancy ? 1.0 : 0.0);
        repl.push_back(a.replacement ? 1.0 : 0.0);
      }
    }
    tb.columns.push_back(Column::integer("ion", "Ion number", ion));
    tb.columns.push_back(Column::integer("cascade", "Cascade (row of the collisions table)", cascade));
    tb.columns.push_back(Column::integer("number", "Number in the cascade", num));
    tb.columns.push_back(Column::strings("atom", "Atom", atom));
    tb.columns.push_back(Column::integer("atom_z", "Z of the atom", zz));
    tb.columns.push_back(Column::real("energy", "Energy", Dimension::energy, "eV", e));
    tb.columns.push_back(Column::real("x", "Depth X", Dimension::length, "A", x));
    tb.columns.push_back(Column::real("y", "Lateral Y", Dimension::length, "A", y));
    tb.columns.push_back(Column::real("z", "Lateral Z", Dimension::length, "A", z));
    tb.columns.push_back(Column::integer("vacancy", "Vacancy", vac));
    tb.columns.push_back(Column::integer("replacement", "Replacement", repl));
    d.tables.push_back(std::move(tb));
  }
  return d;
}

} // namespace trim
} // namespace amcset
