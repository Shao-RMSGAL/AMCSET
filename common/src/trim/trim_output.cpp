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

/*! \file trim_output.cpp
 *  \brief TRIM-format output files and the JSON summary of a SRIM-mode run.
 *
 * The files mirror the layout of SRIM-2013's output (same file names, header
 * lines, column order and units) so that existing post-processing scripts
 * written for TRIM output keep working. Numbers are printed in plain
 * scientific notation.
 */

#include "trim/trim.h"
#include "trim/single_precision.h"
#include "trim/trim_tables.h"
#include "output/text.h"
#include "output/writers.h"

#include "config.h"
#include "element_data.h"
#include "parallel_tasks.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace trim {

using namespace common;

namespace {

// The helpers below forward printf-style formats that are always string
// literals at the call sites.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif
std::string fmt(const char *format, double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, format, v);
  return buf;
}

std::string fmt(const char *format, const char *text) {
  char buf[128];
  std::snprintf(buf, sizeof buf, format, text);
  return buf;
}

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

std::string sci(double v) { return fmt("%11.4E", v); }

std::string ion_symbol(const TrimConfig &c) {
  return std::string(Element::by_z(c.ion.z).symbol());
}

std::string target_label(const TrimConfig &c) {
  std::ostringstream s;
  for (std::size_t i = 0; i < c.layers.size(); ++i) {
    if (i) s << "+";
    s << c.layers[i].name;
  }
  return s.str();
}

std::ofstream open_file(const std::filesystem::path &p) {
  std::ofstream out(p, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write " + p.string());
  return out;
}

class Writer {
public:
  Writer(const TrimSimulation &sim) : sim_(sim), c_(sim.config()), t_(sim.tallies()) {}

  //! One task per file (the files are independent; run_tasks writes them
  //! concurrently, which matters for the long particle lists).
  std::vector<std::function<void()>> tasks(const std::filesystem::path &dir) const {
    std::filesystem::create_directories(dir);
    std::vector<std::function<void()>> t;
    t.push_back([this, dir] { write_range(dir / "RANGE.txt"); });
    t.push_back([this, dir] { write_vacancy(dir / "VACANCY.txt"); });
    t.push_back([this, dir] { write_novac(dir / "NOVAC.txt"); });
    t.push_back([this, dir] {
      write_two_column(dir / "IONIZ.txt", "Ion and Recoil IONIZATION",
                       "Ionization Energy Units are  >>>>  eV /(Angstrom-Ion)  <<<<",
                       {"  TARGET       IONIZ.       IONIZ.  ", "  DEPTH         by           by     ",
                        "  (Ang)        IONS         RECOILS "},
                       t_.ionization_by_ion, t_.ionization_by_recoils);
    });
    t.push_back([this, dir] {
      write_two_column(dir / "PHONON.txt", "Ion and Recoil PHONONS",
                       "Table Units are  >>>> Phonons /(Angstrom-Ion)  <<<<",
                       {"  DEPTH       PHONONS     PHONONS    ", "  (Ang.)      by ION      by RECOILS "},
                       t_.phonons_by_ion, t_.phonons_by_recoils);
    });
    t.push_back([this, dir] { write_e2recoil(dir / "E2RECOIL.txt"); });
    t.push_back([this, dir] { write_lateral(dir / "LATERAL.txt"); });
    t.push_back([this, dir] { write_tdata(dir / "TDATA.txt"); });
    t.push_back([this, dir] { write_exits(dir / "BACKSCAT.txt", 'B', "File of Backscattered Ions"); });
    t.push_back([this, dir] { write_exits(dir / "TRANSMIT.txt", 'T', "File of Transmitted Ions/Recoils"); });
    t.push_back([this, dir] { write_exits(dir / "SPUTTER.txt", 'S', "File of Sputtered Atoms"); });
    t.push_back([this, dir] {
      write_exits(dir / "TRIMOUT.txt", 0, "File of Transmitted / Backscattered / Sputtered Atoms");
    });
    if (c_.exyz_interval > 0.0 * electron_volt) t.push_back([this, dir] { write_exyz(dir / "EXYZ.txt"); });
    if (c_.record_ion_positions) t.push_back([this, dir] { write_range3d(dir / "RANGE_3D.txt"); });
    if (c_.collision_details_ions > 0) t.push_back([this, dir] { write_collisions(dir / "COLLISON.txt"); });
    t.push_back([this, dir] { write_spectrum(dir / "PKA_SPECTRUM.txt"); });
    return t;
  }

private:
  double ions() const { return std::max<double>(1.0, static_cast<double>(t_.ions)); }
  double width() const { return sim_.bin_width(); }
  //! A count-table bin as TRIM.exe prints it: summed in single precision in
  //! steps of 1 / (bin width) (see trim/single_precision.h), unless the
  //! configuration asks for exact tables. \p removed of \p count were
  //! subtracted again (replacements).
  double counted(double count, double removed = 0.0) const {
    return c_.single_precision_tables ? trim_table_count(count, removed, width()) : count;
  }
  double depth_label(std::size_t i) const { return sim_.bin_start(i) + width(); }

  void banner(std::ostream &o, const std::string &what) const {
    o << "====== " << c_.title << " =======\n";
    o << "        AMCSET " << AMCSET_VERSION << " (SRIM mode, TRIM-compatible output)\n";
    o << "==================================================\n";
    o << "      " << what << "\n";
    o << " See TDATA.txt for calculation details\n";
    o << "==================================================\n";
    o << "====== TRIM Calc.= " << ion_symbol(c_) << "("
      << fmt("%g", in_units(c_.ion.energy, kilo_electron_volt)) << " keV) ==> "
      << target_label(c_) << "(" << fmt("%7.0f", in_units(c_.target_thickness(), angstrom))
      << " A) ===\n";
    o << "------------------------------------------------------------------------------\n";
    o << " Ion    = " << ion_symbol(c_) << "   Energy = "
      << fmt("%g", in_units(c_.ion.energy, kilo_electron_volt)) << "  keV\n";
    target(o);
    o << " Total Ions calculated =" << fmt("%.2f", static_cast<double>(t_.ions)) << "\n";
  }

  void target(std::ostream &o) const {
    o << "============= TARGET MATERIAL ======================================\n";
    for (std::size_t i = 0; i < c_.layers.size(); ++i) {
      const auto &l = c_.layers[i];
      double total = 0.0;
      mass_quantity mean_mass = 0.0 * kilogram;
      for (const auto &e : l.elements) total += e.stoichiometry;
      for (const auto &e : l.elements) mean_mass += (e.stoichiometry / total) * e.mass;
      const number_density_quantity n = l.density / mean_mass;
      o << "Layer  " << i + 1 << " : " << l.name << "\n";
      o << "Layer Width = " << fmt("%10.3E", in_units(l.width, angstrom)) << " A   Layer # "
        << i + 1 << "- Density = " << fmt("%.4G", in_units(n, per_cubic_centimeter))
        << " atoms/cm3 = " << fmt("%.4f", in_units(l.density, gram_per_cubic_centimeter))
        << " g/cm3\n";
      for (const auto &e : l.elements) {
        const double at = 100.0 * e.stoichiometry / total;
        const double ms = 100.0 * (e.stoichiometry / total) * (e.mass / mean_mass).value();
        o << "  Layer # " << i + 1 << "- " << Element::by_z(e.z).symbol() << " = "
          << fmt("%.4g", at) << "  Atomic Percent = " << fmt("%.4g", ms) << "  Mass Percent\n";
      }
    }
    o << "====================================================================\n";
  }

  void write_range(const std::filesystem::path &p) const {
    auto o = open_file(p);
    banner(o, "ION and final RECOIL ATOM Distributions");
    o << " Ion Average Range = " << sci(t_.range.mean()) << " A   Straggling = "
      << sci(t_.range.stddev()) << " A\n";
    o << " Ion Lateral Range = " << sci(t_.lateral_abs_y.mean()) << " A   Straggling = "
      << sci(t_.lateral_y.rms()) << " A\n";
    o << " Ion Radial  Range = " << sci(t_.radial.mean()) << " A   Straggling = "
      << sci(t_.radial.stddev()) << " A\n";
    o << "====================================================================\n";
    o << " Transmitted Ions =" << t_.transmitted << ";  Backscattered Ions =" << t_.backscattered
      << "\n";
    o << " (These are not included in Skewne- and Kurtosis below.)\n\n";
    o << " Range Skewne- = " << fmt("%9.4f", t_.range.skewness())
      << "  &= [-(X-Rp)^3]/[N*Straggle^3]\n";
    o << " Range Kurtosis = " << fmt("%9.4f", t_.range.kurtosis())
      << "  &= [-(X-Rp)^4]/[N*Straggle^4]\n";
    o << " Statistical definitions above are those used in VLSI implant modelling.\n";
    o << "====================================================================\n";
    o << "=================================================================\n";
    o << "  Table Distribution Units are >>>  (Atoms/cm3) / (Atoms/cm2)  <<<  \n";
    o << "=================================================================\n";
    o << "   DEPTH     " << fmt("%10s", ion_symbol(c_).c_str());
    for (const auto &n : sim_.element_names()) o << "  " << fmt("%10s", n.c_str());
    o << "\n   (Ang.)       Ions   ";
    for (std::size_t i = 0; i < sim_.element_names().size(); ++i) o << "  Tgt. Atoms";
    o << "\n-----------  ----------";
    for (std::size_t i = 0; i < sim_.element_names().size(); ++i) o << "  ----------";
    o << "\n";
    // count / (ions * bin width in cm) = (atoms/cm3) / (atoms/cm2)
    const double norm = 1.0 / (ions() * in_units(width() * angstrom, centimeter));
    for (std::size_t i = 0; i < c_.depth_bins; ++i) {
      o << sci(depth_label(i)) << "  " << fmt("%10.4E", counted(t_.ion_stopped.values[i]) * norm);
      for (const auto &h : t_.recoil_stopped) o << "  " << fmt("%10.4E", counted(h.values[i]) * norm);
      o << "\n";
    }
  }

  void write_vacancy(const std::filesystem::path &p) const {
    auto o = open_file(p);
    banner(o, "Ion and Target VACANCY production");
    const bool kp = c_.damage == DamageCalculation::kinchin_pease;
    o << " Total Target Vacancies     = " << fmt("%.6g", t_.total_vacancies / ions()) << " /Ion\n";
    if (!kp) {
      o << " Total Target Displacements = " << fmt("%.6g", t_.total_displacements / ions())
        << " /Ion\n";
      o << " Total Target Replacement Collisions = "
        << fmt("%.6g", t_.total_replacements / ions()) << " /Ion\n";
    }
    o << "\n";
    const double norm = 1.0 / (ions() * width());
    if (kp) {
      o << "==========================================================\n";
      o << "  Table Units are  >>>> Vacancies/(Angstrom-Ion)  <<<<  \n";
      o << "==========================================================\n";
      o << "   TARGET     VACANCIES     VACANCIES \n";
      o << "   DEPTH         by            by     \n";
      o << "   (Ang.)       IONS         RECOILS  \n";
      o << "-----------  -----------  ------------\n";
      for (std::size_t i = 0; i < c_.depth_bins; ++i) {
        o << sci(depth_label(i)) << "  " << sci(counted(t_.vacancies_by_ion.values[i]) * norm) << "  "
          << sci(t_.vacancies_by_recoils.values[i] * norm) << "\n";
      }
    } else {
      o << "!!!! NOTE : 2nd Column below is number of Primary Knock-Ons !!!!\n";
      o << "           ( PKO are number of Target Atoms Recoiling from the Ion. )\n";
      o << "==========================================================\n";
      o << "  Table Units are  >>>> Vacancies/(Angstrom-Ion)  <<<<  \n";
      o << "==========================================================\n";
      o << "   TARGET    \n   DEPTH     " << fmt("%11s", ion_symbol(c_).c_str());
      for (const auto &n : sim_.element_names()) o << "  " << fmt("%11s", n.c_str());
      o << "\n   (Ang.)     Knock-Ons  ";
      for (std::size_t i = 0; i < sim_.element_names().size(); ++i) o << "   Vacancies ";
      o << "\n-----------  -----------";
      for (std::size_t i = 0; i < sim_.element_names().size(); ++i) o << "  -----------";
      o << "\n";
      for (std::size_t i = 0; i < c_.depth_bins; ++i) {
        o << sci(depth_label(i)) << "  " << sci(counted(t_.knock_ons.values[i]) * norm);
        for (std::size_t e = 0; e < t_.vacancies.size(); ++e) {
          const double removed = e < t_.vacancies_removed.size() ? t_.vacancies_removed[e].values[i] : 0.0;
          o << "  " << sci(counted(t_.vacancies[e].values[i], removed) * norm);
        }
        o << "\n";
      }
    }
    double eb = 0.0;
    double n_el = 0.0;
    for (const auto &l : c_.layers) {
      for (const auto &e : l.elements) {
        eb += in_units(e.lattice_binding_energy, electron_volt);
        n_el += 1.0;
      }
    }
    o << "\n To convert to Energy Lost - multiply by Average Binding Energy =  "
      << fmt("%g", n_el > 0 ? eb / n_el : 0.0) << "  eV/Vacancy\n";
  }

  void write_novac(const std::filesystem::path &p) const {
    auto o = open_file(p);
    banner(o, "Replacement Collisions");
    o << " Total Target Vacancies     = " << fmt("%.6g", t_.total_vacancies / ions()) << " /Ion\n";
    o << " Total Target Displacements = " << fmt("%.6g", t_.total_displacements / ions()) << " /Ion\n";
    o << " Total Target Replacement Collisions = " << fmt("%.6g", t_.total_replacements / ions())
      << " /Ion\n";
    o << "================================================================\n";
    o << "   Table Vacancy Units are  >>>>  Number /(Angstrom-Ion)  <<<<  \n";
    o << "================================================================\n";
    o << "  DEPTH (A)      Number\n";
    o << "------------   -----------\n";
    const double norm = 1.0 / (ions() * width());
    for (std::size_t i = 0; i < c_.depth_bins; ++i) {
      o << sci(depth_label(i)) << "   " << sci(counted(t_.replacements.values[i]) * norm) << "\n";
    }
  }

  void write_two_column(const std::filesystem::path &p, const std::string &what,
                        const std::string &units, const std::vector<std::string> &header,
                        const DepthHistogram &a, const DepthHistogram &b) const {
    auto o = open_file(p);
    banner(o, what);
    o << "===============================================================\n";
    o << "  " << units << "\n";
    o << "===============================================================\n\n";
    for (const auto &h : header) o << h << "\n";
    o << "-----------  -----------  -----------\n";
    const double norm = 1.0 / (ions() * width());
    for (std::size_t i = 0; i < c_.depth_bins; ++i) {
      o << sci(depth_label(i)) << "  " << sci(a.values[i] * norm) << "  "
        << sci(b.values[i] * norm) << "\n";
    }
  }

  void write_e2recoil(const std::filesystem::path &p) const {
    auto o = open_file(p);
    banner(o, "Energy transferred to RECOILS");
    o << "===========================================================\n";
    o << "  Table Energy Units are  >>>>  eV /(Angstrom-Ion)  <<<<  \n";
    o << "===========================================================\n";
    // Columns as TRIM writes them: energy from the ions, then (full
    // cascades) the energy absorbed by each target atom, (quick
    // calculation) by all recoils, or (surface sputtering) nothing else.
    const bool quick = c_.damage == DamageCalculation::kinchin_pease;
    // TRIM.DAT recoils with monolayer steps (code 7) run in TRIM's ion loop
    // and have the full-cascade layout.
    const bool per_atom = c_.damage == DamageCalculation::full_cascade ||
                          (c_.damage == DamageCalculation::monolayer &&
                           c_.source == TrimConfig::Source::particle_recoils);
    if (quick) {
      o << "   DEPTH       ENERGY       ENERGY    \n";
      o << "   DEPTH        from        Absorbed  \n";
      o << "   (Ang.)       IONS        by Recoils\n";
      o << "----------   ------------  -----------\n";
    } else {
      const std::size_t n = per_atom ? sim_.element_names().size() : 0;
      o << "   DEPTH       ENERGY      ";
      for (std::size_t i = 0; i < n; ++i) o << "Energy       ";
      o << "\n   DEPTH        from       ";
      for (std::size_t i = 0; i < n; ++i) o << "Absorbed     ";
      o << "\n   (Ang.)       IONS       ";
      for (std::size_t i = 0; i < n; ++i) o << "by  " << fmt("%-9s", sim_.element_names()[i].c_str());
      o << "\n----------   ------------  ";
      for (std::size_t i = 0; i < n; ++i) o << "-----------  ";
      o << "\n";
    }
    const double norm = 1.0 / (ions() * width());
    for (std::size_t i = 0; i < c_.depth_bins; ++i) {
      o << sci(depth_label(i)) << "  " << sci(t_.energy_to_recoils.values[i] * norm);
      if (quick) {
        double total = 0.0;
        for (const auto &h : t_.energy_absorbed) total += h.values[i];
        o << "  " << sci(total * norm);
      } else if (per_atom) {
        for (const auto &h : t_.energy_absorbed) o << "  " << sci(h.values[i] * norm);
      }
      o << "\n";
    }
  }

  void write_lateral(const std::filesystem::path &p) const {
    auto o = open_file(p);
    banner(o, "LATERAL SPREAD OF IONS");
    o << "  NOTE : The Projected LATERAL Spread is based on a projection on a plane.\n";
    o << "  The Lateral Projected Range is the average of the absolute values of\n";
    o << "  the projected lateral displacements from the X-axis; the Radial Range is\n";
    o << "  the mean radial displacement from the X-axis.\n";
    o << " Ion Average Range = " << sci(t_.range.mean()) << " A   Straggling = "
      << sci(t_.range.stddev()) << " A\n";
    o << " Ion Lateral Range = " << sci(t_.lateral_abs_y.mean()) << " A   Straggling = "
      << sci(t_.lateral_y.rms()) << " A\n";
    o << " Ion Radial  Range = " << sci(t_.radial.mean()) << " A   Straggling = "
      << sci(t_.radial.stddev()) << " A\n";
    o << "====================================================================\n";
    o << "                LATERAL     Projected      LATERAL        Radial    \n";
    o << "  TARGET       PROJ.RANGE   Straggling      RADIAL      Straggling  \n";
    o << "  DEPTH (A)   (Angstroms)   (Angstroms)   (Angstroms)   (Angstroms) \n";
    o << "------------ ------------- ------------- ------------- -------------\n";
    for (std::size_t i = 0; i < c_.depth_bins; ++i) {
      const auto &b = t_.lateral_bins[i];
      const double n = b[0];
      double lr = 0.0, ls = 0.0, rr = 0.0, rs = 0.0;
      if (n > 0) {
        lr = b[1] / n;
        ls = std::sqrt(b[2] / n);
        rr = b[3] / n;
        rs = std::sqrt(std::max(0.0, b[4] / n - rr * rr));
      }
      o << sci(depth_label(i)) << "   " << sci(lr) << "   " << sci(ls) << "   " << sci(rr)
        << "   " << sci(rs) << "\n";
    }
  }

  void write_tdata(const std::filesystem::path &p) const {
    auto o = open_file(p);
    o << "====== " << c_.title << " =======\n";
    o << "        AMCSET " << AMCSET_VERSION << " (SRIM mode)\n";
    o << "=======================================\n";
    o << "      Data  of  TRIM  Calculation      \n";
    o << "=======================================\n";
    o << "Ion = " << ion_symbol(c_) << " (" << c_.ion.z << ")   Ion Mass= "
      << fmt("%.4f", (c_.ion.mass / atomic_mass_unit).value()) << "\n";
    o << "Energy  = " << fmt("%g", in_units(c_.ion.energy, kilo_electron_volt)) << " keV\n";
    o << "Ion Angle to Surface = " << fmt("%g", in_units(c_.ion.angle, degree)) << " degrees\n";
    target(o);
    // One entry per target atom, in TRIM's order, with the energies the
    // calculation used (TRIM's replacements of zeros included).
    for (std::size_t k = 0; k < sim_.element_names().size(); ++k) {
      const TargetElement *found = nullptr;
      for (const auto &l : sim_.effective_layers()) {
        for (const auto &e : l.elements) {
          const std::string name = e.name.empty() ? std::string(Element::by_z(e.z).symbol()) : e.name;
          if (found == nullptr && e.z == sim_.element_z()[k] && name == sim_.element_names()[k]) found = &e;
        }
      }
      if (found == nullptr) continue;
      const auto &e = *found;
      o << "Target energies for target atom = " << Element::by_z(e.z).symbol() << "\n";
      o << "    Displacement = " << fmt("%g", in_units(e.displacement_energy, electron_volt))
        << " eV, Binding = " << fmt("%g", in_units(e.lattice_binding_energy, electron_volt))
        << " eV, Surface = " << fmt("%g", in_units(e.surface_binding_energy, electron_volt))
        << " eV\n";
    }
    o << "====================================================================\n";
    o << "Type of calculation    = " << to_string(c_.damage) << "\n";
    o << "Electronic stopping    = " << physics::to_string(c_.stopping_model)
      << (c_.stopping_table.empty() ? "" : " (" + c_.stopping_table + ")") << "\n";
    o << "Interatomic potential  = "
      << (c_.srim_data_directory.empty() ? "ZBL universal" : "SRIM-2013 pair potentials") << "\n";
    o << "Random number seed     = " << c_.seed << "\n";
    if (c_.fast_math) {
      o << "Fast approximations    = on (tabulated scattering angle; not number-for-number TRIM)\n";
    }
    o << "Depth Range of Tabulated Data=  " << sci(sim_.bin_start(0)) << "  -  "
      << sci(sim_.bin_start(c_.depth_bins)) << "  Angstroms\n";
    o << "====================================================================\n";
    o << "Total Ions calculated =  " << t_.ions << "\n";
    o << "Average Range         =  " << sci(t_.range.mean()) << " Angstroms\n";
    o << "Average Straggling    =  " << sci(t_.range.stddev()) << " Angstroms\n";
    o << "Average Vacancy/Ion   =  " << sci(t_.total_vacancies / ions()) << "\n";
    o << "====================================================================\n";
    o << "Total Backscattered Ions=  " << t_.backscattered << "\n";
    o << "Total Transmitted Ions  =  " << t_.transmitted << "\n";
    o << "====================================================================\n";
    const double e0 = in_units(c_.ion.energy, electron_volt);
    const auto pct = [&](double v) { return fmt("%7.2f", 100.0 * v / ions() / e0); };
    o << "Energy loss (percent of the incident energy)\n";
    o << "                 Ionization   Vacancies   Phonons\n";
    o << "   Ions:          " << pct(t_.total_ionization_ion) << "     " << pct(0.0) << "     "
      << pct(t_.total_phonons_ion) << "\n";
    o << "   Recoils:       " << pct(t_.total_ionization_recoils) << "     "
      << pct(t_.total_binding_energy_loss) << "     " << pct(t_.total_phonons_recoils) << "\n";
    o << "   Backscattered: " << pct(t_.backscattered_energy) << "   Transmitted: "
      << pct(t_.transmitted_energy) << "\n";
    o << "====================================================================\n";
    o << "Sputtering yields (atoms/ion):";
    for (std::size_t i = 0; i < sim_.element_names().size(); ++i) {
      o << " " << sim_.element_names()[i] << "=" << fmt("%.5g", t_.sputtered[i] / ions());
    }
    o << "\nCalculation time = " << fmt("%.1f", sim_.elapsed_seconds()) << " s\n";
  }

  void write_exits(const std::filesystem::path &p, char kind, const std::string &what) const {
    auto o = open_file(p);
    o << "==============================================================================\n";
    o << "=================  " << p.filename().string() << " : " << what << "  ==================\n";
    o << "=  This file tabulates the kinetics of ions or atoms leaving the target.     =\n";
    o << "=  Column #1: S= Sputtered Atom, B= Backscattered Ion, T= Transmitted Ion.   =\n";
    o << "=  Col.#2: Ion Number, Col.#3: Z of atom leaving, Col.#4: Atom energy (eV).  =\n";
    o << "=  Col.#5-7: Last location:  X= Depth into target, Y,Z= Transverse axes.     =\n";
    o << "=  Col.#8-10: Cosines of final trajectory.                                   =\n";
    o << "= *** The data lines can be read back as a TRIM.DAT particle list.           =\n";
    if (kind == 'S' || kind == 0) {
      if (c_.transport == TrimConfig::Transport::trim2013) {
        o << "=  As in TRIM: the surface binding energy is NOT subtracted from the energy, =\n";
        o << "=  and the cosines are those of the atom inside the surface.                 =\n";
      } else {
        o << "=  Energy after the surface barrier; cosines after refraction at it.        =\n";
      }
    }
    o << "==============================================================================\n";
    o << " Ion  Atom   Energy        Depth       Lateral-Position        Atom Direction      \n";
    o << " Numb Numb    (eV)          X(A)        Y(A)       Z(A)      Cos(X)  Cos(Y) Cos(Z) \n";
    // "%c %6zu %2u %13.6E %11.4E %11.4E %11.4E %9.6f %9.6f %9.6f\n" (output/text.h
    // prints the same characters without printf).
    output::text::BufferedWriter w(o);
    for (const auto &r : t_.exits) {
      if (kind != 0 && r.kind != kind) continue;
      using namespace output::text;
      char *q = w.reserve(10 * max_field);
      *q++ = r.kind;
      *q++ = ' ';
      q = put_unsigned(q, r.ion, 6);
      *q++ = ' ';
      q = put_unsigned(q, r.z, 2);
      *q++ = ' ';
      q = put_sci(q, r.energy_eV, 13, 6);
      for (const double v : {r.x, r.y, r.z_pos}) {
        *q++ = ' ';
        q = put_sci(q, v, 11, 4);
      }
      for (const double v : {r.cx, r.cy, r.cz}) {
        *q++ = ' ';
        q = put_fixed(q, v, 9, 6);
      }
      *q++ = '\n';
      w.commit(q);
    }
    w.finish();
  }

  void write_exyz(const std::filesystem::path &p) const {
    auto o = open_file(p);
    o << " =====================================================================\n";
    o << " =                   Ion Energy vs Position File                     =\n";
    o << " =====================================================================\n";
    o << " =   AXIS DEFINITIONS: X=Depth, Y,Z= Lateral plane of target surface.=\n";
    o << " =   (If beam enters target at an angle, this tilt is in Y direction)=\n";
    o << " =   Shown are: Ion Number, Energy (keV), X, Y, Z Position           =\n";
    o << " ==================  CALCULATION   DATA  =============================\n";
    o << "  Ion Data: Name, Mass,   Energy , Energy Interval\n";
    o << "           " << ion_symbol(c_) << "     " << fmt("%06.2f", (c_.ion.mass / atomic_mass_unit).value())
      << "  " << fmt("%g", in_units(c_.ion.energy, kilo_electron_volt)) << "keV    "
      << fmt("%g", in_units(c_.exyz_interval, electron_volt)) << "eV\n";
    o << " =====================================================================\n";
    o << "Ion       Energy     Depth (X)     Y           Z       Electronic   Energy Lost to\n";
    o << "Number    (keV)    (Angstrom)  (Angstrom)  (Angstrom)  Stop.(eV/A)  Last Recoil(eV)\n";
    o << "------- ----------- ---------- ----------- ----------- -----------  ---------------\n";
    // TRIM's column widths (7, 11, then 12 for each of the rest):
    // "%07zu%11.4E%12.4E%12.4E%12.4E%12.4E%12.4E\n".
    output::text::BufferedWriter w(o);
    for (const auto &r : t_.exyz) {
      using namespace output::text;
      char *q = w.reserve(8 * max_field);
      q = put_unsigned(q, r.ion, 7, '0');
      q = put_sci(q, r.energy_keV, 11, 4);
      for (const double v : {r.x, r.y, r.z, r.se_eV_per_A, r.last_recoil_eV}) q = put_sci(q, v, 12, 4);
      *q++ = '\n';
      w.commit(q);
    }
    w.finish();
  }

  void write_collisions(const std::filesystem::path &p) const {
    auto o = open_file(p);
    banner(o, "COLLISION HISTORY (ion collisions that displaced a target atom)");
    o << "  NOTES: Target DISPlacements = VACancies + REPLACement Collisions.\n";
    o << "         Target VACancies = INTERstitials + Sputtered + Transmitted Atoms.\n";
    o << "         Quick (Kinchin-Pease) mode lists the modified Kinchin-Pease\n";
    o << "         displacements of each primary recoil.\n";
    if (c_.collision_recoil_details) {
      o << "         Each cascade lists its displaced atoms (Vac 1, with the energy they\n";
      o << "         start with) and replacements (Repl 1) in the order they happen.\n";
    }
    const char *header =
        "  Ion    Energy      Depth       Lateral Distance (A)      Se     Atom   Recoil       Target     Target     Target     Target\n"
        "  Numb   (keV)        (A)        Y Axis      Z Axis      (eV/A)   Hit   Energy(eV)    DISP.      VAC.       REPLAC     INTER\n"
        "-------------------------------------------------------------------------------------------------------------------------------\n";
    o << header;
    char buf[400];
    auto cascade_line = [&](const CascadeRecord &r) {
      std::snprintf(buf, sizeof buf,
                    "%7zu %10.4E %11.4E %11.4E %11.4E %9.3f  %-3s %11.4E %10.3f %10.3f %10.3f %10.3f%s\n",
                    r.ion, r.ion_energy_keV, r.x, r.y, r.z, r.se_eV_per_A,
                    std::string(Element::by_z(r.z_hit).symbol()).c_str(), r.recoil_energy_eV,
                    r.displacements, r.vacancies, r.replacements, r.interstitials,
                    r.atoms.empty() ? "" : "  <== Start of New Cascade");
      o << buf;
      if (r.atoms.empty()) return;
      std::snprintf(buf, sizeof buf, "    Recoil  Atom  Energy(eV)      X (A)       Y (A)       Z (A)    Vac  Repl   Ion %zu\n", r.ion);
      o << buf;
      for (const auto &a : r.atoms) {
        std::snprintf(buf, sizeof buf, "    %6u  %4u  %11.4E %11.4E %11.4E %11.4E   %d     %d%s\n", a.number, a.z,
                      a.energy_eV, a.x, a.y, a.z_pos, a.vacancy ? 1 : 0, a.replacement ? 1 : 0,
                      a.number == 1 ? "   <-- Prime Recoil" : "");
        o << buf;
      }
      std::snprintf(buf, sizeof buf, "    Summary of Above Cascade ==> %11.4E %8.0f %8.0f %8.0f %8.0f\n", r.recoil_energy_eV,
                    r.displacements, r.vacancies, r.replacements, r.interstitials);
      o << buf << "\n";
    };
    // Cascades and ion summaries are both in ion order.
    std::size_t next = 0;
    double sd = 0.0, sr = 0.0, sv = 0.0, si = 0.0, ss = 0.0, st = 0.0;
    std::size_t n = 0;
    for (const auto &ion : t_.ion_summaries) {
      while (next < t_.cascades.size() && t_.cascades[next].ion <= ion.ion) cascade_line(t_.cascades[next++]);
      ++n;
      sd += ion.displacements;
      sr += ion.replacements;
      sv += ion.vacancies;
      si += ion.interstitials;
      ss += ion.sputtered;
      st += ion.transmitted;
      const double k = static_cast<double>(n);
      o << "===============================================================================================\n";
      std::snprintf(buf, sizeof buf, "- For Ion %07zu:               For All Ions to date :\n", ion.ion);
      o << buf;
      const std::pair<const char *, std::pair<double, double>> rows[] = {
          {"Displacements    ", {ion.displacements, sd / k}}, {"Replacements     ", {ion.replacements, sr / k}},
          {"Vacancies        ", {ion.vacancies, sv / k}},     {"Interstitials    ", {ion.interstitials, si / k}},
          {"Sputtered Atoms  ", {ion.sputtered, ss / k}},     {"Transmitted Atoms", {ion.transmitted, st / k}}};
      for (const auto &row : rows) {
        std::snprintf(buf, sizeof buf, "- %s = %10.1f   Average per Ion = %12.3f\n", row.first, row.second.first,
                      row.second.second);
        o << buf;
      }
      o << "===============================================================================================\n";
      if (c_.collision_recoil_details) o << header;
    }
    while (next < t_.cascades.size()) cascade_line(t_.cascades[next++]);
  }

  void write_spectrum(const std::filesystem::path &p) const {
    auto o = open_file(p);
    o << "# AMCSET SRIM mode: energy spectrum of displaced target atoms\n";
    o << "# recoil energy = transferred energy - lattice binding energy\n";
    o << "# columns: E_low(eV) E_high(eV) primary_knock_ons_per_ion all_recoils_per_ion\n";
    o << "# (all recoils are only followed in full-cascade and monolayer modes)\n";
    const int ppd = TrimTallies::spectrum_bins_per_decade;
    for (std::size_t i = 0; i < t_.pka_spectrum.size(); ++i) {
      if (t_.pka_spectrum[i] == 0.0 && t_.recoil_spectrum[i] == 0.0) continue;
      char buf[160];
      std::snprintf(buf, sizeof buf, "%12.5E %12.5E %12.5E %12.5E\n",
                    std::pow(10.0, static_cast<double>(i) / ppd),
                    std::pow(10.0, static_cast<double>(i + 1) / ppd), t_.pka_spectrum[i] / ions(),
                    t_.recoil_spectrum[i] / ions());
      o << buf;
    }
  }

  void write_range3d(const std::filesystem::path &p) const {
    auto o = open_file(p);
    o << "====== RANGE_3D.txt : File of Final Ion Positions =====\n\n";
    o << "Ion = " << ion_symbol(c_) << " (" << c_.ion.z << ")   Ion Mass= "
      << fmt("%.4f", (c_.ion.mass / atomic_mass_unit).value()) << "\n";
    o << "Energy  = " << fmt("%g", in_units(c_.ion.energy, kilo_electron_volt)) << " keV\n";
    o << "Ion Angle to Surface = " << fmt("%g", in_units(c_.ion.angle, degree)) << " degrees\n";
    target(o);
    o << "Ion       Depth  X   Lateral Y   Lateral Z  \n";
    o << "Number    (Angstrom) (Angstrom)  (Angstrom) \n";
    o << "-------  ----------- ----------- -----------\n";
    // "%07zu   %11.4E %11.4E %11.4E\n" (one line per ion: the longest file).
    output::text::BufferedWriter w(o);
    for (const auto &r : t_.positions) {
      using namespace output::text;
      char *q = w.reserve(4 * max_field + 8);
      q = put_unsigned(q, r.ion, 7, '0');
      *q++ = ' ';
      *q++ = ' ';
      for (const double v : {r.x, r.y, r.z}) {
        *q++ = ' ';
        q = put_sci(q, v, 11, 4);
      }
      *q++ = '\n';
      w.commit(q);
    }
    w.finish();
  }

  const TrimSimulation &sim_;
  const TrimConfig &c_;
  const TrimTallies &t_;
};

std::string json_escape(const std::string &s) {
  std::string o;
  for (char c : s) {
    if (c == '"' || c == '\\') {
      o += '\\';
      o += c;
    } else if (static_cast<unsigned char>(c) < 0x20) {
      o += ' ';
    } else {
      o += c;
    }
  }
  return o;
}

std::string json_number(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.10g", v);
  return buf;
}

} // namespace

void TrimSimulation::write_outputs(const std::filesystem::path &directory) const {
  const auto &spec = config_.outputs;
  std::filesystem::create_directories(directory);
  const Writer writer(*this);
  std::vector<std::function<void()>> tasks;
  if (spec.formats.count(output::Format::trim)) tasks = writer.tasks(directory);
  tasks.push_back([this, &directory] {
    // AMCSET's own summary (always; tools/trim_compare reads it).
    std::ofstream js = open_file(directory / "summary.json");
    js << summary_json() << "\n";
    if (!js) throw std::runtime_error("Error writing " + (directory / "summary.json").string());
  });
  tasks.push_back([this, &directory, &spec] {
    auto data = make_dataset(*this);
    data.configuration_json = config_.configuration_json;
    output::write_dataset(directory, data, spec);
  });
  common::run_tasks(tasks);
}

std::string TrimSimulation::summary_json() const {
  const auto &t = tallies_;
  const double n = std::max<double>(1.0, static_cast<double>(t.ions));
  const double e0 = in_units(config_.ion.energy, electron_volt);
  std::ostringstream j;
  j << "{\n";
  j << "  \"program\": \"AMCSET " << AMCSET_VERSION << " SRIM mode\",\n";
  j << "  \"title\": \"" << json_escape(config_.title) << "\",\n";
  j << "  \"ion\": {\"z\": " << config_.ion.z << ", \"symbol\": \""
    << Element::by_z(config_.ion.z).symbol() << "\", \"mass_amu\": "
    << json_number((config_.ion.mass / atomic_mass_unit).value()) << ", \"energy_keV\": "
    << json_number(in_units(config_.ion.energy, kilo_electron_volt)) << ", \"angle_deg\": "
    << json_number(in_units(config_.ion.angle, degree)) << "},\n";
  j << "  \"damage\": \"" << to_string(config_.damage) << "\",\n";
  j << "  \"stopping_model\": \"" << physics::to_string(config_.stopping_model) << "\",\n";
  j << "  \"pair_potentials\": " << (config_.srim_data_directory.empty() ? "false" : "true")
    << ",\n";
  if (config_.fast_math) j << "  \"fast_math\": true,\n";
  j << "  \"seed\": " << config_.seed << ",\n";
  if (config_.random_seed) j << "  \"seed_random\": true,\n"; // (drawn for this run; "seed" reruns it)
  j << "  \"ions\": " << t.ions << ",\n";
  j << "  \"elapsed_seconds\": " << json_number(elapsed_seconds_.load()) << ",\n";
  j << "  \"range\": {\"mean_A\": " << json_number(t.range.mean())
    << ", \"straggling_A\": " << json_number(t.range.stddev())
    << ", \"skewness\": " << json_number(t.range.skewness())
    << ", \"kurtosis\": " << json_number(t.range.kurtosis())
    << ", \"lateral_A\": " << json_number(t.lateral_abs_y.mean())
    << ", \"lateral_straggling_A\": " << json_number(t.lateral_y.rms())
    << ", \"radial_A\": " << json_number(t.radial.mean())
    << ", \"radial_straggling_A\": " << json_number(t.radial.stddev())
    << ", \"stopped_ions\": " << json_number(t.range.n) << "},\n";
  j << "  \"backscattered\": {\"count\": " << t.backscattered << ", \"fraction\": "
    << json_number(t.backscattered / n) << ", \"energy_fraction\": "
    << json_number(t.backscattered_energy / n / e0) << "},\n";
  j << "  \"transmitted\": {\"count\": " << t.transmitted << ", \"fraction\": "
    << json_number(t.transmitted / n) << ", \"energy_fraction\": "
    << json_number(t.transmitted_energy / n / e0) << "},\n";
  // The quick calculation has no replacement collisions: every vacancy is a
  // displaced atom. (Its own displacement count is the ion's knock-ons.)
  const bool quick = config_.damage == DamageCalculation::kinchin_pease;
  j << "  \"per_ion\": {\"vacancies\": " << json_number(t.total_vacancies / n)
    << ", \"displacements\": " << json_number((quick ? t.total_vacancies : t.total_displacements) / n)
    << ", \"replacements\": " << json_number(t.total_replacements / n);
  if (quick) {
    j << ", \"knock_ons\": " << json_number(t.total_displacements / n);
  } else {
    j << ", \"interstitials\": " << json_number(t.total_interstitials / n);
  }
  j << ", \"recoils_followed\": " << json_number(t.recoils_followed / n) << "},\n";
  j << "  \"energy_percent\": {\"ionization_ion\": "
    << json_number(100.0 * t.total_ionization_ion / n / e0)
    << ", \"ionization_recoils\": " << json_number(100.0 * t.total_ionization_recoils / n / e0)
    << ", \"phonons_ion\": " << json_number(100.0 * t.total_phonons_ion / n / e0)
    << ", \"phonons_recoils\": " << json_number(100.0 * t.total_phonons_recoils / n / e0)
    << ", \"vacancies\": " << json_number(100.0 * t.total_binding_energy_loss / n / e0)
    << ", \"to_recoils_from_ion\": " << json_number(100.0 * t.total_energy_to_recoils / n / e0)
    << "},\n";
  j << "  \"elements\": [";
  for (std::size_t i = 0; i < element_names_.size(); ++i) {
    if (i) j << ", ";
    j << "{\"name\": \"" << json_escape(element_names_[i]) << "\", \"z\": " << element_z_[i]
      << ", \"sputtering_yield\": " << json_number(t.sputtered[i] / n)
      << ", \"sputtered_mean_energy_eV\": "
      << json_number(t.sputtered[i] > 0 ? t.sputtered_energy[i] / t.sputtered[i] : 0.0)
      << ", \"transmitted_recoils_per_ion\": " << json_number(t.transmitted_recoils[i] / n) << "}";
  }
  j << "],\n";
  if (const auto crystals = crystal_reports(config_, true); !crystals.empty()) {
    // The crystals' orientation and a polycrystal's grains (not TRIM;
    // docs/OUTPUT.md "Crystals").
    const auto vec = [&](const crystal::Vec3 &v) {
      return "[" + json_number(v.x) + ", " + json_number(v.y) + ", " + json_number(v.z) + "]";
    };
    const auto mat = [&](const crystal::Mat3 &m) { return "[" + vec(m[0]) + ", " + vec(m[1]) + ", " + vec(m[2]) + "]"; };
    const auto dir = [&](const std::array<int, 3> &uvw, double off) {
      return "{\"uvw\": \"" + uvw_text(uvw) + "\", \"off_deg\": " + json_number(off) + "}";
    };
    j << "  \"crystals\": [";
    for (std::size_t i = 0; i < crystals.size(); ++i) {
      const auto &r = crystals[i];
      if (i) j << ",";
      j << "\n    {\"layer\": " << r.layer + 1 << ", \"layer_name\": \"" << json_escape(r.layer_name)
        << "\", \"cell\": \"" << json_escape(r.cell_name) << "\", \"followed\": " << (r.followed ? "true" : "false")
        << ", \"top_A\": " << json_number(r.top_A) << ",\n     \"lattice_vectors_A\": " << mat(r.lattice_A)
        << ",\n     \"rotation_crystal_to_trim_axes\": " << mat(r.to_trim)
        << ",\n     \"rotation_crystal_to_output_axes\": " << mat(r.to_output)
        << ",\n     \"euler_bunge_deg\": [" << json_number(r.euler_deg[0]) << ", " << json_number(r.euler_deg[1])
        << ", " << json_number(r.euler_deg[2]) << "], \"depth_direction\": " << dir(r.depth_uvw, r.depth_off_deg)
        << ", \"trim_y_direction\": " << dir(r.y_uvw, r.y_off_deg)
        << ", \"trim_z_direction\": " << dir(r.z_uvw, r.z_off_deg)
        << ",\n     \"beam_direction\": " << dir(r.beam_uvw, r.beam_off_deg)
        << ", \"beam_unit_vector_crystal\": " << vec(r.beam_crystal) << ", \"tilt_deg\": " << json_number(r.tilt_deg)
        << ", \"twist_deg\": " << json_number(r.twist_deg) << ", \"temperature_K\": " << json_number(r.temperature_K)
        << ",\n     \"vibrations\": [";
      for (std::size_t k = 0; k < r.vibrations.size(); ++k) {
        const auto &v = r.vibrations[k];
        j << (k ? ", " : "") << "{\"z\": " << v.z << ", \"rms_A\": " << json_number(v.rms_A)
          << ", \"from_table\": " << (v.from_table ? "true" : "false") << "}";
      }
      j << "],\n     \"axes_near_beam\": [";
      for (std::size_t k = 0; k < r.axes.size(); ++k) {
        const auto &a = r.axes[k];
        j << (k ? ", " : "") << "{\"uvw\": \"" << uvw_text(a.uvw) << "\", \"angle_deg\": " << json_number(a.angle_deg)
          << ", \"spacing_A\": " << json_number(a.spacing_A) << ", \"psi1_deg\": " << json_number(a.psi1_deg)
          << ", \"psi2_deg\": " << json_number(a.psi2_deg) << ", \"critical_deg\": " << json_number(a.critical_deg)
          << "}";
      }
      j << "]";
      if (r.grains) j << ",\n     \"grains\": " << crystal::to_json(*r.grains);
      j << "}";
    }
    j << "],\n";
  }
  if (config_.transport == TrimConfig::Transport::trim2013) {
    // Damage models of the primary knock-on atoms (docs/OUTPUT.md "dpa"), per ion.
    j << "  \"damage_models\": {\"damage_energy_eV_per_ion\": " << json_number(t.total_pka_damage_energy / n)
      << ", \"nrt_per_ion\": " << json_number(t.total_nrt / n)
      << ", \"kinchin_pease_per_ion\": " << json_number(t.total_kp / n);
    if (t.pka_without_efficiency == 0.0) {
      j << ", \"arc_dpa_per_ion\": " << json_number(t.total_arc / n)
        << ", \"rpa_per_ion\": " << json_number(t.total_rpa / n);
    }
    j << ", \"trim_vacancies_per_ion\": " << json_number(t.total_vacancies / n)
      << ", \"dpa_fluence_per_cm2\": " << json_number(config_.dpa_fluence_per_cm2) << "},\n";
  }
  if (config_.md_events.enabled()) {
    // Hybrid BCA-MD events (docs/OUTPUT.md "MD event file"): per E_c, per ion.
    j << "  \"md_events\": {\"mode\": \"" << (config_.md_events.handoff ? "handoff" : "record")
      << "\", \"file\": \"" << json_escape(config_.md_events.file) << "\", \"E_min_eV\": "
      << json_number(in_units(config_.md_events.min_energy, electron_volt)) << ", \"transitions\": [";
    for (std::size_t i = 0; i < t.md.size(); ++i) {
      const auto &m = t.md[i];
      if (i) j << ", ";
      j << "{\"E_c_eV\": " << json_number(m.transition_energy_eV)
        << ", \"events_per_ion\": " << json_number(m.events / n)
        << ", \"birth_events_per_ion\": " << json_number(m.birth_events / n)
        << ", \"crossing_events_per_ion\": " << json_number(m.crossing_events / n)
        << ", \"energy_eV_per_ion\": " << json_number(m.energy_eV / n)
        << ", \"energy_percent\": " << json_number(100.0 * m.energy_eV / n / e0)
        << ", \"below_E_min_per_ion\": " << json_number(m.below_min / n)
        << ", \"below_E_min_energy_eV_per_ion\": " << json_number(m.below_min_energy_eV / n)
        << ", \"ions_handed\": " << json_number(m.ions_handed) << "}";
    }
    j << "]},\n";
  }
  j << "  \"depth_bins\": {\"start_A\": " << json_number(bin_start(0)) << ", \"width_A\": "
    << json_number(bin_width()) << ", \"count\": " << config_.depth_bins << "}\n";
  j << "}";
  return j.str();
}

} // namespace trim
} // namespace amcset
