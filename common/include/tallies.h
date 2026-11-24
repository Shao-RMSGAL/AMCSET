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

/*! \file tallies.h
 *  \brief Accumulators (running moments, fixed-bin histograms) used to
 *  aggregate simulation statistics across bombardments.
 *
 *  These types are not internally synchronized: each Bombardment fills its
 *  own Tallies and Simulation::merge_tallies() merges it under a mutex.
 */

#pragma once

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

namespace amcset {
namespace common {

/**
 * @brief Online algorithm accumulator for mean and variance (Welford's method).
 *
 * Use `add()` to include a new sample and `merge()` to combine two
 * accumulators produced in parallel.
 */
struct Welford {
  size_t n = 0;
  long double mean = 0.0L;
  long double M2 = 0.0L;

  void add(long double x) {
    n++;
    auto delta = x - mean;
    mean += delta / n;
    auto delta2 = x - mean;
    M2 += delta * delta2;
  }

  void merge(const Welford &o) {
    if (o.n == 0) return;
    if (n == 0) { *this = o; return; }
    auto total = n + o.n;
    auto delta = o.mean - mean;
    M2 = M2 + o.M2 + delta * delta * n * o.n / total;
    mean = (mean * n + o.mean * o.n) / total;
    n = total;
  }

  long double variance() const { return (n > 1) ? (M2 / (n - 1)) : 0.0L; }
};

/**
 * @brief Simple fixed-bin histogram.
 *
 * Bins are uniformly spaced between `min` and `max`. Values outside the
 * interval [min,max] are ignored; a value exactly equal to `max` falls in
 * the last bin.
 */
struct Histogram {
  std::vector<double> bins;
  double min = 0.0;
  double max = 1.0;

  Histogram() = default;
  Histogram(size_t bin_count, double min_, double max_)
      : bins(bin_count, 0.0), min(min_), max(max_) {}

  //! Count one event at position x.
  void add(double x) { add_value(x, 1.0); }

  //! Add a weight (e.g. deposited energy) to the bin containing x.
  void add_value(double x, double value) {
    if (bins.empty()) return;
    if (x < min || x > max) return;
    const auto width = max - min;
    if (width <= 0.0) return;
    auto idx = static_cast<size_t>((x - min) / width * bins.size());
    if (idx >= bins.size()) idx = bins.size() - 1;
    bins[idx] += value;
  }

  /*!
   * \brief Add another histogram bin-by-bin.
   *
   * An empty histogram adopts \p o (including its range). Histograms with a
   * different bin count are silently ignored; the ranges are not compared.
   */
  void merge(const Histogram &o) {
    if (o.bins.empty()) return;
    if (bins.empty()) {
      bins = o.bins;
      min = o.min;
      max = o.max;
      return;
    }
    if (bins.size() != o.bins.size()) return; // incompatible
    for (size_t i = 0; i < bins.size(); ++i) bins[i] += o.bins[i];
  }
};

/**
 * @brief Collection of simulation tallies/statistics aggregated per
 *        bombardment or across a simulation.
 *
 * Contains simple counters, running moments for range and energy loss
 * (Welford), and depth-binned histograms. All positions are depths in
 * meters; all histograms span [0, Settings::range_].
 *
 * - range_histogram: final depth of each incident particle.
 * - depth_histogram: legacy combined histogram; receives both incident final
 *   depths and displacement depths. Only exported in tallies.json.
 * - displacement_histogram: depth of each displacement (knock-on).
 * - *_energy_loss_histogram: deposited energy (eV) per depth bin.
 */
struct Tallies {
  // Lightweight counters
  size_t incidents = 0;
  size_t total_displacements = 0;
  size_t sputtered = 0;

  // Moments
  Welford range_stats;
  Welford energy_loss_stats;

  // Depth histograms
  Histogram range_histogram;
  Histogram depth_histogram;
  Histogram displacement_histogram;
  Histogram electronic_energy_loss_histogram;
  Histogram nuclear_energy_loss_histogram;
  Histogram total_energy_loss_histogram;

  /**
   * @brief Default-constructed empty tallies (all histograms have 0 bins and
   * adopt the first non-empty histogram merged into them).
   */
  Tallies() = default;
  /**
   * @brief Tallies with pre-sized histograms over [0, max_depth].
   * @param bin_count Bins of the range and depth histograms.
   */
  Tallies(size_t bin_count, double max_depth,
          size_t displacement_bins = 100,
          size_t electronic_bins = 100,
          size_t nuclear_bins = 100,
          size_t total_loss_bins = 100)
      : incidents(0), total_displacements(0), sputtered(0),
        range_stats(), energy_loss_stats(),
        range_histogram(bin_count, 0.0, max_depth),
        depth_histogram(bin_count, 0.0, max_depth),
        displacement_histogram(displacement_bins, 0.0, max_depth),
        electronic_energy_loss_histogram(electronic_bins, 0.0, max_depth),
        nuclear_energy_loss_histogram(nuclear_bins, 0.0, max_depth),
        total_energy_loss_histogram(total_loss_bins, 0.0, max_depth) {}

  /**
   * @brief Merge another `Tallies` into this one (in-place).
   *
   * Counters are summed and running statistics/histograms are merged.
   */
  void merge(const Tallies &o) {
    incidents += o.incidents;
    total_displacements += o.total_displacements;
    sputtered += o.sputtered;
    range_stats.merge(o.range_stats);
    energy_loss_stats.merge(o.energy_loss_stats);
    // Histogram::merge() adopts the other histogram when this one is empty.
    depth_histogram.merge(o.depth_histogram);
    range_histogram.merge(o.range_histogram);
    displacement_histogram.merge(o.displacement_histogram);
    electronic_energy_loss_histogram.merge(o.electronic_energy_loss_histogram);
    nuclear_energy_loss_histogram.merge(o.nuclear_energy_loss_histogram);
    total_energy_loss_histogram.merge(o.total_energy_loss_histogram);
  }

  /**
   * @brief Render a compact JSON summary of the tallies.
   *
   * The produced JSON includes counters, means (SI units: meters, joules),
   * counts and all histogram bins. Written to output/tallies.json by
   * Simulation::run_simulation() and read back by the GUI summary tab.
   */
  std::string to_json() const {
    const auto write_bins = [](std::ostringstream &out, const char *key,
                               const Histogram &histogram) {
      out << "\"" << key << "\":[";
      for (size_t i = 0; i < histogram.bins.size(); ++i) {
        if (i) out << ",";
        out << histogram.bins[i];
      }
      out << "]";
    };
    std::ostringstream ss;
    ss << "{";
    ss << "\"incidents\":" << incidents << ",";
    ss << "\"total_displacements\":" << total_displacements << ",";
    ss << "\"sputtered\":" << sputtered << ",";
    ss << "\"range_mean\":" << (double)range_stats.mean << ",";
    ss << "\"range_count\":" << range_stats.n << ",";
    ss << "\"range_straggling\":" << (double)std::sqrt(range_stats.variance()) << ",";
    ss << "\"energy_loss_mean\":" << (double)energy_loss_stats.mean << ",";
    ss << "\"energy_loss_count\":" << energy_loss_stats.n << ",";
    write_bins(ss, "depth_histogram_bins", depth_histogram);
    ss << ",";
    write_bins(ss, "range_histogram_bins", range_histogram);
    ss << ",";
    write_bins(ss, "displacement_histogram_bins", displacement_histogram);
    ss << ",";
    write_bins(ss, "electronic_loss_histogram_bins",
               electronic_energy_loss_histogram);
    ss << ",";
    write_bins(ss, "nuclear_loss_histogram_bins", nuclear_energy_loss_histogram);
    ss << ",";
    write_bins(ss, "total_loss_histogram_bins", total_energy_loss_histogram);
    ss << "}";
    return ss.str();
  }
};

/*!
 * \brief Write a depth histogram as CSV with the columns
 * depth_min_angstrom, depth_max_angstrom, <value_label> (the format the GUI
 * reads). Parent directories are created. \throws std::runtime_error
 */
void write_histogram_csv(const std::filesystem::path &file, const Histogram &histogram,
                         const std::string &value_label = "count");

/*!
 * \brief Write the GUI result files (displacement_histogram.csv,
 * range_histogram.csv, electronic/nuclear/total_loss_histogram.csv and
 * tallies.json) of \p tallies into \p directory.
 */
void write_tallies_files(const std::filesystem::path &directory, const Tallies &tallies);

} // namespace common
} // namespace amcset
