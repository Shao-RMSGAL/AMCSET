#pragma once

#include <atomic>
#include <cstddef>
#include <mutex>
#include <vector>
#include <sstream>

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
 * interval [min,max] are ignored.
 */
struct Histogram {
  std::vector<double> bins;
  double min = 0.0;
  double max = 1.0;

  Histogram() = default;
  Histogram(size_t bin_count, double min_, double max_)
      : bins(bin_count, 0.0), min(min_), max(max_) {}

  void add(double x) { add_value(x, 1.0); }

  void add_value(double x, double value) {
    if (bins.empty()) return;
    if (x < min || x > max) return;
    const auto width = max - min;
    if (width <= 0.0) return;
    auto idx = static_cast<size_t>((x - min) / width * bins.size());
    if (idx >= bins.size()) idx = bins.size() - 1;
    bins[idx] += value;
  }

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
 * (Welford), and a depth histogram.
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
  Histogram depth_histogram;
  Histogram displacement_histogram;
  Histogram electronic_energy_loss_histogram;
  Histogram nuclear_energy_loss_histogram;
  Histogram total_energy_loss_histogram;

  /**
   * @brief Default-constructed empty tallies.
   */
  Tallies() = default;
  Tallies(size_t bin_count, double max_depth,
          size_t displacement_bins = 100,
          size_t electronic_bins = 100,
          size_t nuclear_bins = 100,
          size_t total_loss_bins = 100)
      : incidents(0), total_displacements(0), sputtered(0),
        range_stats(), energy_loss_stats(),
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
    if (depth_histogram.bins.empty() && !o.depth_histogram.bins.empty()) {
      depth_histogram = o.depth_histogram;
    } else {
      depth_histogram.merge(o.depth_histogram);
    }
    if (displacement_histogram.bins.empty() && !o.displacement_histogram.bins.empty()) {
      displacement_histogram = o.displacement_histogram;
    } else {
      displacement_histogram.merge(o.displacement_histogram);
    }
    if (electronic_energy_loss_histogram.bins.empty() && !o.electronic_energy_loss_histogram.bins.empty()) {
      electronic_energy_loss_histogram = o.electronic_energy_loss_histogram;
    } else {
      electronic_energy_loss_histogram.merge(o.electronic_energy_loss_histogram);
    }
    if (nuclear_energy_loss_histogram.bins.empty() && !o.nuclear_energy_loss_histogram.bins.empty()) {
      nuclear_energy_loss_histogram = o.nuclear_energy_loss_histogram;
    } else {
      nuclear_energy_loss_histogram.merge(o.nuclear_energy_loss_histogram);
    }
    if (total_energy_loss_histogram.bins.empty() && !o.total_energy_loss_histogram.bins.empty()) {
      total_energy_loss_histogram = o.total_energy_loss_histogram;
    } else {
      total_energy_loss_histogram.merge(o.total_energy_loss_histogram);
    }
  }

  /**
   * @brief Render a compact JSON summary of the tallies.
   *
   * The produced JSON includes counters, means, counts and the depth
   * histogram bins. This is intended for quick diagnostics and may be
   * extended later.
   */
  std::string to_json() const {
    std::ostringstream ss;
    ss << "{";
    ss << "\"incidents\":" << incidents << ",";
    ss << "\"total_displacements\":" << total_displacements << ",";
    ss << "\"sputtered\":" << sputtered << ",";
    ss << "\"range_mean\":" << (double)range_stats.mean << ",";
    ss << "\"range_count\":" << range_stats.n << ",";
    ss << "\"energy_loss_mean\":" << (double)energy_loss_stats.mean << ",";
    ss << "\"energy_loss_count\":" << energy_loss_stats.n << ",";
    ss << "\"depth_histogram_bins\":[";
    for (size_t i = 0; i < depth_histogram.bins.size(); ++i) {
      if (i) ss << ",";
      ss << depth_histogram.bins[i];
    }
    ss << "],";
    ss << "\"displacement_histogram_bins\":[";
    for (size_t i = 0; i < displacement_histogram.bins.size(); ++i) {
      if (i) ss << ",";
      ss << displacement_histogram.bins[i];
    }
    ss << "],";
    ss << "\"electronic_loss_histogram_bins\":[";
    for (size_t i = 0; i < electronic_energy_loss_histogram.bins.size(); ++i) {
      if (i) ss << ",";
      ss << electronic_energy_loss_histogram.bins[i];
    }
    ss << "],";
    ss << "\"nuclear_loss_histogram_bins\":[";
    for (size_t i = 0; i < nuclear_energy_loss_histogram.bins.size(); ++i) {
      if (i) ss << ",";
      ss << nuclear_energy_loss_histogram.bins[i];
    }
    ss << "],";
    ss << "\"total_loss_histogram_bins\":[";
    for (size_t i = 0; i < total_energy_loss_histogram.bins.size(); ++i) {
      if (i) ss << ",";
      ss << total_energy_loss_histogram.bins[i];
    }
    ss << "]";
    ss << "}";
    return ss.str();
  }
};

} // namespace common
} // namespace amcset
