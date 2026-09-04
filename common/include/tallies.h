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
  std::vector<uint64_t> bins;
  double min = 0.0;
  double max = 1.0;

  Histogram() = default;
  Histogram(size_t bin_count, double min_, double max_)
      : bins(bin_count, 0), min(min_), max(max_) {}

  void add(double x) {
    if (x < min || x > max) return;
    auto idx = (size_t)((x - min) / (max - min) * bins.size());
    if (idx >= bins.size()) idx = bins.size() - 1;
    bins[idx]++;
  }

  void merge(const Histogram &o) {
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

  // Depth histogram
  Histogram depth_histogram;

  /**
   * @brief Default-constructed empty tallies.
   */
  Tallies() = default;
  Tallies(size_t bin_count, double max_depth)
      : incidents(0), total_displacements(0), sputtered(0),
        range_stats(), energy_loss_stats(),
        depth_histogram(bin_count, 0.0, max_depth) {}

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
    depth_histogram.merge(o.depth_histogram);
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
    ss << "]";
    ss << "}";
    return ss.str();
  }
};

} // namespace common
} // namespace amcset
