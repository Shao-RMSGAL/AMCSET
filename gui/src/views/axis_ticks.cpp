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

#include "views/axis_ticks.h"

#include "output/units.h"

#include <algorithm>
#include <cmath>

namespace amcset {
namespace gui {

std::vector<double> nice_ticks(double lo, double hi, int max_count) {
  std::vector<double> t;
  if (!(hi > lo) || !std::isfinite(lo) || !std::isfinite(hi)) return t;
  max_count = std::max(max_count, 2);
  const double raw = (hi - lo) / max_count;
  const double mag = std::pow(10.0, std::floor(std::log10(raw)));
  double step = mag;
  for (double m : {1.0, 2.0, 2.5, 5.0, 10.0}) {
    step = m * mag;
    if (step >= raw) break;
  }
  const double first = std::ceil(lo / step) * step;
  for (double v = first; v <= hi + 1e-9 * step; v += step) {
    t.push_back(std::abs(v) < 1e-12 * step ? 0.0 : v);
    if (t.size() > 100) break;
  }
  return t;
}

std::vector<double> log_ticks(double lo, double hi, int max_count) {
  std::vector<double> t;
  if (!(lo > 0) || !(hi > lo)) return t;
  const int d0 = static_cast<int>(std::floor(std::log10(lo)));
  const int d1 = static_cast<int>(std::ceil(std::log10(hi)));
  const int decades = d1 - d0;
  const int every = std::max(1, decades / std::max(max_count, 1) + 1);
  for (int d = d0; d <= d1; ++d) {
    if ((d - d0) % every != 0) continue;
    const double base = std::pow(10.0, d);
    if (base >= lo && base <= hi) t.push_back(base);
    if (decades <= 2) {
      for (double m : {2.0, 5.0}) {
        const double v = m * base;
        if (v >= lo && v <= hi) t.push_back(v);
      }
    }
  }
  std::sort(t.begin(), t.end());
  return t;
}

QString format_tick(double v) {
  if (v == 0.0) return QStringLiteral("0");
  const double a = std::abs(v);
  if (a >= 1e5 || a < 1e-3) return QString::number(v, 'g', 3).replace(QLatin1String("e+0"), QLatin1String("E"))
                                   .replace(QLatin1String("e-0"), QLatin1String("E-"))
                                   .replace(QLatin1String("e+"), QLatin1String("E"))
                                   .replace(QLatin1String("e-"), QLatin1String("E-"));
  return QString::number(v, 'g', 6);
}

AxisUnit length_axis_unit(double m) {
  if (m >= 1e7) return {QStringLiteral("mm"), 1e7};
  if (m >= 1e4) return {QStringLiteral("µm"), 1e4};
  if (m >= 100.0) return {QStringLiteral("nm"), 10.0};
  return {QStringLiteral("Å"), 1.0};
}

QString format_length_value(double a) {
  const auto u = length_axis_unit(std::abs(a));
  return QStringLiteral("%1 %2").arg(QString::number(a / u.factor, 'g', 4), u.name);
}

} // namespace gui
} // namespace amcset
