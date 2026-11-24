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

#include "run/datasets.h"

#include "amcset_units.h"
#include "tallies.h"

#include <QFile>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include <cmath>

namespace amcset {
namespace gui {

using output::Column;
using output::Dimension;
using output::Table;

output::Dataset native_dataset(const common::Tallies &t, const QString &title) {
  output::Dataset d;
  d.title = title.toStdString();
  d.metadata = {{"program", "AMCSET (native engine)"},
                {"histories", std::to_string(t.incidents)}};
  const double n = std::max<double>(1.0, static_cast<double>(t.incidents));
  const auto bins_of = [](const common::Histogram &h) {
    std::vector<double> depth(h.bins.size());
    const double w = h.bins.empty() ? 0.0 : (h.max - h.min) / static_cast<double>(h.bins.size());
    for (std::size_t i = 0; i < h.bins.size(); ++i) {
      // Histogram limits are in metres: report Angstrom like SRIM mode.
      depth[i] = common::in_units((h.min + (static_cast<double>(i) + 1.0) * w) * common::meter,
                                  common::angstrom);
    }
    return std::make_pair(depth, w);
  };
  const auto per_A = [&](const common::Histogram &h) {
    const auto [depth, w_m] = bins_of(h);
    const double w_A = common::in_units(w_m * common::meter, common::angstrom);
    std::vector<double> v(h.bins.size(), 0.0);
    for (std::size_t i = 0; i < v.size(); ++i) v[i] = w_A > 0.0 ? h.bins[i] / (n * w_A) : 0.0;
    return std::make_pair(depth, v);
  };
  {
    Table tb;
    tb.key = "range";
    tb.title = "Range distribution";
    tb.description = "Final depth of the incident particles per unit depth and per particle.";
    auto [depth, v] = per_A(t.range_histogram);
    // per Angstrom -> per cm, like SRIM mode's (atoms/cm3)/(atoms/cm2)
    const double f = output::conversion_factor("1/A", "1/cm");
    for (auto &x : v) x *= f;
    tb.columns.push_back(Column::real("depth", "Depth", Dimension::length, "A", depth));
    tb.columns.push_back(Column::real("ions", "Incident particles", Dimension::inverse_length, "1/cm", v));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "damage";
    tb.title = "Displacements";
    tb.description = "Displaced target atoms per unit depth and per incident particle.";
    auto [depth, v] = per_A(t.displacement_histogram);
    tb.columns.push_back(Column::real("depth", "Depth", Dimension::length, "A", depth));
    tb.columns.push_back(Column::real("displacements", "Displacements", Dimension::inverse_length, "1/A", v));
    d.tables.push_back(std::move(tb));
  }
  {
    Table tb;
    tb.key = "energy_loss";
    tb.title = "Energy loss";
    tb.description = "Electronic, nuclear and total energy loss per unit depth and per incident particle.";
    const auto energy = [&](const common::Histogram &h) {
      auto [depth, v] = per_A(h);
      // Histogram values are joules.
      for (auto &x : v) x = common::in_units(x * common::joule, common::electron_volt);
      return std::make_pair(depth, v);
    };
    auto [depth, el] = energy(t.electronic_energy_loss_histogram);
    tb.columns.push_back(Column::real("depth", "Depth", Dimension::length, "A", depth));
    tb.columns.push_back(Column::real("electronic", "Electronic", Dimension::energy_per_length, "eV/A", el));
    auto nu = energy(t.nuclear_energy_loss_histogram).second;
    auto to = energy(t.total_energy_loss_histogram).second;
    if (nu.size() == el.size()) {
      tb.columns.push_back(Column::real("nuclear", "Nuclear", Dimension::energy_per_length, "eV/A", nu));
    }
    if (to.size() == el.size()) {
      tb.columns.push_back(Column::real("total", "Total", Dimension::energy_per_length, "eV/A", to));
    }
    d.tables.push_back(std::move(tb));
  }
  d.scalars.push_back({"ions", "Incident particles", Dimension::dimensionless, "1",
                       static_cast<double>(t.incidents), ""});
  d.scalars.push_back({"range_mean", "Mean range", Dimension::length, "A",
                       common::in_units(static_cast<double>(t.range_stats.mean) * common::meter,
                                        common::angstrom),
                       ""});
  d.scalars.push_back({"range_straggling", "Range straggling", Dimension::length, "A",
                       common::in_units(std::sqrt(static_cast<double>(t.range_stats.variance())) *
                                            common::meter,
                                        common::angstrom),
                       ""});
  d.scalars.push_back({"displacements_per_ion", "Displacements per particle", Dimension::dimensionless,
                       "1", static_cast<double>(t.total_displacements) / n, ""});
  d.scalars.push_back({"left_surface", "Particles leaving the surface", Dimension::dimensionless, "1",
                       static_cast<double>(t.sputtered), ""});
  return d;
}

bool read_results_json(const QString &path, output::Dataset &out, QString *error) {
  QFile f(path);
  if (!f.open(QIODevice::ReadOnly)) {
    if (error) *error = QStringLiteral("Cannot open %1").arg(path);
    return false;
  }
  QJsonParseError pe;
  const auto doc = QJsonDocument::fromJson(f.readAll(), &pe);
  if (doc.isNull() || !doc.isObject()) {
    if (error) *error = QStringLiteral("%1 is not a JSON object: %2").arg(path, pe.errorString());
    return false;
  }
  const auto root = doc.object();
  if (root.value("format").toString() != "amcset-results") {
    if (error) *error = QStringLiteral("%1 is not an AMCSET results file").arg(path);
    return false;
  }
  output::Dataset d;
  d.title = root.value("title").toString().toStdString();
  const auto meta = root.value("metadata").toObject();
  for (auto it = meta.begin(); it != meta.end(); ++it) {
    d.metadata.emplace_back(it.key().toStdString(), it.value().toString().toStdString());
  }
  if (root.contains("configuration")) {
    d.configuration_json =
        QJsonDocument(root.value("configuration").toObject()).toJson(QJsonDocument::Compact).toStdString();
  }
  try {
    for (const auto &sv : root.value("scalars").toArray()) {
      const auto s = sv.toObject();
      output::Scalar x;
      x.key = s.value("key").toString().toStdString();
      x.label = s.value("label").toString().toStdString();
      x.value = s.value("value").toDouble(std::nan(""));
      x.unit = s.value("unit").toString().toStdString();
      x.dimension = output::parse_dimension(s.value("dimension").toString().toStdString());
      x.description = s.value("description").toString().toStdString();
      d.scalars.push_back(std::move(x));
    }
    for (const auto &tv : root.value("tables").toArray()) {
      const auto to = tv.toObject();
      Table tb;
      tb.key = to.value("key").toString().toStdString();
      tb.title = to.value("title").toString().toStdString();
      tb.description = to.value("description").toString().toStdString();
      const auto data = to.value("data").toObject();
      for (const auto &cv : to.value("columns").toArray()) {
        const auto co = cv.toObject();
        Column c;
        c.key = co.value("key").toString().toStdString();
        c.label = co.value("label").toString().toStdString();
        c.description = co.value("description").toString().toStdString();
        const QString type = co.value("type").toString();
        const auto values = data.value(QString::fromStdString(c.key)).toArray();
        if (type == "text") {
          c.type = output::ColumnType::text;
          c.unit.clear();
          for (const auto &v : values) c.text.push_back(v.toString().toStdString());
        } else {
          c.type = type == "integer" ? output::ColumnType::integer : output::ColumnType::real;
          if (c.type == output::ColumnType::real) {
            c.unit = co.value("unit").toString().toStdString();
            c.dimension = output::parse_dimension(co.value("dimension").toString().toStdString());
            (void)output::find_unit(c.unit, c.dimension);
          }
          for (const auto &v : values) c.values.push_back(v.isNull() ? std::nan("") : v.toDouble());
        }
        tb.columns.push_back(std::move(c));
      }
      d.tables.push_back(std::move(tb));
    }
  } catch (const std::exception &e) {
    if (error) *error = QStringLiteral("%1: %2").arg(path, QString::fromUtf8(e.what()));
    return false;
  }
  out = std::move(d);
  return true;
}

} // namespace gui
} // namespace amcset
