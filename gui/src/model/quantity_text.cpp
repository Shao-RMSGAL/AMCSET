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

#include "model/quantity_text.h"

#include "element_data.h"
#include "run_config.h"

#include <stdexcept>

namespace amcset {
namespace gui {

namespace {

// Append the default unit to a bare number.
std::string with_unit(const QString &text, const char *default_unit) {
  const QString t = text.trimmed();
  bool ok = false;
  t.toDouble(&ok);
  return (ok ? t + QLatin1Char(' ') + QLatin1String(default_unit) : t).toStdString();
}

} // namespace

common::length_quantity parse_length(const QString &text, const char *u) {
  return config::parse_length(with_unit(text, u), "length");
}
common::energy_quantity parse_energy(const QString &text, const char *u) {
  return config::parse_energy(with_unit(text, u), "energy");
}
common::mass_density_quantity parse_density(const QString &text, const char *u) {
  return config::parse_density(with_unit(text, u), "density");
}
common::angle_quantity parse_angle(const QString &text, const char *u) {
  return config::parse_angle(with_unit(text, u), "angle");
}
common::mass_quantity parse_mass(const QString &text, const char *u) {
  return config::parse_mass(with_unit(text, u), "mass");
}

QString format_length(common::length_quantity v) { return QString::fromStdString(config::format_length(v)); }
QString format_energy(common::energy_quantity v) { return QString::fromStdString(config::format_energy(v)); }
QString format_density(common::mass_density_quantity v) {
  return QString::fromStdString(config::format_density(v));
}
QString format_number(double v, int significant) { return QString::number(v, 'g', significant); }

unsigned parse_element(const QString &text) {
  const QString t = text.trimmed();
  bool ok = false;
  const int z = t.toInt(&ok);
  if (ok) {
    if (z < 1 || z > static_cast<int>(common::Element::max_z)) {
      throw std::invalid_argument("atomic number must be 1-" + std::to_string(common::Element::max_z));
    }
    return static_cast<unsigned>(z);
  }
  // "Fe", "fe", "Fe (26)".
  QString sym = t.section(QLatin1Char(' '), 0, 0);
  if (sym.isEmpty()) throw std::invalid_argument("empty element");
  sym = sym.left(1).toUpper() + sym.mid(1).toLower();
  return common::Element::by_symbol(sym.toStdString()).z();
}

QString element_symbol(unsigned z) {
  if (z < 1 || z > common::Element::max_z) return QStringLiteral("?");
  return QString::fromStdString(std::string(common::Element::by_z(z).symbol()));
}

} // namespace gui
} // namespace amcset
