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

#include "app/guard.h"
#include "model/target_models.h"

#include "element_data.h"
#include "model/quantity_text.h"
#include "model/scenario.h"

#include <cmath>

namespace amcset {
namespace gui {

using common::electron_volt;
using common::gram_per_cubic_centimeter;

// ---------------------------------------------------------------------------
// Layers
// ---------------------------------------------------------------------------

LayerTableModel::LayerTableModel(Scenario *scenario, QObject *parent)
    : QAbstractTableModel(parent), scenario_(scenario) {
  rows_ = static_cast<int>(scenario_->config().trim.layers.size());
  connect(scenario_, &Scenario::changed, this, guarded(nullptr, tr("Target table (changed)"), [this] { sync(); }));
}

void LayerTableModel::sync() {
  const int n = static_cast<int>(scenario_->config().trim.layers.size());
  if (n != rows_) {
    beginResetModel();
    rows_ = n;
    endResetModel();
  } else if (n > 0) {
    emit dataChanged(index(0, 0), index(n - 1, column_count - 1));
  }
}

int LayerTableModel::rowCount(const QModelIndex &parent) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    return parent.isValid() ? 0 : rows_;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (rowCount)"), std::current_exception());
    return 0;
  }
}
int LayerTableModel::columnCount(const QModelIndex &parent) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    return parent.isValid() ? 0 : column_count;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (columnCount)"), std::current_exception());
    return 0;
  }
}

QVariant LayerTableModel::headerData(int section, Qt::Orientation o, int role) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (role != Qt::DisplayRole) return {};
    if (o == Qt::Vertical) return section + 1;
    switch (section) {
    case name: return tr("Layer");
    case thickness: return tr("Thickness");
    case density: return tr("Density");
    case gas: return tr("Gas");
    case compound_correction: return tr("Bragg corr.");
    }
    return {};
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (headerData)"), std::current_exception());
    return {};
  }
}

Qt::ItemFlags LayerTableModel::flags(const QModelIndex &index) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (!index.isValid()) return Qt::NoItemFlags;
    if (index.column() == gas) return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsUserCheckable;
    // The density of a crystal layer is its unit cell's (Crystal tab).
    const auto &layers = scenario_->config().trim.layers;
    if (index.column() == density && index.row() < static_cast<int>(layers.size()) &&
        layers[static_cast<std::size_t>(index.row())].crystal) {
      return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
    }
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (flags)"), std::current_exception());
    return Qt::NoItemFlags;
  }
}

QVariant LayerTableModel::data(const QModelIndex &index, int role) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    const auto &layers = scenario_->config().trim.layers;
    if (!index.isValid() || index.row() >= static_cast<int>(layers.size())) return {};
    const auto &l = layers[static_cast<std::size_t>(index.row())];
    if (role == Qt::CheckStateRole && index.column() == gas) return l.gas ? Qt::Checked : Qt::Unchecked;
    if (role == Qt::ToolTipRole) {
      switch (index.column()) {
      case thickness: return tr("Layer thickness, e.g. \"2000\" (Angstrom), \"1.5 um\", \"0.2 mm\"");
      case density:
        return l.crystal ? tr("The density of the layer's unit cell (set in the Crystal tab)")
                         : tr("Mass density, e.g. \"7.87\" (g/cm3)");
      case gas: return tr("Gas phase: TRIM's gas stopping for H, He, N, O, Ne, Ar, Kr, Xe "
                          "(other elements use the solid stopping)");
      case compound_correction: return tr("Multiplies the Bragg's-rule electronic stopping of this layer");
      default: break;
      }
      return {};
    }
    if (role != Qt::DisplayRole && role != Qt::EditRole) return {};
    switch (index.column()) {
    case name: return QString::fromStdString(l.name);
    case thickness: return format_length(l.width);
    case density:
      // A crystal's density is computed from its cell: six digits say it.
      return l.crystal ? QStringLiteral("%1 g/cm3").arg(format_number(
                             common::in_units(l.density, common::gram_per_cubic_centimeter), 6))
                       : format_density(l.density);
    case gas: return {};
    case compound_correction: return l.compound_correction;
    }
    return {};
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (data)"), std::current_exception());
    return {};
  }
}

bool LayerTableModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (!index.isValid()) return false;
    const auto row = static_cast<std::size_t>(index.row());
    const int col = index.column();
    if (col == gas && role != Qt::CheckStateRole) return false;
    if (col != gas && role != Qt::EditRole) return false;
    const QString text = value.toString();
    const QString err = scenario_->modify(
        tr("Edit layer %1").arg(index.row() + 1),
        [&](config::RunConfig &c) {
          auto &l = c.trim.layers.at(row);
          switch (col) {
          case name: l.name = text.trimmed().toStdString(); break;
          case thickness: {
            const auto w = parse_length(text);
            if (!(w > 0.0 * common::angstrom)) throw std::invalid_argument("thickness must be positive");
            l.width = w;
            break;
          }
          case density: {
            const auto d = parse_density(text);
            if (!(d > 0.0 * gram_per_cubic_centimeter)) throw std::invalid_argument("density must be positive");
            l.density = d;
            break;
          }
          case gas: {
            const bool on = value.toInt() == Qt::Checked;
            // An element's default density follows its phase, as in TRIM.
            if (l.elements.size() == 1) {
              const auto &el = common::Element::by_z(l.elements.front().z);
              const auto from = on ? el.solid_density() : el.gas_density();
              const auto to = on ? el.gas_density() : el.solid_density();
              if (std::abs((l.density / from).value() - 1.0) < 1e-6 && to > 0.0 * gram_per_cubic_centimeter) {
                l.density = to;
              }
            }
            l.gas = on;
            break;
          }
          case compound_correction: {
            bool ok = false;
            const double v = text.toDouble(&ok);
            if (!ok || !(v > 0.0)) throw std::invalid_argument("the Bragg correction must be a positive number");
            l.compound_correction = v;
            break;
          }
          }
        },
        1000 + index.row() * 10 + col);
    if (!err.isEmpty()) {
      emit edit_failed(err);
      return false;
    }
    return true;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (setData)"), std::current_exception());
    return false;
  }
}

// ---------------------------------------------------------------------------
// Elements of one layer
// ---------------------------------------------------------------------------

ElementTableModel::ElementTableModel(Scenario *scenario, QObject *parent)
    : QAbstractTableModel(parent), scenario_(scenario) {
  connect(scenario_, &Scenario::changed, this, guarded(nullptr, tr("Target table (changed)"), [this] { sync(); }));
}

void ElementTableModel::set_layer(int layer) {
  beginResetModel();
  layer_ = layer;
  const auto &layers = scenario_->config().trim.layers;
  rows_ = layer_ >= 0 && layer_ < static_cast<int>(layers.size())
              ? static_cast<int>(layers[static_cast<std::size_t>(layer_)].elements.size())
              : 0;
  endResetModel();
}

void ElementTableModel::sync() {
  const auto &layers = scenario_->config().trim.layers;
  if (layer_ >= static_cast<int>(layers.size())) {
    set_layer(layers.empty() ? -1 : static_cast<int>(layers.size()) - 1);
    return;
  }
  const int n = layer_ >= 0 ? static_cast<int>(layers[static_cast<std::size_t>(layer_)].elements.size()) : 0;
  if (n != rows_) {
    beginResetModel();
    rows_ = n;
    endResetModel();
  } else if (n > 0) {
    emit dataChanged(index(0, 0), index(n - 1, column_count - 1));
  }
}

int ElementTableModel::rowCount(const QModelIndex &parent) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    return parent.isValid() ? 0 : rows_;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (rowCount)"), std::current_exception());
    return 0;
  }
}
int ElementTableModel::columnCount(const QModelIndex &parent) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    return parent.isValid() ? 0 : column_count;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (columnCount)"), std::current_exception());
    return 0;
  }
}

QVariant ElementTableModel::headerData(int section, Qt::Orientation o, int role) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (role == Qt::ToolTipRole && o == Qt::Horizontal) {
      switch (section) {
      case stoichiometry: return tr("Relative number of atoms (normalised per layer)");
      case ed: return tr("Displacement energy: energy needed to displace the atom from its site");
      case eb: return tr("Lattice binding energy: lost by every displaced atom");
      case es: return tr("Surface binding energy: barrier for sputtering");
      default: return {};
      }
    }
    if (role != Qt::DisplayRole) return {};
    if (o == Qt::Vertical) return section + 1;
    switch (section) {
    case element: return tr("Element");
    case stoichiometry: return tr("Stoich.");
    case mass: return tr("Mass (amu)");
    case ed: return tr("Ed (eV)");
    case eb: return tr("Eb (eV)");
    case es: return tr("Es (eV)");
    }
    return {};
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (headerData)"), std::current_exception());
    return {};
  }
}

Qt::ItemFlags ElementTableModel::flags(const QModelIndex &index) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (!index.isValid()) return Qt::NoItemFlags;
    // The elements of a crystal layer and their numbers are its unit cell's
    // (Crystal tab); masses and energies stay editable.
    const auto &layers = scenario_->config().trim.layers;
    if ((index.column() == element || index.column() == stoichiometry) && layer_ >= 0 &&
        layer_ < static_cast<int>(layers.size()) && layers[static_cast<std::size_t>(layer_)].crystal) {
      return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
    }
    return Qt::ItemIsEnabled | Qt::ItemIsSelectable | Qt::ItemIsEditable;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (flags)"), std::current_exception());
    return Qt::NoItemFlags;
  }
}

QVariant ElementTableModel::data(const QModelIndex &index, int role) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    const auto &layers = scenario_->config().trim.layers;
    if (!index.isValid() || layer_ < 0 || layer_ >= static_cast<int>(layers.size())) return {};
    const auto &els = layers[static_cast<std::size_t>(layer_)].elements;
    if (index.row() >= static_cast<int>(els.size())) return {};
    const auto &e = els[static_cast<std::size_t>(index.row())];
    if (role == Qt::ToolTipRole && index.column() == element) {
      const auto &el = common::Element::by_z(e.z);
      return tr("%1 (Z = %2)").arg(QString::fromStdString(std::string(el.name()))).arg(e.z);
    }
    if (role != Qt::DisplayRole && role != Qt::EditRole) return {};
    switch (index.column()) {
    case element: return element_symbol(e.z);
    case stoichiometry: return e.stoichiometry;
    case mass: return format_number((e.mass / common::atomic_mass_unit).value());
    case ed: return format_number(common::in_units(e.displacement_energy, electron_volt));
    case eb: return format_number(common::in_units(e.lattice_binding_energy, electron_volt));
    case es: return format_number(common::in_units(e.surface_binding_energy, electron_volt));
    }
    return {};
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (data)"), std::current_exception());
    return {};
  }
}

bool ElementTableModel::setData(const QModelIndex &index, const QVariant &value, int role) {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (!index.isValid() || role != Qt::EditRole || layer_ < 0) return false;
    const auto layer = static_cast<std::size_t>(layer_);
    const auto row = static_cast<std::size_t>(index.row());
    const int col = index.column();
    const QString text = value.toString();
    const QString err = scenario_->modify(
        tr("Edit element"),
        [&](config::RunConfig &c) {
          auto &e = c.trim.layers.at(layer).elements.at(row);
          switch (col) {
          case element: {
            const unsigned z = parse_element(text);
            if (z == e.z) return;
            auto fresh = trim::TargetElement::defaults(z, e.stoichiometry);
            e = fresh; // new element: its own mass and energies
            break;
          }
          case stoichiometry: {
            bool ok = false;
            const double v = text.toDouble(&ok);
            if (!ok || !(v >= 0.0)) throw std::invalid_argument("stoichiometry must be a number >= 0");
            e.stoichiometry = v;
            break;
          }
          case mass: {
            const auto m = parse_mass(text);
            if (!(m > 0.0 * common::kilogram)) throw std::invalid_argument("mass must be positive");
            e.mass = m;
            break;
          }
          case ed:
          case eb:
          case es: {
            const auto v = parse_energy(text);
            if (!(v >= 0.0 * electron_volt)) throw std::invalid_argument("energies must be >= 0");
            (col == ed ? e.displacement_energy : col == eb ? e.lattice_binding_energy : e.surface_binding_energy) = v;
            break;
          }
          }
        },
        2000 + layer_ * 100 + index.row() * 10 + col);
    if (!err.isEmpty()) {
      emit edit_failed(err);
      return false;
    }
    return true;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Target table (setData)"), std::current_exception());
    return false;
  }
}

} // namespace gui
} // namespace amcset
