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

#include "views/stats_panel.h"

#include <QCoreApplication>
#include <QHeaderView>
#include <QTreeWidget>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <map>

namespace amcset {
namespace gui {

namespace {

QString group_of(const std::string &key) {
  static const std::map<std::string, const char *> g = {
      {"ions", "Calculation"},
      {"elapsed_time", "Calculation"},
      {"range_mean", "Ion ranges"},
      {"range_straggling", "Ion ranges"},
      {"range_skewness", "Ion ranges"},
      {"range_kurtosis", "Ion ranges"},
      {"lateral_range", "Ion ranges"},
      {"lateral_straggling", "Ion ranges"},
      {"radial_range", "Ion ranges"},
      {"radial_straggling", "Ion ranges"},
      {"stopped_ions", "Ion ranges"},
      {"backscattered_fraction", "Leaving the target"},
      {"backscattered_energy_fraction", "Leaving the target"},
      {"transmitted_fraction", "Leaving the target"},
      {"transmitted_energy_fraction", "Leaving the target"},
      {"sputtering_yield", "Leaving the target"},
      {"left_surface", "Leaving the target"},
      {"vacancies_per_ion", "Damage"},
      {"displacements_per_ion", "Damage"},
      {"replacements_per_ion", "Damage"},
      {"interstitials_per_ion", "Damage"},
      {"knock_ons_per_ion", "Damage"},
      {"ionization_ion_percent", "Energy loss (% of incident energy)"},
      {"ionization_recoils_percent", "Energy loss (% of incident energy)"},
      {"phonons_ion_percent", "Energy loss (% of incident energy)"},
      {"phonons_recoils_percent", "Energy loss (% of incident energy)"},
      {"vacancies_energy_percent", "Energy loss (% of incident energy)"},
  };
  const auto it = g.find(key);
  return QString::fromLatin1(it != g.end() ? it->second : "Other");
}

//! Short row names, read under their group heading (the full names, as in
//! the output files, are in the tooltips).
QString short_label(const std::string &key) {
  static const std::map<std::string, const char *> l = {
      {"range_mean", QT_TRANSLATE_NOOP("StatsPanel", "Projected range")},
      {"range_straggling", QT_TRANSLATE_NOOP("StatsPanel", "Straggling")},
      {"range_skewness", QT_TRANSLATE_NOOP("StatsPanel", "Skewness")},
      {"range_kurtosis", QT_TRANSLATE_NOOP("StatsPanel", "Kurtosis")},
      {"lateral_range", QT_TRANSLATE_NOOP("StatsPanel", "Lateral range")},
      {"lateral_straggling", QT_TRANSLATE_NOOP("StatsPanel", "Lateral straggling")},
      {"radial_range", QT_TRANSLATE_NOOP("StatsPanel", "Radial range")},
      {"radial_straggling", QT_TRANSLATE_NOOP("StatsPanel", "Radial straggling")},
      {"stopped_ions", QT_TRANSLATE_NOOP("StatsPanel", "Stopped in target")},
      {"backscattered_fraction", QT_TRANSLATE_NOOP("StatsPanel", "Backscattered")},
      {"backscattered_energy_fraction", QT_TRANSLATE_NOOP("StatsPanel", "Backscattered energy")},
      {"transmitted_fraction", QT_TRANSLATE_NOOP("StatsPanel", "Transmitted")},
      {"transmitted_energy_fraction", QT_TRANSLATE_NOOP("StatsPanel", "Transmitted energy")},
      {"sputtering_yield", QT_TRANSLATE_NOOP("StatsPanel", "Sputtering yield")},
      {"vacancies_per_ion", QT_TRANSLATE_NOOP("StatsPanel", "Vacancies / ion")},
      {"displacements_per_ion", QT_TRANSLATE_NOOP("StatsPanel", "Displacements / ion")},
      {"replacements_per_ion", QT_TRANSLATE_NOOP("StatsPanel", "Replacements / ion")},
      {"interstitials_per_ion", QT_TRANSLATE_NOOP("StatsPanel", "Interstitials / ion")},
      {"knock_ons_per_ion", QT_TRANSLATE_NOOP("StatsPanel", "Knock-ons / ion")},
      {"ionization_ion_percent", QT_TRANSLATE_NOOP("StatsPanel", "Ionization, ions")},
      {"ionization_recoils_percent", QT_TRANSLATE_NOOP("StatsPanel", "Ionization, recoils")},
      {"phonons_ion_percent", QT_TRANSLATE_NOOP("StatsPanel", "Phonons, ions")},
      {"phonons_recoils_percent", QT_TRANSLATE_NOOP("StatsPanel", "Phonons, recoils")},
      {"vacancies_energy_percent", QT_TRANSLATE_NOOP("StatsPanel", "Vacancy creation")},
  };
  const auto it = l.find(key);
  return it != l.end() ? QCoreApplication::translate("StatsPanel", it->second) : QString();
}

QString value_text(double v, const std::string &unit) {
  QString s = std::isfinite(v) ? QString::number(v, 'g', 5) : QStringLiteral("-");
  if (!unit.empty() && unit != "1") s += QLatin1Char(' ') + QString::fromStdString(unit);
  return s;
}

QString duration(double s) {
  if (!std::isfinite(s) || s < 0) return QStringLiteral("-");
  const auto t = static_cast<long long>(std::llround(s));
  if (t < 60) return QStringLiteral("%1 s").arg(t);
  if (t < 3600) return QStringLiteral("%1 min %2 s").arg(t / 60).arg(t % 60);
  return QStringLiteral("%1 h %2 min").arg(t / 3600).arg((t % 3600) / 60);
}

} // namespace

StatsPanel::StatsPanel(QWidget *parent) : QWidget(parent) {
  tree_ = new QTreeWidget(this);
  tree_->setColumnCount(2);
  tree_->setHeaderLabels({tr("Quantity"), tr("Value")});
  tree_->header()->setStretchLastSection(false);
  tree_->header()->setSectionResizeMode(0, QHeaderView::Stretch);
  tree_->header()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
  tree_->setRootIsDecorated(false);
  tree_->setUniformRowHeights(true);
  auto *lay = new QVBoxLayout(this);
  lay->setContentsMargins(0, 0, 0, 0);
  lay->addWidget(tree_);
  refresh();
}

QTreeWidgetItem *StatsPanel::group(const QString &name) {
  for (int i = 0; i < tree_->topLevelItemCount(); ++i) {
    if (tree_->topLevelItem(i)->text(0) == name) return tree_->topLevelItem(i);
  }
  auto *g = new QTreeWidgetItem(tree_, {name});
  QFont f = g->font(0);
  f.setBold(true);
  g->setFont(0, f);
  g->setFirstColumnSpanned(true);
  g->setExpanded(true);
  return g;
}

void StatsPanel::set_dataset(std::shared_ptr<const output::Dataset> d) {
  data_ = std::move(d);
  refresh();
}

void StatsPanel::set_output_units(const output::OutputSpec &spec) {
  spec_ = spec;
  refresh();
}

void StatsPanel::set_progress(std::size_t done, std::size_t total, double per_second, double elapsed) {
  if (!progress_) return;
  progress_->setText(1, total > 0 ? tr("%L1 of %L2 (%3 %)").arg(done).arg(total).arg(100.0 * done / total, 0, 'f', 1)
                                  : QString::number(done));
  rate_->setText(1, per_second > 0 ? tr("%1 /s").arg(per_second, 0, 'g', 4) : QStringLiteral("-"));
  const double left = per_second > 0 && total > done ? (total - done) / per_second : 0.0;
  remaining_->setText(1, per_second > 0 ? duration(left) : duration(elapsed));
  remaining_->setText(0, per_second > 0 ? tr("Time remaining") : tr("Calculation time"));
}

void StatsPanel::refresh() {
  // Rows: (key, group, label, value, tooltip).
  struct Row {
    QString key, group, label, value, tip;
  };
  std::vector<Row> rows;
  if (data_) {
    std::vector<output::Scalar> scalars;
    try {
      scalars = output::apply(data_->scalars, spec_);
    } catch (const std::exception &) {
      scalars = data_->scalars;
    }
    for (const auto &s : scalars) {
      if (s.key == "ions" || s.key == "elapsed_time") continue;
      // A label set in the output settings is shown as given; the built-in
      // one is replaced by its short form.
      const auto original = std::find_if(data_->scalars.begin(), data_->scalars.end(),
                                         [&](const output::Scalar &o) { return o.key == s.key; });
      const bool custom = original != data_->scalars.end() && original->label != s.label;
      const QString full = QString::fromStdString(s.label);
      const QString brief = short_label(s.key);
      QString tip = full;
      if (!s.description.empty()) tip += QLatin1Char('\n') + QString::fromStdString(s.description);
      rows.push_back({QString::fromStdString(s.key), group_of(s.key), custom || brief.isEmpty() ? full : brief,
                      value_text(s.value, s.unit), tip});
    }
    if (const auto *sp = data_->find("sputtering")) {
      const auto *el = sp->find("element");
      const auto *y = sp->find("yield");
      if (el && y && el->text.size() > 1) {
        for (std::size_t i = 0; i < el->text.size() && i < y->values.size(); ++i) {
          rows.push_back({QStringLiteral("sputter:%1").arg(i), tr("Leaving the target"),
                          tr("  sputtered %1 per ion").arg(QString::fromStdString(el->text[i])),
                          value_text(y->values[i], "1"), {}});
        }
      }
    }
  }
  QStringList keys;
  for (const auto &r : rows) keys << r.key;
  if (keys == keys_ && progress_ != nullptr) {
    for (std::size_t i = 0; i < rows.size(); ++i) values_[i]->setText(1, rows[i].value);
    return;
  }
  tree_->setUpdatesEnabled(false);
  tree_->clear();
  values_.clear();
  auto *calc = group(tr("Calculation"));
  progress_ = new QTreeWidgetItem(calc, {tr("Particles"), QStringLiteral("-")});
  rate_ = new QTreeWidgetItem(calc, {tr("Rate"), QStringLiteral("-")});
  remaining_ = new QTreeWidgetItem(calc, {tr("Time remaining"), QStringLiteral("-")});
  for (const auto &r : rows) {
    auto *it = new QTreeWidgetItem(group(r.group), {r.label, r.value});
    if (!r.tip.isEmpty()) it->setToolTip(0, r.tip);
    values_.push_back(it);
  }
  keys_ = keys;
  tree_->expandAll();
  tree_->setUpdatesEnabled(true);
}

} // namespace gui
} // namespace amcset
