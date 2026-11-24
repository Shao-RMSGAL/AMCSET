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
#include "views/plots_panel.h"

#include "views/plot_widget.h"

#include <QCheckBox>
#include <QClipboard>
#include <QGuiApplication>
#include <QComboBox>
#include <QCoreApplication>
#include <QDir>
#include <QFileDialog>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QListWidget>
#include <QPushButton>
#include <QSplitter>
#include <QVBoxLayout>

namespace amcset {
namespace gui {

namespace {

const QList<QColor> &palette_colors() {
  static const QList<QColor> c{QColor(213, 94, 0),  QColor(0, 114, 178),   QColor(0, 158, 115),
                               QColor(204, 121, 167), QColor(230, 159, 0), QColor(86, 180, 233),
                               QColor(120, 120, 120), QColor(170, 60, 60)};
  return c;
}

QString with_unit(const output::Column &c) {
  if (c.type != output::ColumnType::real || c.unit.empty() || c.unit == "1") return QString::fromStdString(c.label);
  return QStringLiteral("%1 (%2)").arg(QString::fromStdString(c.label), QString::fromStdString(c.unit));
}

QString unit_text(const output::Column &c) {
  if (c.unit.empty() || c.unit == "1") return {};
  return QString::fromStdString(c.unit);
}

} // namespace

PlotsPanel::PlotsPanel(QWidget *parent) : QWidget(parent) {
  auto *split = new QSplitter(Qt::Horizontal, this);
  auto *left = new QWidget(split);
  auto *ll = new QVBoxLayout(left);
  ll->setContentsMargins(0, 0, 0, 0);
  list_ = new QListWidget(left);
  list_->setMinimumWidth(180);
  ll->addWidget(list_, 1);

  // Custom plot builder.
  custom_box_ = new QWidget(left);
  auto *cf = new QFormLayout(custom_box_);
  cf->setContentsMargins(0, 4, 0, 0);
  custom_table_ = new QComboBox(custom_box_);
  custom_x_ = new QComboBox(custom_box_);
  custom_y_ = new QListWidget(custom_box_);
  custom_y_->setMaximumHeight(140);
  cf->addRow(tr("Table"), custom_table_);
  cf->addRow(tr("X"), custom_x_);
  cf->addRow(tr("Y"), custom_y_);
  ll->addWidget(custom_box_);
  custom_box_->hide();

  auto *opts = new QFormLayout();
  length_unit_ = new QComboBox(left);
  energy_unit_ = new QComboBox(left);
  length_unit_->addItem(tr("as output files"), QString());
  energy_unit_->addItem(tr("as output files"), QString());
  for (const auto *u : output::units_of(output::Dimension::length)) {
    length_unit_->addItem(QString::fromStdString(u->name), QString::fromStdString(u->name));
  }
  for (const auto *u : output::units_of(output::Dimension::energy)) {
    energy_unit_->addItem(QString::fromStdString(u->name), QString::fromStdString(u->name));
  }
  log_y_ = new QCheckBox(tr("Logarithmic"), left);
  opts->addRow(tr("Lengths"), length_unit_);
  opts->addRow(tr("Energies"), energy_unit_);
  opts->addRow(tr("Y axis"), log_y_);
  ll->addLayout(opts);
  auto *buttons = new QHBoxLayout();
  auto *save = new QPushButton(tr("Save image..."), left);
  auto *copy = new QPushButton(tr("Copy data"), left);
  buttons->addWidget(save);
  buttons->addWidget(copy);
  ll->addLayout(buttons);
  note_ = new QLabel(left);
  note_->setWordWrap(true);
  note_->setStyleSheet(QStringLiteral("color: palette(placeholder-text);"));
  ll->addWidget(note_);

  plot_ = new PlotWidget(split);
  split->addWidget(left);
  split->addWidget(plot_);
  split->setStretchFactor(1, 1);
  split->setSizes({230, 800});
  auto *lay = new QVBoxLayout(this);
  lay->setContentsMargins(0, 0, 0, 0);
  lay->addWidget(split);

  custom_.id = QStringLiteral("custom");
  custom_.title = tr("Custom plot");

  connect(list_, &QListWidget::currentRowChanged, this, guarded(this, tr("Plots (currentRowChanged)"), [this](int) {
    const auto *d = current();
    if (!d) return;
    selected_id_ = d->id;
    custom_box_->setVisible(d->id == custom_.id);
    {
      QSignalBlocker b(log_y_);
      log_y_->setChecked(d->log_y);
    }
    plot_->set_log_x(d->log_x);
    plot_->set_log_y(d->log_y);
    refresh(false);
  }));
  connect(log_y_, &QCheckBox::toggled, plot_, guarded(this, tr("Plots (log scale)"), [p = plot_](bool on) { p->set_log_y(on); }));
  connect(plot_, &PlotWidget::log_scale_changed, this, guarded(this, tr("Plots (log scale changed)"), [this](bool, bool ly) {
    QSignalBlocker b(log_y_);
    log_y_->setChecked(ly);
  }));
  connect(length_unit_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Plots (currentIndexChanged)"), [this] { refresh(false); }));
  connect(energy_unit_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Plots (currentIndexChanged)"), [this] { refresh(false); }));
  connect(plot_, &PlotWidget::hover_text, this, guarded(this, tr("Plots (hover text)"), [this](auto &&...a) { hover_text(std::forward<decltype(a)>(a)...); }));
  connect(save, &QPushButton::clicked, this, guarded(this, tr("Plots (clicked)"), [this] {
    const auto *d = current();
    const QString f = QFileDialog::getSaveFileName(this, tr("Save plot"),
                                                   (d ? d->id : QStringLiteral("plot")) + QStringLiteral(".png"),
                                                   tr("Images (*.png *.jpg *.bmp)"));
    if (!f.isEmpty()) plot_->save_image(f);
  }));
  connect(copy, &QPushButton::clicked, this, guarded(this, tr("Plots (clicked)"), [this] {
    QGuiApplication::clipboard()->setText(plot_->data_text());
  }));
  connect(custom_table_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Plots (currentIndexChanged)"), [this] {
    update_custom_columns();
    refresh(false);
  }));
  connect(custom_x_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Plots (currentIndexChanged)"), [this] { refresh(false); }));
  connect(custom_y_, &QListWidget::itemChanged, this, guarded(this, tr("Plots (itemChanged)"), [this] { refresh(false); }));
}

void PlotsPanel::set_output_units(const output::OutputSpec &spec) {
  spec_ = spec;
  refresh(true);
}

void PlotsPanel::set_surface_binding_energy(double eV) {
  sbe_ = eV;
  refresh(true);
}

void PlotsPanel::set_dataset(std::shared_ptr<const output::Dataset> d) {
  const bool tables_changed = [&] {
    if (!data_ || !d) return true;
    if (data_->tables.size() != d->tables.size()) return true;
    for (std::size_t i = 0; i < d->tables.size(); ++i) {
      if (data_->tables[i].key != d->tables[i].key ||
          data_->tables[i].columns.size() != d->tables[i].columns.size()) {
        return true;
      }
    }
    return false;
  }();
  data_ = std::move(d);
  if (tables_changed) rebuild_list();
  refresh(true);
}

void PlotsPanel::save_plot_images(const QString &dir) {
  const int keep = list_->currentRow();
  for (int i = 0; i < static_cast<int>(defs_.size()); ++i) {
    list_->setCurrentRow(i);
    QCoreApplication::processEvents();
    grab().save(QDir(dir).filePath(QStringLiteral("plot_%1.png").arg(defs_[static_cast<std::size_t>(i)].id)));
  }
  list_->setCurrentRow(keep);
}

void PlotsPanel::rebuild_list() {
  defs_.clear();
  const auto has = [&](const char *t) { return data_ && data_->find(t) != nullptr; };
  const auto add = [&](PlotDef d) { defs_.push_back(std::move(d)); };
  using S = PlotSeries::Style;
  if (has("range")) {
    add({"ions", tr("Ion distribution"), "range", "depth", {"ions"}, false, false, int(S::steps), true});
    if (data_->find("range")->columns.size() > 2) {
      add({"ions_recoils", tr("Ion and recoil distributions"), "range", "depth", {"ions", "recoils_*"},
           false, false, int(S::steps), false});
    }
  }
  if (has("lateral")) {
    add({"lateral", tr("Lateral spread of ions"), "lateral", "depth",
         {"lateral_range", "lateral_straggling", "radial_range", "radial_straggling"}, false, false,
         int(S::steps), false});
  }
  if (has("ionization")) {
    add({"ionization", tr("Ionization (energy to electrons)"), "ionization", "depth", {"by_ions", "by_recoils"},
         false, false, int(S::steps), true});
  }
  if (has("phonons")) {
    add({"phonons", tr("Phonons"), "phonons", "depth", {"by_ions", "by_recoils"}, false, false, int(S::steps), true});
  }
  if (has("energy_to_recoils")) {
    add({"e2recoil", tr("Energy transferred to recoils"), "energy_to_recoils", "depth", {"from_ions"}, false,
         false, int(S::steps), true});
    add({"absorbed", tr("Energy absorbed by target atoms"), "energy_to_recoils", "depth", {"absorbed_*"}, false,
         false, int(S::steps), false});
  }
  if (has("vacancies")) {
    add({"vacancies", tr("Vacancies"), "vacancies", "depth",
         {"by_ions", "by_recoils", "knock_ons", "vacancies_*", "total"}, false, false, int(S::steps), false});
  }
  if (has("dpa")) {
    add({"dpa", tr("Displacements per atom (dpa) by damage model"), "dpa", "depth",
         {"trim_vacancies", "nrt", "kinchin_pease", "arc_dpa"}, false, false, int(S::steps), false});
    add({"rpa", tr("Replacements per atom (rpa)"), "dpa", "depth", {"rpa"}, false, false, int(S::steps), false});
  }
  if (has("damage")) {
    add({"damage", tr("Collision events (damage)"), "damage", "depth",
         {"displacements", "vacancies", "replacements", "interstitials"}, false, false, int(S::steps), false});
  }
  if (has("replacements")) {
    add({"replacements", tr("Replacement collisions"), "replacements", "depth", {"replacements"}, false, false,
         int(S::steps), true});
  }
  if (has("energy_loss")) { // native engine
    add({"energy_loss", tr("Energy loss"), "energy_loss", "depth", {"electronic", "nuclear", "total"}, false,
         false, int(S::steps), false});
  }
  if (has("surface_energy")) {
    add({"sputter_integral", tr("Sputtering yield vs surface binding energy"), "surface_energy", "energy",
         {"yield_above_*"}, false, true, int(S::steps), false});
    add({"sputter_differential", tr("Atoms reaching the surface (differential)"), "surface_energy", "energy",
         {"arrivals_*"}, false, false, int(S::steps), true});
  }
  if (has("recoil_spectrum")) {
    add({"recoil_spectrum", tr("Energy spectrum of displaced atoms"), "recoil_spectrum", "energy_high",
         {"primary_knock_ons", "all_recoils"}, true, true, int(S::steps), false});
  }
  if (has("exyz")) {
    add({"exyz", tr("Ion energy vs depth (EXYZ)"), "exyz", "x", {"energy"}, false, false, int(S::points), false});
  }
  if (has("ion_positions")) {
    add({"positions", tr("Final ion positions (depth vs Y)"), "ion_positions", "x", {"y"}, false, false,
         int(S::points), false});
  }
  QSignalBlocker block(list_);
  list_->clear();
  int select = 0;
  for (std::size_t i = 0; i < defs_.size(); ++i) {
    list_->addItem(defs_[i].title);
    if (defs_[i].id == selected_id_) select = static_cast<int>(i);
  }
  list_->addItem(custom_.title);
  if (selected_id_ == custom_.id) select = static_cast<int>(defs_.size());
  // Custom tables.
  {
    QSignalBlocker b2(custom_table_);
    const QString keep = custom_table_->currentData().toString();
    custom_table_->clear();
    if (data_) {
      for (const auto &t : data_->tables) {
        custom_table_->addItem(QString::fromStdString(t.title.empty() ? t.key : t.title),
                               QString::fromStdString(t.key));
      }
    }
    const int k = custom_table_->findData(keep);
    if (k >= 0) custom_table_->setCurrentIndex(k);
  }
  update_custom_columns();
  list_->setCurrentRow(select);
  const auto *d = current();
  custom_box_->setVisible(d && d->id == custom_.id);
  if (d) {
    QSignalBlocker b(log_y_);
    log_y_->setChecked(d->log_y);
    plot_->set_log_x(d->log_x);
    plot_->set_log_y(d->log_y);
  }
}

void PlotsPanel::update_custom_columns() {
  QSignalBlocker bx(custom_x_);
  QSignalBlocker by(custom_y_);
  const QString keep_x = custom_x_->currentData().toString();
  QStringList keep_y;
  for (int i = 0; i < custom_y_->count(); ++i) {
    if (custom_y_->item(i)->checkState() == Qt::Checked) keep_y << custom_y_->item(i)->data(Qt::UserRole).toString();
  }
  custom_x_->clear();
  custom_y_->clear();
  if (!data_) return;
  const auto *t = data_->find(custom_table_->currentData().toString().toStdString());
  if (!t) return;
  for (const auto &c : t->columns) {
    if (c.type == output::ColumnType::text) continue;
    const QString key = QString::fromStdString(c.key);
    custom_x_->addItem(QString::fromStdString(c.label), key);
    auto *item = new QListWidgetItem(QString::fromStdString(c.label), custom_y_);
    item->setData(Qt::UserRole, key);
    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
    const bool on = keep_y.isEmpty() ? (custom_y_->count() == 2) : keep_y.contains(key);
    item->setCheckState(on ? Qt::Checked : Qt::Unchecked);
  }
  const int k = custom_x_->findData(keep_x);
  custom_x_->setCurrentIndex(k >= 0 ? k : 0);
}

const PlotsPanel::PlotDef *PlotsPanel::current() const {
  const int row = list_->currentRow();
  if (row < 0) return nullptr;
  if (row < static_cast<int>(defs_.size())) return &defs_[static_cast<std::size_t>(row)];
  return &custom_;
}

output::OutputSpec PlotsPanel::unit_spec() const {
  output::OutputSpec s;
  s.units = spec_.units;
  const QString l = length_unit_->currentData().toString();
  const QString e = energy_unit_->currentData().toString();
  if (!l.isEmpty()) s.units[output::Dimension::length] = l.toStdString();
  if (!e.isEmpty()) s.units[output::Dimension::energy] = e.toStdString();
  return s;
}

void PlotsPanel::refresh(bool keep_zoom) {
  const PlotDef *d = current();
  if (!d || !data_) {
    plot_->set_series({}, keep_zoom);
    plot_->set_title({});
    note_->clear();
    return;
  }
  PlotDef def = *d;
  if (def.id == custom_.id) {
    def.table = custom_table_->currentData().toString();
    def.x = custom_x_->currentData().toString();
    def.y.clear();
    for (int i = 0; i < custom_y_->count(); ++i) {
      if (custom_y_->item(i)->checkState() == Qt::Checked) def.y << custom_y_->item(i)->data(Qt::UserRole).toString();
    }
    def.style = static_cast<int>(PlotSeries::Style::line);
    def.title = custom_table_->currentText();
  }
  const auto *raw = data_->find(def.table.toStdString());
  if (!raw) {
    plot_->set_series({}, keep_zoom);
    return;
  }
  output::Table t;
  try {
    t = output::apply(*raw, unit_spec());
  } catch (const std::exception &) {
    t = *raw;
  }
  const auto *xc = t.find(def.x.toStdString());
  if (!xc) {
    plot_->set_series({}, keep_zoom);
    return;
  }
  QVector<PlotSeries> series;
  QString yunit;
  int color = 0;
  const auto add_column = [&](const output::Column &c) {
    PlotSeries s;
    s.name = QString::fromStdString(c.label);
    s.x = QVector<double>(xc->values.begin(), xc->values.end());
    s.y = QVector<double>(c.values.begin(), c.values.end());
    s.color = palette_colors()[color++ % palette_colors().size()];
    s.style = static_cast<PlotSeries::Style>(def.style);
    s.fill = def.fill;
    if (yunit.isEmpty()) yunit = unit_text(c);
    series.push_back(std::move(s));
  };
  for (const auto &pattern : def.y) {
    if (pattern.endsWith(QLatin1Char('*'))) {
      const std::string prefix = pattern.chopped(1).toStdString();
      for (const auto &c : t.columns) {
        if (c.type != output::ColumnType::text && c.key.rfind(prefix, 0) == 0) add_column(c);
      }
    } else if (const auto *c = t.find(pattern.toStdString()); c && c->type != output::ColumnType::text) {
      add_column(*c);
    }
  }
  plot_->set_title(def.title);
  QString ylabel = series.size() == 1 ? series.front().name : QString();
  if (!yunit.isEmpty()) ylabel += (ylabel.isEmpty() ? QString() : QStringLiteral(" ")) + QStringLiteral("(%1)").arg(yunit);
  plot_->set_axis_labels(with_unit(*xc), ylabel);
  QVector<QPair<double, QString>> markers;
  if ((def.id == "sputter_integral" || def.id == "sputter_differential") && sbe_ > 0) {
    double f = 1.0;
    try {
      f = output::conversion_factor(output::find_unit("eV"), output::find_unit(xc->unit));
    } catch (const std::exception &) {
    }
    markers.push_back({sbe_ * f, tr("Es = %1 eV").arg(sbe_)});
  }
  plot_->set_markers(markers);
  plot_->set_series(std::move(series), keep_zoom);
  if (def.id == "ions" || def.id == "ions_recoils") {
    note_->setText(tr("Multiply by the fluence (atoms/cm2) to get the concentration (atoms/cm3)."));
  } else if (def.id == "sputter_integral") {
    note_->setText(tr("Atoms per ion whose normal energy at the surface is at least the value on the x axis: "
                      "the sputtering yield the target would have with that surface binding energy."));
  } else {
    note_->setText(QString::fromStdString(raw->description));
  }
}

} // namespace gui
} // namespace amcset
