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

#include "setup/target_panel.h"

#include "app/guard.h"

#include "element_data.h"
#include "model/quantity_text.h"
#include "model/scenario.h"
#include "model/target_models.h"
#include "setup/periodic_table.h"

#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPushButton>
#include <QTableView>
#include <QToolButton>
#include <QVBoxLayout>

#include <map>

namespace amcset {
namespace gui {

using common::Element;

namespace {

// Chemical formula ("Fe2O3", "Ca5(PO4)3OH", "Fe0.91Cr0.09") -> atoms.
std::vector<std::pair<unsigned, double>> parse_formula(const QString &f) {
  std::vector<std::map<unsigned, double>> stack(1);
  int i = 0;
  const auto number = [&]() {
    const int s = i;
    while (i < f.size() && (f[i].isDigit() || f[i] == QLatin1Char('.'))) ++i;
    if (s == i) return 1.0;
    bool ok = false;
    const double v = f.mid(s, i - s).toDouble(&ok);
    if (!ok) throw std::invalid_argument("bad number in the formula");
    return v;
  };
  while (i < f.size()) {
    const QChar c = f[i];
    if (c.isSpace()) {
      ++i;
    } else if (c == QLatin1Char('(')) {
      stack.emplace_back();
      ++i;
    } else if (c == QLatin1Char(')')) {
      ++i;
      if (stack.size() < 2) throw std::invalid_argument("unbalanced ')'");
      const double n = number();
      auto top = stack.back();
      stack.pop_back();
      for (const auto &kv : top) stack.back()[kv.first] += kv.second * n;
    } else if (c.isUpper()) {
      int s = i++;
      while (i < f.size() && f[i].isLower()) ++i;
      const auto z = Element::by_symbol(f.mid(s, i - s).toStdString()).z();
      stack.back()[z] += number();
    } else {
      throw std::invalid_argument(std::string("unexpected '") + c.toLatin1() + "' in the formula");
    }
  }
  if (stack.size() != 1) throw std::invalid_argument("unbalanced '('");
  std::vector<std::pair<unsigned, double>> out(stack.front().begin(), stack.front().end());
  if (out.empty()) throw std::invalid_argument("empty formula");
  return out;
}

struct Compound {
  const char *name;
  const char *formula;
  double density; // g/cm3
  bool gas;
};

// A few common compounds (densities from the CRC Handbook).
const std::vector<Compound> &compounds() {
  static const std::vector<Compound> c{
      {"Silicon dioxide (quartz/glass)", "SiO2", 2.32, false},
      {"Aluminium oxide (sapphire)", "Al2O3", 3.97, false},
      {"Uranium dioxide", "UO2", 10.97, false},
      {"Silicon carbide", "SiC", 3.21, false},
      {"Silicon nitride", "Si3N4", 3.17, false},
      {"Gallium arsenide", "GaAs", 5.32, false},
      {"Gallium nitride", "GaN", 6.15, false},
      {"Indium phosphide", "InP", 4.81, false},
      {"Titanium dioxide (rutile)", "TiO2", 4.23, false},
      {"Zirconium dioxide", "ZrO2", 5.68, false},
      {"Magnesium oxide", "MgO", 3.58, false},
      {"Iron(III) oxide (hematite)", "Fe2O3", 5.24, false},
      {"Lithium fluoride", "LiF", 2.64, false},
      {"Boron carbide", "B4C", 2.52, false},
      {"Tungsten carbide", "WC", 15.63, false},
      {"Water", "H2O", 1.0, false},
      {"Polyethylene", "CH2", 0.93, false},
      {"PMMA (acrylic)", "C5H8O2", 1.18, false},
      {"Kapton (polyimide)", "C22H10N2O5", 1.42, false},
      {"Mylar (PET)", "C10H8O4", 1.40, false},
      {"FLiNaK (LiF-NaF-KF eutectic)", "Li0.465Na0.115K0.42F", 2.02, false},
      {"Carbon dioxide (gas)", "CO2", 1.977e-3, true},
      {"Air (gas)", "N0.78O0.21Ar0.01", 1.205e-3, true},
  };
  return c;
}

} // namespace

TargetPanel::TargetPanel(Scenario *scenario, QWidget *parent)
    : QWidget(parent), scenario_(scenario), layers_(new LayerTableModel(scenario, this)),
      elements_(new ElementTableModel(scenario, this)) {
  layer_view_ = new QTableView(this);
  layer_view_->setModel(layers_);
  layer_view_->setSelectionBehavior(QAbstractItemView::SelectRows);
  layer_view_->setSelectionMode(QAbstractItemView::SingleSelection);
  layer_view_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  layer_view_->horizontalHeader()->setSectionResizeMode(LayerTableModel::name, QHeaderView::Stretch);
  layer_view_->verticalHeader()->setDefaultSectionSize(layer_view_->fontMetrics().height() + 8);
  layer_view_->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed |
                               QAbstractItemView::AnyKeyPressed);
  element_view_ = new QTableView(this);
  element_view_->setModel(elements_);
  element_view_->setSelectionBehavior(QAbstractItemView::SelectRows);
  element_view_->setSelectionMode(QAbstractItemView::SingleSelection);
  element_view_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  element_view_->horizontalHeader()->setSectionResizeMode(ElementTableModel::element, QHeaderView::Stretch);
  element_view_->verticalHeader()->setDefaultSectionSize(element_view_->fontMetrics().height() + 8);
  element_view_->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed |
                                 QAbstractItemView::AnyKeyPressed);

  const auto button = [this](const QString &text, const QString &tip) {
    auto *b = new QToolButton(this);
    b->setText(text);
    b->setToolTip(tip);
    return b;
  };
  auto *add_l = button(tr("Add layer"), tr("Add a layer below the selected one"));
  auto *add_c = button(tr("Compound..."), tr("Add a layer from a chemical formula or the compound list"));
  auto *rem_l = button(tr("Remove"), tr("Remove the selected layer"));
  auto *up_l = button(tr("Up"), tr("Move the selected layer towards the surface"));
  auto *down_l = button(tr("Down"), tr("Move the selected layer deeper"));
  auto *add_e = button(tr("Add element"), tr("Add an element to the selected layer"));
  auto *rem_e = button(tr("Remove element"), tr("Remove the selected element"));

  auto *lrow = new QHBoxLayout();
  lrow->addWidget(add_l);
  lrow->addWidget(add_c);
  lrow->addWidget(rem_l);
  lrow->addWidget(up_l);
  lrow->addWidget(down_l);
  lrow->addStretch(1);
  auto *erow = new QHBoxLayout();
  erow->addWidget(add_e);
  erow->addWidget(rem_e);
  erow->addStretch(1);
  summary_ = new QLabel(this);
  summary_->setWordWrap(true);
  element_title_ = new QLabel(this);

  auto *lay = new QVBoxLayout(this);
  auto *heading = new QLabel(tr("<b>Layers</b> from the surface (double-click a cell to edit; values accept units)"), this);
  heading->setWordWrap(true);
  lay->addWidget(heading);
  lay->addWidget(layer_view_, 2);
  lay->addLayout(lrow);
  lay->addWidget(summary_);
  lay->addWidget(element_title_);
  lay->addWidget(element_view_, 2);
  lay->addLayout(erow);
  auto *note = new QLabel(tr("Ed: displacement, Eb: lattice binding, Es: surface binding energy. "
                             "Elements listed in several layers keep separate tallies."),
                          this);
  note->setWordWrap(true);
  lay->addWidget(note);

  connect(layers_, &LayerTableModel::edit_failed, this, &TargetPanel::message);
  connect(elements_, &ElementTableModel::edit_failed, this, &TargetPanel::message);
  // Every slot is guarded: no exception may reach Qt (app/guard.h).
  connect(layer_view_->selectionModel(), &QItemSelectionModel::currentRowChanged, this,
          guarded(this, tr("Selecting a layer"), [this](const QModelIndex &cur) { select_layer(cur.row()); }));
  connect(layers_, &QAbstractItemModel::modelReset, this, guarded(this, tr("Updating the layers"), [this] {
    const int n = layers_->rowCount();
    // (std::clamp needs lo <= hi: no layers gives n - 1 = -1.)
    select_layer(n > 0 ? std::clamp(elements_->layer(), 0, n - 1) : -1);
  }));
  connect(scenario_, &Scenario::changed, this, guarded_quiet(tr("Target summary"), [this] { update_summary(); }));
  connect(add_l, &QToolButton::clicked, this, guarded(this, tr("Add layer"), [this] { add_layer(); }));
  connect(add_c, &QToolButton::clicked, this, guarded(this, tr("Add compound layer"), [this] { add_compound(); }));
  connect(rem_l, &QToolButton::clicked, this, guarded(this, tr("Remove layer"), [this] { remove_layer(); }));
  connect(up_l, &QToolButton::clicked, this, guarded(this, tr("Move layer"), [this] { move_layer(-1); }));
  connect(down_l, &QToolButton::clicked, this, guarded(this, tr("Move layer"), [this] { move_layer(1); }));
  connect(add_e, &QToolButton::clicked, this, guarded(this, tr("Add element"), [this] { add_element(); }));
  connect(rem_e, &QToolButton::clicked, this, guarded(this, tr("Remove element"), [this] { remove_element(); }));
  select_layer(0);
  update_summary();
}

int TargetPanel::current_layer() const { return elements_->layer(); }

void TargetPanel::select_layer(int row) {
  const int n = static_cast<int>(scenario_->config().trim.layers.size());
  if (n == 0) {
    elements_->set_layer(-1);
    element_title_->setText(tr("<b>Elements</b>"));
    return;
  }
  row = std::clamp(row, 0, n - 1);
  if (row != elements_->layer()) elements_->set_layer(row);
  if (layer_view_->currentIndex().row() != row) {
    QSignalBlocker b(layer_view_->selectionModel());
    layer_view_->selectRow(row);
  }
  element_title_->setText(tr("<b>Elements of layer %1</b> (%2)")
                              .arg(row + 1)
                              .arg(QString::fromStdString(scenario_->config().trim.layers[static_cast<std::size_t>(row)].name)));
}

void TargetPanel::update_summary() {
  const auto &t = scenario_->config().trim;
  common::length_quantity total = 0.0 * common::angstrom;
  for (const auto &l : t.layers) total += l.width;
  summary_->setText(tr("%n layer(s), total thickness %1", nullptr, static_cast<int>(t.layers.size()))
                        .arg(format_length(total)));
  if (elements_->layer() >= 0 && elements_->layer() < static_cast<int>(t.layers.size())) {
    element_title_->setText(tr("<b>Elements of layer %1</b> (%2)")
                                .arg(elements_->layer() + 1)
                                .arg(QString::fromStdString(t.layers[static_cast<std::size_t>(elements_->layer())].name)));
  }
}

void TargetPanel::add_layer() {
  const unsigned z = PeriodicTableDialog::pick(this, 26);
  if (z == 0) return;
  const int at = current_layer() + 1;
  const QString err = scenario_->modify(tr("Add layer"), [&](config::RunConfig &c) {
    trim::TargetLayer l;
    const auto &el = Element::by_z(z);
    l.name = std::string(el.symbol());
    l.width = 1.0 * common::micrometer;
    l.gas = el.is_gas_at_stp();
    l.density = l.gas ? el.gas_density() : el.solid_density();
    l.elements.push_back(trim::TargetElement::defaults(z));
    const auto pos = static_cast<std::size_t>(std::clamp(at, 0, static_cast<int>(c.trim.layers.size())));
    c.trim.layers.insert(c.trim.layers.begin() + static_cast<std::ptrdiff_t>(pos), l);
  });
  if (!err.isEmpty()) emit message(err);
  select_layer(at);
}

void TargetPanel::add_compound() {
  QDialog d(this);
  d.setWindowTitle(tr("Add a compound layer"));
  auto *list = new QComboBox(&d);
  list->addItem(tr("(chemical formula)"));
  for (const auto &c : compounds()) list->addItem(QStringLiteral("%1  -  %2").arg(QString::fromLatin1(c.name), QString::fromLatin1(c.formula)));
  auto *formula = new QLineEdit(&d);
  formula->setPlaceholderText(tr("e.g. Fe2O3, Ca5(PO4)3OH, Fe0.91Cr0.09"));
  auto *density = new QLineEdit(&d);
  density->setPlaceholderText(tr("g/cm3"));
  auto *thickness = new QLineEdit(QStringLiteral("1 um"), &d);
  auto *form = new QFormLayout();
  form->addRow(tr("Compound"), list);
  form->addRow(tr("Formula"), formula);
  form->addRow(tr("Density"), density);
  form->addRow(tr("Thickness"), thickness);
  auto *bb = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &d);
  auto *lay = new QVBoxLayout(&d);
  lay->addLayout(form);
  lay->addWidget(new QLabel(tr("Stoichiometries come from the formula; set the energies of each element "
                               "and the Bragg correction afterwards if needed."),
                            &d));
  lay->addWidget(bb);
  bool gas = false;
  connect(list, &QComboBox::currentIndexChanged, &d, guarded(&d, tr("Choose a compound"), [&](int i) {
    if (i <= 0 || i > static_cast<int>(compounds().size())) return;
    const auto &c = compounds()[static_cast<std::size_t>(i - 1)];
    formula->setText(QString::fromLatin1(c.formula));
    density->setText(QString::number(c.density, 'g', 6));
    gas = c.gas;
  }));
  connect(bb, &QDialogButtonBox::accepted, &d, &QDialog::accept);
  connect(bb, &QDialogButtonBox::rejected, &d, &QDialog::reject);
  if (d.exec() != QDialog::Accepted) return;
  const int at = current_layer() + 1;
  const QString err = scenario_->modify(tr("Add compound layer"), [&](config::RunConfig &c) {
    const auto atoms = parse_formula(formula->text().trimmed());
    trim::TargetLayer l;
    l.name = formula->text().trimmed().toStdString();
    l.density = parse_density(density->text());
    if (!(l.density > 0.0 * common::gram_per_cubic_centimeter)) throw std::invalid_argument("enter the density");
    l.width = parse_length(thickness->text());
    l.gas = gas;
    for (const auto &[z, n] : atoms) l.elements.push_back(trim::TargetElement::defaults(z, n));
    const auto pos = static_cast<std::size_t>(std::clamp(at, 0, static_cast<int>(c.trim.layers.size())));
    c.trim.layers.insert(c.trim.layers.begin() + static_cast<std::ptrdiff_t>(pos), l);
  });
  if (!err.isEmpty()) {
    QMessageBox::warning(this, tr("Compound"), err);
    return;
  }
  select_layer(at);
}

void TargetPanel::remove_layer() {
  const int row = current_layer();
  if (row < 0) return;
  if (scenario_->config().trim.layers.size() <= 1) {
    emit message(tr("The target needs at least one layer."));
    return;
  }
  const QString err = scenario_->modify(tr("Remove layer"), [&](config::RunConfig &c) {
    if (row >= static_cast<int>(c.trim.layers.size())) throw std::out_of_range("no such layer");
    c.trim.layers.erase(c.trim.layers.begin() + row);
  });
  if (!err.isEmpty()) {
    emit message(err);
    return;
  }
  select_layer(row - 1);
}

void TargetPanel::move_layer(int delta) {
  const int row = current_layer();
  const int to = row + delta;
  const int n = static_cast<int>(scenario_->config().trim.layers.size());
  if (row < 0 || to < 0 || to >= n) return;
  const QString err = scenario_->modify(tr("Move layer"), [&](config::RunConfig &c) {
    std::swap(c.trim.layers.at(static_cast<std::size_t>(row)), c.trim.layers.at(static_cast<std::size_t>(to)));
  });
  if (!err.isEmpty()) {
    emit message(err);
    return;
  }
  select_layer(to);
}

void TargetPanel::add_element() {
  const int row = current_layer();
  if (row < 0) return;
  const unsigned z = PeriodicTableDialog::pick(this, 0);
  if (z == 0) return;
  const QString err = scenario_->modify(tr("Add element"), [&](config::RunConfig &c) {
    c.trim.layers.at(static_cast<std::size_t>(row)).elements.push_back(trim::TargetElement::defaults(z));
  });
  if (!err.isEmpty()) emit message(err);
}

void TargetPanel::remove_element() {
  const int row = current_layer();
  const int e = element_view_->currentIndex().row();
  const auto &layers = scenario_->config().trim.layers;
  if (row < 0 || e < 0 || row >= static_cast<int>(layers.size())) return;
  if (layers[static_cast<std::size_t>(row)].elements.size() <= 1) {
    emit message(tr("A layer needs at least one element."));
    return;
  }
  const QString err = scenario_->modify(tr("Remove element"), [&](config::RunConfig &c) {
    auto &els = c.trim.layers.at(static_cast<std::size_t>(row)).elements;
    if (e >= static_cast<int>(els.size())) throw std::out_of_range("no such element");
    els.erase(els.begin() + e);
  });
  if (!err.isEmpty()) emit message(err);
}

} // namespace gui
} // namespace amcset
