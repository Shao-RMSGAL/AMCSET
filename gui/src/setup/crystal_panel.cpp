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

#include "setup/crystal_panel.h"

#include "app/guard.h"
#include "crystal/cell_io.h"
#include "crystal/crystal_database.h"
#include "crystal/lattice_search.h"
#include "crystal/materials.h"
#include "element_data.h"
#include "model/quantity_text.h"
#include "model/scenario.h"
#include "setup/form_binder.h"
#include "setup/periodic_table.h"
#include "views/cell_view.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QAbstractItemView>
#include <QHeaderView>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMenu>
#include <QMessageBox>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QStyle>
#include <QTableWidget>
#include <QToolButton>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

namespace amcset {
namespace gui {

using common::Element;
using common::in_units;
using crystal::Mat3;
using crystal::Vec3;

namespace {

enum SiteColumn { site_label, site_x, site_y, site_z, site_kind, site_occupants, site_columns };
enum VibrationColumn { vib_element, vib_debye, vib_rms, vib_used, vib_columns };

QString view_of(std::string_view s) { return QString::fromUtf8(s.data(), static_cast<qsizetype>(s.size())); }

QString sym(unsigned z) { return view_of(Element::by_z(z).symbol()); }

QString angstrom_text(common::length_quantity v) {
  return QStringLiteral("%1 A").arg(format_number(in_units(v, common::angstrom), 6));
}

QString degree_text(common::angle_quantity v) {
  return QStringLiteral("%1 deg").arg(format_number(in_units(v, common::degree), 6));
}

//! Numbers separated by blanks, commas or semicolons.
std::vector<double> parse_numbers(const QString &text, const char *hint) {
  std::vector<double> v;
  QString part;
  const auto flush = [&] {
    if (part.isEmpty()) return;
    bool ok = false;
    const double x = part.toDouble(&ok);
    if (!ok || !std::isfinite(x)) throw std::invalid_argument("'" + part.toStdString() + "' is not a number (" + hint + ")");
    v.push_back(x);
    part.clear();
  };
  for (const QChar ch : text) {
    if (ch.isSpace() || ch == QLatin1Char(',') || ch == QLatin1Char(';')) flush();
    else part += ch;
  }
  flush();
  return v;
}

/*!
 * Miller indices as typed: "1 1 0", "[1,1,0]", "(0 0 0 1)", "<110>" or the
 * compact "110" and "1-10" (one digit per index). What follows the closing
 * bracket is a remark and ignored. Empty text gives no indices.
 */
std::vector<double> parse_indices(const QString &text) {
  QString t;
  for (const QChar ch : text.trimmed()) {
    if (ch == QLatin1Char(']') || ch == QLatin1Char(')') || ch == QLatin1Char('>') || ch == QLatin1Char('}')) break;
    if (ch == QLatin1Char('[') || ch == QLatin1Char('(') || ch == QLatin1Char('<') || ch == QLatin1Char('{')) continue;
    t += ch;
  }
  t = t.trimmed();
  std::vector<double> v;
  if (t.isEmpty()) return v;
  bool compact = true;
  int digits = 0;
  for (const QChar ch : t) {
    if (ch.isDigit()) ++digits;
    else if (ch != QLatin1Char('-')) compact = false;
  }
  if (compact && digits >= 3) {
    double sign = 1.0;
    for (const QChar ch : t) {
      if (ch == QLatin1Char('-')) {
        sign = -1.0;
      } else {
        v.push_back(sign * ch.digitValue());
        sign = 1.0;
      }
    }
    return v;
  }
  return parse_numbers(t, "write indices as 1 1 0");
}

//! A fractional coordinate: "0.25" or "1/3".
double parse_coordinate(const QString &text) {
  const QString t = text.trimmed();
  bool ok = false;
  const int slash = t.indexOf(QLatin1Char('/'));
  if (slash > 0) {
    bool ok2 = false;
    const double p = t.left(slash).toDouble(&ok);
    const double q = t.mid(slash + 1).toDouble(&ok2);
    if (ok && ok2 && q != 0.0) return p / q;
    throw std::invalid_argument("'" + t.toStdString() + "' is not a fraction (write 1/3 or 0.3333)");
  }
  const double v = t.toDouble(&ok);
  if (!ok || !std::isfinite(v)) throw std::invalid_argument("'" + t.toStdString() + "' is not a number");
  return v;
}

Vec3 parse_position(const QString &text) {
  QString t = text;
  t.replace(QLatin1Char(','), QLatin1Char(' ')).replace(QLatin1Char(';'), QLatin1Char(' '));
  const QStringList parts = t.split(QLatin1Char(' '), Qt::SkipEmptyParts);
  if (parts.size() != 3) throw std::invalid_argument("give three fractional coordinates, e.g. 1/2 1/2 0");
  return crystal::wrap_fractional(Vec3{parse_coordinate(parts[0]), parse_coordinate(parts[1]), parse_coordinate(parts[2])});
}

//! "Fe", "Fe 0.9, Cr 0.1", "Ga 50 %"; empty: a vacant site.
std::vector<crystal::Occupant> parse_occupants(const QString &text) {
  std::vector<crystal::Occupant> out;
  QString t = text;
  t.replace(QLatin1Char(';'), QLatin1Char(','));
  for (const QString &entry : t.split(QLatin1Char(','), Qt::SkipEmptyParts)) {
    QString e = entry.trimmed();
    if (e.isEmpty()) continue;
    bool percent = false;
    if (e.endsWith(QLatin1Char('%'))) {
      percent = true;
      e.chop(1);
      e = e.trimmed();
    }
    const QStringList parts = e.split(QLatin1Char(' '), Qt::SkipEmptyParts);
    if (parts.isEmpty() || parts.size() > 2) {
      throw std::invalid_argument("write occupants as \"Fe 0.9, Cr 0.1\" (element and fraction)");
    }
    crystal::Occupant o;
    o.z = parse_element(parts[0]);
    if (parts.size() == 2) {
      bool ok = false;
      o.fraction = parts[1].toDouble(&ok) * (percent ? 0.01 : 1.0);
      if (!ok) throw std::invalid_argument("'" + parts[1].toStdString() + "' is not a fraction");
    } else if (percent) {
      throw std::invalid_argument("a fraction is missing before '%'");
    }
    const auto same = std::find_if(out.begin(), out.end(), [&](const crystal::Occupant &x) { return x.z == o.z; });
    if (same != out.end()) same->fraction += o.fraction;
    else out.push_back(o);
  }
  return out;
}

QString occupants_text(const std::vector<crystal::Occupant> &occupants) {
  if (occupants.size() == 1 && occupants.front().fraction == 1.0) return sym(occupants.front().z);
  QString s;
  for (const auto &o : occupants) {
    s += (s.isEmpty() ? QString() : QStringLiteral(", ")) + QStringLiteral("%1 %2").arg(sym(o.z), format_number(o.fraction, 6));
  }
  return s;
}

//! The low-index direction nearest to \p v: "[1 1 0]", or "[1 1 8] + 0.35 deg".
QString direction_text(const crystal::Lattice &lattice, const Vec3 &v) {
  const auto n = crystal::nearest_direction(lattice, v, 12);
  const double off = in_units(n.off, common::degree);
  QString s = QStringLiteral("[%1 %2 %3]").arg(n.uvw[0]).arg(n.uvw[1]).arg(n.uvw[2]);
  if (off >= 0.005) s += QStringLiteral(" + %1 deg").arg(QString::number(off, 'f', 2));
  return s;
}

//! Masses for the density of \p cell: the layer's own, natural masses for new elements.
std::map<unsigned, common::mass_quantity> masses_for(const trim::TargetLayer &layer, const crystal::UnitCell &cell) {
  std::map<unsigned, common::mass_quantity> m;
  for (const auto &[z, atoms] : cell.composition()) {
    (void)atoms;
    const auto it = std::find_if(layer.elements.begin(), layer.elements.end(),
                                 [z = z](const trim::TargetElement &e) { return e.z == z; });
    m[z] = it != layer.elements.end() ? it->mass : trim::TargetElement::defaults(z).mass;
  }
  return m;
}

/*!
 * A name for a layer with this unit cell: the elements in the order of the
 * sites, with whole numbers of atoms where the cell has them ("GaAs", "UO2",
 * "BaTiO3"), or the symbols alone for a mixture ("FeCrC").
 */
std::string formula_of(const crystal::UnitCell &cell) {
  std::vector<unsigned> order;
  for (const auto &s : cell.sites) {
    for (const auto &o : s.occupants) {
      if (std::find(order.begin(), order.end(), o.z) == order.end()) order.push_back(o.z);
    }
  }
  const auto composition = cell.composition();
  double smallest = 0.0;
  for (const auto &[z, atoms] : composition) {
    (void)z;
    if (atoms > 0.0 && (smallest == 0.0 || atoms < smallest)) smallest = atoms;
  }
  bool whole = smallest > 0.0;
  for (const auto &[z, atoms] : composition) {
    (void)z;
    const double ratio = whole ? atoms / smallest : 0.0;
    if (std::abs(ratio - std::round(ratio)) > 1e-6 || ratio > 99.0) whole = false;
  }
  std::string name;
  for (const unsigned z : order) {
    name += std::string(Element::by_z(z).symbol());
    const auto it = composition.find(z);
    if (whole && it != composition.end()) {
      const long long n = std::llround(it->second / smallest);
      if (n > 1) name += std::to_string(n);
    }
  }
  return name;
}

//! The elements of a layer, sorted (to see whether an edit changed them).
std::vector<unsigned> elements_of(const trim::TargetLayer &layer) {
  std::vector<unsigned> z;
  for (const auto &e : layer.elements) z.push_back(e.z);
  std::sort(z.begin(), z.end());
  return z;
}

//! Let the list of a combo box show its entries in full, however narrow the box is.
void widen_popup(QComboBox *box) {
  auto *view = box->view();
  view->setTextElideMode(Qt::ElideNone);
  view->setMinimumWidth(view->sizeHintForColumn(0) + 2 * view->frameWidth() +
                        box->style()->pixelMetric(QStyle::PM_ScrollBarExtent) + 8);
}

//! A label made by AMCSET ("Si3") follows its element; any other label stays.
void relabel(crystal::Site &site, unsigned old_z, unsigned new_z) {
  const std::string old_symbol(Element::by_z(old_z).symbol());
  if (site.label.rfind(old_symbol, 0) != 0) return;
  const std::string rest = site.label.substr(old_symbol.size());
  if (!std::all_of(rest.begin(), rest.end(), [](char ch) { return ch >= '0' && ch <= '9'; })) return;
  site.label = std::string(Element::by_z(new_z).symbol()) + rest;
}

QTableWidgetItem *set_item(QTableWidget *table, int row, int column, const QString &text, bool editable = true) {
  auto *item = table->item(row, column);
  if (item == nullptr) {
    item = new QTableWidgetItem();
    table->setItem(row, column, item);
  }
  item->setText(text);
  item->setFlags(editable ? (Qt::ItemIsSelectable | Qt::ItemIsEnabled | Qt::ItemIsEditable)
                          : (Qt::ItemIsSelectable | Qt::ItemIsEnabled));
  return item;
}

//! What the "Structure..." dialog returns.
struct StructureChoice {
  crystal::StructureSpec spec;
  std::optional<common::mass_density_quantity> density; //!< instead of the lattice constant
};

/*!
 * Ask for a built-in structure: the elements of its sublattices, the lattice
 * constant (or the density), c/a and u where the structure has them, and the
 * occupants of its interstitial sites.
 */
std::optional<StructureChoice> ask_structure(QWidget *parent, const trim::TargetLayer *layer) {
  QDialog d(parent);
  d.setWindowTitle(QObject::tr("Build a crystal structure"));
  d.setMinimumWidth(480);
  const auto &structures = crystal::builtin_structures();
  auto *structure = new QComboBox(&d);
  for (const auto &s : structures) {
    structure->addItem(QStringLiteral("%1  -  %2  (%3)").arg(view_of(s.name), view_of(s.title), view_of(s.examples)));
  }
  constexpr int max_sublattices = 3;
  constexpr int max_interstitials = 2;
  QLabel *sub_label[max_sublattices];
  QLineEdit *sub_edit[max_sublattices];
  QToolButton *sub_pick[max_sublattices];
  QWidget *sub_row[max_sublattices];
  auto *form = new QFormLayout();
  form->addRow(QObject::tr("Structure"), structure);
  for (int i = 0; i < max_sublattices; ++i) {
    sub_label[i] = new QLabel(&d);
    sub_edit[i] = new QLineEdit(&d);
    sub_edit[i]->setPlaceholderText(QObject::tr("element, or a mixture: Fe 0.9, Cr 0.1"));
    sub_edit[i]->setToolTip(QObject::tr("The element on this sublattice. Several elements share its sites at random "
                                        "(a substitutional alloy): \"Fe 0.9, Cr 0.1\". Fractions below 1 in total "
                                        "leave vacancies."));
    sub_pick[i] = new QToolButton(&d);
    sub_pick[i]->setText(QObject::tr("PT"));
    sub_pick[i]->setToolTip(QObject::tr("Choose from the periodic table"));
    sub_row[i] = new QWidget(&d);
    auto *h = new QHBoxLayout(sub_row[i]);
    h->setContentsMargins(0, 0, 0, 0);
    h->addWidget(sub_edit[i], 1);
    h->addWidget(sub_pick[i]);
    form->addRow(sub_label[i], sub_row[i]);
    QObject::connect(sub_pick[i], &QToolButton::clicked, &d,
                     guarded(&d, QObject::tr("Choose an element"), [&d, edit = sub_edit[i]] {
                       const unsigned z = PeriodicTableDialog::pick(&d, 0);
                       if (z != 0) edit->setText(sym(z));
                     }));
  }
  auto *a = new QLineEdit(&d);
  a->setPlaceholderText(QObject::tr("e.g. 2.867 A"));
  a->setToolTip(QObject::tr("The lattice constant a. Leave it empty to give the density instead."));
  auto *density = new QLineEdit(&d);
  density->setPlaceholderText(QObject::tr("used when a is empty, e.g. 7.87 g/cm3"));
  density->setToolTip(QObject::tr("Without a lattice constant the cell takes the size that gives this density."));
  auto *c_over_a = new QLineEdit(&d);
  c_over_a->setPlaceholderText(QObject::tr("empty: ideal"));
  auto *c_label = new QLabel(QObject::tr("c/a"), &d);
  auto *u = new QLineEdit(&d);
  u->setPlaceholderText(QStringLiteral("0.375"));
  auto *u_label = new QLabel(QObject::tr("Internal parameter u"), &d);
  form->addRow(QObject::tr("Lattice constant a"), a);
  form->addRow(QObject::tr("or density"), density);
  form->addRow(c_label, c_over_a);
  form->addRow(u_label, u);
  QLabel *int_label[max_interstitials];
  QLineEdit *int_edit[max_interstitials];
  for (int i = 0; i < max_interstitials; ++i) {
    int_label[i] = new QLabel(&d);
    int_edit[i] = new QLineEdit(&d);
    int_edit[i]->setPlaceholderText(QObject::tr("empty, or e.g. C 0.01"));
    int_edit[i]->setToolTip(QObject::tr("Atoms on these interstitial sites, as element and the fraction of the "
                                        "sites that hold one: \"C 0.01\" puts a carbon atom on 1 % of them."));
    form->addRow(int_label[i], int_edit[i]);
  }
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &d);
  auto *lay = new QVBoxLayout(&d);
  lay->addLayout(form);
  lay->addWidget(buttons);
  QObject::connect(buttons, &QDialogButtonBox::accepted, &d, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &d, &QDialog::reject);
  const auto show_structure = [&](int index) {
    if (index < 0 || index >= static_cast<int>(structures.size())) return;
    const auto &s = structures[static_cast<std::size_t>(index)];
    for (int i = 0; i < max_sublattices; ++i) {
      const bool used = i < static_cast<int>(s.sublattices.size());
      sub_label[i]->setVisible(used);
      sub_row[i]->setVisible(used);
      if (used) {
        sub_label[i]->setText(s.sublattices.size() == 1 ? QObject::tr("Atoms")
                                                        : view_of(s.sublattices[static_cast<std::size_t>(i)]));
      }
    }
    const bool has_c = s.hexagonal || s.tetragonal;
    c_label->setVisible(has_c);
    c_over_a->setVisible(has_c);
    u_label->setVisible(s.has_u);
    u->setVisible(s.has_u);
    for (int i = 0; i < max_interstitials; ++i) {
      const bool used = i < static_cast<int>(s.interstitials.size());
      int_label[i]->setVisible(used);
      int_edit[i]->setVisible(used);
      if (used) {
        int_label[i]->setText(QObject::tr("%1 interstitial sites").arg(view_of(s.interstitials[static_cast<std::size_t>(i)])));
      }
    }
    d.adjustSize();
  };
  QObject::connect(structure, &QComboBox::currentIndexChanged, &d,
                   guarded(&d, QObject::tr("Choose a structure"), [&](int index) { show_structure(index); }));
  // Start from the layer: its elements on the sublattices, its density.
  if (layer != nullptr) {
    for (std::size_t i = 0; i < layer->elements.size() && i < static_cast<std::size_t>(max_sublattices); ++i) {
      sub_edit[i]->setText(sym(layer->elements[i].z));
    }
    if (layer->density > 0.0 * common::gram_per_cubic_centimeter) {
      density->setText(QStringLiteral("%1 g/cm3").arg(format_number(in_units(layer->density, common::gram_per_cubic_centimeter), 6)));
    }
  }
  show_structure(structure->currentIndex());
  while (d.exec() == QDialog::Accepted) {
    try {
      const auto &s = structures[static_cast<std::size_t>(structure->currentIndex())];
      StructureChoice choice;
      choice.spec.name = std::string(s.name);
      for (std::size_t i = 0; i < s.sublattices.size(); ++i) {
        auto occupants = parse_occupants(sub_edit[i]->text());
        if (occupants.empty()) throw std::invalid_argument("give the element of every sublattice");
        choice.spec.sublattices.push_back(std::move(occupants));
      }
      if (!a->text().trimmed().isEmpty()) {
        choice.spec.a_A = in_units(parse_length(a->text(), "A"), common::angstrom);
        if (!(choice.spec.a_A > 0.0)) throw std::invalid_argument("the lattice constant must be positive");
      } else if (!density->text().trimmed().isEmpty()) {
        choice.density = parse_density(density->text());
        if (!(*choice.density > 0.0 * common::gram_per_cubic_centimeter)) {
          throw std::invalid_argument("the density must be positive");
        }
      } else {
        throw std::invalid_argument("give the lattice constant a or the density");
      }
      if ((s.hexagonal || s.tetragonal) && !c_over_a->text().trimmed().isEmpty()) {
        choice.spec.c_over_a = parse_coordinate(c_over_a->text());
      }
      if (s.has_u && !u->text().trimmed().isEmpty()) choice.spec.u = parse_coordinate(u->text());
      for (std::size_t i = 0; i < s.interstitials.size() && i < static_cast<std::size_t>(max_interstitials); ++i) {
        auto occupants = parse_occupants(int_edit[i]->text());
        if (!occupants.empty()) choice.spec.interstitials.push_back({std::string(s.interstitials[i]), std::move(occupants)});
      }
      (void)crystal::build_structure(choice.spec); // refuses what cannot be built
      return choice;
    } catch (const std::exception &e) {
      QMessageBox::warning(&d, QObject::tr("Crystal structure"), QString::fromUtf8(e.what()));
    }
  }
  return std::nullopt;
}

} // namespace

CrystalPanel::CrystalPanel(Scenario *scenario, QWidget *parent)
    : QWidget(parent), scenario_(scenario), binder_(new FormBinder(scenario, this)) {
  const auto button = [](QWidget *owner, const QString &text, const QString &tip) {
    auto *b = new QToolButton(owner);
    b->setText(text);
    b->setToolTip(tip);
    return b;
  };

  // --- the layer --------------------------------------------------------------
  layer_ = new QComboBox(this);
  layer_->setToolTip(tr("The target layer whose crystal is shown (layers are defined in the Target tab)"));
  layer_->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  layer_->setMinimumContentsLength(8);
  kind_ = new QComboBox(this);
  kind_->addItems({tr("Amorphous (as TRIM)"), tr("Single crystal")});
  kind_->setToolTip(tr("Amorphous: atoms at random, as in every TRIM target. Single crystal: the atoms sit on the "
                       "sites of the unit cell below, and the layer's elements, stoichiometry and density are those "
                       "of the cell. Not TRIM: only the Crystal transport follows the lattice (channeling). "
                       "(Scenario key layers.N.crystal.)"));
  needs_layer_.push_back(kind_);
  auto *top = new QHBoxLayout();
  top->addWidget(new QLabel(tr("Layer"), this));
  top->addWidget(layer_, 1);
  top->addWidget(kind_);
  follow_ = new QCheckBox(tr("Crystal transport (follows the lattice)"), this);
  follow_->setToolTip(tr("On: the calculation follows the lattice of crystal layers, so that ions can channel. "
                         "Off: every layer is treated as amorphous, as in TRIM, with the unit cell's composition "
                         "and density. The same switch is in the Calculation tab. Not TRIM: it takes none of "
                         "TRIM.exe's bugs, and its division of the electronic stopping is not calibrated (see "
                         "below). Switching it on also switches the automatic target depth off, which this "
                         "transport does not offer."));
  transport_note_ = new QLabel(tr("The crystal transport is off: this calculation treats every layer as amorphous "
                                  "(as TRIM does), with the unit cell's composition and density."),
                               this);
  transport_note_->setWordWrap(true);
  match_ = new QLabel(this);
  match_->setWordWrap(true);
  match_->setTextFormat(Qt::RichText);
  match_->setToolTip(tr("The layer's elements and stoichiometry are those of a crystal of the database (the form "
                        "most common at room temperature). \"Give it this crystal\" makes the layer that single "
                        "crystal ((0 0 1) surface, 293 K); its density becomes the crystal's."));
  auto *note_row = new QVBoxLayout();
  note_row->addWidget(follow_);
  note_row->addWidget(transport_note_);
  note_row->addWidget(match_);

  // --- unit cell --------------------------------------------------------------
  auto *cell_group = new QGroupBox(tr("Unit cell"), this);
  material_ = new QComboBox(cell_group);
  material_->addItem(tr("Common crystals..."));
  for (const auto &r : crystal::reference_crystals()) {
    material_->addItem(QStringLiteral("%1  (%2, a = %3 A)").arg(view_of(r.name), view_of(r.structure), format_number(r.lattice_constant_A, 5)));
  }
  material_->setMaxVisibleItems(24);
  material_->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  material_->setMinimumContentsLength(10);
  material_->setToolTip(tr("Give the layer one of 43 common cubic crystals: structure and lattice constant at room "
                           "temperature and, for the 19 elements, the Debye temperature, from Gemmell, Rev. Mod. "
                           "Phys. 46, 129 (1974), Tables IV and V. The layer's elements and density are replaced."));
  auto *database = button(cell_group, tr("Database..."),
                          tr("Search the crystal database (%1 elements and compounds in the built-in structures, "
                             "with lattice constants, Debye temperatures where known, and the source of every "
                             "value) and give the layer one of them. The layer's elements and density are "
                             "replaced.")
                              .arg(crystal::crystal_database().size()));
  needs_layer_.push_back(database);
  connect(database, &QToolButton::clicked, this,
          guarded(this, tr("Crystal database"), [this] { choose_from_database({}); }));
  auto *structure = button(cell_group, tr("Structure..."),
                           tr("Build one of the built-in structures (bcc, fcc, hcp, diamond, zinc blende, wurtzite, "
                              "rock salt, CsCl, fluorite, perovskite, L1_2, L1_0, sc) from elements and a lattice "
                              "constant, with substitutional mixtures and interstitial atoms"));
  auto *import_button = button(cell_group, tr("Import..."),
                        tr("Read the unit cell from a CIF, VASP POSCAR/CONTCAR, LAMMPS data or extended XYZ file"));
  auto *export_button = button(cell_group, tr("Export..."),
                               tr("Write the unit cell as CIF, VASP POSCAR, LAMMPS data or extended XYZ (only CIF "
                                  "keeps partial occupancy)"));
  needs_layer_.insert(needs_layer_.end(), {material_, structure, import_button});
  needs_crystal_.push_back(export_button);
  widen_popup(material_);
  auto *source_row = new QHBoxLayout();
  source_row->addWidget(database);
  source_row->addWidget(structure);
  source_row->addWidget(import_button);
  source_row->addWidget(export_button);
  source_row->addStretch(1);

  QLineEdit *lattice[6];
  static const char *const lattice_name[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
  auto *grid = new QGridLayout();
  for (int i = 0; i < 6; ++i) {
    lattice[i] = new QLineEdit(cell_group);
    lattice[i]->setMinimumWidth(46);
    lattice[i]->setToolTip(i < 3 ? tr("Length of the lattice vector %1 (Angstrom unless a unit is given)").arg(QLatin1String(lattice_name[i]))
                                 : tr("Angle between the two other lattice vectors (alpha: b and c; beta: a and c; "
                                      "gamma: a and b), in degrees"));
    grid->addWidget(new QLabel(QLatin1String(lattice_name[i]), cell_group), i / 3, 2 * (i % 3));
    grid->addWidget(lattice[i], i / 3, 2 * (i % 3) + 1);
    needs_crystal_.push_back(lattice[i]);
  }
  keep_shape_ = new QCheckBox(tr("Keep the shape when a length changes"), cell_group);
  keep_shape_->setChecked(true);
  keep_shape_->setToolTip(tr("On: changing a, b or c scales the whole cell (a cubic cell stays cubic). Off: each "
                             "length and angle is set by itself."));
  needs_crystal_.push_back(keep_shape_);
  summary_ = new QLabel(cell_group);
  summary_->setWordWrap(true);
  summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);

  view_ = new CellView(cell_group);
  view_->setMinimumHeight(280);
  repeat_ = new QSpinBox(cell_group);
  repeat_->setRange(1, 4);
  repeat_->setToolTip(tr("Unit cells drawn along each lattice vector"));
  look_ = new QLineEdit(cell_group);
  look_->setPlaceholderText(tr("[u v w]"));
  look_->setMaximumWidth(90);
  look_->setToolTip(tr("Type a crystal direction, e.g. 1 1 0, and press Enter: the view turns to look along it "
                       "(the rows of atoms then line up and the channels show)"));
  auto *look_beam = button(cell_group, tr("Beam"), tr("Look along the beam: the crystal as the ions see it"));
  auto *look_normal = button(cell_group, tr("Normal"), tr("Look along the surface normal, into the target"));
  auto *look_snap = button(cell_group, tr("Snap"), tr("Turn the view to the nearest low-index direction"));
  auto *look_reset = button(cell_group, tr("Reset"), tr("The standard oblique view"));
  auto *view_row = new QHBoxLayout();
  view_row->addWidget(new QLabel(tr("Cells"), cell_group));
  view_row->addWidget(repeat_);
  view_row->addSpacing(8);
  view_row->addWidget(new QLabel(tr("View along"), cell_group));
  view_row->addWidget(look_);
  view_row->addStretch(1);
  auto *look_row = new QHBoxLayout();
  look_row->addWidget(new QLabel(tr("Look along the"), cell_group));
  look_row->addWidget(look_beam);
  look_row->addWidget(look_normal);
  look_row->addWidget(look_snap);
  look_row->addWidget(look_reset);
  look_row->addStretch(1);
  auto *beam_view = new QPushButton(tr("Point the beam along this view"), cell_group);
  beam_view->setToolTip(tr("Sets the beam's angle of incidence (Ion tab) and the crystal's twist so that the ions "
                           "travel the way you are looking at the cell. The view must look into the surface."));
  auto *orient_view = new QPushButton(tr("Turn the crystal to this view"), cell_group);
  orient_view->setToolTip(tr("Sets the crystal's orientation to the view: the direction you look along becomes the "
                             "surface normal (you look into the target), and the top of the view the lateral Y axis, "
                             "towards which the beam tilts. A view along the normal of a low-index plane (type a "
                             "direction, or use Snap) gives that plane as the surface, in Miller indices; any other "
                             "view gives Euler angles. The beam's angle of incidence stays as it is."));
  auto *hint = new QLabel(tr("Drag: turn. Shift+drag: roll. Wheel: zoom. Click an atom: replace its element, share "
                             "the site, set its occupancy, make it interstitial, or remove it."),
                          cell_group);
  hint->setWordWrap(true);
  for (QWidget *w : std::initializer_list<QWidget *>{view_, repeat_, look_, look_beam, look_normal, look_snap, look_reset, beam_view, orient_view}) {
    needs_crystal_.push_back(w);
  }

  sites_ = new QTableWidget(0, site_columns, cell_group);
  sites_->setHorizontalHeaderLabels({tr("Site"), tr("x"), tr("y"), tr("z"), tr("Kind"), tr("Occupants")});
  sites_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  sites_->horizontalHeader()->setSectionResizeMode(site_occupants, QHeaderView::Stretch);
  sites_->verticalHeader()->setDefaultSectionSize(sites_->fontMetrics().height() + 8);
  sites_->setSelectionBehavior(QAbstractItemView::SelectRows);
  sites_->setSelectionMode(QAbstractItemView::SingleSelection);
  sites_->setMinimumHeight(150);
  sites_->setToolTip(tr("The sites of the unit cell. x, y, z: fractional coordinates (0.25 or 1/4). Kind: lattice "
                        "(a site of the host lattice; elements sharing it are substitutional) or interstitial. "
                        "Occupants: \"Fe\", or \"Fe 0.9, Cr 0.1\" for elements sharing the site at random; "
                        "fractions below 1 in total leave vacancies; empty: a vacant site."));
  auto *add_site_button = button(cell_group, tr("Add site..."), tr("Add an atom at given fractional coordinates"));
  auto *remove_site_button = button(cell_group, tr("Remove site"), tr("Remove the selected site"));
  auto *site_row = new QHBoxLayout();
  site_row->addWidget(add_site_button);
  site_row->addWidget(remove_site_button);
  site_row->addStretch(1);
  needs_crystal_.insert(needs_crystal_.end(), {sites_, add_site_button, remove_site_button});

  auto *cell_lay = new QVBoxLayout(cell_group);
  cell_lay->addWidget(material_);
  cell_lay->addLayout(source_row);
  cell_lay->addLayout(grid);
  cell_lay->addWidget(keep_shape_);
  cell_lay->addWidget(summary_);
  cell_lay->addWidget(view_, 1);
  cell_lay->addLayout(view_row);
  cell_lay->addLayout(look_row);
  auto *apply_row = new QHBoxLayout();
  apply_row->addWidget(beam_view);
  apply_row->addWidget(orient_view);
  cell_lay->addLayout(apply_row);
  cell_lay->addWidget(hint);
  cell_lay->addWidget(sites_);
  cell_lay->addLayout(site_row);

  // --- orientation ------------------------------------------------------------
  auto *orientation_group = new QGroupBox(tr("Orientation and beam"), this);
  auto *given_as = new QComboBox(orientation_group);
  given_as->addItems({tr("Miller indices"), tr("Euler angles (Bunge)")});
  given_as->setToolTip(tr("How the crystal sits in the target: by the plane of the surface and a direction in it, "
                          "or by Bunge Euler angles relative to the sample (texture data)"));
  auto *surface = new QLineEdit(orientation_group);
  surface->setToolTip(tr("The crystal plane of the surface, (h k l): its normal is the depth axis. Hexagonal cells "
                         "take four indices (h k i l). Example: 0 0 1 for a (100) wafer. "
                         "(layers.N.crystal.orientation.surface)"));
  auto *reference = new QLineEdit(orientation_group);
  reference->setPlaceholderText(tr("automatic"));
  reference->setToolTip(tr("A crystal direction [u v w] whose part in the surface lies along the lateral Y axis, "
                           "the axis the beam tilts towards. Empty: the lattice vector a (or b). "
                           "(layers.N.crystal.orientation.reference)"));
  auto *euler = new QLineEdit(orientation_group);
  euler->setToolTip(tr("Bunge Euler angles phi1, Phi, phi2 in degrees (layers.N.crystal.orientation.euler)"));
  auto *twist = new QLineEdit(orientation_group);
  twist->setToolTip(tr("Turns the crystal about the surface normal, as when a wafer is rotated on its holder "
                       "(right-handed about the depth axis). With the tilt of the beam it decides which rows and "
                       "planes the ions meet. (layers.N.crystal.orientation.twist)"));
  beam_ = new QLineEdit(orientation_group);
  beam_->setToolTip(tr("The crystal direction the ions travel along. Type one, e.g. 1 1 0: the beam's angle of "
                       "incidence (Ion tab) and the twist are set to it, and the view turns to it. A beam along a "
                       "low-index direction is channeled; implanters avoid that with about 7 deg tilt and 22 deg "
                       "twist."));
  directions_ = new QLabel(orientation_group);
  directions_->setWordWrap(true);
  directions_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  auto *of = new QFormLayout(orientation_group);
  of->addRow(tr("Given as"), given_as);
  of->addRow(tr("Surface (h k l)"), surface);
  of->addRow(tr("In the surface [u v w]"), reference);
  of->addRow(tr("Euler angles"), euler);
  of->addRow(tr("Twist"), twist);
  of->addRow(tr("Beam along [u v w]"), beam_);
  of->addRow(directions_);
  needs_crystal_.insert(needs_crystal_.end(), {given_as, surface, reference, euler, twist, beam_});

  // --- temperature and vibrations ---------------------------------------------
  auto *thermal_group = new QGroupBox(tr("Temperature and thermal vibrations"), this);
  auto *temperature = new QLineEdit(thermal_group);
  temperature->setToolTip(tr("Temperature of the crystal, e.g. 300 K or 25 C. With the Debye temperatures it gives "
                             "the thermal displacements of the atoms, which limit channeling. "
                             "(layers.N.crystal.temperature)"));
  vibrations_ = new QTableWidget(0, vib_columns, thermal_group);
  vibrations_->setHorizontalHeaderLabels({tr("Element"), tr("Debye temp."), tr("rms"), tr("In use")});
  vibrations_->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
  vibrations_->horizontalHeader()->setMinimumSectionSize(64);
  vibrations_->horizontalHeader()->setSectionResizeMode(vib_used, QHeaderView::Stretch);
  vibrations_->verticalHeader()->setVisible(false);
  vibrations_->verticalHeader()->setDefaultSectionSize(vibrations_->fontMetrics().height() + 8);
  vibrations_->setToolTip(tr("Thermal vibrations of each element: a Debye temperature (the amplitude then follows "
                             "from the temperature and the mass, Debye model) or the rms amplitude along one axis "
                             "itself, which wins. Empty: a pure element of Gemmell's table (Rev. Mod. Phys. 46, 129 "
                             "(1974)) takes its tabulated Debye temperature; any other atom stands still, which "
                             "makes channeling too strong. \"In use\" is what the crystal transport takes. "
                             "(layers.N.crystal.vibrations)"));
  auto *tf = new QFormLayout(thermal_group);
  tf->addRow(tr("Temperature"), temperature);
  tf->addRow(vibrations_);
  needs_crystal_.insert(needs_crystal_.end(), {temperature, vibrations_});

  // --- the crystal transport ----------------------------------------------------
  auto *transport_group = new QGroupBox(tr("Crystal transport (not TRIM)"), this);
  auto *fraction = new QLineEdit(transport_group);
  fraction->setToolTip(tr("The fraction of the electronic stopping that is lost in the collisions, falling off with "
                          "the distance of closest approach; the rest is lost along the path. Channeled ions avoid "
                          "close collisions and so lose less: this decides how deep they go. 0: no dependence on "
                          "the path; 1: everything in the collisions. Not calibrated: measured boron profiles in "
                          "silicon are compatible with 0.25 to 1 (docs/CLI.md \"Crystals\"). Default 0.5. "
                          "(physics.crystal.local_stopping_fraction)"));
  auto *decay = new QLineEdit(transport_group);
  decay->setToolTip(tr("How fast that part falls off: exp(-decay R/a) of the distance of closest approach R, with "
                       "the screening length a. Default 0.3, the slowest exponential of the Moliere screening "
                       "function. (physics.crystal.local_stopping_decay)"));
  auto *impact = new QLineEdit(transport_group);
  impact->setToolTip(tr("Largest impact parameter of a collision. \"auto\": for each crystal the radius of the "
                        "largest sphere that fits between its atoms, the smallest value that leaves no part of "
                        "the lattice free of collisions. (physics.crystal.max_impact_parameter)"));
  auto *simultaneous = new QLineEdit(transport_group);
  simultaneous->setToolTip(tr("Atoms whose closest approaches lie within this distance along the path are struck "
                              "together (as MARLOWE does). Default 0.5 A. (physics.crystal.simultaneous_distance)"));
  auto *vibrate = new QCheckBox(tr("Thermal vibrations"), transport_group);
  vibrate->setToolTip(tr("Off: every atom stands exactly on its site whatever the crystals say (a static lattice). "
                         "Channeling is then much stronger than in a real crystal. "
                         "(physics.crystal.thermal_vibrations)"));
  auto *spot = new QComboBox(transport_group);
  spot->addItems({tr("Point"), tr("Rectangle"), tr("Circle"), tr("Gaussian")});
  spot->setToolTip(tr("Where the ions meet the surface. In a crystal the entry point matters (on a row of atoms or "
                      "between two): each ion enters at a random point of the spot. Point: every ion at the same "
                      "point, the origin of the lattice. Results stay relative to each ion's entry point. "
                      "(ion.spot.shape; TRIM's beam is a point on an amorphous target)"));
  auto *spot_width = new QLineEdit(transport_group);
  spot_width->setToolTip(tr("Width of the spot along Y (rectangle), its diameter (circle) or its full width at "
                            "half maximum (Gaussian). (ion.spot.width)"));
  auto *spot_height = new QLineEdit(transport_group);
  spot_height->setToolTip(tr("Height of the spot along Z (rectangle, Gaussian). (ion.spot.height)"));
  auto *pf = new QFormLayout(transport_group);
  pf->addRow(tr("Stopping in collisions"), fraction);
  pf->addRow(tr("Its decay"), decay);
  pf->addRow(tr("Impact parameter up to"), impact);
  pf->addRow(tr("Simultaneous within"), simultaneous);
  pf->addRow(QString(), vibrate);
  auto *auto_crystal = new QCheckBox(tr("Amorphous layers of a database crystal's composition become that crystal"),
                                     transport_group);
  auto_crystal->setToolTip(tr("When the calculation starts, every layer without a crystal whose elements and "
                              "stoichiometry are those of a crystal of the database (Fe: bcc iron; Ga + As 1:1: "
                              "zinc-blende GaAs; the form most common at room temperature) is followed as that "
                              "single crystal, (0 0 1) surface, 293 K; the log says which, and the layer's density "
                              "becomes the crystal's. Off: such layers stay amorphous. "
                              "(physics.crystal.auto_crystal)"));
  pf->addRow(QString(), auto_crystal);
  binder_->check(auto_crystal, tr("Automatic crystals"), [](const config::RunConfig &c) { return c.trim.crystal.auto_crystal; },
                 [](config::RunConfig &c, bool on) { c.trim.crystal.auto_crystal = on; });
  pf->addRow(tr("Beam spot"), spot);
  pf->addRow(tr("Spot width"), spot_width);
  pf->addRow(tr("Spot height"), spot_height);

  auto *lay = new QVBoxLayout(this);
  lay->addLayout(top);
  lay->addLayout(note_row);
  lay->addWidget(cell_group);
  lay->addWidget(orientation_group);
  lay->addWidget(thermal_group);
  lay->addWidget(transport_group);
  lay->addStretch(1);

  // --- bindings ---------------------------------------------------------------
  connect(binder_, &FormBinder::error, this,
          guarded(this, tr("Crystal settings (error)"), [this](const QString &m) { emit message(m); }));
  // The layer list comes first: every field below shows the crystal of the
  // layer it selects (the binder refreshes in the order of registration).
  binder_->custom([this](const config::RunConfig &c) { refresh_layers(c); });
  for (int i = 0; i < 6; ++i) {
    binder_->text(lattice[i], tr("Lattice parameter %1").arg(QLatin1String(lattice_name[i])),
                  [this, i](const config::RunConfig &) {
                    const auto *cr = current_crystal();
                    if (cr == nullptr) return QString();
                    const auto &l = cr->cell.lattice;
                    switch (i) {
                    case 0: return angstrom_text(l.a());
                    case 1: return angstrom_text(l.b());
                    case 2: return angstrom_text(l.c());
                    case 3: return degree_text(l.alpha());
                    case 4: return degree_text(l.beta());
                    default: return degree_text(l.gamma());
                    }
                  },
                  [this, i](config::RunConfig &c, const QString &t) {
                    auto &layer = layer_in(c);
                    auto &cr = crystal_in(c);
                    const auto &l = cr.cell.lattice;
                    common::length_quantity len[3] = {l.a(), l.b(), l.c()};
                    common::angle_quantity ang[3] = {l.alpha(), l.beta(), l.gamma()};
                    if (i < 3) {
                      const auto v = parse_length(t, "A");
                      if (!(v > 0.0 * common::angstrom)) throw std::invalid_argument("must be positive");
                      if (keep_shape_->isChecked()) {
                        cr.cell.lattice = l.scaled(in_units(v, common::angstrom) / in_units(len[i], common::angstrom));
                      } else {
                        len[i] = v;
                        cr.cell.lattice = crystal::Lattice::from_parameters(len[0], len[1], len[2], ang[0], ang[1], ang[2]);
                      }
                    } else {
                      ang[i - 3] = parse_angle(t);
                      cr.cell.lattice = crystal::Lattice::from_parameters(len[0], len[1], len[2], ang[0], ang[1], ang[2]);
                    }
                    trim::apply_crystal(layer);
                  });
  }
  binder_->combo(given_as, tr("Orientation"),
                 [this](const config::RunConfig &) {
                   const auto *cr = current_crystal();
                   return cr != nullptr && cr->orientation.kind == crystal::OrientationKind::euler ? 1 : 0;
                 },
                 [this](config::RunConfig &c, int i) {
                   if (current_crystal() == nullptr) return;
                   crystal_in(c).orientation.kind = i == 1 ? crystal::OrientationKind::euler : crystal::OrientationKind::miller;
                 });
  binder_->text(surface, tr("Surface plane"),
                [this](const config::RunConfig &) {
                  const auto *cr = current_crystal();
                  return cr != nullptr ? QString::fromStdString(crystal::format_indices(cr->orientation.surface, '(', ')')) : QString();
                },
                [this](config::RunConfig &c, const QString &t) {
                  auto v = parse_indices(t);
                  if (v.empty()) throw std::invalid_argument("give the plane of the surface, e.g. 0 0 1");
                  crystal_in(c).orientation.surface = std::move(v);
                });
  binder_->text(reference, tr("Direction in the surface"),
                [this](const config::RunConfig &) {
                  const auto *cr = current_crystal();
                  return cr != nullptr && !cr->orientation.reference.empty()
                             ? QString::fromStdString(crystal::format_indices(cr->orientation.reference))
                             : QString();
                },
                [this](config::RunConfig &c, const QString &t) { crystal_in(c).orientation.reference = parse_indices(t); });
  binder_->text(euler, tr("Euler angles"),
                [this](const config::RunConfig &) {
                  const auto *cr = current_crystal();
                  if (cr == nullptr) return QString();
                  const auto &e = cr->orientation.euler_deg;
                  return QStringLiteral("%1, %2, %3").arg(format_number(e[0]), format_number(e[1]), format_number(e[2]));
                },
                [this](config::RunConfig &c, const QString &t) {
                  QString s = t;
                  s.remove(QStringLiteral("deg"));
                  const auto v = parse_numbers(s, "write the angles as 30, 45, 0");
                  if (v.size() != 3) throw std::invalid_argument("give three angles in degrees: phi1, Phi, phi2");
                  crystal_in(c).orientation.euler_deg = {v[0], v[1], v[2]};
                });
  binder_->text(twist, tr("Twist"),
                [this](const config::RunConfig &) {
                  const auto *cr = current_crystal();
                  return cr != nullptr ? degree_text(cr->orientation.twist) : QString();
                },
                [this](config::RunConfig &c, const QString &t) { crystal_in(c).orientation.twist = parse_angle(t); });
  binder_->text(temperature, tr("Crystal temperature"),
                [this](const config::RunConfig &) {
                  const auto *cr = current_crystal();
                  return cr != nullptr ? QString::fromStdString(config::format_temperature(cr->temperature)) : QString();
                },
                [this](config::RunConfig &c, const QString &t) {
                  crystal_in(c).temperature = config::parse_temperature(t.toStdString(), "temperature");
                });
  binder_->text(fraction, tr("Stopping lost in collisions"),
                [](const config::RunConfig &c) { return format_number(c.trim.crystal.local_stopping_fraction); },
                [](config::RunConfig &c, const QString &t) { c.trim.crystal.local_stopping_fraction = parse_coordinate(t); });
  binder_->text(decay, tr("Decay of the stopping in collisions"),
                [](const config::RunConfig &c) { return format_number(c.trim.crystal.local_stopping_decay); },
                [](config::RunConfig &c, const QString &t) { c.trim.crystal.local_stopping_decay = parse_coordinate(t); });
  binder_->text(impact, tr("Largest impact parameter"),
                [](const config::RunConfig &c) {
                  const auto p = c.trim.crystal.max_impact_parameter;
                  return p > 0.0 * common::angstrom ? angstrom_text(p) : QStringLiteral("auto");
                },
                [](config::RunConfig &c, const QString &t) {
                  const QString s = t.trimmed().toLower();
                  c.trim.crystal.max_impact_parameter =
                      s.isEmpty() || s.startsWith(QLatin1String("auto")) ? 0.0 * common::angstrom : parse_length(t, "A");
                });
  binder_->text(simultaneous, tr("Simultaneous collisions"),
                [](const config::RunConfig &c) { return angstrom_text(c.trim.crystal.simultaneous_distance); },
                [](config::RunConfig &c, const QString &t) { c.trim.crystal.simultaneous_distance = parse_length(t, "A"); });
  binder_->check(vibrate, tr("Thermal vibrations"),
                 [](const config::RunConfig &c) { return c.trim.crystal.thermal_vibrations; },
                 [](config::RunConfig &c, bool on) { c.trim.crystal.thermal_vibrations = on; });
  binder_->combo(spot, tr("Beam spot"), [](const config::RunConfig &c) { return static_cast<int>(c.trim.ion.spot.shape); },
                 [](config::RunConfig &c, int i) { c.trim.ion.spot.shape = static_cast<trim::BeamSpot::Shape>(i); });
  binder_->text(spot_width, tr("Spot width"), [](const config::RunConfig &c) { return format_length(c.trim.ion.spot.width); },
                [](config::RunConfig &c, const QString &t) { c.trim.ion.spot.width = parse_length(t, "um"); });
  binder_->text(spot_height, tr("Spot height"), [](const config::RunConfig &c) { return format_length(c.trim.ion.spot.height); },
                [](config::RunConfig &c, const QString &t) { c.trim.ion.spot.height = parse_length(t, "um"); });
  binder_->custom([this, spot_width, spot_height, transport_group, surface, reference, euler](const config::RunConfig &c) {
    const auto shape = c.trim.ion.spot.shape;
    spot_width->setEnabled(shape != trim::BeamSpot::Shape::point);
    spot_height->setEnabled(shape == trim::BeamSpot::Shape::rectangle || shape == trim::BeamSpot::Shape::gaussian);
    transport_group->setEnabled(c.trim.transport == trim::TrimConfig::Transport::crystal);
    refresh(c);
    // Miller indices or Euler angles: only the fields that apply.
    if (const auto *cr = current_crystal()) {
      const bool by_euler = cr->orientation.kind == crystal::OrientationKind::euler;
      surface->setEnabled(!by_euler);
      reference->setEnabled(!by_euler);
      euler->setEnabled(by_euler);
    }
  });

  connect(layer_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Choosing a layer"), [this] {
    view_->set_selected_site(-1);
    binder_->refresh();
  }));
  connect(kind_, &QComboBox::activated, this, guarded(this, tr("Crystal or amorphous"), [this](int i) { set_kind(i); }));
  connect(follow_, &QCheckBox::clicked, this,
          guarded(this, tr("Crystal transport"), [this](bool on) { set_crystal_transport(on); }));
  connect(material_, &QComboBox::activated, this,
          guarded(this, tr("Common crystal"), [this](int i) { choose_material(i); }));
  connect(match_, &QLabel::linkActivated, this, guarded(this, tr("Database crystal"), [this](const QString &) {
    const auto *l = current_layer();
    if (l == nullptr || l->crystal) return;
    std::vector<std::pair<unsigned, double>> composition;
    for (const auto &e : l->elements) composition.emplace_back(e.z, e.stoichiometry);
    if (const auto *row = crystal::match_composition(composition)) {
      const auto made = crystal::make_database_crystal(*row);
      const std::string layer_name = l->name; // (l does not survive the change)
      set_cell(tr("Crystal %1").arg(view_of(row->name)), made.cell, &made.vibrations, layer_name, std::nullopt);
    }
  }));
  connect(structure, &QToolButton::clicked, this, guarded(this, tr("Crystal structure"), [this] { build_structure(); }));
  connect(import_button, &QToolButton::clicked, this, guarded(this, tr("Import a unit cell"), [this] { import_cell(); }));
  connect(export_button, &QToolButton::clicked, this, guarded(this, tr("Export the unit cell"), [this] { export_cell(); }));
  connect(repeat_, &QSpinBox::valueChanged, this,
          guarded(this, tr("Unit cell view"), [this](int n) { view_->set_repeat(n); }));
  connect(look_, &QLineEdit::returnPressed, this, guarded(this, tr("View direction"), [this] { look_from_text(); }));
  connect(look_beam, &QToolButton::clicked, this, guarded(this, tr("View along the beam"), [this] {
    if (const auto *cr = current_crystal()) {
      view_->look_along(crystal::beam_in_crystal(cr->cell.lattice, cr->orientation, scenario_->config().trim.ion.angle));
    }
  }));
  connect(look_normal, &QToolButton::clicked, this, guarded(this, tr("View along the normal"), [this] {
    if (const auto *cr = current_crystal()) view_->look_along(cr->rotation()[0]);
  }));
  connect(look_snap, &QToolButton::clicked, this, guarded(this, tr("Snap the view"), [this] {
    if (const auto *cr = current_crystal()) {
      const auto n = crystal::nearest_direction(cr->cell.lattice, view_->view_direction(), 6);
      view_->look_along(crystal::direction_vector(cr->cell.lattice, Vec3{double(n.uvw[0]), double(n.uvw[1]), double(n.uvw[2])}));
    }
  }));
  connect(look_reset, &QToolButton::clicked, this, guarded(this, tr("Reset the view"), [this] { view_->reset_view(); }));
  connect(beam_view, &QPushButton::clicked, this, guarded(this, tr("Beam along the view"), [this] { beam_along_view(); }));
  connect(orient_view, &QPushButton::clicked, this,
          guarded(this, tr("Crystal orientation from the view"), [this] { orient_to_view(); }));
  connect(beam_, &QLineEdit::editingFinished, this, guarded(this, tr("Beam direction"), [this] { beam_from_text(); }));
  connect(view_, &CellView::site_clicked, this,
          guarded(this, tr("Editing a site"), [this](int site, const QPoint &pos) { site_menu(site, pos); }));
  connect(sites_, &QTableWidget::cellChanged, this,
          guarded(this, tr("Editing a site"), [this](int row, int column) { site_edited(row, column); }));
  connect(sites_, &QTableWidget::currentCellChanged, this,
          guarded_quiet(tr("Selecting a site"), [this](int row) { view_->set_selected_site(row); }));
  connect(add_site_button, &QToolButton::clicked, this, guarded(this, tr("Add a site"), [this] { add_site(); }));
  connect(remove_site_button, &QToolButton::clicked, this, guarded(this, tr("Remove a site"), [this] { remove_site(); }));
  connect(vibrations_, &QTableWidget::cellChanged, this,
          guarded(this, tr("Thermal vibrations"), [this](int row, int column) { vibration_edited(row, column); }));
  binder_->refresh();
}

// --- the selected layer ---------------------------------------------------------

int CrystalPanel::layer() const {
  const int n = static_cast<int>(scenario_->config().trim.layers.size());
  const int i = layer_->currentIndex();
  return n == 0 ? -1 : std::clamp(i, 0, n - 1);
}

const trim::TargetLayer *CrystalPanel::current_layer() const {
  const int i = layer();
  return i < 0 ? nullptr : &scenario_->config().trim.layers[static_cast<std::size_t>(i)];
}

const crystal::Crystal *CrystalPanel::current_crystal() const {
  const auto *l = current_layer();
  return l != nullptr && l->crystal ? &*l->crystal : nullptr;
}

trim::TargetLayer &CrystalPanel::layer_in(config::RunConfig &c) const {
  const int i = layer();
  if (i < 0 || i >= static_cast<int>(c.trim.layers.size())) throw std::invalid_argument("the target has no such layer");
  return c.trim.layers[static_cast<std::size_t>(i)];
}

crystal::Crystal &CrystalPanel::crystal_in(config::RunConfig &c) const {
  auto &l = layer_in(c);
  if (!l.crystal) throw std::invalid_argument("the layer is amorphous: make it a single crystal first");
  return *l.crystal;
}

bool CrystalPanel::edit(const QString &what, const Edit &f) {
  const QString err = scenario_->modify(what, [&](config::RunConfig &c) {
    auto &l = layer_in(c);
    const auto before = elements_of(l);
    f(crystal_in(c), l, c);
    trim::apply_crystal(l);
    // Another material: the layer is named after it (the name can be
    // changed in the Target tab).
    if (elements_of(l) != before) l.name = formula_of(l.crystal->cell);
  });
  if (err.isEmpty()) return true;
  emit message(tr("%1: %2").arg(what, err));
  binder_->refresh();
  return false;
}

bool CrystalPanel::set_cell(const QString &what, const crystal::UnitCell &cell,
                            const std::vector<crystal::Vibration> *vibrations, const std::string &layer_name,
                            std::optional<common::mass_density_quantity> density) {
  const QString err = scenario_->modify(what, [&](config::RunConfig &c) {
    auto &l = layer_in(c);
    if (l.gas) throw std::invalid_argument("a gas layer cannot be a crystal (clear \"Gas\" in the Target tab first)");
    crystal::Crystal cr;
    cr.cell = cell;
    if (density) cr.cell.scale_to_density(*density, masses_for(l, cr.cell));
    const auto composition = cr.cell.composition();
    if (l.crystal) {
      // The same wafer with another cell: orientation and temperature stay,
      // and so do the vibrations of the elements that are still there.
      cr.orientation = l.crystal->orientation;
      cr.temperature = l.crystal->temperature;
      cr.polycrystal = l.crystal->polycrystal; // (no editor for the grains yet: they are kept)
      if (vibrations == nullptr) {
        for (const auto &v : l.crystal->vibrations) {
          if (composition.count(v.z) != 0) cr.vibrations.push_back(v);
        }
      }
      // An orientation the new lattice cannot take (four indices on a cubic
      // cell, ...) starts afresh.
      try {
        (void)crystal::rotation(cr.cell.lattice, cr.orientation);
      } catch (const std::exception &) {
        cr.orientation = {};
      }
    }
    if (vibrations != nullptr) cr.vibrations = *vibrations;
    const auto before = elements_of(l);
    l.crystal = std::move(cr);
    trim::apply_crystal(l);
    // Another material: the layer is named after it.
    if (!layer_name.empty()) l.name = layer_name;
    else if (elements_of(l) != before) l.name = formula_of(l.crystal->cell);
  });
  if (err.isEmpty()) return true;
  QMessageBox::warning(this, what, err);
  binder_->refresh();
  return false;
}

// --- showing the scenario -------------------------------------------------------

void CrystalPanel::refresh_layers(const config::RunConfig &c) {
  QStringList names;
  for (std::size_t i = 0; i < c.trim.layers.size(); ++i) {
    const auto &l = c.trim.layers[i];
    names << QStringLiteral("%1: %2 (%3)")
                 .arg(i + 1)
                 .arg(QString::fromStdString(l.name),
                      !l.crystal ? tr("amorphous") : l.crystal->polycrystal ? tr("polycrystal") : tr("crystal"));
  }
  QStringList shown;
  for (int i = 0; i < layer_->count(); ++i) shown << layer_->itemText(i);
  if (names == shown) return;
  const QSignalBlocker block(layer_);
  const bool other_target = layer_->count() != names.size();
  int current = layer_->currentIndex();
  layer_->clear();
  layer_->addItems(names);
  if (names.isEmpty()) return;
  widen_popup(layer_);
  current = std::clamp(current, 0, static_cast<int>(names.size()) - 1);
  // Another target (a file was opened, a layer added or removed): show its
  // first crystal.
  if (other_target) {
    for (std::size_t i = 0; i < c.trim.layers.size(); ++i) {
      if (c.trim.layers[i].crystal) {
        current = static_cast<int>(i);
        break;
      }
    }
  }
  layer_->setCurrentIndex(current);
}

void CrystalPanel::refresh(const config::RunConfig &c) {
  refresh_layers(c);
  const int row = layer();
  const trim::TargetLayer *l = row >= 0 ? &c.trim.layers[static_cast<std::size_t>(row)] : nullptr;
  const crystal::Crystal *cr = l != nullptr && l->crystal ? &*l->crystal : nullptr;
  {
    const QSignalBlocker block(kind_);
    kind_->setCurrentIndex(cr != nullptr ? 1 : 0);
  }
  {
    const QSignalBlocker block(material_);
    material_->setCurrentIndex(0);
  }
  for (QWidget *w : needs_layer_) w->setEnabled(l != nullptr);
  for (QWidget *w : needs_crystal_) w->setEnabled(cr != nullptr);
  bool any = false;
  for (const auto &layer : c.trim.layers) any = any || layer.crystal.has_value();
  const bool following = c.trim.transport == trim::TrimConfig::Transport::crystal;
  follow_->setChecked(following); // (clicked() is not emitted for this)
  transport_note_->setVisible(any && !following);
  // An amorphous layer whose composition is a crystal of the database.
  {
    const crystal::DatabaseCrystal *match = nullptr;
    if (l != nullptr && cr == nullptr && !l->gas) {
      std::vector<std::pair<unsigned, double>> composition;
      for (const auto &e : l->elements) composition.emplace_back(e.z, e.stoichiometry);
      match = crystal::match_composition(composition);
    }
    match_->setVisible(match != nullptr);
    if (match != nullptr) {
      match_->setText(tr("This layer's composition is the crystal %1 (%2, a = %3 A) of the database. "
                         "<a href=\"use\">Give it this crystal</a>%4")
                          .arg(view_of(match->name), view_of(match->structure), format_number(match->a_A, 5),
                               c.trim.crystal.auto_crystal && following
                                   ? tr(" (with automatic crystals on, the calculation does so by itself).")
                                   : QStringLiteral(".")));
    }
  }

  if (cr == nullptr) {
    summary_->setText(l == nullptr ? tr("The target has no layer.")
                                   : tr("This layer is amorphous. Choose \"Single crystal\", a common crystal, a "
                                        "structure or a cell file to give it a lattice."));
    view_->clear(tr("Amorphous layer"));
  } else {
    QString composition;
    for (const auto &[z, atoms] : cr->cell.composition()) {
      composition += (composition.isEmpty() ? QString() : QStringLiteral(", ")) +
                     QStringLiteral("%1 %2").arg(sym(z), format_number(atoms, 6));
    }
    QString text = tr("%1: %2 sites, atoms per cell: %3; density %4")
                       .arg(cr->cell.name.empty() ? tr("Unit cell") : QString::fromStdString(cr->cell.name))
                       .arg(cr->cell.sites.size())
                       .arg(composition, QStringLiteral("%1 g/cm3").arg(format_number(
                                             in_units(l->density, common::gram_per_cubic_centimeter), 5)));
    try {
      if (!cr->cell.sites.empty() && cr->cell.sites.size() <= 64) {
        const double reach = 1.01 * std::min({in_units(cr->cell.lattice.a(), common::angstrom),
                                              in_units(cr->cell.lattice.b(), common::angstrom),
                                              in_units(cr->cell.lattice.c(), common::angstrom)});
        const auto shells = cr->cell.neighbour_shells(cr->cell.sites.front().position, reach);
        if (!shells.empty()) {
          text += tr("; nearest neighbours of %1: %2 at %3 A")
                      .arg(QString::fromStdString(cr->cell.sites.front().label), format_number(shells.front().count, 4),
                           format_number(shells.front().distance_A, 5));
        }
        text += tr("; largest gap between the atoms: radius %1 A").arg(format_number(crystal::largest_hole_radius_A(cr->cell), 4));
      }
    } catch (const std::exception &) {
      // A cell that is being edited may not allow these; the rest is still shown.
    }
    summary_->setText(text);
    view_->set_crystal(*cr);
  }
  refresh_sites(cr);
  refresh_vibrations(c, l);
  refresh_directions(c, cr);
}

void CrystalPanel::refresh_sites(const crystal::Crystal *crystal) {
  const QSignalBlocker block(sites_);
  const int n = crystal != nullptr ? static_cast<int>(crystal->cell.sites.size()) : 0;
  sites_->setRowCount(n);
  for (int i = 0; i < n; ++i) {
    const auto &s = crystal->cell.sites[static_cast<std::size_t>(i)];
    set_item(sites_, i, site_label, QString::fromStdString(s.label));
    set_item(sites_, i, site_x, format_number(s.position.x, 6));
    set_item(sites_, i, site_y, format_number(s.position.y, 6));
    set_item(sites_, i, site_z, format_number(s.position.z, 6));
    set_item(sites_, i, site_kind, QString::fromStdString(crystal::to_string(s.kind)));
    set_item(sites_, i, site_occupants, occupants_text(s.occupants));
  }
}

void CrystalPanel::refresh_vibrations(const config::RunConfig &c, const trim::TargetLayer *layer) {
  const QSignalBlocker block(vibrations_);
  // As tall as its rows (one per element).
  const auto fit = [this](int rows) {
    const int row_height = vibrations_->verticalHeader()->defaultSectionSize();
    vibrations_->setFixedHeight(vibrations_->horizontalHeader()->sizeHint().height() + std::max(rows, 1) * row_height + 6);
  };
  if (layer == nullptr || !layer->crystal) {
    vibrations_->setRowCount(0);
    fit(0);
    return;
  }
  const auto &cr = *layer->crystal;
  const auto composition = cr.cell.composition();
  vibrations_->setRowCount(static_cast<int>(composition.size()));
  fit(static_cast<int>(composition.size()));
  int row = 0;
  for (const auto &[z, atoms] : composition) {
    (void)atoms;
    const auto given = std::find_if(cr.vibrations.begin(), cr.vibrations.end(),
                                    [z = z](const crystal::Vibration &v) { return v.z == z; });
    auto *item = set_item(vibrations_, row, vib_element, sym(z), false);
    item->setData(Qt::UserRole, z);
    set_item(vibrations_, row, vib_debye,
             given != cr.vibrations.end() && given->debye_temperature
                 ? QString::fromStdString(config::format_temperature(*given->debye_temperature))
                 : QString());
    set_item(vibrations_, row, vib_rms,
             given != cr.vibrations.end() && given->rms ? angstrom_text(*given->rms) : QString());
    QString used;
    QString why;
    try {
      const auto element = std::find_if(layer->elements.begin(), layer->elements.end(),
                                        [z = z](const trim::TargetElement &e) { return e.z == z; });
      const auto mass = element != layer->elements.end() ? element->mass : trim::TargetElement::defaults(z).mass;
      bool from_table = false;
      const auto u = crystal::thermal_amplitude(cr, z, mass, &from_table);
      if (!c.trim.crystal.thermal_vibrations) {
        used = tr("off");
        why = tr("Thermal vibrations are switched off in the settings of the crystal transport below.");
      } else if (!u) {
        used = tr("static");
        why = tr("No Debye temperature or amplitude is given and none is tabulated for this crystal: the atoms "
                 "stand exactly on their sites, which makes channeling stronger than in a real crystal.");
      } else if (from_table) {
        used = tr("%1 A (table: %2 K)")
                   .arg(format_number(in_units(*u, common::angstrom), 3),
                        format_number(in_units(*crystal::reference_debye_temperature(z), common::kelvin), 4));
        why = tr("The rms displacement along one axis, from the Debye temperature of Gemmell's table (Rev. Mod. "
                 "Phys. 46, 129 (1974)), the crystal's temperature and the atom's mass.");
      } else {
        used = tr("%1 A").arg(format_number(in_units(*u, common::angstrom), 3));
        why = tr("The rms displacement along one axis that the crystal transport uses.");
      }
    } catch (const std::exception &e) {
      used = tr("?");
      why = QString::fromUtf8(e.what());
    }
    set_item(vibrations_, row, vib_used, used, false)->setToolTip(why);
    ++row;
  }
}

void CrystalPanel::refresh_directions(const config::RunConfig &c, const crystal::Crystal *crystal) {
  beam_text_.clear();
  if (crystal == nullptr) {
    directions_->clear();
    beam_->clear();
    view_->set_directions(std::nullopt, std::nullopt);
    return;
  }
  try {
    const auto &lattice = crystal->cell.lattice;
    const Mat3 r = crystal->rotation(); // rows: the crystal directions along depth, Y and Z
    const Vec3 beam = crystal::beam_in_crystal(lattice, crystal->orientation, c.trim.ion.angle);
    beam_text_ = direction_text(lattice, beam);
    directions_->setText(tr("Depth (surface normal) along %1; in the surface, Y along %2 and Z along %3. The beam "
                            "(angle of incidence %4, Ion tab; twist %5) travels along %6.")
                             .arg(direction_text(lattice, r[0]), direction_text(lattice, r[1]), direction_text(lattice, r[2]),
                                  degree_text(c.trim.ion.angle), degree_text(crystal->orientation.twist), beam_text_));
    view_->set_directions(beam, r[0]);
  } catch (const std::exception &e) {
    directions_->setText(tr("The orientation cannot be used: %1").arg(QString::fromUtf8(e.what())));
    view_->set_directions(std::nullopt, std::nullopt);
  }
  if (!beam_->hasFocus()) beam_->setText(beam_text_);
}

// --- actions ----------------------------------------------------------------------

void CrystalPanel::set_kind(int index) {
  const auto *l = current_layer();
  if (l == nullptr) return;
  if (index == 0) {
    if (!l->crystal) return;
    const QString err = scenario_->modify(tr("Amorphous layer"), [&](config::RunConfig &c) { layer_in(c).crystal.reset(); });
    if (!err.isEmpty()) emit message(err);
    return;
  }
  if (l->crystal) return;
  // A pure element of the table of common crystals becomes that crystal;
  // anything else is asked for.
  if (l->elements.size() == 1) {
    for (const auto &r : crystal::reference_crystals()) {
      if (r.sublattices.size() == 1 && r.sublattices.front() == l->elements.front().z) {
        const auto cr = crystal::make_reference_crystal(r);
        set_cell(tr("Single crystal"), cr.cell, &cr.vibrations, {}, std::nullopt);
        return;
      }
    }
  }
  build_structure();
  binder_->refresh(); // back to "Amorphous" if nothing was built
}

void CrystalPanel::choose_from_database(const QString &query) {
  const auto *l = current_layer();
  if (l == nullptr) return;
  QDialog d(this);
  d.setWindowTitle(tr("Crystal database"));
  d.resize(860, 520);
  auto *search = new QLineEdit(query, &d);
  search->setPlaceholderText(tr("Search: a name, elements, a structure or a word of the note (\"GaAs\", \"Ga As\", "
                                "\"iron\", \"nitride\", \"wurtzite\", \"hcp\")"));
  search->setClearButtonEnabled(true);
  auto *table = new QTableWidget(&d);
  const QStringList headers{tr("Name"), tr("Structure"), tr("a (A)"), tr("c (A)"), tr("Debye (K)"), tr("Source"),
                            tr("Note")};
  table->setColumnCount(static_cast<int>(headers.size()));
  table->setHorizontalHeaderLabels(headers);
  table->setSelectionBehavior(QAbstractItemView::SelectRows);
  table->setSelectionMode(QAbstractItemView::SingleSelection);
  table->setEditTriggers(QAbstractItemView::NoEditTriggers);
  table->verticalHeader()->setVisible(false);
  table->horizontalHeader()->setStretchLastSection(true);
  auto *citation = new QLabel(&d);
  citation->setWordWrap(true);
  citation->setTextInteractionFlags(Qt::TextSelectableByMouse);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &d);
  buttons->button(QDialogButtonBox::Ok)->setText(tr("Give the layer this crystal"));
  auto *lay = new QVBoxLayout(&d);
  lay->addWidget(search);
  lay->addWidget(table, 1);
  lay->addWidget(citation);
  lay->addWidget(buttons);
  std::vector<const crystal::DatabaseCrystal *> shown;
  const auto fill = [&] {
    shown = crystal::search_crystal_database(search->text().toStdString());
    table->setRowCount(static_cast<int>(shown.size()));
    for (int i = 0; i < static_cast<int>(shown.size()); ++i) {
      const auto &r = *shown[static_cast<std::size_t>(i)];
      QString source = view_of(r.lattice_source);
      if (!r.debye_source.empty() && r.debye_source != r.lattice_source) source += QStringLiteral(" / ") + view_of(r.debye_source);
      QString note = view_of(r.note);
      if (!r.preferred) note += (note.isEmpty() ? QString() : QStringLiteral("; ")) + tr("not chosen for its composition");
      const QStringList cells{view_of(r.name), view_of(r.structure), format_number(r.a_A, 6),
                              r.c_A > 0.0 ? format_number(r.c_A, 6) : QString(),
                              r.debye_temperature_K > 0.0 ? QString::number(r.debye_temperature_K) : QString(), source,
                              note};
      for (int k = 0; k < static_cast<int>(cells.size()); ++k) table->setItem(i, k, new QTableWidgetItem(cells[k]));
    }
    table->resizeColumnsToContents();
    if (!shown.empty()) table->selectRow(0);
    buttons->button(QDialogButtonBox::Ok)->setEnabled(!shown.empty());
  };
  const auto cite = [&] {
    const int i = table->currentRow();
    if (i < 0 || i >= static_cast<int>(shown.size())) {
      citation->setText(tr("%1 crystals match.").arg(shown.size()));
      return;
    }
    const auto &r = *shown[static_cast<std::size_t>(i)];
    QString text = tr("Lattice: %1").arg(view_of(crystal::database_citation(r.lattice_source)));
    if (!r.debye_source.empty()) text += QStringLiteral("\n") + tr("Debye temperature: %1").arg(view_of(crystal::database_citation(r.debye_source)));
    else text += QStringLiteral("\n") + tr("No Debye temperature: its atoms stand still unless you give one (Thermal vibrations).");
    citation->setText(text);
  };
  QObject::connect(search, &QLineEdit::textChanged, &d, [&] {
    fill();
    cite();
  });
  QObject::connect(table, &QTableWidget::currentCellChanged, &d, [&] { cite(); });
  QObject::connect(table, &QTableWidget::cellDoubleClicked, &d, [&] { d.accept(); });
  QObject::connect(buttons, &QDialogButtonBox::accepted, &d, &QDialog::accept);
  QObject::connect(buttons, &QDialogButtonBox::rejected, &d, &QDialog::reject);
  fill();
  cite();
  if (d.exec() != QDialog::Accepted) return;
  const int i = table->currentRow();
  if (i < 0 || i >= static_cast<int>(shown.size())) return;
  const auto &row = *shown[static_cast<std::size_t>(i)];
  const auto made = crystal::make_database_crystal(row);
  if (set_cell(tr("Crystal %1").arg(view_of(row.name)), made.cell, &made.vibrations, std::string(row.name),
               std::nullopt) &&
      !row.note.empty()) {
    emit message(tr("%1: %2").arg(view_of(row.name), view_of(row.note)));
  }
}

void CrystalPanel::choose_material(int index) {
  const auto &rows = crystal::reference_crystals();
  if (index <= 0 || index > static_cast<int>(rows.size())) return;
  const auto &row = rows[static_cast<std::size_t>(index - 1)];
  const auto cr = crystal::make_reference_crystal(row);
  if (set_cell(tr("Crystal %1").arg(view_of(row.name)), cr.cell, &cr.vibrations, std::string(row.name), std::nullopt) &&
      !row.note.empty()) {
    emit message(tr("%1: %2").arg(view_of(row.name), view_of(row.note)));
  }
}

void CrystalPanel::build_structure() {
  const auto *l = current_layer();
  if (l == nullptr) return;
  const auto choice = ask_structure(this, l);
  if (!choice) return;
  set_cell(tr("Crystal structure"), crystal::build_structure(choice->spec), nullptr, {}, choice->density);
}

void CrystalPanel::import_cell() {
  if (current_layer() == nullptr) return;
  const QString path = QFileDialog::getOpenFileName(
      this, tr("Import a unit cell"), {},
      tr("Crystal structures (*.cif *.vasp *.poscar POSCAR* CONTCAR* *.data *.lmp *.lammps *.xyz *.extxyz);;All files (*)"));
  if (path.isEmpty()) return;
  const std::filesystem::path file(path.toStdWString());
  auto format = crystal::format_from_path(file);
  if (!format) {
    bool ok = false;
    const QStringList formats{QStringLiteral("cif"), QStringLiteral("poscar"), QStringLiteral("lammps"), QStringLiteral("xyz")};
    const QString chosen = QInputDialog::getItem(this, tr("Import a unit cell"),
                                                 tr("The file name does not tell the format of %1:").arg(QFileInfo(path).fileName()),
                                                 formats, 0, false, &ok);
    if (!ok) return;
    format = crystal::parse_cell_format(chosen.toStdString());
  }
  crystal::CellReadOptions options;
  crystal::UnitCell cell;
  try {
    cell = crystal::read_cell(file, format, options);
  } catch (const std::invalid_argument &e) {
    // A VASP 4 POSCAR has no species line and LAMMPS types may not say their
    // elements: ask for them once.
    if (*format != crystal::CellFormat::poscar && *format != crystal::CellFormat::lammps_data) throw;
    bool ok = false;
    const QString elements = QInputDialog::getText(
        this, tr("Import a unit cell"),
        tr("%1\n\nElements in the order of the file (species of a POSCAR, atom types 1, 2, ... of LAMMPS), "
           "e.g. Ga As:")
            .arg(QString::fromUtf8(e.what())),
        QLineEdit::Normal, QString(), &ok);
    if (!ok || elements.trimmed().isEmpty()) return;
    QString list = elements;
    list.replace(QLatin1Char(','), QLatin1Char(' '));
    int type = 1;
    for (const QString &s : list.split(QLatin1Char(' '), Qt::SkipEmptyParts)) {
      const unsigned z = parse_element(s);
      options.elements.push_back(z);
      options.lammps_types[type++] = z;
    }
    cell = crystal::read_cell(file, format, options);
  }
  if (cell.name.empty()) cell.name = QFileInfo(path).completeBaseName().toStdString();
  if (set_cell(tr("Import a unit cell"), cell, nullptr, {}, std::nullopt)) {
    emit message(tr("%1: %2 sites read").arg(QFileInfo(path).fileName()).arg(cell.sites.size()));
  }
}

void CrystalPanel::export_cell() {
  const auto *cr = current_crystal();
  if (cr == nullptr) return;
  const QStringList filters{tr("CIF (*.cif)"), tr("VASP POSCAR (*.vasp)"), tr("LAMMPS data (*.data)"), tr("Extended XYZ (*.xyz)")};
  QString chosen = filters.front();
  const QString path = QFileDialog::getSaveFileName(this, tr("Export the unit cell"), {}, filters.join(QStringLiteral(";;")), &chosen);
  if (path.isEmpty()) return;
  const crystal::CellFormat format = chosen == filters[1]   ? crystal::CellFormat::poscar
                                     : chosen == filters[2] ? crystal::CellFormat::lammps_data
                                     : chosen == filters[3] ? crystal::CellFormat::extxyz
                                                            : crystal::CellFormat::cif;
  crystal::save_cell(cr->cell, std::filesystem::path(path.toStdWString()), format);
  emit message(tr("Unit cell written to %1").arg(path));
}

void CrystalPanel::add_site() {
  const auto *cr = current_crystal();
  if (cr == nullptr) return;
  QDialog d(this);
  d.setWindowTitle(tr("Add a site"));
  d.setMinimumWidth(400);
  auto *position = new QLineEdit(&d);
  position->setPlaceholderText(tr("fractional x y z, e.g. 1/2 1/2 0"));
  auto *element = new QLineEdit(&d);
  element->setPlaceholderText(tr("e.g. C, or a mixture: Fe 0.9, Cr 0.1"));
  auto *pick = new QToolButton(&d);
  pick->setText(tr("PT"));
  auto *element_row = new QHBoxLayout();
  element_row->addWidget(element, 1);
  element_row->addWidget(pick);
  auto *kind = new QComboBox(&d);
  kind->addItems({tr("Lattice site"), tr("Interstitial site")});
  auto *form = new QFormLayout();
  form->addRow(tr("Position"), position);
  form->addRow(tr("Occupants"), element_row);
  form->addRow(tr("Kind"), kind);
  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &d);
  auto *lay = new QVBoxLayout(&d);
  lay->addLayout(form);
  lay->addWidget(buttons);
  connect(pick, &QToolButton::clicked, &d, guarded(&d, tr("Choose an element"), [&d, element] {
    const unsigned z = PeriodicTableDialog::pick(&d, 0);
    if (z != 0) element->setText(sym(z));
  }));
  connect(buttons, &QDialogButtonBox::accepted, &d, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, &d, &QDialog::reject);
  if (d.exec() != QDialog::Accepted) return;
  int added = -1;
  if (edit(tr("Add a site"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
        crystal::Site s;
        s.position = parse_position(position->text());
        s.occupants = parse_occupants(element->text());
        s.kind = kind->currentIndex() == 1 ? crystal::SiteKind::interstitial : crystal::SiteKind::lattice;
        // Labelled like the built-in structures: element and a running number.
        const std::string symbol = s.occupants.empty() ? std::string("X") : std::string(Element::by_z(s.occupants.front().z).symbol());
        int number = 1;
        const auto taken = [&](const std::string &label) {
          return std::any_of(c.cell.sites.begin(), c.cell.sites.end(), [&](const crystal::Site &o) { return o.label == label; });
        };
        while (taken(symbol + std::to_string(number))) ++number;
        s.label = symbol + std::to_string(number);
        c.cell.sites.push_back(std::move(s));
        added = static_cast<int>(c.cell.sites.size()) - 1;
      })) {
    select_site(added);
  }
}

void CrystalPanel::remove_site() {
  const int row = sites_->currentRow();
  if (row < 0 || current_crystal() == nullptr) return;
  if (edit(tr("Remove a site"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
        if (row >= static_cast<int>(c.cell.sites.size())) throw std::out_of_range("no such site");
        c.cell.sites.erase(c.cell.sites.begin() + row);
      })) {
    select_site(-1);
  }
}

void CrystalPanel::site_edited(int row, int column) {
  const auto *item = sites_->item(row, column);
  if (item == nullptr) return;
  const QString text = item->text();
  edit(tr("Edit a site"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
    auto &s = c.cell.sites.at(static_cast<std::size_t>(row));
    switch (column) {
    case site_label: s.label = text.trimmed().toStdString(); break;
    case site_x: s.position.x = crystal::wrap_fraction(parse_coordinate(text)); break;
    case site_y: s.position.y = crystal::wrap_fraction(parse_coordinate(text)); break;
    case site_z: s.position.z = crystal::wrap_fraction(parse_coordinate(text)); break;
    case site_kind: s.kind = crystal::parse_site_kind(text.trimmed().toLower().toStdString()); break;
    case site_occupants: s.occupants = parse_occupants(text); break;
    default: break;
    }
  });
  // (A refused value: edit() has put the scenario's back.)
}

void CrystalPanel::vibration_edited(int row, int column) {
  const auto *element = vibrations_->item(row, vib_element);
  const auto *item = vibrations_->item(row, column);
  if (element == nullptr || item == nullptr || (column != vib_debye && column != vib_rms)) return;
  const unsigned z = element->data(Qt::UserRole).toUInt();
  const QString text = item->text().trimmed();
  edit(tr("Thermal vibrations of %1").arg(sym(z)), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
    auto it = std::find_if(c.vibrations.begin(), c.vibrations.end(), [z](const crystal::Vibration &v) { return v.z == z; });
    crystal::Vibration v = it != c.vibrations.end() ? *it : crystal::Vibration{z, std::nullopt, std::nullopt};
    if (column == vib_debye) {
      v.debye_temperature = text.isEmpty() ? std::nullopt
                                           : std::optional(config::parse_temperature(text.toStdString(), "Debye temperature"));
    } else {
      v.rms = text.isEmpty() ? std::nullopt : std::optional(parse_length(text, "A"));
    }
    const bool empty = !v.debye_temperature && !v.rms;
    if (it != c.vibrations.end()) {
      if (empty) c.vibrations.erase(it);
      else *it = v;
    } else if (!empty) {
      c.vibrations.push_back(v);
    }
  });
}

void CrystalPanel::select_site(int site) {
  view_->set_selected_site(site);
  const QSignalBlocker block(sites_);
  if (site >= 0 && site < sites_->rowCount()) {
    sites_->selectRow(site);
    sites_->scrollToItem(sites_->item(site, site_label));
  } else {
    sites_->clearSelection();
    sites_->setCurrentCell(-1, -1);
  }
}

void CrystalPanel::site_menu(int site, const QPoint &global_position) {
  const auto *cr = current_crystal();
  if (cr == nullptr || site < 0 || site >= static_cast<int>(cr->cell.sites.size())) return;
  select_site(site);
  const crystal::Site s = cr->cell.sites[static_cast<std::size_t>(site)];
  // The sites that look the same: same kind, same occupants.
  std::vector<std::size_t> alike;
  for (std::size_t i = 0; i < cr->cell.sites.size(); ++i) {
    if (cr->cell.sites[i].kind == s.kind && cr->cell.sites[i].occupants == s.occupants) alike.push_back(i);
  }
  enum Action { act_replace, act_share, act_occupancy, act_lattice, act_interstitial, act_remove, act_all = 16 };
  QMenu menu(this);
  auto *title = menu.addAction(tr("%1 at (%2, %3, %4): %5")
                                   .arg(QString::fromStdString(s.label), format_number(s.position.x, 4),
                                        format_number(s.position.y, 4), format_number(s.position.z, 4),
                                        s.occupants.empty() ? tr("empty") : occupants_text(s.occupants)));
  title->setEnabled(false);
  menu.addSeparator();
  const auto add = [&](QMenu *m, int flag) {
    m->addAction(tr("Replace the element..."))->setData(act_replace | flag);
    m->addAction(tr("Share the site with another element..."))->setData(act_share | flag);
    m->addAction(tr("Set the occupancy..."))->setData(act_occupancy | flag);
  };
  add(&menu, 0);
  if (alike.size() > 1) add(menu.addMenu(tr("All %1 sites like this one").arg(alike.size())), act_all);
  menu.addSeparator();
  auto *lattice_action = menu.addAction(tr("Lattice site"));
  lattice_action->setCheckable(true);
  lattice_action->setChecked(s.kind == crystal::SiteKind::lattice);
  lattice_action->setData(act_lattice);
  auto *interstitial_action = menu.addAction(tr("Interstitial site"));
  interstitial_action->setCheckable(true);
  interstitial_action->setChecked(s.kind == crystal::SiteKind::interstitial);
  interstitial_action->setData(act_interstitial);
  menu.addSeparator();
  menu.addAction(tr("Remove this site"))->setData(act_remove);
  const QAction *chosen = menu.exec(global_position);
  if (chosen == nullptr || !chosen->data().isValid()) return;
  const int code = chosen->data().toInt();
  const std::vector<std::size_t> targets = (code & act_all) != 0 ? alike : std::vector<std::size_t>{static_cast<std::size_t>(site)};
  const double held = s.occupancy();
  const int action = code & ~static_cast<int>(act_all);
  switch (action) {
  case act_replace: {
    const unsigned z = PeriodicTableDialog::pick(this, s.occupants.empty() ? 0 : s.occupants.front().z);
    if (z == 0) return;
    edit(tr("Replace the element"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
      for (const std::size_t i : targets) {
        auto &t = c.cell.sites.at(i);
        const double occupancy_now = t.occupancy();
        if (t.occupants.size() == 1) relabel(t, t.occupants.front().z, z);
        t.occupants = {crystal::Occupant{z, occupancy_now > 0.0 ? std::min(occupancy_now, 1.0) : 1.0}};
      }
    });
    break;
  }
  case act_share: {
    const unsigned z = PeriodicTableDialog::pick(this, 0);
    if (z == 0) return;
    bool ok = false;
    const double percent = QInputDialog::getDouble(this, tr("Share the site"),
                                                   tr("Part of the site that %1 takes (%):").arg(sym(z)), 10.0, 0.0001,
                                                   100.0, 4, &ok);
    if (!ok) return;
    edit(tr("Share the site with %1").arg(sym(z)), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
      const double f = 0.01 * percent;
      for (const std::size_t i : targets) {
        auto &t = c.cell.sites.at(i);
        const double occupancy_now = t.occupancy();
        // The new element takes its part from those already there; an
        // empty or partly vacant site is filled first.
        const double free_part = std::max(0.0, 1.0 - occupancy_now);
        const double from_others = std::max(0.0, f - free_part);
        if (from_others > 0.0 && occupancy_now > 0.0) {
          for (auto &o : t.occupants) o.fraction *= std::max(0.0, occupancy_now - from_others) / occupancy_now;
        }
        t.occupants.erase(std::remove_if(t.occupants.begin(), t.occupants.end(),
                                         [](const crystal::Occupant &o) { return !(o.fraction > 1e-12); }),
                          t.occupants.end());
        const auto same = std::find_if(t.occupants.begin(), t.occupants.end(), [z](const crystal::Occupant &o) { return o.z == z; });
        if (same != t.occupants.end()) same->fraction = std::min(1.0, same->fraction + f);
        else t.occupants.push_back({z, f});
      }
    });
    break;
  }
  case act_occupancy: {
    if (s.occupants.empty()) {
      emit message(tr("The site is empty: replace its element first."));
      return;
    }
    bool ok = false;
    const double percent = QInputDialog::getDouble(this, tr("Occupancy"),
                                                   tr("Probability that the site holds an atom (%):"), 100.0 * held, 0.0,
                                                   100.0, 4, &ok);
    if (!ok) return;
    edit(tr("Set the occupancy"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
      for (const std::size_t i : targets) {
        auto &t = c.cell.sites.at(i);
        const double occupancy_now = t.occupancy();
        if (!(occupancy_now > 0.0)) continue;
        if (percent <= 0.0) {
          t.occupants.clear();
        } else {
          for (auto &o : t.occupants) o.fraction *= 0.01 * percent / occupancy_now;
        }
      }
    });
    break;
  }
  case act_lattice:
  case act_interstitial:
    if ((action == act_interstitial) == (s.kind == crystal::SiteKind::interstitial)) return; // as it is
    edit(tr("Kind of site"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
      c.cell.sites.at(static_cast<std::size_t>(site)).kind =
          action == act_interstitial ? crystal::SiteKind::interstitial : crystal::SiteKind::lattice;
    });
    break;
  case act_remove:
    if (edit(tr("Remove a site"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
          c.cell.sites.erase(c.cell.sites.begin() + site);
        })) {
      select_site(-1);
    }
    break;
  default: break;
  }
}

void CrystalPanel::beam_from_text() {
  const auto *cr = current_crystal();
  if (cr == nullptr) return;
  const QString text = beam_->text();
  if (text.trimmed() == beam_text_ || text.trimmed().isEmpty()) {
    beam_->setText(beam_text_);
    return;
  }
  Vec3 direction;
  const bool ok = edit(tr("Beam direction"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &config) {
    const auto &lattice = c.cell.lattice;
    direction = crystal::direction_vector(lattice, crystal::direction_indices(lattice, parse_indices(text), "beam direction"));
    const auto tt = crystal::tilt_twist_for_beam(lattice, c.orientation, direction, "beam direction");
    c.orientation.twist = tt.twist;
    config.trim.ion.angle = tt.tilt;
  });
  beam_->setText(beam_text_);
  if (ok) view_->look_along(direction);
}

void CrystalPanel::beam_along_view() {
  if (current_crystal() == nullptr) return;
  const Vec3 direction = view_->view_direction();
  edit(tr("Beam along the view"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &config) {
    const auto tt = crystal::tilt_twist_for_beam(c.cell.lattice, c.orientation, direction, "view direction");
    c.orientation.twist = tt.twist;
    config.trim.ion.angle = tt.tilt;
  });
  beam_->setText(beam_text_);
}

void CrystalPanel::orient_to_view() {
  if (current_crystal() == nullptr) return;
  // The view's frame: into the screen is the depth axis, its top lateral Y.
  const Vec3 normal = view_->view_direction();
  const Vec3 up = view_->view_up();
  edit(tr("Crystal orientation from the view"), [&](crystal::Crystal &c, trim::TargetLayer &, config::RunConfig &) {
    c.orientation = crystal::orientation_for_frame(c.cell.lattice, normal, up, 6, 0.05 * common::degree, "view direction");
  });
  beam_->setText(beam_text_);
}

void CrystalPanel::look_from_text() {
  const auto *cr = current_crystal();
  if (cr == nullptr) return;
  try {
    const auto &lattice = cr->cell.lattice;
    view_->look_along(crystal::direction_vector(lattice, crystal::direction_indices(lattice, parse_indices(look_->text()), "view direction")));
  } catch (const std::exception &e) {
    emit message(tr("View direction: %1").arg(QString::fromUtf8(e.what())));
  }
}

void CrystalPanel::set_crystal_transport(bool on) {
  const QString err = scenario_->set_crystal_transport(on);
  if (!err.isEmpty()) {
    QMessageBox::warning(this, tr("Crystal transport"), err);
    binder_->refresh(); // the box goes back to what the scenario says
  }
}

} // namespace gui
} // namespace amcset
