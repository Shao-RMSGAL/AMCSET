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
#include "setup/output_spec_dialog.h"

#include "views/tables_panel.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDialogButtonBox>
#include <QFile>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QMessageBox>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QTabWidget>
#include <QTableView>
#include <QTableWidget>
#include <QVBoxLayout>

namespace amcset {
namespace gui {

using output::Dimension;
using output::NumberFormat;

namespace {

const std::vector<std::pair<Dimension, const char *>> &dimension_names() {
  static const std::vector<std::pair<Dimension, const char *>> d{
      {Dimension::length, "Lengths"},
      {Dimension::energy, "Energies"},
      {Dimension::energy_per_length, "Energy per length"},
      {Dimension::inverse_length, "Per-depth densities"},
      {Dimension::inverse_energy, "Per-energy spectra"},
      {Dimension::dimensionless, "Pure numbers"},
      {Dimension::angle, "Angles"},
      {Dimension::time, "Times"},
  };
  return d;
}

void fill_notation(QComboBox *c, bool with_default) {
  if (with_default) c->addItem(QObject::tr("(default)"), -1);
  c->addItem(QObject::tr("General (%g)"), static_cast<int>(NumberFormat::Notation::general));
  c->addItem(QObject::tr("Scientific (1.234E+05)"), static_cast<int>(NumberFormat::Notation::scientific));
  c->addItem(QObject::tr("Engineering (123.4E+03)"), static_cast<int>(NumberFormat::Notation::engineering));
  c->addItem(QObject::tr("Fixed decimals"), static_cast<int>(NumberFormat::Notation::fixed));
  c->addItem(QObject::tr("Shortest exact"), static_cast<int>(NumberFormat::Notation::shortest));
}

} // namespace

OutputSpecDialog::OutputSpecDialog(const output::OutputSpec &spec, std::shared_ptr<const output::Dataset> sample,
                                   QWidget *parent)
    : QDialog(parent), spec_(spec), sample_(std::move(sample)) {
  setWindowTitle(tr("Output tables, units and number format"));
  resize(1100, 720);
  auto *tabs = new QTabWidget(this);

  // --- General ------------------------------------------------------------------
  auto *general = new QWidget(tabs);
  auto *gl = new QHBoxLayout(general);
  auto *ubox = new QGroupBox(tr("Units (all tables)"), general);
  auto *uf = new QFormLayout(ubox);
  for (const auto &[d, name] : dimension_names()) {
    auto *c = new QComboBox(ubox);
    c->addItem(tr("as TRIM (each column's own)"), QString());
    for (const auto *u : output::units_of(d)) c->addItem(QString::fromStdString(u->name), QString::fromStdString(u->name));
    uf->addRow(tr(name), c);
    units_[d] = c;
    connect(c, &QComboBox::currentIndexChanged, this, guarded(this, tr("Output tables dialog (currentIndexChanged)"), [this] { store_general(); }));
  }
  gl->addWidget(ubox);
  auto *nbox = new QGroupBox(tr("Numbers and files"), general);
  auto *nf = new QFormLayout(nbox);
  notation_ = new QComboBox(nbox);
  fill_notation(notation_, false);
  digits_ = new QSpinBox(nbox);
  digits_->setRange(0, 17);
  digits_->setToolTip(tr("Significant digits (general), decimals (fixed, scientific, engineering)"));
  uppercase_ = new QCheckBox(tr("Upper-case exponent (E)"), nbox);
  delimiter_ = new QComboBox(nbox);
  delimiter_->addItem(tr("Comma"), QStringLiteral(","));
  delimiter_->addItem(tr("Semicolon"), QStringLiteral(";"));
  delimiter_->addItem(tr("Tab"), QStringLiteral("\t"));
  comments_ = new QCheckBox(tr("Description and run data as # comment lines"), nbox);
  pretty_ = new QCheckBox(tr("Indented (readable) JSON"), nbox);
  pq_units_ = new QCheckBox(tr("Units in Parquet column names, e.g. \"depth [nm]\""), nbox);
  nf->addRow(tr("Notation"), notation_);
  nf->addRow(tr("Digits"), digits_);
  nf->addRow(QString(), uppercase_);
  nf->addRow(tr("CSV separator"), delimiter_);
  nf->addRow(QString(), comments_);
  nf->addRow(QString(), pretty_);
  nf->addRow(QString(), pq_units_);
  gl->addWidget(nbox);
  for (auto *w : {static_cast<QWidget *>(notation_), static_cast<QWidget *>(delimiter_)}) {
    connect(static_cast<QComboBox *>(w), &QComboBox::currentIndexChanged, this, guarded(this, tr("Output tables dialog (currentIndexChanged)"), [this] { store_general(); }));
  }
  connect(digits_, &QSpinBox::valueChanged, this, guarded(this, tr("Output tables dialog (valueChanged)"), [this] { store_general(); }));
  for (auto *c : {uppercase_, comments_, pretty_, pq_units_}) connect(c, &QCheckBox::toggled, this, guarded(this, tr("Output tables dialog (toggled)"), [this] { store_general(); }));
  tabs->addTab(general, tr("General"));

  // --- Tables ----------------------------------------------------------------------
  auto *tw = new QWidget(tabs);
  auto *tl = new QHBoxLayout(tw);
  tables_ = new QListWidget(tw);
  tables_->setMaximumWidth(230);
  tl->addWidget(tables_);
  auto *right = new QVBoxLayout();
  include_ = new QCheckBox(tr("Write this table"), tw);
  right->addWidget(include_);
  auto *mid = new QHBoxLayout();
  auto *cbox = new QGroupBox(tr("Columns (check to write, drag to reorder)"), tw);
  auto *cl = new QVBoxLayout(cbox);
  columns_ = new QListWidget(cbox);
  columns_->setDragDropMode(QAbstractItemView::InternalMove);
  cl->addWidget(columns_);
  auto *cf = new QFormLayout();
  col_unit_ = new QComboBox(cbox);
  col_label_ = new QLineEdit(cbox);
  col_label_->setPlaceholderText(tr("(default heading)"));
  col_notation_ = new QComboBox(cbox);
  fill_notation(col_notation_, true);
  col_digits_ = new QSpinBox(cbox);
  col_digits_->setRange(0, 17);
  cf->addRow(tr("Unit"), col_unit_);
  cf->addRow(tr("Heading"), col_label_);
  cf->addRow(tr("Notation"), col_notation_);
  cf->addRow(tr("Digits"), col_digits_);
  cl->addLayout(cf);
  mid->addWidget(cbox, 1);
  auto *rbox = new QGroupBox(tr("Rows and file"), tw);
  auto *rl = new QVBoxLayout(rbox);
  rl->addWidget(new QLabel(tr("Keep rows within these ranges (values may have units, e.g. \"500 nm\"):"), rbox));
  filters_ = new QTableWidget(0, 3, rbox);
  filters_->setHorizontalHeaderLabels({tr("Column"), tr("Minimum"), tr("Maximum")});
  filters_->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
  rl->addWidget(filters_);
  auto *frow = new QHBoxLayout();
  auto *add_filter = new QPushButton(tr("Add range"), rbox);
  auto *remove_filter = new QPushButton(tr("Remove"), rbox);
  frow->addWidget(add_filter);
  frow->addWidget(remove_filter);
  frow->addStretch(1);
  rl->addLayout(frow);
  auto *rf = new QFormLayout();
  drop_zero_ = new QCheckBox(tr("Leave out rows that are all zero"), rbox);
  every_ = new QSpinBox(rbox);
  every_->setRange(1, 1000000);
  every_->setPrefix(tr("every "));
  max_rows_ = new QSpinBox(rbox);
  max_rows_->setRange(0, 2'000'000'000);
  max_rows_->setSpecialValueText(tr("all"));
  file_name_ = new QLineEdit(rbox);
  table_notation_ = new QComboBox(rbox);
  fill_notation(table_notation_, true);
  table_digits_ = new QSpinBox(rbox);
  table_digits_->setRange(0, 17);
  rf->addRow(QString(), drop_zero_);
  rf->addRow(tr("Rows"), every_);
  rf->addRow(tr("At most"), max_rows_);
  rf->addRow(tr("File name"), file_name_);
  rf->addRow(tr("Notation"), table_notation_);
  rf->addRow(tr("Digits"), table_digits_);
  rl->addLayout(rf);
  mid->addWidget(rbox, 1);
  right->addLayout(mid, 1);
  tl->addLayout(right, 1);
  tabs->addTab(tw, tr("Tables"));

  // --- Preview and buttons ------------------------------------------------------
  auto *split = new QSplitter(Qt::Vertical, this);
  split->addWidget(tabs);
  auto *pw = new QWidget(split);
  auto *pl = new QVBoxLayout(pw);
  pl->setContentsMargins(0, 0, 0, 0);
  status_ = new QLabel(pw);
  status_->setWordWrap(true);
  preview_ = new QTableView(pw);
  preview_model_ = new TableModel(this);
  preview_->setModel(preview_model_);
  preview_->verticalHeader()->setDefaultSectionSize(preview_->fontMetrics().height() + 4);
  pl->addWidget(status_);
  pl->addWidget(preview_);
  split->addWidget(pw);
  split->setSizes({430, 260});
  auto *bb = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel | QDialogButtonBox::Reset, this);
  auto *load = bb->addButton(tr("Load..."), QDialogButtonBox::ActionRole);
  auto *save = bb->addButton(tr("Save..."), QDialogButtonBox::ActionRole);
  auto *lay = new QVBoxLayout(this);
  lay->addWidget(split, 1);
  lay->addWidget(bb);

  connect(bb, &QDialogButtonBox::accepted, this, guarded(this, tr("Output tables dialog (accepted)"), [this] {
    if (sample_) {
      try {
        spec_.validate(*sample_);
      } catch (const std::exception &e) {
        QMessageBox::warning(this, tr("Output settings"), QString::fromUtf8(e.what()));
        return;
      }
    }
    accept();
  }));
  connect(bb, &QDialogButtonBox::rejected, this, &QDialog::reject);
  connect(bb->button(QDialogButtonBox::Reset), &QPushButton::clicked, this, guarded(this, tr("Output tables dialog (clicked)"), [this] {
    const auto formats = spec_.formats;
    spec_ = output::OutputSpec{};
    spec_.formats = formats;
    load_general();
    select_table(tables_->currentRow());
  }));
  connect(load, &QPushButton::clicked, this, guarded(this, tr("Output tables dialog (clicked)"), [this] {
    const QString f = QFileDialog::getOpenFileName(this, tr("Load output settings"), {}, tr("JSON (*.json)"));
    if (f.isEmpty()) return;
    QFile file(f);
    if (!file.open(QIODevice::ReadOnly)) return;
    try {
      spec_ = output::parse_output_spec(file.readAll().toStdString());
    } catch (const std::exception &e) {
      QMessageBox::warning(this, tr("Output settings"), QString::fromUtf8(e.what()));
      return;
    }
    load_general();
    select_table(tables_->currentRow());
  }));
  connect(save, &QPushButton::clicked, this, guarded(this, tr("Output tables dialog (clicked)"), [this] {
    const QString f = QFileDialog::getSaveFileName(this, tr("Save output settings"), QStringLiteral("output_spec.json"),
                                                   tr("JSON (*.json)"));
    if (f.isEmpty()) return;
    QFile file(f);
    if (file.open(QIODevice::WriteOnly)) file.write(QByteArray::fromStdString(output::to_json(spec_)));
  }));

  // Tables list.
  if (sample_) {
    for (const auto &t : sample_->tables) {
      auto *it = new QListWidgetItem(QString::fromStdString(t.title.empty() ? t.key : t.title), tables_);
      it->setData(Qt::UserRole, QString::fromStdString(t.key));
      it->setToolTip(QString::fromStdString(t.description));
    }
  }
  connect(tables_, &QListWidget::currentRowChanged, this, guarded(this, tr("Output tables dialog (currentRowChanged)"), [this](auto &&...a) { select_table(std::forward<decltype(a)>(a)...); }));
  connect(include_, &QCheckBox::toggled, this, guarded(this, tr("Output tables dialog (toggled)"), [this] { store_table(); }));
  connect(columns_, &QListWidget::itemChanged, this, guarded(this, tr("Output tables dialog (itemChanged)"), [this] { store_table(); }));
  connect(columns_->model(), &QAbstractItemModel::rowsMoved, this, guarded(this, tr("Output tables dialog (rowsMoved)"), [this] { store_table(); }));
  connect(columns_, &QListWidget::currentRowChanged, this, guarded(this, tr("Output tables dialog (currentRowChanged)"), [this](auto &&...a) { select_column(std::forward<decltype(a)>(a)...); }));
  connect(col_unit_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Output tables dialog (currentIndexChanged)"), [this] { store_column(); }));
  connect(col_label_, &QLineEdit::editingFinished, this, guarded(this, tr("Output tables dialog (editingFinished)"), [this] { store_column(); }));
  connect(col_notation_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Output tables dialog (currentIndexChanged)"), [this] { store_column(); }));
  connect(col_digits_, &QSpinBox::valueChanged, this, guarded(this, tr("Output tables dialog (valueChanged)"), [this] { store_column(); }));
  connect(filters_, &QTableWidget::itemChanged, this, guarded(this, tr("Output tables dialog (itemChanged)"), [this] { store_table(); }));
  connect(add_filter, &QPushButton::clicked, this, guarded(this, tr("Output tables dialog (clicked)"), [this] {
    const auto *t = current_table();
    if (!t || t->columns.empty()) return;
    loading_ = true;
    const int r = filters_->rowCount();
    filters_->insertRow(r);
    filters_->setItem(r, 0, new QTableWidgetItem(QString::fromStdString(t->columns.front().key)));
    filters_->setItem(r, 1, new QTableWidgetItem());
    filters_->setItem(r, 2, new QTableWidgetItem());
    loading_ = false;
    store_table();
  }));
  connect(remove_filter, &QPushButton::clicked, this, guarded(this, tr("Output tables dialog (clicked)"), [this] {
    const int r = filters_->currentRow();
    if (r >= 0) filters_->removeRow(r);
    store_table();
  }));
  connect(drop_zero_, &QCheckBox::toggled, this, guarded(this, tr("Output tables dialog (toggled)"), [this] { store_table(); }));
  connect(every_, &QSpinBox::valueChanged, this, guarded(this, tr("Output tables dialog (valueChanged)"), [this] { store_table(); }));
  connect(max_rows_, &QSpinBox::valueChanged, this, guarded(this, tr("Output tables dialog (valueChanged)"), [this] { store_table(); }));
  connect(file_name_, &QLineEdit::editingFinished, this, guarded(this, tr("Output tables dialog (editingFinished)"), [this] { store_table(); }));
  connect(table_notation_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Output tables dialog (currentIndexChanged)"), [this] { store_table(); }));
  connect(table_digits_, &QSpinBox::valueChanged, this, guarded(this, tr("Output tables dialog (valueChanged)"), [this] { store_table(); }));

  load_general();
  if (tables_->count() > 0) tables_->setCurrentRow(0);
  else status_->setText(tr("Run a calculation (or check the target) to list the tables."));
}

void OutputSpecDialog::load_general() {
  loading_ = true;
  for (auto &[d, c] : units_) {
    const auto it = spec_.units.find(d);
    const int k = it == spec_.units.end() ? 0 : c->findData(QString::fromStdString(it->second));
    c->setCurrentIndex(std::max(k, 0));
  }
  notation_->setCurrentIndex(std::max(0, notation_->findData(static_cast<int>(spec_.format.notation))));
  digits_->setValue(spec_.format.digits);
  uppercase_->setChecked(spec_.format.uppercase);
  delimiter_->setCurrentIndex(std::max(0, delimiter_->findData(QString(QChar(spec_.csv_delimiter)))));
  comments_->setChecked(spec_.csv_comments);
  pretty_->setChecked(spec_.json_pretty);
  pq_units_->setChecked(spec_.parquet_units_in_names);
  loading_ = false;
  update_preview();
}

void OutputSpecDialog::store_general() {
  if (loading_) return;
  spec_.units.clear();
  for (const auto &[d, c] : units_) {
    const QString u = c->currentData().toString();
    if (!u.isEmpty()) spec_.units[d] = u.toStdString();
  }
  spec_.format.notation = static_cast<NumberFormat::Notation>(notation_->currentData().toInt());
  spec_.format.digits = digits_->value();
  spec_.format.uppercase = uppercase_->isChecked();
  const QString del = delimiter_->currentData().toString();
  spec_.csv_delimiter = del.isEmpty() ? ',' : del[0].toLatin1();
  spec_.csv_comments = comments_->isChecked();
  spec_.json_pretty = pretty_->isChecked();
  spec_.parquet_units_in_names = pq_units_->isChecked();
  update_preview();
}

const output::Table *OutputSpecDialog::current_table() const {
  if (!sample_ || !tables_->currentItem()) return nullptr;
  return sample_->find(tables_->currentItem()->data(Qt::UserRole).toString().toStdString());
}

output::TableSpec &OutputSpecDialog::table_spec() {
  const auto *t = current_table();
  return spec_.tables[t ? t->key : std::string("*")];
}

void OutputSpecDialog::select_table(int) { load_table(); }

void OutputSpecDialog::load_table() {
  const auto *t = current_table();
  if (!t) return;
  loading_ = true;
  const output::TableSpec ts = spec_.table(t->key);
  include_->setChecked(ts.include);
  columns_->clear();
  std::vector<std::string> order = ts.columns;
  for (const auto &c : t->columns) {
    if (std::find(order.begin(), order.end(), c.key) == order.end()) order.push_back(c.key);
  }
  for (const auto &k : order) {
    const auto *c = t->find(k);
    if (!c) continue;
    auto *it = new QListWidgetItem(QStringLiteral("%1   [%2]").arg(QString::fromStdString(c->label), QString::fromStdString(c->key)),
                                   columns_);
    it->setData(Qt::UserRole, QString::fromStdString(c->key));
    it->setFlags(it->flags() | Qt::ItemIsUserCheckable | Qt::ItemIsDragEnabled);
    const bool listed = ts.columns.empty() || std::find(ts.columns.begin(), ts.columns.end(), k) != ts.columns.end();
    const bool excluded = std::find(ts.exclude.begin(), ts.exclude.end(), k) != ts.exclude.end();
    it->setCheckState(listed && !excluded ? Qt::Checked : Qt::Unchecked);
  }
  filters_->setRowCount(0);
  for (const auto &f : ts.filters) {
    const int r = filters_->rowCount();
    filters_->insertRow(r);
    const auto txt = [](const std::optional<double> &v, const std::string &u) {
      if (!v) return QString();
      return QString::number(*v, 'g', 10) + (u.empty() ? QString() : QStringLiteral(" ") + QString::fromStdString(u));
    };
    filters_->setItem(r, 0, new QTableWidgetItem(QString::fromStdString(f.column)));
    filters_->setItem(r, 1, new QTableWidgetItem(txt(f.min, f.min_unit)));
    filters_->setItem(r, 2, new QTableWidgetItem(txt(f.max, f.max_unit)));
  }
  drop_zero_->setChecked(ts.drop_zero_rows);
  every_->setValue(static_cast<int>(ts.every));
  max_rows_->setValue(static_cast<int>(std::min<std::size_t>(ts.max_rows, 2'000'000'000)));
  file_name_->setText(QString::fromStdString(ts.file_name));
  file_name_->setPlaceholderText(QString::fromStdString(t->key));
  table_notation_->setCurrentIndex(ts.format ? std::max(0, table_notation_->findData(static_cast<int>(ts.format->notation))) : 0);
  table_digits_->setValue(ts.format ? ts.format->digits : spec_.format.digits);
  loading_ = false;
  columns_->setCurrentRow(0);
  select_column(0);
  update_preview();
}

void OutputSpecDialog::store_table() {
  if (loading_) return;
  const auto *t = current_table();
  if (!t) return;
  output::TableSpec &ts = table_spec();
  ts.include = include_->isChecked();
  std::vector<std::string> checked;
  bool default_order = true;
  for (int i = 0; i < columns_->count(); ++i) {
    const auto key = columns_->item(i)->data(Qt::UserRole).toString().toStdString();
    if (i >= static_cast<int>(t->columns.size()) || t->columns[static_cast<std::size_t>(i)].key != key) default_order = false;
    if (columns_->item(i)->checkState() == Qt::Checked) checked.push_back(key);
  }
  ts.exclude.clear();
  if (default_order && checked.size() == t->columns.size()) ts.columns.clear();
  else ts.columns = checked;
  ts.filters.clear();
  QString problems;
  for (int r = 0; r < filters_->rowCount(); ++r) {
    output::RowFilter f;
    f.column = filters_->item(r, 0) ? filters_->item(r, 0)->text().trimmed().toStdString() : std::string();
    const auto parse = [&](int col, std::optional<double> &v, std::string &unit) {
      const QString s = filters_->item(r, col) ? filters_->item(r, col)->text().trimmed() : QString();
      if (s.isEmpty()) return;
      const int sp = s.indexOf(QLatin1Char(' '));
      bool ok = false;
      const double x = (sp < 0 ? s : s.left(sp)).toDouble(&ok);
      if (!ok) {
        problems += tr("'%1' is not a number. ").arg(s);
        return;
      }
      v = x;
      if (sp >= 0) unit = s.mid(sp + 1).trimmed().toStdString();
    };
    parse(1, f.min, f.min_unit);
    parse(2, f.max, f.max_unit);
    if (!f.column.empty() && (f.min || f.max)) ts.filters.push_back(std::move(f));
  }
  ts.drop_zero_rows = drop_zero_->isChecked();
  ts.every = static_cast<std::size_t>(every_->value());
  ts.max_rows = static_cast<std::size_t>(max_rows_->value());
  ts.file_name = file_name_->text().trimmed().toStdString();
  const int tn = table_notation_->currentData().toInt();
  if (tn < 0) ts.format.reset();
  else ts.format = NumberFormat{static_cast<NumberFormat::Notation>(tn), table_digits_->value(), spec_.format.uppercase};
  update_preview();
  if (!problems.isEmpty()) status_->setText(problems);
}

void OutputSpecDialog::select_column(int row) {
  const auto *t = current_table();
  if (!t || row < 0 || row >= columns_->count()) return;
  const auto key = columns_->item(row)->data(Qt::UserRole).toString().toStdString();
  const auto *c = t->find(key);
  if (!c) return;
  loading_ = true;
  const auto &ts = spec_.table(t->key);
  const auto cs = ts.column_specs.find(key);
  col_unit_->clear();
  col_unit_->addItem(tr("(default)"), QString());
  const bool real = c->type == output::ColumnType::real;
  if (real) {
    for (const auto *u : output::units_of(c->dimension)) col_unit_->addItem(QString::fromStdString(u->name), QString::fromStdString(u->name));
  }
  col_unit_->setEnabled(real);
  col_notation_->setEnabled(real);
  col_digits_->setEnabled(real);
  if (cs != ts.column_specs.end()) {
    col_unit_->setCurrentIndex(std::max(0, col_unit_->findData(QString::fromStdString(cs->second.unit))));
    col_label_->setText(QString::fromStdString(cs->second.label));
    col_notation_->setCurrentIndex(cs->second.format ? std::max(0, col_notation_->findData(static_cast<int>(cs->second.format->notation))) : 0);
    col_digits_->setValue(cs->second.format ? cs->second.format->digits : spec_.format.digits);
  } else {
    col_unit_->setCurrentIndex(0);
    col_label_->clear();
    col_notation_->setCurrentIndex(0);
    col_digits_->setValue(spec_.format.digits);
  }
  col_label_->setPlaceholderText(QString::fromStdString(c->label));
  loading_ = false;
}

void OutputSpecDialog::store_column() {
  if (loading_) return;
  const auto *t = current_table();
  const int row = columns_->currentRow();
  if (!t || row < 0) return;
  const auto key = columns_->item(row)->data(Qt::UserRole).toString().toStdString();
  auto &ts = table_spec();
  output::ColumnSpec cs;
  cs.unit = col_unit_->currentData().toString().toStdString();
  cs.label = col_label_->text().trimmed().toStdString();
  const int n = col_notation_->currentData().toInt();
  if (n >= 0) cs.format = NumberFormat{static_cast<NumberFormat::Notation>(n), col_digits_->value(), spec_.format.uppercase};
  if (cs.unit.empty() && cs.label.empty() && !cs.format) ts.column_specs.erase(key);
  else ts.column_specs[key] = cs;
  update_preview();
}

void OutputSpecDialog::update_preview() {
  const auto *t = current_table();
  if (!t) {
    preview_model_->set_table({}, spec_);
    return;
  }
  try {
    auto shaped = output::apply(*t, spec_);
    const std::size_t n = shaped.rows();
    constexpr std::size_t limit = 200;
    if (n > limit) {
      for (auto &c : shaped.columns) {
        if (c.type == output::ColumnType::text) c.text.resize(limit);
        else c.values.resize(limit);
      }
    }
    preview_model_->set_table(std::move(shaped), spec_);
    status_->setText(tr("Preview of %1: %L2 rows%3. %4")
                         .arg(QString::fromStdString(t->key))
                         .arg(n)
                         .arg(n > limit ? tr(" (first %1 shown)").arg(limit) : QString())
                         .arg(spec_.table(t->key).include ? QString() : tr("This table will not be written.")));
  } catch (const std::exception &e) {
    status_->setText(tr("Problem: %1").arg(QString::fromUtf8(e.what())));
  }
}

} // namespace gui
} // namespace amcset
