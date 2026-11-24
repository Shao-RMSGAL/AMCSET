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
#include "views/tables_panel.h"

#include "output/writers.h"

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QComboBox>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QLabel>
#include <QMessageBox>
#include <QPushButton>
#include <QTableView>
#include <QVBoxLayout>

#include <filesystem>
#include <fstream>

namespace amcset {
namespace gui {

void TableModel::set_table(output::Table table, output::OutputSpec spec) {
  beginResetModel();
  table_ = std::move(table);
  spec_ = std::move(spec);
  formats_.clear();
  for (const auto &c : table_.columns) formats_.push_back(output::column_format(table_, c, spec_));
  endResetModel();
}

int TableModel::rowCount(const QModelIndex &parent) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    return parent.isValid() ? 0 : static_cast<int>(table_.rows());
  } catch (...) {
    report_exception_quietly(QStringLiteral("Tables (rowCount)"), std::current_exception());
    return 0;
  }
}

int TableModel::columnCount(const QModelIndex &parent) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    return parent.isValid() ? 0 : static_cast<int>(table_.columns.size());
  } catch (...) {
    report_exception_quietly(QStringLiteral("Tables (columnCount)"), std::current_exception());
    return 0;
  }
}

QVariant TableModel::data(const QModelIndex &index, int role) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (!index.isValid()) return {};
    const auto &c = table_.columns[static_cast<std::size_t>(index.column())];
    const auto r = static_cast<std::size_t>(index.row());
    if (role == Qt::TextAlignmentRole) {
      return c.type == output::ColumnType::text ? QVariant(Qt::AlignLeft | Qt::AlignVCenter)
                                                : QVariant(Qt::AlignRight | Qt::AlignVCenter);
    }
    if (role != Qt::DisplayRole) return {};
    switch (c.type) {
    case output::ColumnType::text: return r < c.text.size() ? QString::fromStdString(c.text[r]) : QString();
    case output::ColumnType::integer: return r < c.values.size() ? QString::number(c.values[r], 'f', 0) : QString();
    case output::ColumnType::real:
      return r < c.values.size()
                 ? QString::fromStdString(output::format_number(c.values[r], formats_[static_cast<std::size_t>(index.column())]))
                 : QString();
    }
    return {};
  } catch (...) {
    report_exception_quietly(QStringLiteral("Tables (data)"), std::current_exception());
    return {};
  }
}

QVariant TableModel::headerData(int section, Qt::Orientation o, int role) const {
  // Called by Qt: no exception may leave (app/guard.h).
  try {
    if (o == Qt::Vertical) return role == Qt::DisplayRole ? QVariant(section + 1) : QVariant();
    if (section < 0 || section >= static_cast<int>(table_.columns.size())) return {};
    const auto &c = table_.columns[static_cast<std::size_t>(section)];
    if (role == Qt::ToolTipRole) {
      return QStringLiteral("%1\n%2").arg(QString::fromStdString(c.key), QString::fromStdString(c.description));
    }
    if (role != Qt::DisplayRole) return {};
    QString h = QString::fromStdString(c.label);
    if (c.type == output::ColumnType::real && !c.unit.empty() && c.unit != "1") {
      h += QStringLiteral("\n(%1)").arg(QString::fromStdString(c.unit));
    }
    return h;
  } catch (...) {
    report_exception_quietly(QStringLiteral("Tables (headerData)"), std::current_exception());
    return {};
  }
}

TablesPanel::TablesPanel(QWidget *parent) : QWidget(parent) {
  tables_ = new QComboBox(this);
  tables_->setMinimumContentsLength(24);
  tables_->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  shaped_ = new QCheckBox(tr("As in the output files"), this);
  shaped_->setToolTip(tr("Show the table as the output settings shape it: units, columns, rows and number format"));
  shaped_->setChecked(true);
  auto *copy = new QPushButton(tr("Copy"), this);
  auto *csv = new QPushButton(tr("Save CSV..."), this);
  auto *pq = new QPushButton(tr("Save Parquet..."), this);
  info_ = new QLabel(this);
  info_->setWordWrap(true);
  view_ = new QTableView(this);
  model_ = new TableModel(this);
  view_->setModel(model_);
  view_->horizontalHeader()->setSectionResizeMode(QHeaderView::Interactive);
  view_->horizontalHeader()->setDefaultSectionSize(110);
  view_->verticalHeader()->setDefaultSectionSize(view_->fontMetrics().height() + 4);
  view_->setAlternatingRowColors(true);
  auto *top = new QHBoxLayout();
  top->addWidget(new QLabel(tr("Table"), this));
  top->addWidget(tables_);
  top->addWidget(shaped_, 1);
  top->addWidget(copy);
  top->addWidget(csv);
  top->addWidget(pq);
  auto *lay = new QVBoxLayout(this);
  lay->addLayout(top);
  lay->addWidget(info_);
  lay->addWidget(view_, 1);
  connect(tables_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Tables (currentIndexChanged)"), [this] { refresh(); }));
  connect(shaped_, &QCheckBox::toggled, this, guarded(this, tr("Tables (toggled)"), [this] { refresh(); }));
  connect(copy, &QPushButton::clicked, this, guarded(this, tr("Tables (clicked)"), [this] {
    output::Dataset ctx;
    std::ostringstream o;
    output::OutputSpec s = shaped_->isChecked() ? spec_ : output::OutputSpec{};
    s.csv_delimiter = '\t';
    s.csv_comments = false;
    output::write_csv(o, model_->table(), ctx, s);
    QApplication::clipboard()->setText(QString::fromStdString(o.str()));
  }));
  connect(csv, &QPushButton::clicked, this, guarded(this, tr("Tables (clicked)"), [this] { export_table(output::Format::csv); }));
  connect(pq, &QPushButton::clicked, this, guarded(this, tr("Tables (clicked)"), [this] { export_table(output::Format::parquet); }));
}

void TablesPanel::set_output_spec(const output::OutputSpec &spec) {
  spec_ = spec;
  refresh();
}

void TablesPanel::set_dataset(std::shared_ptr<const output::Dataset> d) {
  data_ = std::move(d);
  const QString keep = tables_->currentData().toString();
  {
    QSignalBlocker b(tables_);
    tables_->clear();
    if (data_) {
      for (const auto &t : data_->tables) {
        tables_->addItem(QStringLiteral("%1  [%2]").arg(QString::fromStdString(t.title), QString::fromStdString(t.key)),
                         QString::fromStdString(t.key));
      }
    }
    const int k = tables_->findData(keep);
    tables_->setCurrentIndex(k >= 0 ? k : 0);
  }
  refresh();
}

void TablesPanel::refresh() {
  if (!data_) {
    model_->set_table({}, {});
    info_->clear();
    return;
  }
  const auto *t = data_->find(tables_->currentData().toString().toStdString());
  if (!t) {
    model_->set_table({}, {});
    return;
  }
  QString note = QString::fromStdString(t->description);
  try {
    if (shaped_->isChecked()) {
      const bool included = spec_.table(t->key).include;
      model_->set_table(output::apply(*t, spec_), spec_);
      if (!included) note = tr("(This table is excluded from the output files.) ") + note;
    } else {
      output::OutputSpec plain;
      plain.format.notation = output::NumberFormat::Notation::general;
      plain.format.digits = 8;
      model_->set_table(*t, plain);
    }
  } catch (const std::exception &e) {
    note = tr("Output settings error: %1").arg(QString::fromUtf8(e.what()));
    model_->set_table(*t, {});
  }
  info_->setText(note);
}

void TablesPanel::export_table(output::Format f) {
  const auto &t = model_->table();
  if (t.columns.empty()) return;
  const bool csv = f == output::Format::csv;
  const QString path = QFileDialog::getSaveFileName(
      this, csv ? tr("Save table as CSV") : tr("Save table as Parquet"),
      QString::fromStdString(t.key) + (csv ? QStringLiteral(".csv") : QStringLiteral(".parquet")),
      csv ? tr("CSV files (*.csv)") : tr("Parquet files (*.parquet)"));
  if (path.isEmpty()) return;
  try {
    const output::OutputSpec s = shaped_->isChecked() ? spec_ : output::OutputSpec{};
    output::Dataset ctx;
    if (data_) {
      ctx.title = data_->title;
      ctx.metadata = data_->metadata;
    }
    if (csv) {
      std::ofstream o(std::filesystem::path(path.toStdWString()), std::ios::binary);
      if (!o) throw std::runtime_error("cannot write the file");
      output::write_csv(o, t, ctx, s);
    } else {
      output::write_parquet(std::filesystem::path(path.toStdWString()), t, ctx, s);
    }
  } catch (const std::exception &e) {
    QMessageBox::warning(this, tr("Export failed"), QString::fromUtf8(e.what()));
  }
}

} // namespace gui
} // namespace amcset
