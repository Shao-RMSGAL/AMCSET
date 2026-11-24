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
#include "setup/periodic_table.h"

#include "element_data.h"

#include <QGridLayout>
#include <QLabel>
#include <QToolButton>
#include <QVBoxLayout>

namespace amcset {
namespace gui {

namespace {

// Row and column (0-based) of element z in the standard layout; the
// lanthanides and actinides go to rows 8 and 9.
std::pair<int, int> position(unsigned z) {
  if (z == 1) return {0, 0};
  if (z == 2) return {0, 17};
  const auto block = [](unsigned z0, int row) -> std::pair<int, int> { return {row, static_cast<int>(z0)}; };
  if (z <= 10) return z <= 4 ? block(z - 3, 1) : block(z - 5 + 12, 1);
  if (z <= 18) return z <= 12 ? block(z - 11, 2) : block(z - 13 + 12, 2);
  if (z <= 36) return block(z - 19, 3);
  if (z <= 54) return block(z - 37, 4);
  if (z <= 56) return block(z - 55, 5);
  if (z <= 70) return {8, static_cast<int>(z - 57) + 2};
  if (z <= 86) return block(z - 71 + 2, 5);
  if (z <= 88) return block(z - 87, 6);
  return {9, static_cast<int>(z - 89) + 2};
}

} // namespace

PeriodicTableDialog::PeriodicTableDialog(QWidget *parent, unsigned current) : QDialog(parent) {
  setWindowTitle(tr("Periodic table"));
  auto *grid = new QGridLayout();
  grid->setSpacing(2);
  for (unsigned z = 1; z <= common::Element::max_z; ++z) {
    const auto &el = common::Element::by_z(z);
    auto *b = new QToolButton(this);
    b->setText(QStringLiteral("%1\n%2").arg(z).arg(QString::fromStdString(std::string(el.symbol()))));
    b->setToolTip(QStringLiteral("%1  (%2 amu, %3 g/cm3)")
                      .arg(QString::fromStdString(std::string(el.name())))
                      .arg((el.natural_mass() / common::atomic_mass_unit).value(), 0, 'g', 6)
                      .arg(common::in_units(el.solid_density(), common::gram_per_cubic_centimeter), 0, 'g', 4));
    b->setMinimumSize(40, 40);
    b->setCheckable(true);
    b->setChecked(z == current);
    b->setAutoRaise(false);
    const auto [r, c] = position(z);
    grid->addWidget(b, r, c);
    connect(b, &QToolButton::clicked, this, guarded(this, tr("Periodic table (clicked)"), [this, z] {
      selected_ = z;
      accept();
    }));
  }
  grid->setRowMinimumHeight(7, 10);
  auto *lay = new QVBoxLayout(this);
  lay->addLayout(grid);
  lay->addWidget(new QLabel(tr("Click an element. Elements up to uranium (Z = 92) are available."), this));
}

unsigned PeriodicTableDialog::pick(QWidget *parent, unsigned current) {
  PeriodicTableDialog d(parent, current);
  return d.exec() == QDialog::Accepted ? d.selected() : 0;
}

} // namespace gui
} // namespace amcset
