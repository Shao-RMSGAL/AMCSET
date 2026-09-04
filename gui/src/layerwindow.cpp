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

// Standard includes
#include <iostream>

// Qt includes
#include <QBoxLayout>
#include <QLabel>
#include <QStandardItem>
#include <QStandardItemModel>
#include <QTableView>

// Local includes
#include "layerwindow.h"

LayerWindow::LayerWindow(QStandardItemModel *compositionModel,
                         const QString &layerName, QWidget *parent)
    : QWidget(parent) {

  auto *compositionTable = new QTableView(this);
  compositionTable->setModel(compositionModel);
  compositionTable->resizeColumnsToContents();
  compositionTable->resizeRowsToContents();

  QPushButton *okButton = new QPushButton("OK", this);
  connect(okButton, &QPushButton::clicked, this, &LayerWindow::onOkClicked);

  QPushButton *cancelButton = new QPushButton("Cancel", this);
  connect(cancelButton, &QPushButton::clicked, this,
          &LayerWindow::onCancelClicked);

  // Layout
  auto *centralWidget = new QWidget(this);
  auto *hLayout = new QHBoxLayout();
  hLayout->addWidget(okButton);
  hLayout->addWidget(cancelButton);

  auto *vLayout = new QVBoxLayout(centralWidget);
  vLayout->addWidget(compositionTable);
  vLayout->addLayout(hLayout);

  setLayout(vLayout);
  setWindowTitle("Edit " + layerName);
  setMinimumSize(480, 480);
  resize(480, 480);
}
