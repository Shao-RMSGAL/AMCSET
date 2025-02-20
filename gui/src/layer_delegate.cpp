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

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QPersistentModelIndex>
#include <QPushButton>
#include <QStandardItemModel>
#include <QTableView>
#include <qnamespace.h>

#include "layer_delegate.h"

QWidget *LayerWidgetDelegate::createEditor(QWidget *parent,
                                           const QStyleOptionViewItem &,
                                           const QModelIndex &index) const {
  if (!index.isValid()) {
    return nullptr;
  }

  LayerWidgetType type =
      static_cast<LayerWidgetType>(index.data(Qt::UserRole).toInt());

  switch (type) {
  case StartingDepthLabel: {
    auto *startingDepthLabel = new QLabel(parent);
    startingDepthLabel->setText(index.data(Qt::DisplayRole).toString());
    return startingDepthLabel;
    ;
  }
  case DepthSpinBox: {
    auto *depthSpinBox = new QDoubleSpinBox(parent);
    const auto minDepth = index.data(Qt::DisplayRole).toString().toDouble();
    const auto maxDepth = 1e10;
    depthSpinBox->setRange(minDepth, maxDepth);
    depthSpinBox->setSingleStep(0.01);
    depthSpinBox->setValue(minDepth);
    connect(depthSpinBox, &QDoubleSpinBox::valueChanged, this,
            [this]() { emit updateLayerDepths(); });
    index.data(Qt::DisplayRole);
    return depthSpinBox;
    // connect(depthSpinBox, &QDoubleSpinBox::valueChanged, this,
    //         &MainWindow::updateLayerDepths);
  }
  case LayerLineEdit: {
    auto *layerLineEdit = new QLineEdit(parent);
    layerLineEdit->setText(index.data(Qt::DisplayRole).toString());
    return layerLineEdit;
  }
  case DensitySpinBox: {
    auto *densityField = new QDoubleSpinBox(parent);
    const auto minDensity = index.data(Qt::DisplayRole).toDouble();
    const auto maxDensity = 1e10;
    densityField->setRange(minDensity, maxDensity);
    densityField->setSingleStep(0.01);
    return densityField;
  }
  case PhaseDropDown: {
    auto *phaseDropDown = new QComboBox(parent);
    auto phaseList = QStringList{"Solid", "Liquid", "Gas"};
    phaseDropDown->addItems(phaseList);
    phaseDropDown->setCurrentIndex(0);
    return phaseDropDown;
  }
  default:
    return nullptr;
  }
}
