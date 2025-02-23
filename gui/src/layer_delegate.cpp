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

#include <iostream>

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
  // case DepthSpinBox: {
  //   auto *depthSpinBox = new QDoubleSpinBox(parent);
  //   const auto minDepth = index.data(Qt::DisplayRole).toString().toDouble();
  //   const auto maxDepth = 1e10;
  //   depthSpinBox->setRange(minDepth, maxDepth);
  //   depthSpinBox->setSingleStep(0.01);
  //   depthSpinBox->setValue(minDepth);
  //   connect(depthSpinBox, &QDoubleSpinBox::valueChanged, this,
  //           [this]() { emit updateLayerDepths(); });
  //   return depthSpinBox;
  // }
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
  case ElementLineEdit: {
    auto *elementLineEntry = new QLineEdit(parent);
    elementLineEntry->setText(index.data(Qt::DisplayRole).toString());
    return elementLineEntry;
  }
  case FractionSpinBox: {
    auto *fractionSpinBox = new QDoubleSpinBox(parent);
    // fractionSpinBox->setValue(
    //     index.data(Qt::DisplayRole).toString().toDouble());
    return fractionSpinBox;
  }
  case PercentLabel: {
    auto *percentLabel = new QLabel(parent);
    percentLabel->setText(index.data(Qt::DisplayRole).toString());
    return percentLabel;
  }
  case ElementDeleteButton: {
    auto *elementDeleteButton = new QPushButton(parent);
    elementDeleteButton->setIcon(
        QApplication::style()->standardIcon(QStyle::SP_DialogCloseButton));
    connect(elementDeleteButton, &QPushButton::clicked, this,
            [this, index]() { emit deleteElement(index); });
    return elementDeleteButton;
  }
  case DisplacementEnergySpinBox: {
    auto *displacementEnergySpinBox = new QDoubleSpinBox(parent);
    auto minValue = 0.01;
    auto maxValue = 1e6;
    displacementEnergySpinBox->setRange(minValue, maxValue);
    return displacementEnergySpinBox;
  }
  default:
    return nullptr;
  }
}
