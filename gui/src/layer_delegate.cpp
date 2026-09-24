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
#include <QStylePainter>
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
  case LayerLineEdit: {
    auto *layerLineEdit = new QLineEdit(parent);
    layerLineEdit->setText(index.data(Qt::DisplayRole).toString());
    return layerLineEdit;
  }
  case ElementLineEdit: {
    auto *elementLineEntry = new QLineEdit(parent);
    elementLineEntry->setText(index.data(Qt::DisplayRole).toString());
    return elementLineEntry;
  }
  case FractionSpinBox: {
    auto *fractionSpinBox = new QDoubleSpinBox(parent);
    fractionSpinBox->setDecimals(6);
    fractionSpinBox->setRange(0.0, 1.0e9);
    fractionSpinBox->setValue(index.data(Qt::DisplayRole).toDouble());
    return fractionSpinBox;
  }
  case PercentLabel: {
    auto *percentLabel = new QLabel(parent);
    percentLabel->setText(index.data(Qt::DisplayRole).toString());
    return percentLabel;
  }
  case ElementDeleteButton: {
    auto *elementDeleteButton = new QPushButton(parent);
    elementDeleteButton->setText(QString());
    elementDeleteButton->setToolTip("Remove element");
    elementDeleteButton->setIcon(
        QApplication::style()->standardIcon(QStyle::SP_TitleBarCloseButton));
    elementDeleteButton->setIconSize(QSize(12, 12));
    connect(elementDeleteButton, &QPushButton::clicked, this,
            [this, index]() { emit deleteElement(index); });
    return elementDeleteButton;
  }
  case DisplacementEnergySpinBox: {
    auto *displacementEnergySpinBox = new QDoubleSpinBox(parent);
    auto minValue = 0.01;
    auto maxValue = 1e6;
    displacementEnergySpinBox->setRange(minValue, maxValue);
    displacementEnergySpinBox->setValue(index.data(Qt::DisplayRole).toDouble());
    return displacementEnergySpinBox;
  }
  default:
    return nullptr;
  }
}

void LayerWidgetDelegate::paint(QPainter *painter,
                               const QStyleOptionViewItem &option,
                               const QModelIndex &index) const {
  if (index.data(Qt::UserRole).toInt() == ElementDeleteButton) {
    QStyleOptionButton buttonOption;
    buttonOption.rect = option.rect;
    buttonOption.text.clear();
    buttonOption.icon = QApplication::style()->standardIcon(
        QStyle::SP_TitleBarCloseButton);
    buttonOption.iconSize = QSize(12, 12);
    buttonOption.state = option.state;
    buttonOption.features = QStyleOptionButton::Flat;
    if (option.state & QStyle::State_Selected) {
      buttonOption.state |= QStyle::State_Sunken;
    }
    QApplication::style()->drawControl(QStyle::CE_PushButton, &buttonOption,
                                       painter, option.widget);
    return;
  }

  QStyledItemDelegate::paint(painter, option, index);
}

bool LayerWidgetDelegate::editorEvent(QEvent *event,
                                     QAbstractItemModel *model,
                                     const QStyleOptionViewItem &option,
                                     const QModelIndex &index) {
  if (!index.isValid()) {
    return false;
  }

  if (index.data(Qt::UserRole).toInt() == ElementDeleteButton &&
      (event->type() == QEvent::MouseButtonPress ||
       event->type() == QEvent::MouseButtonRelease)) {
    auto *mouseEvent = static_cast<QMouseEvent *>(event);
    if (mouseEvent && mouseEvent->button() == Qt::LeftButton &&
        option.rect.contains(mouseEvent->pos())) {
      emit deleteElement(index);
      return true;
    }
  }

  return QStyledItemDelegate::editorEvent(event, model, option, index);
}
