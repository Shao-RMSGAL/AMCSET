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

#pragma once

#include <QStyledItemDelegate>
#include <qcoreevent.h>

enum LayerWidgetType {
  // DepthSpinBox,
  LayerLineEdit,
  DensitySpinBox,
  ElementLineEdit,
  FractionSpinBox,
  PercentLabel,
  ElementDeleteButton,
  DisplacementEnergySpinBox,
};

class LayerWidgetDelegate : public QStyledItemDelegate {
  Q_OBJECT

public:
  QWidget *createEditor(QWidget *parent, const QStyleOptionViewItem &,
                        const QModelIndex &index) const override;

signals:

  void deleteElement(QModelIndex index) const;
  //   void updateLayerDepths() const;

protected:
  using QStyledItemDelegate::QStyledItemDelegate;
};
