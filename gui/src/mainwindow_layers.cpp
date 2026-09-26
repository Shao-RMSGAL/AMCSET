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

#include "mainwindow.h"

#include <QApplication>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QHeaderView>
#include <QLabel>
#include <QStandardItem>
#include <QStandardItemModel>
#include <QTableView>
#include <QVBoxLayout>

#include "layer_delegate.h"
#include "layerwindow.h"

namespace amcset {
namespace gui {

QGroupBox *MainWindow::createIncidentParticleForms(QWidget *parent) {
  auto formGroupBox = new QGroupBox("Incident particle properties", parent);
  incidentFormLayout = new QFormLayout(formGroupBox);

  incidentElectron = new QCheckBox(this);
  incidentElectron->setChecked(false);
  incidentFormLayout->addRow("Incident electron?", incidentElectron);

  damageCascade = new QCheckBox(this);
  damageCascade->setChecked(true);
  incidentFormLayout->addRow("Simulate damage cascade?", damageCascade);

  incidentEnergy = new QDoubleSpinBox(this);
  incidentEnergy->setRange(0.0, 10000.0);
  incidentEnergy->setDecimals(1);
  const auto defaultEnergy = 10;
  incidentEnergy->setValue(defaultEnergy);
  incidentFormLayout->addRow("Incident Energy (keV)", incidentEnergy);

  incidentElement = new QLineEdit("H", this);
  incidentElement->setPlaceholderText("Element symbol");
  incidentElement->setMaxLength(2);
  incidentElementFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Element", incidentElement);

  incidentMass = new QSpinBox(this);
  incidentMass->setRange(1, 295);
  incidentMass->setValue(1);
  incidentMassFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Isotope (A)", incidentMass);

  incidentCharge = new QSpinBox(this);
  incidentCharge->setRange(1, 118);
  incidentCharge->setValue(1);
  incidentCharge->setReadOnly(true);
  incidentCharge->setButtonSymbols(QAbstractSpinBox::NoButtons);
  incidentCharge->setFocusPolicy(Qt::NoFocus);
  incidentChargeFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Z-number", incidentCharge);

  incidentAbundanceLabel = new QLabel("Relative abundance: n/a", this);
  incidentAbundanceFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Relative abundance", incidentAbundanceLabel);

  incidentAtomicMassLabel = new QLabel("Relative atomic mass: n/a", this);
  incidentAtomicMassFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Relative atomic mass", incidentAtomicMassLabel);

  connect(incidentElement, &QLineEdit::editingFinished, this, [this]() {
    const auto rawElement = incidentElement->text().trimmed();
    const auto symbol = rawElement.isEmpty() ? std::string("H") : rawElement.toStdString();
    const auto atomicNumber = common::IsotopeData::element_to_atomic_number(symbol);
    if (atomicNumber == 0) {
      incidentElement->setText("H");
    }
    updateIncidentIonSelection();
  });
  connect(incidentMass, &QSpinBox::editingFinished, this, [this]() {
    if (!incidentElement || !incidentMass || !incidentCharge || !incidentAbundanceLabel || !incidentAtomicMassLabel) {
      return;
    }

    const auto rawElement = incidentElement->text().trimmed();
    const auto symbol = rawElement.isEmpty() ? std::string("H") : rawElement.toStdString();
    const auto atomicNumber = common::IsotopeData::element_to_atomic_number(symbol);
    const size_t effectiveZ = atomicNumber != 0 ? atomicNumber : 1;
    auto isotopeEntries = common::IsotopeData::get_relative_abundance(effectiveZ);

    const auto selectedMass = static_cast<size_t>(incidentMass->value());
    const auto it = std::find_if(isotopeEntries.begin(), isotopeEntries.end(),
                                [selectedMass](const auto &entry) {
                                  return entry.first == selectedMass;
                                });

    if (it == isotopeEntries.end()) {
      size_t mostCommonMass = 1;
      double mostCommonAbundance = 0.0;
      for (const auto &[mass, abundance] : isotopeEntries) {
        if (abundance > mostCommonAbundance) {
          mostCommonAbundance = abundance;
          mostCommonMass = mass;
        }
      }
      incidentMass->setValue(static_cast<int>(mostCommonMass));
    }

    updateIncidentIonSelection();
  });
  connect(incidentElectron, &QCheckBox::checkStateChanged, this,
          &MainWindow::toggleIncidentIonFormVisibility);

  incidentAngle = new QDoubleSpinBox(this);
  incidentAngle->setRange(0.0, 90.0);
  incidentAngle->setDecimals(2);
  incidentFormLayout->addRow("Incident Angle (degree)", incidentAngle);

  bombardmentCountSpinBox = new QSpinBox(this);
  bombardmentCountSpinBox->setRange(1, std::numeric_limits<int>::max());
  const auto defaultBombardments = 1;
  bombardmentCountSpinBox->setValue(defaultBombardments);
  incidentFormLayout->addRow("Number of Bombardments", bombardmentCountSpinBox);

  threadCountSpinBox = new QSpinBox(this);
  const auto numberOfThreads = std::thread::hardware_concurrency();
  threadCountSpinBox->setRange(1, numberOfThreads);
  threadCountSpinBox->setValue(numberOfThreads);
  incidentFormLayout->addRow("Number of Threads", threadCountSpinBox);

  rangeSpinBox = new QDoubleSpinBox(this);
  const auto minRange = 0.1;
  const auto maxRange = std::numeric_limits<decltype(minRange)>::max();
  const auto defaultDepth = 1000;
  rangeSpinBox->setRange(minRange, maxRange);
  rangeSpinBox->setValue(defaultDepth);
  connect(rangeSpinBox, &QDoubleSpinBox::valueChanged, view3D,
          &Widget3DView::setRange);
    connect(rangeSpinBox, &QDoubleSpinBox::editingFinished, this,
      &MainWindow::saveState);
      connect(rangeSpinBox, qOverload<double>(&QDoubleSpinBox::valueChanged),
        this, [this](double) { saveState(); });
  rangeSpinBox->setRange(minRange, maxRange);
  rangeSpinBox->setValue(defaultDepth);
  incidentFormLayout->addRow("Range (Angstrom)", rangeSpinBox);

  use3DView = new QCheckBox(this);
  use3DView->setChecked(false);
  connect(use3DView, &QCheckBox::checkStateChanged, this,
          &MainWindow::update3DViewMode);
  use3DView->hide();

  use2DView = new QCheckBox(this);
  use2DView->setChecked(true);
  connect(use2DView, &QCheckBox::checkStateChanged, this,
          &MainWindow::update2DViewMode);
  use2DView->hide();

  counterLabel = new QLabel(this);
  counterLabel->setText(
      QString("Bombardment %1 out of %2").arg(0).arg(defaultBombardments));
  counter = new Counter(this);

  return formGroupBox;
}

void MainWindow::toggleIncidentIonFormVisibility(int state) {
  const bool isIon = state != Qt::Checked;
  incidentFormLayout->setRowVisible(incidentElementFormIndex, isIon);
  incidentFormLayout->setRowVisible(incidentMassFormIndex, isIon);
  incidentFormLayout->setRowVisible(incidentChargeFormIndex, isIon);
  incidentFormLayout->setRowVisible(incidentAbundanceFormIndex, isIon);
  incidentFormLayout->setRowVisible(incidentAtomicMassFormIndex, isIon);
  incidentElement->setVisible(isIon);
  incidentMass->setVisible(isIon);
  incidentCharge->setVisible(isIon);
  incidentAbundanceLabel->setVisible(isIon);
  incidentAtomicMassLabel->setVisible(isIon);
  if (isIon) {
    updateIncidentIonSelection();
  }
}

void MainWindow::updateIncidentIonSelection() {
  if (!incidentElement || !incidentMass || !incidentCharge || !incidentAbundanceLabel || !incidentAtomicMassLabel) {
    return;
  }

  const auto rawElement = incidentElement->text().trimmed();
  const std::string symbol = rawElement.isEmpty() ? "H" : rawElement.toStdString();
  const auto atomicNumber = common::IsotopeData::element_to_atomic_number(symbol);
  const size_t effectiveZ = atomicNumber != 0 ? atomicNumber : 1;
  const std::string effectiveSymbol = atomicNumber != 0 ? symbol : "H";

  if (rawElement.isEmpty()) {
    incidentElement->setText(QString::fromStdString(effectiveSymbol));
  }

  const auto isotopeEntries = common::IsotopeData::get_relative_abundance(effectiveZ);
  size_t mostCommonMass = 1;
  double mostCommonAbundance = 0.0;
  for (const auto &[mass, abundance] : isotopeEntries) {
    if (abundance > mostCommonAbundance) {
      mostCommonAbundance = abundance;
      mostCommonMass = mass;
    }
  }

  size_t selectedMass = static_cast<size_t>(incidentMass->value());
  if (selectedMass == 0 || selectedMass > 295) {
    selectedMass = mostCommonMass;
  }

  bool validIsotope = false;
  for (const auto &[mass, abundance] : isotopeEntries) {
    if (mass == selectedMass) {
      validIsotope = true;
      break;
    }
  }

  if (!validIsotope) {
    selectedMass = mostCommonMass;
  }

  incidentCharge->setValue(static_cast<int>(effectiveZ));
  incidentMass->setValue(static_cast<int>(selectedMass));

  double isotopeAbundance = 0.0;
  double relativeAtomicMass = 0.0;
  for (const auto &entry : common::IsotopeData::database) {
    if (entry.atomic_number == effectiveZ && entry.mass_number == selectedMass) {
      isotopeAbundance = entry.isotopic_composition;
      relativeAtomicMass = entry.relative_atomic_mass;
      break;
    }
  }
  if (isotopeAbundance <= 0.0) {
    for (const auto &entry : common::IsotopeData::database) {
      if (entry.atomic_number == effectiveZ && entry.mass_number == mostCommonMass) {
        isotopeAbundance = entry.isotopic_composition;
        relativeAtomicMass = entry.relative_atomic_mass;
        break;
      }
    }
  }

  if (effectiveZ == 1 && effectiveSymbol == "H") {
    incidentMass->setRange(1, 7);
  } else {
    incidentMass->setRange(1, 295);
  }

  incidentAbundanceLabel->setText(QString("Relative abundance: %1%")
                                      .arg(isotopeAbundance > 0.0 ? QString::number(isotopeAbundance * 100.0, 'f', 4) : "0.0000"));
  incidentAtomicMassLabel->setText(QString("Relative atomic mass: %1")
                                      .arg(relativeAtomicMass > 0.0 ? QString::number(relativeAtomicMass, 'f', 6) : "n/a"));

  incidentElement->setText(QString::fromStdString(effectiveSymbol));
}

void MainWindow::handleLayerButtonClicked(int row) { openLayerWindow(row); }

void MainWindow::updateLayerDepths() {
  auto totalDepth = 0.0;
  auto range = rangeSpinBox->value();
  for (int row = 1; row < layerModel->rowCount(); row++) {
    auto prevLabelIndex = layerModel->index(row - 1, 0);
    auto prevEntryIndex = layerModel->index(row - 1, 1);
    auto labelIndex = layerModel->index(row, 0);
    auto prevDepth =
        qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(prevEntryIndex))
            ->value();
    auto prevStartingDepth =
        qobject_cast<QLabel *>(layerTable->indexWidget(prevLabelIndex))
            ->text()
            .toDouble();
    auto *startingDepthLabel =
        qobject_cast<QLabel *>(layerTable->indexWidget(labelIndex));
    totalDepth = prevStartingDepth + prevDepth;
    startingDepthLabel->setText(QString::number(totalDepth));
  }

  auto depthEntryIndex = layerModel->index(layerModel->rowCount() - 1, 1);
  auto lastDepth =
      qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(depthEntryIndex))
          ->value();
  totalDepth += lastDepth;

  if (totalDepth > range) {
    rangeSpinBox->setValue(totalDepth);
  }
}

double calculateLayerDensity(size_t row) { return row; }

void MainWindow::createLayerRow() {
  const int row = layerModel->rowCount();
  layerModel->insertRow(row);

  const int startingDepthColumn = 0;
  const int depthSpinBoxColumn = 1;
  const int layerLineEditColumn = 2;
  const int densitySpinBoxColumn = 3;
  const int compositionViewColumn = 4;
  const int phaseDropDownColumn = 5;
  const int removeButtonColumn = 6;

  const auto startingDepthLabelIndex = layerModel->index(row, startingDepthColumn);
  const auto depthSpinBoxIndex = layerModel->index(row, depthSpinBoxColumn);
  const auto densitySpinBoxIndex = layerModel->index(row, densitySpinBoxColumn);
  const auto compositionViewIndex = layerModel->index(row, compositionViewColumn);
  const auto phaseDropDownIndex = layerModel->index(row, phaseDropDownColumn);
  const auto removeButtonIndex = layerModel->index(row, removeButtonColumn);

  auto *startingDepthLabel = new QLabel(this);
  startingDepthLabel->setText(QString::number(row));
  const auto minDepth = 0.01;
  const auto maxDepth = std::numeric_limits<decltype(minDepth)>::max();
  const auto defaultDepth = 1000;
  auto *depthSpinBox = new QDoubleSpinBox();
  depthSpinBox->setValue(minDepth);
  depthSpinBox->setSingleStep(minDepth);
  depthSpinBox->setRange(minDepth, maxDepth);
  depthSpinBox->setValue(defaultDepth);
  connect(depthSpinBox, &QDoubleSpinBox::valueChanged, this,
          &MainWindow::updateLayerDepths);
  auto *layerLineEdit = new QStandardItem(QString("Layer %1").arg(row + 1));
  layerLineEdit->setData(LayerLineEdit, Qt::UserRole);

  const auto minDensity = 0.01;
  const auto maxDensity = std::numeric_limits<decltype(minDensity)>::max();
  const auto defaultDensity = 1000.0;
  auto *densitySpinBox = new QDoubleSpinBox();
  densitySpinBox->setRange(minDensity, maxDensity);
  densitySpinBox->setValue(defaultDensity);

  const auto compositionColumns = 5;
  auto *compositionModel = new QStandardItemModel(0, compositionColumns, this);
  auto *compositionTable = new QTableView(this);
  compositionTable->horizontalHeader()->setSectionResizeMode(
      QHeaderView::ResizeToContents);
  auto *compositionDelegate = new LayerWidgetDelegate(compositionTable);
  compositionTable->setItemDelegate(compositionDelegate);
  connect(compositionDelegate, &LayerWidgetDelegate::deleteElement, this,
          [compositionModel](QModelIndex index) {
            compositionModel->removeRow(index.row());
          });
  compositionModel->setHorizontalHeaderLabels(QStringList{
      "Element",
      "Atomic\nProportion",
      "Atomic %",
      "Displacement\nEnergy (eV)",
      "Remove",
  });

  compositionTable->setModel(compositionModel);
  compositionTable->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  auto *compositionLayoutWidget = new QWidget(this);
  compositionLayoutWidget->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  auto *compositionLayout = new QVBoxLayout();
  compositionLayout->setContentsMargins(0, 0, 0, 0);
  compositionLayoutWidget->setLayout(compositionLayout);
  const auto elementColumn = 0;
  const auto fractionColumn = 1;
  const auto percentColumn = 2;
  const auto displacementEnergyColumn = 3;
  const auto elementDeleteButtonColumn = 4;
  auto recalculateAtomicPercentages = [compositionModel, elementColumn,
                                      fractionColumn, percentColumn,
                                      displacementEnergyColumn,
                                      elementDeleteButtonColumn]() {
    double totalProportion = 0.0;
    for (int row = 0; row < compositionModel->rowCount(); ++row) {
      const auto *fractionItem = compositionModel->item(row, fractionColumn);
      if (fractionItem) {
        totalProportion += fractionItem->data(Qt::DisplayRole).toDouble();
      }
    }

    for (int row = 0; row < compositionModel->rowCount(); ++row) {
      auto *percentItem = compositionModel->item(row, percentColumn);
      auto *fractionItem = compositionModel->item(row, fractionColumn);
      if (!percentItem || !fractionItem) {
        continue;
      }

      const auto proportion = fractionItem->data(Qt::DisplayRole).toDouble();
      const auto percent = totalProportion > 0.0 ? (proportion / totalProportion) * 100.0 : 0.0;
      percentItem->setData(percent, Qt::DisplayRole);
    }
  };

  auto compositionAddElementFunction = [compositionModel, compositionTable,
                                        elementColumn, fractionColumn,
                                        percentColumn,
                                        displacementEnergyColumn,
                                        elementDeleteButtonColumn,
                                        recalculateAtomicPercentages]() {
    const auto defaultElement = "H";
    const auto defaultFraction = 1.0;
    auto index = compositionModel->rowCount();
    auto *elementLabel = new QStandardItem();
    elementLabel->setData(ElementLineEdit, Qt::UserRole);
    elementLabel->setData(defaultElement, Qt::DisplayRole);

    auto *fractionSpinBox = new QStandardItem();
    fractionSpinBox->setData(FractionSpinBox, Qt::UserRole);
    fractionSpinBox->setData(defaultFraction, Qt::DisplayRole);

    auto *percentLabel = new QStandardItem();
    percentLabel->setData(PercentLabel, Qt::UserRole);
    percentLabel->setData(defaultFraction, Qt::DisplayRole);

    const auto displacementEnergyMinValue = 0.01;
    const auto displacementEnergyDefaultValue = 40.0;
    auto *displacementEnergySpinBox = new QStandardItem(displacementEnergyMinValue);
    displacementEnergySpinBox->setData(
        QString::number(displacementEnergyDefaultValue), Qt::DisplayRole);
    displacementEnergySpinBox->setData(DisplacementEnergySpinBox, Qt::UserRole);

    auto *elementDeleteButton = new QStandardItem();
    elementDeleteButton->setData(ElementDeleteButton, Qt::UserRole);
    elementDeleteButton->setData(
        QApplication::style()->standardIcon(QStyle::SP_TitleBarCloseButton),
        Qt::DecorationRole);
    elementDeleteButton->setText(QString());

    compositionModel->insertRow(index);
    compositionModel->setItem(index, elementColumn, elementLabel);
    compositionModel->setItem(index, fractionColumn, fractionSpinBox);
    compositionModel->setItem(index, percentColumn, percentLabel);
    compositionModel->setItem(index, displacementEnergyColumn,
                              displacementEnergySpinBox);
    compositionModel->setItem(index, elementDeleteButtonColumn,
                              elementDeleteButton);
    recalculateAtomicPercentages();
    compositionTable->resizeColumnsToContents();
    compositionTable->resizeRowsToContents();
  };
  connect(compositionModel, &QStandardItemModel::itemChanged, this,
          [compositionModel, fractionColumn, percentColumn,
           recalculateAtomicPercentages](QStandardItem *item) {
            if (!item || item->column() != fractionColumn) {
              return;
            }
            recalculateAtomicPercentages();
          });
  compositionAddElementFunction();

  auto *compositionAddElementButton = new QPushButton("Add Element");
  connect(compositionAddElementButton, &QPushButton::clicked, this,
          compositionAddElementFunction);

  compositionLayout->addWidget(compositionTable, 1);
  compositionLayout->addWidget(compositionAddElementButton);

  for (int compositionRow = 0; compositionRow < compositionModel->rowCount();
       ++compositionRow) {
    QModelIndex index = compositionModel->index(compositionRow, 0);
    QWidget *editor = compositionDelegate->createEditor(
        compositionTable, QStyleOptionViewItem(), index);
    compositionDelegate->updateEditorGeometry(editor, QStyleOptionViewItem(),
                                              index);
  }

  auto *phaseDropDown = new QComboBox(this);
  phaseDropDown->addItems(QStringList{"Solid", "Liquid", "Gas"});
  phaseDropDown->setCurrentIndex(0);

  auto *removeButton = new QPushButton(this);
  QPersistentModelIndex removeButtonPersistentIndex(removeButtonIndex);
  removeButton->setIcon(
      QApplication::style()->standardIcon(QStyle::SP_DialogCloseButton));
  removeButton->setToolTip("Delete");
  connect(removeButton, &QPushButton::clicked, this,
          [this, removeButtonPersistentIndex]() {
            if (removeButtonPersistentIndex.isValid()) {
              layerModel->removeRow(removeButtonPersistentIndex.row());
            }
          });

  layerModel->setItem(row, layerLineEditColumn, layerLineEdit);
  layerTable->setIndexWidget(startingDepthLabelIndex, startingDepthLabel);
  layerTable->setIndexWidget(depthSpinBoxIndex, depthSpinBox);
  layerTable->setIndexWidget(densitySpinBoxIndex, densitySpinBox);
  layerTable->setIndexWidget(compositionViewIndex, compositionLayoutWidget);
  layerTable->setIndexWidget(phaseDropDownIndex, phaseDropDown);
  layerTable->setIndexWidget(removeButtonPersistentIndex, removeButton);

  layerTable->resizeColumnsToContents();
  layerTable->resizeRowsToContents();
  layerTable->show();

  updateLayerDepths();
}

QGroupBox *MainWindow::createLayerForms(QWidget *parent) {
  auto layerGroupBox = new QGroupBox("Layers", parent);
  auto layerVBox = new QVBoxLayout(layerGroupBox);
  auto addLayerButton = new QPushButton("Add Layer");
  connect(addLayerButton, &QPushButton::clicked, this,
          &MainWindow::createLayerRow);

  layerHeaders = QStringList{"Starting Depth\n(Angstrom)",
                             "Depth\n(Angstrom)",
                             "Layer",
                             QString("Density\n(kg/m%1)").arg(QChar(0x00B3)),
                             "Composition",
                             "Phase",
                             "Remove"};
  layerModel = new QStandardItemModel(0, layerHeaders.count());
  layerModel->setHorizontalHeaderLabels(layerHeaders);

  layerTable = new QTableView();
  auto *layerDelegate = new LayerWidgetDelegate(layerTable);
  layerTable->setItemDelegate(layerDelegate);

  layerTable->horizontalHeader()->setStyleSheet(
      "QHeaderView::section { padding: 4px; }");
  layerTable->setModel(layerModel);
  layerTable->horizontalHeader()->setSectionResizeMode(
      QHeaderView::ResizeToContents);
  const auto compositionColumn = 4;
  layerTable->horizontalHeader()->setSectionResizeMode(compositionColumn,
                                                       QHeaderView::Stretch);
  layerTable->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

  layerVBox->addWidget(layerTable);
  layerVBox->addWidget(addLayerButton);

  return layerGroupBox;
}

void MainWindow::openLayerWindow(int row) {
  const int compositionColumn = 4;
  const auto compositionIndex = layerModel->index(row, compositionColumn);
  const int layerNameColumn = 2;

  auto *compositionViewWidget =
      qobject_cast<QWidget *>(layerTable->indexWidget(compositionIndex));
  auto *compositionViewLayout =
      compositionViewWidget->findChild<QVBoxLayout *>();
  auto *compositionView =
      qobject_cast<QTableView *>(compositionViewLayout->itemAt(0)->widget());
  auto *compositionModel =
      qobject_cast<QStandardItemModel *>(compositionView->model());

  compositionView = new QTableView(this);
  auto *compositionDelegate = new LayerWidgetDelegate(compositionView);
  compositionView->setItemDelegate(compositionDelegate);
  compositionView->setModel(compositionModel);
  auto layerName =
      layerModel->index(row, layerNameColumn).data(Qt::DisplayRole).toString();

  LayerWindow *layerWindow = new LayerWindow(compositionModel, layerName,
                                             this);
  layerWindow->setWindowFlags(Qt::Window);
  connect(layerWindow, &LayerWindow::okClicked, this, &MainWindow::onOkClicked);
  connect(layerWindow, &LayerWindow::cancelClicked, this,
          &MainWindow::onCancelClicked);
  layerWindow->resize(400, 400);
  layerWindow->show();
}

} // namespace gui
} // namespace amcset
