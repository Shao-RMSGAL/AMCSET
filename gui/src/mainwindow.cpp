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

// Standard library includes
#include <exception>
#include <limits>

// Qt includes
#include <QCheckBox>
#include <QComboBox>
#include <QDebug>
#include <QFormLayout>
#include <QHeaderView>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QSizePolicy>
#include <QSplitter>
#include <QStatusBar>
#include <QTableView>
#include <QToolBar>
#include <QVBoxLayout>
#include <QtContainerFwd>
#include <QtWidgets>
#include <qabstractitemmodel.h>
#include <qboxlayout.h>
#include <qnamespace.h>
#include <qobject.h>
#include <qstandarditemmodel.h>
#include <sys/wait.h>
#include <unistd.h>

// Local includes
#include "amcset_common.h"
#include "amcset_utilities.h"
#include "application.h"
#include "canvas.h"
#include "layer_delegate.h"
#include "layerwindow.h"
#include "mainwindow.h"

namespace amcset {
namespace gui {

void MainWindow::runSimulation() {
  // Create AMCSET settings (TODO: Read this from a file and allow a user to
  // configure)

  // TODO: Figure out how to set these parameters
  common::energy_quantity electron_stopping_energy =
      40.0 * common::electron_volt; // Energy at which electrons stop.
  common::energy_quantity ion_stopping_energy =
      40.0 * common::electron_volt; // Typical for iron.
  common::energy_quantity ion_displacement_energy =
      40.0 * common::electron_volt; // Typical for iron
  bool log_single_displacement =
      false; // TEST: Test later, will consume a LOT of memory/disk space.
  size_t divisor_angle_number = 1000; // Break angle range into 1000 segements.
                                      // TEST: Test with different numbers
  size_t flying_distance_number =
      1000; // Break flying distnaces into groups of 1000. TEST: Test with
            // different numbers

  bool is_electron =
      incidentElectron
          ->checkState(); // Start with ions. TEST: Test with electrons later.
  bool enable_damage_cascade =
      damageCascade->checkState(); // TEST: Update to test cascades
  size_t mass_number = incidentMass->value();
  size_t z_number = incidentCharge->value();
  // double incident_angle =
  //     incidentAngle->value(); // TODO: Add settings for this.
  common::length_quantity range =
      rangeSpinBox->value() *
      common::angstrom; // Target range. Should be unused. HACK: Remove
                        // once confirmed to be unused.
  size_t bombardment_count =
      bombardmentCountSpinBox->value(); // TEST: Test with more later.
  common::energy_quantity incident_energy =
      incidentEnergy->value() *
      common::kilo_electron_volt; // TEST: Test with other energies
  size_t thread_count =
      threadCountSpinBox->value(); // Try with single threading, PERF:
                                   // Try with multithreading

  // PERF: Add an option to use GPU (Use Kokkos for cross-platform
  // compatability. Note that it is not available in conan. Either make a
  // conan recipe for it or use another installation method.

  auto settings = common::Simulation::Settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count);

  // Create AMCSET volume (TODO: Read this from a file and allow a user to
  // configure)

  const size_t layerCount = layerModel->rowCount();

  if (layerCount == 0) {
    throw std::logic_error("Layer count cannot be zero.");
  }

  std::vector<common::Layer> layers;
  for (size_t layerRow = 0; layerRow < layerCount; layerRow++) {
    const auto depthColumn = 1;
    const auto densityColumn = 3;
    const auto compositionColumn = 4;

    auto depthIndex = layerModel->index(layerRow, depthColumn);
    auto depth =
        qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(depthIndex))
            ->value() *
        common::angstrom;

    std::cout << "Depth: " << depth << std::endl;

    // auto density = layerModel->item(layerRow, densityColumn)
    //                    ->data(Qt::DisplayRole)
    //                    .toDouble() *
    //                common::kg_per_cubic_meter;

    auto densityIndex = layerModel->index(layerRow, densityColumn);
    auto density =
        qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(densityIndex))
            ->value() *
        common::kg_per_cubic_meter;
    std::cout << "Density: " << density << std::endl;

    auto materialIndex = layerModel->index(layerRow, compositionColumn);
    auto *materialWidget =
        qobject_cast<QWidget *>(layerTable->indexWidget(materialIndex));
    auto *materialLayout = materialWidget->layout();
    auto *materialTable =
        qobject_cast<QTableView *>(materialLayout->itemAt(0)->widget());
    auto *materialModel = materialTable->model();
    // TODO: Check for duplicate elements
    size_t elementCount = materialModel->rowCount();
    common::Layer::material_vector material;
    auto proportionSum = 0.0;
    for (size_t elementRow = 0; elementRow < elementCount; elementRow++) {
      const auto atomicSymbolColumn = 0;
      const auto proportionColumn = 1;
      const auto displacementEnergyColumn = 3;

      auto atomicSymbol = materialModel->index(elementRow, atomicSymbolColumn)
                              .data(Qt::DisplayRole)
                              .toString();
      auto atomicProportion = materialModel->index(elementRow, proportionColumn)
                                  .data(Qt::DisplayRole)
                                  .toDouble();
      auto displacementEnergy =
          materialModel->index(elementRow, displacementEnergyColumn)
              .data(Qt::DisplayRole)
              .toDouble();
      // Handle elements
      auto atomicNumber = common::IsotopeData::element_to_atomic_number(
          atomicSymbol.toStdString());
      std::vector<std::pair<size_t, double>> isotopicCompositions;
      if (atomicNumber != 0) {
        isotopicCompositions =
            common::IsotopeData::get_relative_abundance(atomicNumber);
        for (auto abundance : isotopicCompositions) {
          const auto proportion = abundance.second * atomicProportion;
          proportionSum += proportion;
          material.emplace_back(std::pair(
              proportion,
              common::Particle::Properties(atomicNumber, abundance.first)));
        }
      } else {
        // TODO: Handle specific isotopes (When user enters something like
        // 'H-2')
      }
    }

    // Normalize compositions
    for (auto &entry : material) {
      entry.first /= proportionSum;
      std::cout << "Proportion: " << entry.first
                << " Charge: " << entry.second.charge_ << " ~Mass Number: "
                << entry.second.mass_ / common::atomic_mass_unit << std::endl;
    }
    layers.push_back(common::Layer(std::move(material), depth, density));
  }

  std::cout << "Layer count: " << layers.size() << std::endl;

  // TODO: Fix segfault
  auto volume = common::Volume(std::move(layers));

  // Create the AMCSET simulation object
  auto simulation = common::Simulation(settings, std::move(volume));

  // sleep(3);
  try {
    auto exitCode =
        simulation.run_simulation(); // TODO: Check return value. Run in a
                                     // separate thread to prevent blocking.

    std::cout << "Exit status: " << exitCode << std::endl;
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
  }
}

MainWindow::MainWindow(Application *app) : app_(app) {
  // // Load window size and geometry from settings
  // auto size = app_->loadWindowSize(); // TODO: Try to get this to work
  // auto geometry = app_->loadWindowGeometry();
  // if (!geometry.isEmpty()) {
  //   restoreGeometry(geometry);
  // } else {
  //   resize(size);
  // }

  // // Load fullscreen state from settings
  // if (app_->loadFullscreenState()) {
  //   showFullScreen();
  // } else {
  //   showNormal();
  // }

  resize(1920, 1080);

  // Create a central widget
  auto *topWidget = new QWidget(this); // TODO: Check if this is OK
  setCentralWidget(topWidget);

  // Create form widgets
  infoLabel = new QLabel(tr("<i>Choose a menu option, or right-click to "
                            "invoke a context menu</i>"));
  infoLabel->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
  infoLabel->setAlignment(Qt::AlignCenter);

  auto *runCancelButtonLayoutWidget = new QWidget(this);
  auto *runCancelButtonLayout = new QHBoxLayout(runCancelButtonLayoutWidget);
  auto *runButton = new QPushButton("Run Simulation");
  runButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  connect(runButton, &QPushButton::clicked, this, &MainWindow::runSimulation);
  auto *cancelButton = new QPushButton("Cancel Run");
  cancelButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  // TODO: Figure out how to implement cancellation.

  runCancelButtonLayout->addWidget(runButton);
  runCancelButtonLayout->addWidget(cancelButton);

  auto *canvas = new Canvas;
  canvas->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);

  // Layouts
  auto *topLayout = new QVBoxLayout(); // The outermost layout.
  topLayout->setContentsMargins(5, 5, 5, 5);

  auto *formBox = new QGroupBox(this); // Contains forms for user input.
  auto *formLayout = new QVBoxLayout(formBox);

  formLayout->addWidget(infoLabel);
  auto layerSplitter = new QSplitter(Qt::Vertical);
  formLayout->addWidget(
      createIncidentParticleForms(formBox)); // Layout within the widget. //
  layerSplitter->addWidget(createLayerForms());
  layerSplitter->addWidget(runCancelButtonLayoutWidget);
  formLayout->addWidget(layerSplitter);

  auto *hPane = new QSplitter(Qt::Horizontal); // Split the forms and canvas
  hPane->addWidget(formBox);
  hPane->addWidget(canvas);
  auto formRatio = 0.6;
  const auto formWidth = (int)(this->width() * formRatio);
  const auto canvasWidth = (int)(this->width() * (1 - formRatio));
  hPane->setSizes(QList<int>{formWidth, canvasWidth});
  topLayout->addWidget(hPane);
  topWidget->setLayout(topLayout);
  createActions();
  createMenus();
  auto message = tr("A context menu is available by right-clicking");
  statusBar()->showMessage(message);
  setWindowTitle(tr("AMCSET"));
  setMinimumSize(480, 480);
}

QGroupBox *MainWindow::createIncidentParticleForms(QWidget *parent) {
  auto formGroupBox = new QGroupBox("Incident particle properties", parent);
  incidentFormLayout = new QFormLayout(formGroupBox);

  incidentElectron = new QCheckBox(this);
  incidentElectron->setChecked(false); // Default to ion simulation
  incidentFormLayout->addRow("Incident electron?", incidentElectron);

  damageCascade = new QCheckBox(this);
  damageCascade->setChecked(true); // Default to ion simulation
  incidentFormLayout->addRow("Simulate damage cascade?", damageCascade);

  incidentEnergy = new QDoubleSpinBox(this);
  incidentEnergy->setRange(0.0,
                           10000.0); // TODO: Ensure that range is reasonable;
  incidentEnergy->setDecimals(1);
  const auto defaultEnergy = 10;
  incidentEnergy->setValue(defaultEnergy);
  incidentFormLayout->addRow("Incident Energy (keV)", incidentEnergy);

  incidentMass = new QSpinBox(this);
  incidentMass->setRange(1, 295);
  const auto defaultMass = 1; // Hydrogen
  incidentEnergy->setValue(defaultMass);
  incidentMassFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Mass (amu)", incidentMass);

  incidentCharge = new QSpinBox(this);
  const auto defaultCharge = 1; // Hydrogen
  incidentEnergy->setValue(defaultCharge);
  incidentCharge->setRange(1, 118);
  incidentChargeFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Z-number", incidentCharge);

  connect(incidentElectron, &QCheckBox::stateChanged, this,
          &MainWindow::toggleIncidentIonFormVisibility);

  incidentAngle = new QDoubleSpinBox(this);
  incidentAngle->setRange(0.0, 90.0);
  incidentAngle->setDecimals(2);
  incidentFormLayout->addRow("Incident Angle (degree)", incidentAngle);

  bombardmentCountSpinBox = new QSpinBox(this);
  bombardmentCountSpinBox->setRange(
      1,
      std::numeric_limits<int>::max()); // TODO: add a cancel button to cancel
                                        // exceedingly long simulations
  const auto defaultBombardments =
      1; // TODO: Set to something more realistic for release
  bombardmentCountSpinBox->setValue(defaultBombardments);
  incidentFormLayout->addRow("Number of Bombardments", bombardmentCountSpinBox);

  threadCountSpinBox = new QSpinBox(this);
  const auto numberOfThreads = std::thread::hardware_concurrency();
  threadCountSpinBox->setRange(1, numberOfThreads);
  threadCountSpinBox->setValue(numberOfThreads);
  incidentFormLayout->addRow("Number of Threads", threadCountSpinBox);

  rangeSpinBox = new QDoubleSpinBox(this);
  const auto minRange = 0.1;
  const auto maxRange = std::numeric_limits<typeof(minRange)>::max();
  const auto defaultDepth = 1000;
  rangeSpinBox->setRange(minRange, maxRange);
  rangeSpinBox->setValue(defaultDepth);
  incidentFormLayout->addRow("Range (Angstrom)", rangeSpinBox);

  return formGroupBox;
}

void MainWindow::toggleIncidentIonFormVisibility(int state) {
  incidentFormLayout->setRowVisible(incidentMassFormIndex,
                                    state != Qt::Checked);
  incidentFormLayout->setRowVisible(incidentChargeFormIndex,
                                    state != Qt::Checked);
  incidentMass->setVisible(state != Qt::Checked);
}

void MainWindow::handleLayerButtonClicked(size_t row) { openLayerWindow(row); }

// TODO: Finish this. Return false if depths are invalid.
void MainWindow::updateLayerDepths() {
  // Go through every column and update the absolute depths. Verify correct
  // entries
  auto totalDepth = 0.0;
  auto range = rangeSpinBox->value();
  for (size_t row = 1; row < (size_t)layerModel->rowCount(); row++) {
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

  if ((totalDepth) > range) { // Only set if the range is too small.
    rangeSpinBox->setValue(totalDepth);
  }
}

// TODO: Calculate layer density based on the elemental composition of the
// layer.
double calculateLayerDensity(size_t row) { return row; }

void MainWindow::createLayerRow() {
  const size_t row = layerModel->rowCount();
  // const size_t columns = layerModel->columnCount();
  layerModel->insertRow(row);

  // TODO: Figure out how to do this dynamically (optional)
  const size_t startingDepthColumn = 0;
  const size_t depthSpinBoxColumn = 1;
  const size_t layerLineEditColumn = 2;
  const size_t densitySpinBoxColumn = 3;
  const size_t compositionViewColumn = 4;
  const size_t phaseDropDownColumn = 5;
  const size_t editButtonColumn = 6;
  const size_t removeButtonColumn = 7;

  const auto startingDepthLabelIndex =
      layerModel->index(row, startingDepthColumn);
  const auto depthSpinBoxIndex = layerModel->index(row, depthSpinBoxColumn);
  // const auto layerLineEditIndex = layerModel->index(row, layerColumn);
  const auto densitySpinBoxIndex = layerModel->index(row, densitySpinBoxColumn);
  const auto compositionViewIndex =
      layerModel->index(row, compositionViewColumn);
  const auto phaseDropDownIndex = layerModel->index(row, phaseDropDownColumn);
  const auto editButtonIndex = layerModel->index(row, editButtonColumn);
  const auto removeButtonIndex = layerModel->index(row, removeButtonColumn);

  // Depth field and label
  auto *startingDepthLabel = new QLabel(this);
  startingDepthLabel->setText(QString::number(row));
  // auto *depthField = new QDoubleSpinBox(this);
  // const auto minDepth = 0.01;
  // const auto maxDepth = 1e10;
  // depthField->setRange(minDepth, maxDepth);
  // depthField->setSingleStep(0.01);
  // connect(depthField, &QDoubleSpinBox::valueChanged, this,
  //         &MainWindow::updateLayerDepths);

  // auto *startingDepthLabel = new QStandardItem();
  // startingDepthLabel->setData(StartingDepthLabel, Qt::UserRole);
  // startingDepthLabel->setData(QString::number(0), Qt::DisplayRole);

  const auto minDepth = 0.01; // TODO: More elegant solution (optional)
  const auto maxDepth = std::numeric_limits<typeof(minDepth)>::max();
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
  const auto maxDensity = std::numeric_limits<typeof(minDensity)>::max();
  const auto defaultDensity =
      1.0; // TODO: Pick some proper density for a default material
  auto *densitySpinBox = new QDoubleSpinBox();
  densitySpinBox->setRange(minDensity, maxDensity);
  densitySpinBox->setValue(defaultDensity);

  // auto *phaseDropDown = new QStandardItem();
  // phaseDropDown->setData(PhaseDropDown, Qt::UserRole);

  // // Layer  LineEdit
  // auto *layerLineEdit = new QLineEdit(this);
  // layerLineEdit->setText(QString("Layer ") + QString::number(row + 1));

  // // Density QDoubleSpinBox
  // auto *densityField = new QDoubleSpinBox(this);
  // densityField->setRange(0.01, 1000.0);
  // densityField->setSingleStep(0.01);

  // Composition
  const auto compositionColumns = 5;
  auto *compositionModel = new QStandardItemModel(0, compositionColumns, this);
  auto *compositionTable = new QTableView(this);

  // Fit inside the inner column
  // compositionTable->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
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
  }); // TODO: Support input of weight or atomic percent

  compositionTable->setModel(compositionModel);

  // auto *compositionSpinBox = new QDoubleSpinBox(this);
  // compositionSpinBox->setRange(0.0, 100.0);
  // auto *compositionLabel = new QLabel("100");
  // compositionTable->setIndexWidget(compositionModel->index(0, 0),
  // new QLabel("H"));
  // compositionTable->setIndexWidget(compositionModel->index(0, 1),
  // compositionSpinBox);
  // compositionTable->setIndexWidget(compositionModel->index(0, 2),
  // compositionLabel);

  auto *compositionLayoutWidget = new QWidget(this);
  auto *compositionLayout = new QVBoxLayout();
  compositionLayoutWidget->setLayout(compositionLayout);

  auto compositionAddElementFunction = [compositionModel, compositionTable]() {
    const auto defaultElement = "H";
    const auto defaultFraction = 1.0;
    const auto elementColumn = 0;
    const auto fractionColumn = 1;
    const auto percentColumn = 2;
    const auto displacementEnergyColumn = 3;
    const auto elementDeleteButtonColumn = 4;
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
    auto *displacementEnergySpinBox =
        new QStandardItem(displacementEnergyMinValue);
    displacementEnergySpinBox->setData(
        QString::number(displacementEnergyDefaultValue), Qt::DisplayRole);
    displacementEnergySpinBox->setData(DisplacementEnergySpinBox, Qt::UserRole);
    percentLabel->setData(defaultFraction, Qt::DisplayRole);

    auto *elementDeleteButton = new QStandardItem("Dbl Clk");
    elementDeleteButton->setData(ElementDeleteButton, Qt::UserRole);

    compositionModel->insertRow(index);
    compositionModel->setItem(index, elementColumn, elementLabel);
    compositionModel->setItem(index, fractionColumn, fractionSpinBox);
    compositionModel->setItem(index, percentColumn, percentLabel);
    compositionModel->setItem(index, displacementEnergyColumn,
                              displacementEnergySpinBox);
    compositionModel->setItem(index, elementDeleteButtonColumn,
                              elementDeleteButton);
    compositionTable->resizeColumnsToContents();
    compositionTable->resizeRowsToContents();
  };
  compositionAddElementFunction(); // Start with one row

  auto *compositionAddElementButton = new QPushButton("Add Element");
  connect(compositionAddElementButton, &QPushButton::clicked, this,
          compositionAddElementFunction);

  compositionLayout->addWidget(compositionTable);
  compositionLayout->addWidget(compositionAddElementButton);

  for (int compositionRow = 0; compositionRow < compositionModel->rowCount();
       ++compositionRow) {
    QModelIndex index = compositionModel->index(compositionRow, 0);
    QWidget *editor = compositionDelegate->createEditor(
        compositionTable, QStyleOptionViewItem(), index);
    compositionDelegate->updateEditorGeometry(editor, QStyleOptionViewItem(),
                                              index);
  }

  // Phase drop down
  auto *phaseDropDown = new QComboBox(this);
  phaseDropDown->addItems(QStringList{"Solid", "Liquid", "Gas"});
  phaseDropDown->setCurrentIndex(0);

  // Setup edit button
  auto *editButton = new QPushButton(this);
  QPersistentModelIndex editButtonPersistentIndex(editButtonIndex);
  editButton->setIcon(QApplication::style()->standardIcon(QStyle::SP_FileIcon));
  editButton->setToolTip("Edit");
  connect(editButton, &QPushButton::clicked, this,
          [this, editButtonPersistentIndex]() {
            if (editButtonPersistentIndex.isValid()) {
              handleLayerButtonClicked(editButtonPersistentIndex.row());
            }
          });

  // // Setup remove button
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

  // Insert widgets
  // layerModel->setItem(row, startingDepthColumn, startingDepthLabel);
  // layerModel->setItem(row, depthSpinBoxColumn, depthSpinBox);
  layerModel->setItem(row, layerLineEditColumn, layerLineEdit);
  // layerModel->setItem(row, densitySpinBoxColumn, densitySpinBox);
  // layerModel->setItem(row, phaseDropDownColumn, phaseDropDown);
  // layerModel->setItem(row, removeButtonColumn, removeButton);
  // layerModel->setItem(row, editButtonColumn, editButton);

  layerTable->setIndexWidget(startingDepthLabelIndex, startingDepthLabel);
  layerTable->setIndexWidget(depthSpinBoxIndex, depthSpinBox);
  // layerTable->setIndexWidget(layerLineEditIndex, layerLineEdit);
  layerTable->setIndexWidget(densitySpinBoxIndex, densitySpinBox);
  // layerTable->setIndexWidget(compositionViewIndex, compositionTable);
  layerTable->setIndexWidget(compositionViewIndex, compositionLayoutWidget);
  layerTable->setIndexWidget(phaseDropDownIndex, phaseDropDown);
  layerTable->setIndexWidget(editButtonPersistentIndex, editButton);
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

  // TODO: Change to allow arbitrary units
  layerHeaders = QStringList{"Starting Depth\n(Angstrom)",
                             "Depth\n(Angstrom)",
                             "Layer",
                             QString("Density\n(kg/m%1)").arg(QChar(0x00B3)),
                             "Composition",
                             "Phase",
                             "Edit\n(WIP)",
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

void MainWindow::openLayerWindow(size_t row) {
  const size_t compositionColumn = 4;
  const auto compositionIndex = layerModel->index(row, compositionColumn);
  const size_t layerNameColumn = 2;
  // const auto layerNameIndex = layerModel->index(row, layerNameColumn);

  // The QTableView containing the element data is inside of the
  // Layer QTableView, in a QWidget with a QVboxLayout at index 0.
  // We need to do some pointer magic to obtain this pointer, before
  // finally getting the model we need to open the data in a separate window.
  auto *compositionViewWidget =
      qobject_cast<QWidget *>(layerTable->indexWidget(compositionIndex));
  auto *compositionViewLayout =
      compositionViewWidget->findChild<QVBoxLayout *>();
  auto *compositionView =
      qobject_cast<QTableView *>(compositionViewLayout->itemAt(0)->widget());
  auto *compositionModel =
      qobject_cast<QStandardItemModel *>(compositionView->model()); // Get model

  compositionView =
      new QTableView(this); // Reuse this pointer for a new TableView
  auto *compositionDelegate = new LayerWidgetDelegate(compositionView);
  compositionView->setItemDelegate(compositionDelegate);
  compositionView->setModel(compositionModel);
  auto layerName =
      layerModel->index(row, layerNameColumn).data(Qt::DisplayRole).toString();

  LayerWindow *layerWindow = new LayerWindow(compositionModel, layerName,
                                             this); // Pass 'this' as the parent
  layerWindow->setWindowFlags(Qt::Window);
  connect(layerWindow, &LayerWindow::okClicked, this, &MainWindow::onOkClicked);
  connect(layerWindow, &LayerWindow::cancelClicked, this,
          &MainWindow::onCancelClicked);
  layerWindow->resize(400, 400);
  layerWindow->show();
}

//! [2]

//! [3]
#ifndef QT_NO_CONTEXTMENU
void MainWindow::contextMenuEvent(QContextMenuEvent *event) {
  QMenu menu(this);
  menu.addAction(cutAct);
  menu.addAction(copyAct);
  menu.addAction(pasteAct);
  menu.exec(event->globalPos());
}
#endif // QT_NO_CONTEXTMENU
//! [3]

void MainWindow::newFile() { /*{{{*/
  infoLabel->setText(tr("Invoked <b>File|New</b>"));
}

void MainWindow::open() { infoLabel->setText(tr("Invoked <b>File|Open</b>")); }

void MainWindow::save() { infoLabel->setText(tr("Invoked <b>File|Save</b>")); }

void MainWindow::print() {
  infoLabel->setText(tr("Invoked <b>File|Print</b>"));
}

void MainWindow::undo() { infoLabel->setText(tr("Invoked <b>Edit|Undo</b>")); }

void MainWindow::redo() { infoLabel->setText(tr("Invoked <b>Edit|Redo</b>")); }

void MainWindow::cut() { infoLabel->setText(tr("Invoked <b>Edit|Cut</b>")); }

void MainWindow::copy() { infoLabel->setText(tr("Invoked <b>Edit|Copy</b>")); }

void MainWindow::paste() {
  infoLabel->setText(tr("Invoked <b>Edit|Paste</b>"));
}

void MainWindow::bold() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Bold</b>"));
}

void MainWindow::italic() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Italic</b>"));
}

void MainWindow::leftAlign() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Left Align</b>"));
}

void MainWindow::rightAlign() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Right Align</b>"));
}

void MainWindow::justify() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Justify</b>"));
}

void MainWindow::center() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Center</b>"));
}

void MainWindow::setLineSpacing() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Set Line Spacing</b>"));
}

void MainWindow::setParagraphSpacing() {
  infoLabel->setText(tr("Invoked <b>Edit|Format|Set Paragraph Spacing</b>"));
}

void MainWindow::about() {
  infoLabel->setText(tr("Invoked <b>Help|About</b>"));
  QMessageBox::about(this, tr("About Menu"),
                     tr("The <b>Menu</b> example shows how to create "
                        "menu-bar menus and context menus."));
}

void MainWindow::aboutQt() {
  infoLabel->setText(tr("Invoked <b>Help|About Qt</b>"));
} /*}}}*/

//! [4]
void MainWindow::createActions() { /*{{{*/
  //! [5]
  newAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::DocumentNew),
                       tr("&New"), this);
  newAct->setShortcuts(QKeySequence::New);
  newAct->setStatusTip(tr("Create a new file"));
  connect(newAct, &QAction::triggered, this, &MainWindow::newFile);
  //! [4]

  openAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::DocumentOpen),
                        tr("&Open..."), this);
  openAct->setShortcuts(QKeySequence::Open);
  openAct->setStatusTip(tr("Open an existing file"));
  connect(openAct, &QAction::triggered, this, &MainWindow::open);
  //! [5]

  saveAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::DocumentSave),
                        tr("&Save"), this);
  saveAct->setShortcuts(QKeySequence::Save);
  saveAct->setStatusTip(tr("Save the document to disk"));
  connect(saveAct, &QAction::triggered, this, &MainWindow::save);

  printAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::DocumentPrint),
                         tr("&Print..."), this);
  printAct->setShortcuts(QKeySequence::Print);
  printAct->setStatusTip(tr("Print the document"));
  connect(printAct, &QAction::triggered, this, &MainWindow::print);

  exitAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::ApplicationExit),
                        tr("E&xit"), this);
  exitAct->setShortcuts(QKeySequence::Quit);
  exitAct->setStatusTip(tr("Exit the application"));
  connect(exitAct, &QAction::triggered, this, &QWidget::close);

  undoAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditUndo),
                        tr("&Undo"), this);
  undoAct->setShortcuts(QKeySequence::Undo);
  undoAct->setStatusTip(tr("Undo the last operation"));
  connect(undoAct, &QAction::triggered, this, &MainWindow::undo);

  redoAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditRedo),
                        tr("&Redo"), this);
  redoAct->setShortcuts(QKeySequence::Redo);
  redoAct->setStatusTip(tr("Redo the last operation"));
  connect(redoAct, &QAction::triggered, this, &MainWindow::redo);

  cutAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditCut), tr("Cu&t"),
                       this);
  cutAct->setShortcuts(QKeySequence::Cut);
  cutAct->setStatusTip(tr("Cut the current selection's contents to the "
                          "clipboard"));
  connect(cutAct, &QAction::triggered, this, &MainWindow::cut);

  copyAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditCopy),
                        tr("&Copy"), this);
  copyAct->setShortcuts(QKeySequence::Copy);
  copyAct->setStatusTip(tr("Copy the current selection's contents to the "
                           "clipboard"));
  connect(copyAct, &QAction::triggered, this, &MainWindow::copy);

  pasteAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditPaste),
                         tr("&Paste"), this);
  pasteAct->setShortcuts(QKeySequence::Paste);
  pasteAct->setStatusTip(tr("Paste the clipboard's contents into the current "
                            "selection"));
  connect(pasteAct, &QAction::triggered, this, &MainWindow::paste);

  boldAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatTextBold),
                        tr("&Bold"), this);
  boldAct->setCheckable(true);
  boldAct->setShortcut(QKeySequence::Bold);
  boldAct->setStatusTip(tr("Make the text bold"));
  connect(boldAct, &QAction::triggered, this, &MainWindow::bold);

  QFont boldFont = boldAct->font();
  boldFont.setBold(true);
  boldAct->setFont(boldFont);

  italicAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatTextItalic),
                          tr("&Italic"), this);
  italicAct->setCheckable(true);
  italicAct->setShortcut(QKeySequence::Italic);
  italicAct->setStatusTip(tr("Make the text italic"));
  connect(italicAct, &QAction::triggered, this, &MainWindow::italic);

  QFont italicFont = italicAct->font();
  italicFont.setItalic(true);
  italicAct->setFont(italicFont);

  setLineSpacingAct = new QAction(tr("Set &Line Spacing..."), this);
  setLineSpacingAct->setStatusTip(tr("Change the gap between the lines of a "
                                     "paragraph"));
  connect(setLineSpacingAct, &QAction::triggered, this,
          &MainWindow::setLineSpacing);

  setParagraphSpacingAct = new QAction(tr("Set &Paragraph Spacing..."), this);
  setParagraphSpacingAct->setStatusTip(tr("Change the gap between paragraphs"));
  connect(setParagraphSpacingAct, &QAction::triggered, this,
          &MainWindow::setParagraphSpacing);

  aboutAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::HelpAbout),
                         tr("&About"), this);
  aboutAct->setStatusTip(tr("Show the application's About box"));
  connect(aboutAct, &QAction::triggered, this, &MainWindow::about);

  aboutQtAct = new QAction(tr("About &Qt"), this);
  aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
  connect(aboutQtAct, &QAction::triggered, qApp, &QApplication::aboutQt);
  connect(aboutQtAct, &QAction::triggered, this, &MainWindow::aboutQt);

  leftAlignAct =
      new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyLeft),
                  tr("&Left Align"), this);
  leftAlignAct->setCheckable(true);
  leftAlignAct->setShortcut(tr("Ctrl+L"));
  leftAlignAct->setStatusTip(tr("Left align the selected text"));
  connect(leftAlignAct, &QAction::triggered, this, &MainWindow::leftAlign);

  rightAlignAct =
      new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyRight),
                  tr("&Right Align"), this);
  rightAlignAct->setCheckable(true);
  rightAlignAct->setShortcut(tr("Ctrl+R"));
  rightAlignAct->setStatusTip(tr("Right align the selected text"));
  connect(rightAlignAct, &QAction::triggered, this, &MainWindow::rightAlign);

  justifyAct =
      new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyFill),
                  tr("&Justify"), this);
  justifyAct->setCheckable(true);
  justifyAct->setShortcut(tr("Ctrl+J"));
  justifyAct->setStatusTip(tr("Justify the selected text"));
  connect(justifyAct, &QAction::triggered, this, &MainWindow::justify);

  centerAct =
      new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyCenter),
                  tr("&Center"), this);
  centerAct->setCheckable(true);
  centerAct->setShortcut(tr("Ctrl+E"));
  centerAct->setStatusTip(tr("Center the selected text"));
  connect(centerAct, &QAction::triggered, this, &MainWindow::center);

  //! [6] //! [7]
  alignmentGroup = new QActionGroup(this);
  alignmentGroup->addAction(leftAlignAct);
  alignmentGroup->addAction(rightAlignAct);
  alignmentGroup->addAction(justifyAct);
  alignmentGroup->addAction(centerAct);
  leftAlignAct->setChecked(true);
  //! [6]
} /*}}}*/
//! [7]

//! [8]
void MainWindow::createMenus() { /*{{{*/
  //! [9] //! [10]
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newAct);
  //! [9]
  fileMenu->addAction(openAct);
  //! [10]
  fileMenu->addAction(saveAct);
  fileMenu->addAction(printAct);
  //! [11]
  fileMenu->addSeparator();
  //! [11]
  fileMenu->addAction(exitAct);

  editMenu = menuBar()->addMenu(tr("&Edit"));
  editMenu->addAction(undoAct);
  editMenu->addAction(redoAct);
  editMenu->addSeparator();
  editMenu->addAction(cutAct);
  editMenu->addAction(copyAct);
  editMenu->addAction(pasteAct);
  editMenu->addSeparator();

  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(aboutAct);
  helpMenu->addAction(aboutQtAct);
  //! [8]

  //! [12]
  formatMenu = editMenu->addMenu(tr("&Format"));
  formatMenu->addAction(boldAct);
  formatMenu->addAction(italicAct);
  formatMenu->addSeparator()->setText(tr("Alignment"));
  formatMenu->addAction(leftAlignAct);
  formatMenu->addAction(rightAlignAct);
  formatMenu->addAction(justifyAct);
  formatMenu->addAction(centerAct);
  formatMenu->addSeparator();
  formatMenu->addAction(setLineSpacingAct);
  formatMenu->addAction(setParagraphSpacingAct);
} /*}}}*/
//! [12]

} // namespace gui
} // namespace amcset
