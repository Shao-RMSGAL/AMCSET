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

#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QSplitter>
#include <QTabWidget>
#include <QTextBrowser>
#include <QVBoxLayout>

#include "layer_delegate.h"

#include "layer_delegate.h"

namespace amcset {
namespace gui {

void MainWindow::closeEvent(QCloseEvent *event) {
  saveState();
  QMainWindow::closeEvent(event);
}

MainWindow::~MainWindow() {
  saveState();

  if (simulationThread) {
    simulationThread->quit();
    simulationThread->wait();
    simulationThread->deleteLater();
    simulationThread = nullptr;
  }

  if (settingsDialog) {
    settingsDialog->reject();
    settingsDialog->deleteLater();
    settingsDialog = nullptr;
  }

  if (view3D) {
    view3D->deleteLater();
    view3D = nullptr;
  }
}

void MainWindow::saveState() {
  if (!app_) {
    return;
  }

  QVariantMap state;
  state["incidentElectron"] = incidentElectron ? incidentElectron->isChecked() : false;
  state["damageCascade"] = damageCascade ? damageCascade->isChecked() : true;
  state["incidentEnergy"] = incidentEnergy ? incidentEnergy->value() : 10.0;
  state["incidentElement"] = incidentElement ? incidentElement->text() : QString("H");
  state["incidentMass"] = incidentMass ? incidentMass->value() : 1;
  state["incidentCharge"] = incidentCharge ? incidentCharge->value() : 1;
  state["incidentAngle"] = incidentAngle ? incidentAngle->value() : 0.0;
  state["bombardmentCount"] = bombardmentCountSpinBox ? bombardmentCountSpinBox->value() : 1;
  state["threadCount"] = threadCountSpinBox ? threadCountSpinBox->value() : 1;
  state["range"] = rangeSpinBox ? rangeSpinBox->value() : 1000.0;

  QVariantList layers;
  if (layerModel && layerTable) {
    for (int row = 0; row < layerModel->rowCount(); ++row) {
      QVariantMap layer;
      const auto startingDepthIndex = layerModel->index(row, 0);
      const auto depthIndex = layerModel->index(row, 1);
      const auto nameIndex = layerModel->index(row, 2);
      const auto densityIndex = layerModel->index(row, 3);
      const auto compositionIndex = layerModel->index(row, 4);
      const auto phaseIndex = layerModel->index(row, 5);

      layer["startingDepth"] = qobject_cast<QLabel *>(layerTable->indexWidget(startingDepthIndex))
                                   ? qobject_cast<QLabel *>(layerTable->indexWidget(startingDepthIndex))->text().toDouble()
                                   : 0.0;
      layer["depth"] = qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(depthIndex))
                           ? qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(depthIndex))->value()
                           : 0.0;
      layer["name"] = nameIndex.data(Qt::DisplayRole).toString();
      layer["density"] = qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(densityIndex))
                            ? qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(densityIndex))->value()
                            : 0.0;
      layer["phase"] = phaseIndex.data(Qt::DisplayRole).toString();

      auto *compositionViewWidget = qobject_cast<QWidget *>(layerTable->indexWidget(compositionIndex));
      QVariantList composition;
      if (compositionViewWidget) {
        auto *layout = compositionViewWidget->layout();
        if (layout) {
          auto *table = qobject_cast<QTableView *>(layout->itemAt(0)->widget());
          if (table && table->model()) {
            auto *model = table->model();
            for (int compositionRow = 0; compositionRow < model->rowCount(); ++compositionRow) {
              QVariantMap element;
              element["element"] = model->index(compositionRow, 0).data(Qt::DisplayRole).toString();
              element["atomicProportion"] = model->index(compositionRow, 1).data(Qt::DisplayRole).toDouble();
              element["atomicPercent"] = model->index(compositionRow, 2).data(Qt::DisplayRole).toDouble();
              element["displacementEnergy"] = model->index(compositionRow, 3).data(Qt::DisplayRole).toDouble();
              composition.append(element);
            }
          }
        }
      }
      layer["composition"] = composition;
      layers.append(layer);
    }
  }
  state["layers"] = layers;

  if (outputDirectoryEdit) state["outputDirectory"] = outputDirectoryEdit->text();
  if (logFileEdit) state["logFile"] = logFileEdit->text();
  if (displacementOutputEdit) state["displacementOutput"] = displacementOutputEdit->text();
  if (electronicOutputEdit) state["electronicOutput"] = electronicOutputEdit->text();
  if (nuclearOutputEdit) state["nuclearOutput"] = nuclearOutputEdit->text();
  if (totalOutputEdit) state["totalOutput"] = totalOutputEdit->text();
  if (particleTreeEdit) state["particleTreeOutput"] = particleTreeEdit->text();
  if (enableLoggingCheck) state["enableLogging"] = enableLoggingCheck->isChecked();
  if (logSingleDisplacementCheck) state["logSingleDisplacement"] = logSingleDisplacementCheck->isChecked();
  if (logEndOfFlyingDistanceOnlyCheck) state["logEndOfFlyingDistanceOnly"] = logEndOfFlyingDistanceOnlyCheck->isChecked();
  if (logStoppingPointOnlyCheck) state["logStoppingPointOnly"] = logStoppingPointOnlyCheck->isChecked();
  if (angleDivisorSpinBox) state["angleDivisor"] = angleDivisorSpinBox->value();
  if (flyingDistanceSpinBox) state["flyingDistance"] = flyingDistanceSpinBox->value();
  if (electronScreeningParametersEdit) state["electronScreeningParameters"] = electronScreeningParametersEdit->text();
  if (mottScatteringParametersEdit) state["mottScatteringParameters"] = mottScatteringParametersEdit->text();
  if (screeningModelCombo) state["screeningModel"] = screeningModelCombo->currentIndex();

  app_->saveMainWindowState(state);
}

void MainWindow::restoreState() {
  if (!app_) {
    return;
  }

  const auto state = app_->loadMainWindowState();
  if (state.isEmpty()) {
    return;
  }

  if (incidentElectron) incidentElectron->setChecked(state.value("incidentElectron", false).toBool());
  if (damageCascade) damageCascade->setChecked(state.value("damageCascade", true).toBool());
  if (incidentEnergy) incidentEnergy->setValue(state.value("incidentEnergy", 10.0).toDouble());
  if (incidentElement) incidentElement->setText(state.value("incidentElement", "H").toString());
  if (incidentMass) incidentMass->setValue(state.value("incidentMass", 1).toInt());
  if (incidentCharge) incidentCharge->setValue(state.value("incidentCharge", 1).toInt());
  if (incidentAngle) incidentAngle->setValue(state.value("incidentAngle", 0.0).toDouble());
  if (bombardmentCountSpinBox) bombardmentCountSpinBox->setValue(state.value("bombardmentCount", 1).toInt());
  if (threadCountSpinBox) threadCountSpinBox->setValue(state.value("threadCount", 1).toInt());
  if (rangeSpinBox) rangeSpinBox->setValue(state.value("range", 1000.0).toDouble());

  if (outputDirectoryEdit) outputDirectoryEdit->setText(state.value("outputDirectory", "output").toString());
  if (logFileEdit) logFileEdit->setText(state.value("logFile", "amcset.log").toString());
  if (displacementOutputEdit) displacementOutputEdit->setText(state.value("displacementOutput", "displacement_histogram.csv").toString());
  if (electronicOutputEdit) electronicOutputEdit->setText(state.value("electronicOutput", "electronic_loss_histogram.csv").toString());
  if (nuclearOutputEdit) nuclearOutputEdit->setText(state.value("nuclearOutput", "nuclear_loss_histogram.csv").toString());
  if (totalOutputEdit) totalOutputEdit->setText(state.value("totalOutput", "total_loss_histogram.csv").toString());
  if (particleTreeEdit) particleTreeEdit->setText(state.value("particleTreeOutput", "particle_tree.json").toString());
  if (enableLoggingCheck) enableLoggingCheck->setChecked(state.value("enableLogging", false).toBool());
  if (logSingleDisplacementCheck) logSingleDisplacementCheck->setChecked(state.value("logSingleDisplacement", false).toBool());
  if (logEndOfFlyingDistanceOnlyCheck) logEndOfFlyingDistanceOnlyCheck->setChecked(state.value("logEndOfFlyingDistanceOnly", false).toBool());
  if (logStoppingPointOnlyCheck) logStoppingPointOnlyCheck->setChecked(state.value("logStoppingPointOnly", false).toBool());
  if (angleDivisorSpinBox) angleDivisorSpinBox->setValue(state.value("angleDivisor", 1000).toInt());
  if (flyingDistanceSpinBox) flyingDistanceSpinBox->setValue(state.value("flyingDistance", 1000).toInt());
  if (electronScreeningParametersEdit) electronScreeningParametersEdit->setText(state.value("electronScreeningParameters", "common/assets/electron_screening_parameters.csv").toString());
  if (mottScatteringParametersEdit) mottScatteringParametersEdit->setText(state.value("mottScatteringParameters", "common/assets/mott_scattering_parameters.csv").toString());
  if (screeningModelCombo) screeningModelCombo->setCurrentIndex(state.value("screeningModel", 0).toInt());

  const auto layers = state.value("layers", QVariantList()).toList();
  for (const auto &layerVariant : layers) {
    createLayerRow();
    const auto &layerMap = layerVariant.toMap();
    const int row = layerModel->rowCount() - 1;
    const auto startingDepthIndex = layerModel->index(row, 0);
    const auto depthIndex = layerModel->index(row, 1);
    const auto densityIndex = layerModel->index(row, 3);
    const auto phaseIndex = layerModel->index(row, 5);
    const auto layerNameIndex = layerModel->index(row, 2);

    auto *startingDepthLabel = qobject_cast<QLabel *>(layerTable->indexWidget(startingDepthIndex));
    auto *depthSpinBox = qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(depthIndex));
    auto *densitySpinBox = qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(densityIndex));
    auto *phaseCombo = qobject_cast<QComboBox *>(layerTable->indexWidget(phaseIndex));

    if (startingDepthLabel) startingDepthLabel->setText(QString::number(layerMap.value("startingDepth", 0.0).toDouble()));
    if (depthSpinBox) depthSpinBox->setValue(layerMap.value("depth", 0.0).toDouble());
    if (densitySpinBox) densitySpinBox->setValue(layerMap.value("density", 1000.0).toDouble());
    if (phaseCombo) phaseCombo->setCurrentText(layerMap.value("phase", "Solid").toString());
    if (layerNameIndex.isValid()) layerModel->setData(layerNameIndex, layerMap.value("name", QString("Layer %1").arg(row + 1)), Qt::DisplayRole);

    auto *compositionViewWidget = qobject_cast<QWidget *>(layerTable->indexWidget(layerModel->index(row, 4)));
    if (compositionViewWidget) {
      auto *layout = compositionViewWidget->layout();
      if (layout) {
        auto *table = qobject_cast<QTableView *>(layout->itemAt(0)->widget());
        if (table && table->model()) {
          auto *model = qobject_cast<QStandardItemModel *>(table->model());
          if (!model) {
            continue;
          }
          model->removeRows(0, model->rowCount());
          const auto compositionList = layerMap.value("composition", QVariantList()).toList();
          for (const auto &compositionVariant : compositionList) {
            const auto compositionMap = compositionVariant.toMap();
            const int newRow = model->rowCount();
            model->insertRow(newRow);
            model->setData(model->index(newRow, 0), compositionMap.value("element", "H"), Qt::DisplayRole);
            model->setData(model->index(newRow, 1), compositionMap.value("atomicProportion", 1.0), Qt::DisplayRole);
            model->setData(model->index(newRow, 2), compositionMap.value("atomicPercent", 100.0), Qt::DisplayRole);
            model->setData(model->index(newRow, 3), compositionMap.value("displacementEnergy", 40.0), Qt::DisplayRole);

            auto *deleteButtonItem = new QStandardItem();
            deleteButtonItem->setData(ElementDeleteButton, Qt::UserRole);
            deleteButtonItem->setText(QString());
            deleteButtonItem->setData(QApplication::style()->standardIcon(QStyle::SP_TitleBarCloseButton), Qt::DecorationRole);
            model->setItem(newRow, 4, deleteButtonItem);
          }
        }
      }
    }
  }
  if (layerModel && layerModel->rowCount() > 0) {
    updateLayerDepths();
  }
}

MainWindow::MainWindow(Application *app) : app_(app) {
  resize(1920, 1080);

  auto *topWidget = new QWidget(this);
  setCentralWidget(topWidget);

  infoLabel = new QLabel(tr("<i>Choose a menu option, or right-click to "
                            "invoke a context menu</i>"));
  infoLabel->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
  infoLabel->setAlignment(Qt::AlignCenter);

  auto *runCancelButtonLayoutWidget = new QWidget(this);
  auto *runCancelButtonLayout = new QHBoxLayout(runCancelButtonLayoutWidget);
  auto *runButton = new QPushButton("Run Simulation");
  runButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  connect(runButton, &QPushButton::clicked, this, &MainWindow::runSimulation);
  auto *cancelButton = new QPushButton("Cancel Run ");
  cancelButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);

  runCancelButtonLayout->addWidget(runButton);
  runCancelButtonLayout->addWidget(cancelButton);

  canvas = new Canvas;
  canvas->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
  canvasBox = new QGroupBox("2D View");
  auto *canvasVBox = new QVBoxLayout();
  canvasVBox->addWidget(canvas);
  canvasBox->setLayout(canvasVBox);

  view3D = new Widget3DView(this);
  view3DBox = new QGroupBox("3D View");
  view3DBox->setVisible(false);
  auto *view3DVBox = new QVBoxLayout();
  view3DVBox->addWidget(view3D);
  view3DBox->setLayout(view3DVBox);

  auto *topLayout = new QVBoxLayout();
  topLayout->setContentsMargins(5, 5, 5, 5);

  auto *formBox = new QGroupBox(this);
  auto *formLayout = new QVBoxLayout(formBox);

  formLayout->addWidget(infoLabel);
  auto layerSplitter = new QSplitter(Qt::Vertical);
  formLayout->addWidget(createIncidentParticleForms(formBox));
  layerSplitter->addWidget(createLayerForms());
  layerSplitter->addWidget(runCancelButtonLayoutWidget);
  formLayout->addWidget(layerSplitter);

  auto *hPane = new QSplitter(Qt::Horizontal);
  auto *visualizationContainer = new QWidget(this);
  auto *visualizationLayout = new QVBoxLayout(visualizationContainer);
  visualizationLayout->setContentsMargins(0, 0, 0, 0);

  auto *viewToggleContainer = new QWidget(this);
  auto *viewToggleLayout = new QHBoxLayout(viewToggleContainer);
  viewToggleLayout->setContentsMargins(8, 8, 8, 8);
  auto *viewToggleLabel = new QLabel("View mode", this);
  viewToggleLayout->addWidget(viewToggleLabel);

  viewModeToggle = new QCheckBox("3D", this);
  viewModeToggle->setChecked(false);
  viewModeToggle->setToolTip("Toggle between 2D and 3D views");
  viewModeToggle->setStyleSheet(
      "QCheckBox { spacing: 8px; }"
      "QCheckBox::indicator { width: 42px; height: 22px; border-radius: 11px; background: #cccccc; border: 1px solid #9a9a9a; }"
      "QCheckBox::indicator:checked { background: #3d7de8; border-color: #2d66c8; }");
  connect(viewModeToggle, &QCheckBox::toggled, this, [this](bool checked) {
    use2DView->setChecked(!checked);
    use3DView->setChecked(checked);
    update2DViewMode(!checked ? Qt::Checked : Qt::Unchecked);
    update3DViewMode(checked ? Qt::Checked : Qt::Unchecked);
  });
  viewToggleLayout->addWidget(viewModeToggle);
  visualizationLayout->addWidget(viewToggleContainer);

  auto *viewStack = new QWidget(this);
  auto *viewStackLayout = new QVBoxLayout(viewStack);
  viewStackLayout->setContentsMargins(0, 0, 0, 0);
  viewStackLayout->addWidget(canvasBox);
  viewStackLayout->addWidget(view3DBox);
  visualizationLayout->addWidget(viewStack);

  hPane->addWidget(formBox);
  hPane->addWidget(visualizationContainer);
  auto formRatio = 0.6;
  const auto formWidth = (int)(this->width() * formRatio);
  const auto canvasWidth = (int)(this->width() * (1 - formRatio));
  hPane->setSizes(QList<int>{formWidth, canvasWidth});

  mainTabs = new QTabWidget(this);
  mainTabs->addTab(formBox, "Setup");
  mainTabs->addTab(hPane, "Visualization");

  histogramTabs = new QTabWidget(this);
  histogramTabs->setTabPosition(QTabWidget::North);
  QTimer::singleShot(0, this, &MainWindow::populateResultCharts);
  mainTabs->addTab(histogramTabs, "Histograms");

  simulationSummaryView = new QTextBrowser(this);
  simulationSummaryView->setOpenExternalLinks(false);
  simulationSummaryView->setReadOnly(true);
  simulationSummaryView->setStyleSheet("QTextBrowser { background: #f8f8f8; }");
  mainTabs->addTab(simulationSummaryView, "Summary");

  topLayout->addWidget(mainTabs);
  topWidget->setLayout(topLayout);
  createActions();
  createMenus();
  createSettingsDialog();

  auto message = tr("A context menu is available by right-clicking");
  statusBar()->showMessage(message);
  setWindowTitle(tr("AMCSET"));
  setMinimumSize(480, 480);

  const auto savedState = app_ ? app_->loadMainWindowState() : QVariantMap();
  const auto savedLayers = savedState.value("layers", QVariantList()).toList();
  if (savedLayers.isEmpty()) {
    createLayerRow();
  }
  restoreState();
  show();
}

} // namespace gui
} // namespace amcset
