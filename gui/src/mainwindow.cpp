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
#include <QEvent>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QSplitter>
#include <QTabWidget>
#include <QTextBrowser>
#include <QVBoxLayout>

#include "layer_delegate.h"

namespace amcset {
namespace gui {

void MainWindow::closeEvent(QCloseEvent *event) {
  saveState();
  QMainWindow::closeEvent(event);
}

void MainWindow::changeEvent(QEvent *event) {
  QMainWindow::changeEvent(event);
  if (event->type() == QEvent::ApplicationPaletteChange ||
      event->type() == QEvent::PaletteChange) {
    updateStatusBarColors();
  }
}

void MainWindow::updateStatusBarColors() {
  if (!statusBar()) {
    return;
  }
  auto statusPalette = statusBar()->palette();
  const auto windowColor = qApp->palette().color(QPalette::Window);
  const auto foreground = windowColor.lightness() > 128
                              ? QColor(32, 32, 32)
                              : QColor(238, 238, 238);
  const auto foregroundCss = QStringLiteral("rgb(%1, %2, %3)")
                                 .arg(foreground.red())
                                 .arg(foreground.green())
                                 .arg(foreground.blue());
  statusBar()->setStyleSheet(
      QStringLiteral("QStatusBar, QStatusBar QLabel { color: %1; }")
          .arg(foregroundCss));
  statusPalette.setColor(QPalette::WindowText, foreground);
  statusPalette.setColor(QPalette::Text, foreground);
  statusPalette.setColor(QPalette::ButtonText, foreground);
  statusBar()->setPalette(statusPalette);
  for (auto *label : statusBar()->findChildren<QLabel *>()) {
    auto labelPalette = label->palette();
    labelPalette.setColor(QPalette::WindowText, foreground);
    labelPalette.setColor(QPalette::Text, foreground);
    label->setPalette(labelPalette);
  }
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
  if (!app_ || restoringState) {
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
  state["showXYProjection"] = xyProjectionCheck ? xyProjectionCheck->isChecked() : true;
  state["showYZProjection"] = yzProjectionCheck ? yzProjectionCheck->isChecked() : true;
  state["showXZProjection"] = xzProjectionCheck ? xzProjectionCheck->isChecked() : true;
  state["generationFilter"] = generationFilterCombo ? generationFilterCombo->currentIndex() : 3;
  state["customGenerations"] = customGenerationEdit ? customGenerationEdit->text() : QString();

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
  if (rangeOutputEdit) state["rangeOutput"] = rangeOutputEdit->text();
  if (perBombardmentOutputDirEdit) state["perBombardmentOutputDir"] = perBombardmentOutputDirEdit->text();
  if (savePerBombardmentCheck) state["savePerBombardment"] = savePerBombardmentCheck->isChecked();
  if (saveAllGraphsCheck) state["saveAllGraphs"] = saveAllGraphsCheck->isChecked();
  if (enableLoggingCheck) state["enableLogging"] = enableLoggingCheck->isChecked();
  if (angleDivisorSpinBox) state["angleDivisor"] = angleDivisorSpinBox->value();
  if (flyingDistanceSpinBox) state["flyingDistance"] = flyingDistanceSpinBox->value();
  if (electronScreeningParametersEdit) state["electronScreeningParameters"] = electronScreeningParametersEdit->text();
  if (mottScatteringParametersEdit) state["mottScatteringParameters"] = mottScatteringParametersEdit->text();
  if (screeningModelCombo) state["screeningModel"] = screeningModelCombo->currentIndex();
  if (displacementBinCountSpinBox) state["displacementBinCount"] = displacementBinCountSpinBox->value();
  if (electronicLossBinCountSpinBox) state["electronicLossBinCount"] = electronicLossBinCountSpinBox->value();
  if (nuclearLossBinCountSpinBox) state["nuclearLossBinCount"] = nuclearLossBinCountSpinBox->value();
  if (totalLossBinCountSpinBox) state["totalLossBinCount"] = totalLossBinCountSpinBox->value();
  if (rangeBinCountSpinBox) state["rangeBinCount"] = rangeBinCountSpinBox->value();

  app_->saveMainWindowState(state);
}

void MainWindow::restoreState() {
  if (!app_) {
    return;
  }
  restoringState = true;

  const auto state = app_->loadMainWindowState();
  if (state.isEmpty()) {
    restoringState = false;
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
  if (xyProjectionCheck) xyProjectionCheck->setChecked(state.value("showXYProjection", true).toBool());
  if (yzProjectionCheck) yzProjectionCheck->setChecked(state.value("showYZProjection", true).toBool());
  if (xzProjectionCheck) xzProjectionCheck->setChecked(state.value("showXZProjection", true).toBool());
  if (generationFilterCombo) generationFilterCombo->setCurrentIndex(state.value("generationFilter", 3).toInt());
  if (customGenerationEdit) customGenerationEdit->setText(state.value("customGenerations").toString());
  updateProjectionVisibility();
  updateGenerationFilter();

  if (outputDirectoryEdit) outputDirectoryEdit->setText(state.value("outputDirectory", "output").toString());
  if (logFileEdit) logFileEdit->setText(state.value("logFile", "amcset.log").toString());
  if (displacementOutputEdit) displacementOutputEdit->setText(state.value("displacementOutput", "displacement_histogram.csv").toString());
  if (electronicOutputEdit) electronicOutputEdit->setText(state.value("electronicOutput", "electronic_loss_histogram.csv").toString());
  if (nuclearOutputEdit) nuclearOutputEdit->setText(state.value("nuclearOutput", "nuclear_loss_histogram.csv").toString());
  if (totalOutputEdit) totalOutputEdit->setText(state.value("totalOutput", "total_loss_histogram.csv").toString());
  if (rangeOutputEdit) rangeOutputEdit->setText(state.value("rangeOutput", "range_histogram.csv").toString());
  if (perBombardmentOutputDirEdit) {
    auto perBombardmentDir =
        state.value("perBombardmentOutputDir", "bombardments/").toString();
    perBombardmentDir = QDir::fromNativeSeparators(perBombardmentDir);
    if (perBombardmentDir.startsWith(QStringLiteral("output/"))) {
      perBombardmentDir.remove(0, QStringLiteral("output/").size());
    }
    perBombardmentOutputDirEdit->setText(perBombardmentDir);
  }
  if (savePerBombardmentCheck) savePerBombardmentCheck->setChecked(state.value("savePerBombardment", false).toBool());
  if (saveAllGraphsCheck) saveAllGraphsCheck->setChecked(state.value("saveAllGraphs", false).toBool());
  if (enableLoggingCheck) enableLoggingCheck->setChecked(state.value("enableLogging", false).toBool());
  if (angleDivisorSpinBox) angleDivisorSpinBox->setValue(state.value("angleDivisor", 1000).toInt());
  if (flyingDistanceSpinBox) flyingDistanceSpinBox->setValue(state.value("flyingDistance", 1000).toInt());
  if (electronScreeningParametersEdit) electronScreeningParametersEdit->setText(state.value("electronScreeningParameters", "common/assets/electron_screening_parameters.csv").toString());
  if (mottScatteringParametersEdit) mottScatteringParametersEdit->setText(state.value("mottScatteringParameters", "common/assets/mott_scattering_parameters.csv").toString());
  if (screeningModelCombo) screeningModelCombo->setCurrentIndex(state.value("screeningModel", 0).toInt());
  if (displacementBinCountSpinBox) displacementBinCountSpinBox->setValue(state.value("displacementBinCount", 100).toInt());
  if (electronicLossBinCountSpinBox) electronicLossBinCountSpinBox->setValue(state.value("electronicLossBinCount", 100).toInt());
  if (nuclearLossBinCountSpinBox) nuclearLossBinCountSpinBox->setValue(state.value("nuclearLossBinCount", 100).toInt());
  if (totalLossBinCountSpinBox) totalLossBinCountSpinBox->setValue(state.value("totalLossBinCount", 100).toInt());
  if (rangeBinCountSpinBox) rangeBinCountSpinBox->setValue(state.value("rangeBinCount", 100).toInt());

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
  if (rangeSpinBox) rangeSpinBox->setValue(state.value("range", 1000.0).toDouble());
  restoringState = false;
}

MainWindow::MainWindow(Application *app) : app_(app) {
  resize(1920, 1080);

  auto *topWidget = new QWidget(this);
  setCentralWidget(topWidget);

  infoLabel = new QLabel(topWidget);
  infoLabel->hide();
  infoLabel->setFrameStyle(QFrame::StyledPanel | QFrame::Sunken);
  infoLabel->setAlignment(Qt::AlignCenter);

  auto *runCancelButtonLayoutWidget = new QWidget(this);
  auto *runCancelButtonLayout = new QHBoxLayout(runCancelButtonLayoutWidget);
  runSimulationButton = new QPushButton("Run Simulation");
  runSimulationButton->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
  connect(runSimulationButton, &QPushButton::clicked, this, &MainWindow::runSimulation);
  cancelSimulationButton = new QPushButton("Cancel Run");
  cancelSimulationButton->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
  cancelSimulationButton->setEnabled(false);
  connect(cancelSimulationButton, &QPushButton::clicked, this, [this]() {
    if (!activeSimulation) {
      return;
    }
    activeSimulation->request_cancel();
    cancelSimulationButton->setEnabled(false);
    statusBar()->showMessage(tr("Cancelling after active bombardments finish..."));
  });

  runCancelButtonLayout->addWidget(runSimulationButton);
  runCancelButtonLayout->addWidget(cancelSimulationButton);
  const auto controlButtonHeight =
      std::max(runSimulationButton->sizeHint().height(),
               cancelSimulationButton->sizeHint().height()) + 8;
  runSimulationButton->setFixedHeight(controlButtonHeight);
  cancelSimulationButton->setFixedHeight(controlButtonHeight);

  xyCanvas = new Canvas;
  xyCanvas->setProjection(Canvas::Projection::XY);
  yzCanvas = new Canvas;
  yzCanvas->setProjection(Canvas::Projection::YZ);
  xzCanvas = new Canvas;
  xzCanvas->setProjection(Canvas::Projection::XZ);
  canvas = xzCanvas;
  for (auto *view : {xyCanvas, yzCanvas, xzCanvas}) {
    view->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);
  }
  canvasBox = new QGroupBox("2D View");
  projectionGrid = new QGridLayout();
  canvasBox->setLayout(projectionGrid);

  auto *projectionSidebar = new QWidget(this);
  auto *projectionLayout = new QVBoxLayout(projectionSidebar);
  projectionLayout->addWidget(new QLabel(tr("2D projections"), projectionSidebar));
  xyProjectionCheck = new QCheckBox(tr("XY"), projectionSidebar);
  yzProjectionCheck = new QCheckBox(tr("YZ"), projectionSidebar);
  xzProjectionCheck = new QCheckBox(tr("XZ"), projectionSidebar);
  xyProjectionCheck->setChecked(true);
  yzProjectionCheck->setChecked(true);
  xzProjectionCheck->setChecked(true);
  projectionLayout->addWidget(xyProjectionCheck);
  projectionLayout->addWidget(yzProjectionCheck);
  projectionLayout->addWidget(xzProjectionCheck);
    generationFilterCombo = new QComboBox(projectionSidebar);
    generationFilterCombo->addItems(
        {tr("Incident particle only"), tr("Primary knock-ons only"),
         tr("Incident and primary knock-ons"), tr("All"), tr("Custom")});
    generationFilterCombo->setCurrentIndex(3);
    projectionLayout->addWidget(new QLabel(tr("Generations"), projectionSidebar));
    projectionLayout->addWidget(generationFilterCombo);
    customGenerationEdit = new QLineEdit(projectionSidebar);
    customGenerationEdit->setPlaceholderText(tr("e.g. 0-3,5,6"));
    customGenerationEdit->hide();
    projectionLayout->addWidget(customGenerationEdit);
    connect(generationFilterCombo, qOverload<int>(&QComboBox::currentIndexChanged),
      this, [this](int index) {
        customGenerationEdit->setVisible(index == 4);
        updateGenerationFilter();
      });
    connect(customGenerationEdit, &QLineEdit::textChanged, this,
      &MainWindow::updateGenerationFilter);
      connect(customGenerationEdit, &QLineEdit::textChanged, this,
        [this](const QString &) { saveState(); });
  projectionLayout->addStretch();
  for (auto *check : {xyProjectionCheck, yzProjectionCheck, xzProjectionCheck}) {
    connect(check, &QCheckBox::toggled, this,
            &MainWindow::updateProjectionVisibility);
      connect(check, &QCheckBox::toggled, this, &MainWindow::saveState);
  }
    connect(generationFilterCombo, qOverload<int>(&QComboBox::currentIndexChanged),
      this, &MainWindow::saveState);
    connect(customGenerationEdit, &QLineEdit::editingFinished, this,
      &MainWindow::saveState);

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

  auto layerSplitter = new QSplitter(Qt::Vertical);
  formLayout->addWidget(createIncidentParticleForms(formBox));
  layerSplitter->addWidget(createLayerForms());
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
  auto *viewStackLayout = new QHBoxLayout(viewStack);
  viewStackLayout->setContentsMargins(0, 0, 0, 0);
  viewStackLayout->addWidget(projectionSidebar);
  auto *viewPages = new QWidget(viewStack);
  auto *viewPagesLayout = new QVBoxLayout(viewPages);
  viewPagesLayout->setContentsMargins(0, 0, 0, 0);
  viewPagesLayout->addWidget(canvasBox);
  viewPagesLayout->addWidget(view3DBox);
  viewStackLayout->addWidget(viewPages, 1);
  visualizationLayout->addWidget(viewStack);
  updateProjectionVisibility();

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
  mainTabs->addTab(simulationSummaryView, "Summary");

  auto *statisticsBanner = new QFrame(this);
  statisticsBanner->setFrameShape(QFrame::NoFrame);
  auto *statisticsLayout = new QHBoxLayout(statisticsBanner);
  statisticsLayout->setContentsMargins(8, 4, 8, 4);
  auto addStatistic = [statisticsBanner, statisticsLayout](const QString &title,
                                                           QLabel *value = nullptr) {
    auto *tile = new QFrame(statisticsBanner);
    tile->setFrameShape(QFrame::StyledPanel);
    tile->setStyleSheet(QStringLiteral("QFrame { border: 1px solid palette(mid); border-radius: 4px; }"));
    auto *tileLayout = new QVBoxLayout(tile);
    tileLayout->setContentsMargins(8, 3, 8, 3);
    auto *caption = new QLabel(title, tile);
    caption->setForegroundRole(QPalette::WindowText);
    caption->setStyleSheet(QStringLiteral("border: 0; font-size: 9pt;"));
    tileLayout->addWidget(caption);
    if (!value) value = new QLabel(QStringLiteral("0"), tile);
    value->setParent(tile);
    value->setForegroundRole(QPalette::Text);
    value->setStyleSheet(QStringLiteral("border: 0; font-weight: 600;"));
    tileLayout->addWidget(value);
    statisticsLayout->addWidget(tile);
    return value;
  };
  addStatistic(tr("Bombardments"), counterLabel);
  for (const auto &title : {tr("Displacements"), tr("Sputtered"),
                            tr("Mean range (Å)"), tr("Range straggling (Å)"),
                            tr("Mean loss (eV)"), tr("Loss straggling (eV)")}) {
    liveStatisticValues.append(addStatistic(title));
  }

  auto *controlBanner = new QFrame(this);
  controlBanner->setFrameShape(QFrame::NoFrame);
  auto *controlLayout = new QHBoxLayout(controlBanner);
  controlLayout->setContentsMargins(8, 4, 8, 4);
  controlLayout->addStretch();
  controlLayout->addWidget(runCancelButtonLayoutWidget);
  controlLayout->addStretch();

  topLayout->insertWidget(0, statisticsBanner);
  topLayout->addWidget(mainTabs);
  topLayout->addWidget(controlBanner);
  topWidget->setLayout(topLayout);
  updateStatusBarColors();
  createActions();
  createMenus();
  createSettingsDialog();

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
