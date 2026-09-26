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

#include <QSet>

#include <QDir>
#include <QMessageBox>
#include <QThread>
#include <QTimer>
#include <glog/logging.h>

namespace amcset {
namespace gui {

QString MainWindow::formatSimulationError(const std::exception &e,
                                          std::size_t depth) {
  std::ostringstream oss;
  oss << std::string(depth * 2, ' ') << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception &nested) {
    oss << formatSimulationError(nested, depth + 1).toStdString();
  } catch (...) {
  }
  return QString::fromStdString(oss.str());
}

void MainWindow::handleSimulationError(const QString &errorMessage) {
  activeSimulation = nullptr;
  if (runSimulationButton) runSimulationButton->setEnabled(true);
  if (cancelSimulationButton) cancelSimulationButton->setEnabled(false);
  if (bombardmentPollTimer) {
    bombardmentPollTimer->stop();
  }
  if (simulationThread) {
    simulationThread->deleteLater();
    simulationThread = nullptr;
  }

  QMessageBox::critical(this, "Simulation error", errorMessage,
                       QMessageBox::Ok);
  openSettingsDialog();
}

void MainWindow::handleSimulationResult(int result) {
  simulationResult = result;
  simulationCancelled = activeSimulation && activeSimulation->is_cancelled();
  simulationCompletedBombardments =
      activeSimulation ? activeSimulation->get_completed_bombardment_count() : 0;
  if (activeSimulation) {
    updateLiveStatistics(activeSimulation->snapshot_tallies(),
                         simulationCompletedBombardments);
  }
  if (counter) {
    counter->setValue(static_cast<int>(simulationCompletedBombardments));
  }
  if (counterLabel) {
    counterLabel->setText(
        tr("Bombardment %1 out of %2")
            .arg(simulationCompletedBombardments)
            .arg(bombardmentCountSpinBox ? bombardmentCountSpinBox->value() : 0));
  }
  if (cancelSimulationButton) cancelSimulationButton->setEnabled(false);
  if (simulationThread) {
    simulationThread->wait();
  }
  if (activeSimulation) {
    activeSimulation->clear_bombardment_results();
  }
  finishSimulation();
}

void MainWindow::finishSimulation() {
  activeSimulation = nullptr;
  if (runSimulationButton) runSimulationButton->setEnabled(true);
  if (bombardmentPollTimer) bombardmentPollTimer->stop();
  if (simulationThread) {
    simulationThread->deleteLater();
    simulationThread = nullptr;
  }
  if (simulationResult != 0) {
    handleSimulationError(QString("Simulation exited with code %1.")
                             .arg(simulationResult));
    return;
  }

  refreshResultsPane();
  if (simulationCancelled) {
    statusBar()->showMessage(
        tr("Simulation cancelled after %1 bombardments.")
            .arg(simulationCompletedBombardments),
        10000);
    updateResultPane(
        tr("Simulation cancelled after %1 completed bombardments; partial results were saved.")
            .arg(simulationCompletedBombardments));
  } else {
    statusBar()->showMessage(tr("Simulation complete."), 5000);
  }
}

void MainWindow::launchSimulationThread(
    std::unique_ptr<common::Simulation> simulation) {
  simulationResult = 0;
  simulationCancelled = false;
  simulationCompletedBombardments = 0;
  lastLiveTallyBombardments = 0;
  simulationThread = new SimulationThread(std::move(simulation));
  connect(simulationThread, &SimulationThread::resultReady, this,
          &MainWindow::handleSimulationResult);
  connect(simulationThread, &SimulationThread::errorReady, this,
          &MainWindow::handleSimulationError);
    if (runSimulationButton) runSimulationButton->setEnabled(false);
    if (cancelSimulationButton) cancelSimulationButton->setEnabled(true);
  simulationThread->start();
}

void MainWindow::runSimulation() {
  if (enableLoggingCheck && enableLoggingCheck->isChecked()) {
    const auto outputDir = outputDirectoryEdit->text().trimmed();
    const auto logFile = QDir(outputDir).absoluteFilePath(logFileEdit->text().trimmed());
    QFileInfo logInfo(logFile);
    if (!logInfo.absoluteDir().exists()) {
      QDir().mkpath(logInfo.absoluteDir().absolutePath());
    }
    const auto stdPath = logFile.toStdString();
    google::SetLogDestination(google::GLOG_INFO, stdPath.c_str());
    google::SetLogDestination(google::GLOG_WARNING, stdPath.c_str());
    google::SetLogDestination(google::GLOG_ERROR, stdPath.c_str());
    google::SetLogDestination(google::GLOG_FATAL, stdPath.c_str());
    updateResultPane(QString("Logging enabled to %1").arg(logFile));
  } else {
    google::SetLogDestination(google::GLOG_INFO, "");
    google::SetLogDestination(google::GLOG_WARNING, "");
    google::SetLogDestination(google::GLOG_ERROR, "");
    google::SetLogDestination(google::GLOG_FATAL, "");
  }

  const auto outputDir = outputDirectoryEdit->text().trimmed();
  if (!outputDir.isEmpty()) {
    QDir().mkpath(outputDir);
  }

  counter->reset();
  if (counterLabel) {
    counterLabel->setText(tr("Bombardment %1 out of %2")
                              .arg(0)
                              .arg(bombardmentCountSpinBox->value()));
  }
  updateLiveStatistics(common::Tallies{}, 0);

  if (use2DView->isChecked()) {
    for (auto *view : {xyCanvas, yzCanvas, xzCanvas}) {
      if (view) {
        view->clear();
        view->setRangeAngstrom(rangeSpinBox->value());
      }
    }
  }
  if (use3DView->isChecked()) {
    view3D->clearScene();
  }

  const auto displacementOutput = buildOutputPath(displacementOutputEdit->text().trimmed());
  const auto electronicOutput = buildOutputPath(electronicOutputEdit->text().trimmed());
  const auto nuclearOutput = buildOutputPath(nuclearOutputEdit->text().trimmed());
  const auto totalOutput = buildOutputPath(totalOutputEdit->text().trimmed());
  const auto rangeOutput =
      buildOutputPath(rangeOutputEdit ? rangeOutputEdit->text().trimmed()
                                      : QStringLiteral("range_histogram.csv"));
  const auto perBombardmentOutputDir = buildOutputPath(
      perBombardmentOutputDirEdit ? perBombardmentOutputDirEdit->text().trimmed()
                                  : QStringLiteral("bombardments/"))
                                       .toStdString();

  common::energy_quantity electron_stopping_energy = 40.0 * common::electron_volt;
  common::energy_quantity ion_stopping_energy = 40.0 * common::electron_volt;
  common::energy_quantity ion_displacement_energy = 40.0 * common::electron_volt;
  const size_t divisor_angle_number =
      angleDivisorSpinBox ? static_cast<size_t>(angleDivisorSpinBox->value()) : 1000;
  const size_t flying_distance_number =
      flyingDistanceSpinBox ? static_cast<size_t>(flyingDistanceSpinBox->value()) : 1000;

  const std::string electron_screening_path =
      electronScreeningParametersEdit ? electronScreeningParametersEdit->text().trimmed().toStdString()
                                      : "common/assets/electron_screening_parameters.csv";
  const std::string mott_path =
      mottScatteringParametersEdit ? mottScatteringParametersEdit->text().trimmed().toStdString()
                                  : "common/assets/mott_scattering_parameters.csv";

  bool is_electron = (incidentElectron->checkState() == Qt::Checked);
  bool enable_damage_cascade = (damageCascade->checkState() == Qt::Checked);
  size_t mass_number = incidentMass->value();
  size_t z_number = incidentCharge->value();
  common::angle_quantity incident_angle = common::degree * incidentAngle->value();
  common::length_quantity range = rangeSpinBox->value() * common::angstrom;
  size_t bombardment_count = static_cast<size_t>(bombardmentCountSpinBox->value());
  common::energy_quantity incident_energy = incidentEnergy->value() * common::kilo_electron_volt;
  size_t thread_count = threadCountSpinBox->value();

  common::Simulation::Settings::ScreeningModel selected_model =
      (screeningModelCombo && screeningModelCombo->currentIndex() == 1)
          ? common::Simulation::Settings::ScreeningModel::QuickBasic
          : common::Simulation::Settings::ScreeningModel::ThomasFermi;

  auto settings = common::Simulation::Settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle,
      displacementBinCountSpinBox ? static_cast<size_t>(displacementBinCountSpinBox->value()) : 100,
      electronicLossBinCountSpinBox ? static_cast<size_t>(electronicLossBinCountSpinBox->value()) : 100,
      nuclearLossBinCountSpinBox ? static_cast<size_t>(nuclearLossBinCountSpinBox->value()) : 100,
      totalLossBinCountSpinBox ? static_cast<size_t>(totalLossBinCountSpinBox->value()) : 100,
      rangeBinCountSpinBox ? static_cast<size_t>(rangeBinCountSpinBox->value()) : 100,
      displacementOutput.toStdString(), electronicOutput.toStdString(),
      nuclearOutput.toStdString(), totalOutput.toStdString(),
      savePerBombardmentCheck && savePerBombardmentCheck->isChecked(),
      perBombardmentOutputDir,
      electron_screening_path, mott_path, selected_model,
      rangeOutput.toStdString());

  // Perform safety confirmations BEFORE constructing or launching the simulation.
  if (savePerBombardmentCheck && savePerBombardmentCheck->isChecked()) {
    // Warn only when bombardment count is large
    if (bombardment_count > 100 && (!app_ || app_->showLargeCsvWarning())) {
      QMessageBox warning(QMessageBox::Warning, "Large Output Warning",
        "Saving per-bombardment CSVs can generate very large files. It is not recommended for more than 100 bombardments. Continue?",
        QMessageBox::Ok | QMessageBox::Cancel, this);
      auto *dontTellAgain = new QCheckBox(tr("Do not tell me again"), &warning);
      warning.setCheckBox(dontTellAgain);
      const auto ret = warning.exec();
      if (app_ && dontTellAgain->isChecked()) {
        app_->setShowLargeCsvWarning(false);
        if (showLargeCsvWarningCheck) showLargeCsvWarningCheck->setChecked(false);
      }
      if (ret == QMessageBox::Cancel) return;
    }

    // Check whether the per-bombardment output directory already contains files.
    const QString perDir = QString::fromStdString(perBombardmentOutputDir);
    if (!perDir.trimmed().isEmpty()) {
      const QDir dir(perDir);
      if (dir.exists()) {
        const auto entries = dir.entryList(QDir::AllEntries | QDir::NoDotAndDotDot);
        if (!entries.isEmpty() &&
            (!app_ || app_->showOutputOverwriteWarning())) {
          QMessageBox warning(QMessageBox::Warning, "Output Directory Not Empty",
            QStringLiteral("The per-bombardment output directory '%1' is not empty. Running another simulation may overwrite existing files. Continue?").arg(perDir),
            QMessageBox::Ok | QMessageBox::Cancel, this);
          auto *dontTellAgain = new QCheckBox(tr("Do not tell me again"), &warning);
          warning.setCheckBox(dontTellAgain);
          const auto overwriteRet = warning.exec();
          if (app_ && dontTellAgain->isChecked()) {
            app_->setShowOutputOverwriteWarning(false);
            if (showOutputOverwriteWarningCheck) {
              showOutputOverwriteWarningCheck->setChecked(false);
            }
          }
          if (overwriteRet == QMessageBox::Cancel) return;
        }
      }
    }
  }

  try {
    const int layerCount = layerModel->rowCount();
    if (layerCount == 0) {
      throw std::logic_error("Layer count cannot be zero.");
    }

    std::vector<common::Layer> layers;
    for (int layerRow = 0; layerRow < layerCount; layerRow++) {
      const auto startingDepthColumn = 0;
      const auto depthColumn = 1;
      const auto densityColumn = 3;
      const auto compositionColumn = 4;

      auto startingDepthIndex = layerModel->index(layerRow, startingDepthColumn);
      auto depthIndex = layerModel->index(layerRow, depthColumn);
      auto startingDepth =
          qobject_cast<QLabel *>(layerTable->indexWidget(startingDepthIndex))
              ->text()
              .toDouble();
      auto depth = (startingDepth + qobject_cast<QDoubleSpinBox *>(
                                        layerTable->indexWidget(depthIndex))
                                        ->value()) *
                   common::angstrom;

      auto densityIndex = layerModel->index(layerRow, densityColumn);
      auto density =
          qobject_cast<QDoubleSpinBox *>(layerTable->indexWidget(densityIndex))
              ->value() *
          common::kg_per_cubic_meter;

      auto materialIndex = layerModel->index(layerRow, compositionColumn);
      auto *materialWidget =
          qobject_cast<QWidget *>(layerTable->indexWidget(materialIndex));
      auto *materialLayout = materialWidget->layout();
      auto *materialTable =
          qobject_cast<QTableView *>(materialLayout->itemAt(0)->widget());
      auto *materialModel = materialTable->model();
      const int elementCount = materialModel->rowCount();
      common::Layer::material_vector material;
      auto proportionSum = 0.0;
      for (int elementRow = 0; elementRow < elementCount; elementRow++) {
        const auto atomicSymbolColumn = 0;
        const auto proportionColumn = 1;

        auto atomicSymbol = materialModel->index(elementRow, atomicSymbolColumn)
                                .data(Qt::DisplayRole)
                                .toString();
        auto atomicProportion = materialModel->index(elementRow, proportionColumn)
                                    .data(Qt::DisplayRole)
                                    .toDouble();
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
        }
      }

      for (auto &entry : material) {
        entry.first /= proportionSum;
      }
      layers.push_back(common::Layer(std::move(material), depth, density));
    }

    auto volume = common::Volume(std::move(layers));
    auto simulation = std::make_unique<common::Simulation>(
        settings, std::move(volume),
        [](const std::vector<common::Coordinate> &, size_t) {});

    activeSimulation = simulation.get();
    launchSimulationThread(std::move(simulation));
    if (!bombardmentPollTimer) {
      bombardmentPollTimer = new QTimer(this);
      connect(bombardmentPollTimer, &QTimer::timeout, this,
              &MainWindow::pollBombardmentQueue);
    }
    bombardmentPollTimer->start(1);
  } catch (const std::exception &e) {
    common::print_exception(e);
    handleSimulationError(formatSimulationError(e));
  }
}

void MainWindow::update2DViewMode(Qt::CheckState state) {
  if (state == Qt::Checked) {
    canvasBox->setVisible(true);
  } else {
    canvasBox->setVisible(false);
  }
}

void MainWindow::update3DViewMode(Qt::CheckState state) {
  if (state == Qt::Checked) {
    view3DBox->setVisible(true);
  } else {
    view3DBox->setVisible(false);
  }
}

void MainWindow::updateProjectionVisibility() {
  if (!projectionGrid) {
    return;
  }

  while (auto *item = projectionGrid->takeAt(0)) {
    delete item;
  }

  const bool showXY = xyCanvas && xyProjectionCheck && xyProjectionCheck->isChecked();
  const bool showYZ = yzCanvas && yzProjectionCheck && yzProjectionCheck->isChecked();
  const bool showXZ = xzCanvas && xzProjectionCheck && xzProjectionCheck->isChecked();
  if (xyCanvas) xyCanvas->setVisible(showXY);
  if (yzCanvas) yzCanvas->setVisible(showYZ);
  if (xzCanvas) xzCanvas->setVisible(showXZ);
  const int visibleCount = static_cast<int>(showXY) + static_cast<int>(showYZ) +
                           static_cast<int>(showXZ);

  projectionGrid->setRowStretch(0, 0);
  projectionGrid->setRowStretch(1, 0);
  projectionGrid->setColumnStretch(0, 0);
  projectionGrid->setColumnStretch(1, 0);

  if (visibleCount == 1) {
    auto *visibleView = showXY ? xyCanvas : showYZ ? yzCanvas : xzCanvas;
    projectionGrid->addWidget(visibleView, 0, 0);
    projectionGrid->setRowStretch(0, 1);
    projectionGrid->setColumnStretch(0, 1);
  } else if (visibleCount == 2) {
    if (showXY) {
      int column = 0;
      for (auto *view : {xzCanvas, yzCanvas, xyCanvas}) {
        if (view == xyCanvas ? showXY : view == yzCanvas ? showYZ : showXZ) {
          projectionGrid->addWidget(view, 0, column++);
        }
      }
      projectionGrid->setRowStretch(0, 1);
      projectionGrid->setColumnStretch(0, 1);
      projectionGrid->setColumnStretch(1, 1);
    } else {
      projectionGrid->addWidget(xzCanvas, 0, 0);
      projectionGrid->addWidget(yzCanvas, 1, 0);
      projectionGrid->setRowStretch(0, 1);
      projectionGrid->setRowStretch(1, 1);
      projectionGrid->setColumnStretch(0, 1);
    }
  } else if (visibleCount == 3) {
    projectionGrid->addWidget(xzCanvas, 0, 0);
    projectionGrid->addWidget(yzCanvas, 1, 0);
    projectionGrid->addWidget(xyCanvas, 0, 1, 2, 1);
    projectionGrid->setRowStretch(0, 1);
    projectionGrid->setRowStretch(1, 1);
    projectionGrid->setColumnStretch(0, 1);
    projectionGrid->setColumnStretch(1, 1);
  }
}

void MainWindow::updateGenerationFilter() {
  if (!generationFilterCombo) {
    return;
  }

  const auto index = generationFilterCombo->currentIndex();
  QSet<std::size_t> generations;
  bool showAllGenerations = false;
  if (index == 0) {
    generations.insert(0);
  } else if (index == 1) {
    generations.insert(1);
  } else if (index == 2) {
    generations.insert(0);
    generations.insert(1);
  } else if (index == 3) {
    showAllGenerations = true;
  } else if (customGenerationEdit) {
    for (const auto &part : customGenerationEdit->text().split(',', Qt::SkipEmptyParts)) {
      const auto bounds = part.trimmed().split('-', Qt::SkipEmptyParts);
      bool firstOk = false;
      const auto first = bounds.value(0).toULongLong(&firstOk);
      if (!firstOk) {
        continue;
      }
      if (bounds.size() == 1) {
        generations.insert(static_cast<std::size_t>(first));
        continue;
      }
      bool secondOk = false;
      const auto second = bounds.value(1).toULongLong(&secondOk);
      if (!secondOk || second < first) {
        continue;
      }
      for (auto generation = first; generation <= second; ++generation) {
        generations.insert(static_cast<std::size_t>(generation));
      }
    }
  }

  for (auto *view : {xyCanvas, yzCanvas, xzCanvas}) {
    if (view) {
      view->setGenerationFilter(generations, showAllGenerations);
    }
  }
}

void MainWindow::pollBombardmentQueue() {
  if (!activeSimulation) return;

  const auto viewPtr = QPointer<Widget3DView>(view3D);
  const auto counterPtr = QPointer<Counter>(counter);

  const auto range = activeSimulation->get_settings(&common::Simulation::Settings::range_);
  QVector<Canvas::Point> canvasPoints;
  std::vector<common::CollisionPoint> vec;
  if (activeSimulation->try_pop_bombardment_results(vec)) {
    if (use2DView && use2DView->checkState() == Qt::Checked) {
      for (const auto &point : vec) {
        canvasPoints.append({
            (point.coordinate.x_ / common::angstrom).value(),
          (point.coordinate.y_ / common::angstrom).value(),
          (point.coordinate.z_ / common::angstrom).value(),
            point.generation});
      }
    }

    if (viewPtr && use3DView && use3DView->checkState() == Qt::Checked) {
      for (size_t i = 0; i + 1 < vec.size(); ++i) {
        const auto &p1 = vec[i].coordinate;
        const auto &p2 = vec[i + 1].coordinate;
        QVector3D point1((p1.x_ / common::angstrom).value(),
                         (p1.y_ / common::angstrom).value(),
                         (p1.z_ / common::angstrom).value());
        QVector3D point2((p2.x_ / common::angstrom).value(),
                         (p2.y_ / common::angstrom).value(),
                         (p2.z_ / common::angstrom).value());
        viewPtr->requestDrawLine(point1, point2, Qt::green);
      }
    }
  }

  if (!vec.empty() && use2DView && use2DView->checkState() == Qt::Checked) {
    const auto rangeAngstrom = (range / common::angstrom).value();
    for (auto *view : {xyCanvas, yzCanvas, xzCanvas}) {
      if (view) {
        view->setRangeAngstrom(rangeAngstrom);
        view->addPoints(canvasPoints);
      }
    }
  }

  if (counterPtr && activeSimulation != nullptr) {
    const auto completed = activeSimulation->get_completed_bombardment_count();
    if (static_cast<std::size_t>(counterPtr->getValue()) != completed) {
      counterPtr->setValue(static_cast<int>(completed));
    }
    if (counterPtr->hasChanged() && counterLabel) {
      counterLabel->setText(
          tr("Bombardment %1 out of %2")
              .arg(completed)
              .arg(bombardmentCountSpinBox ? bombardmentCountSpinBox->value() : 0));
    }
    if (completed != lastLiveTallyBombardments) {
      updateLiveStatistics(activeSimulation->snapshot_tallies(), completed);
      lastLiveTallyBombardments = completed;
    }
  }

}

void MainWindow::updateLiveStatistics(const common::Tallies &tallies,
                                      std::size_t bombardmentCount) {
  latestLiveTallies = tallies;
  latestLiveTalliesBombardmentCount = bombardmentCount;
  refreshLiveHistograms(tallies, bombardmentCount);
  const auto rangeAngstrom =
      (common::length_quantity(static_cast<double>(tallies.range_stats.mean) *
                   common::si::meter) /
       common::angstrom)
        .value();
  const auto stragglingAngstrom =
      (common::length_quantity(
         std::sqrt(static_cast<double>(tallies.range_stats.variance())) *
         common::si::meter) /
       common::angstrom)
        .value();
  const auto energyLossEv =
        (common::energy_quantity(static_cast<double>(tallies.energy_loss_stats.mean) *
                     common::si::joule) /
         common::electron_volt)
          .value();
  const auto energyLossStragglingEv =
        (common::energy_quantity(
           std::sqrt(static_cast<double>(tallies.energy_loss_stats.variance())) *
           common::si::joule) /
         common::electron_volt)
          .value();
  const QStringList values{
      QString::number(tallies.total_displacements),
      QString::number(tallies.sputtered),
      tallies.range_stats.n == 0 ? QStringLiteral("--")
                                 : QString::number(rangeAngstrom, 'f', 2),
      tallies.range_stats.n < 2 ? QStringLiteral("--")
                                : QString::number(stragglingAngstrom, 'f', 2),
      tallies.energy_loss_stats.n == 0 ? QStringLiteral("--")
                                       : QString::number(energyLossEv, 'f', 2),
      tallies.energy_loss_stats.n < 2 ? QStringLiteral("--")
                                      : QString::number(energyLossStragglingEv, 'f', 2)};
  const auto count = std::min(liveStatisticValues.size(), values.size());
  for (qsizetype i = 0; i < count; ++i) {
    if (liveStatisticValues.at(i)) {
      liveStatisticValues.at(i)->setText(values.at(i));
    }
  }
}

} // namespace gui
} // namespace amcset
