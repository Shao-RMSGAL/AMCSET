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

// Boost includes
// #include <boost/thread/thread.hpp>

// Standard library includes
#include <exception>
#include <limits>

// Qt includes
#include <QCheckBox>
#include <QComboBox>
#include <QDebug>
#include <QDir>
#include <QFileDialog>
#include <QFormLayout>
#include <QHeaderView>
#include <QMenu>
#include <QMenuBar>
#include <QObject>
#include <QSizePolicy>
#include <QSplitter>
#include <QStatusBar>
#include <QTableView>
#include <QTextBrowser>
#include <QToolBar>
#include <QVBoxLayout>
#include <glog/logging.h>
#include <QtContainerFwd>
#include <QtWidgets>
#include <QtGraphs/QBarSet>
#include <QtGraphs/QLineSeries>
#include <QtGraphs/QValueAxis>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <qabstractitemmodel.h>
#include <qboxlayout.h>
#include <qnamespace.h>
#include <qobject.h>
#include <qstandarditemmodel.h>
// #include <sys/wait.h>
// #include <unistd.h>

// Local includes
#include "3dview.h"
#include "amcset_common.h"
#include "amcset_utilities.h"
#include "application.h"
#include "canvas.h"
#include "counter.h"
#include "layer_delegate.h"
#include "layerwindow.h"
#include "mainwindow.h"
#include "simulation.h"

namespace amcset {
namespace gui {

namespace {
enum class HistogramDataKind {
  Generic,
  Displacement,
  ElectronicLoss,
  NuclearLoss,
  TotalLoss
};

QVector<QString> collectCsvFilesForDisplay(const QString &outputDirPath,
                                          const QStringList &fallbackPaths) {
  QSet<QString> seen;
  QVector<QString> csvFiles;

  auto addIfValid = [&](const QString &path) {
    if (path.isEmpty()) {
      return;
    }

    const auto absPath = QFileInfo(path).absoluteFilePath();
    if (absPath.isEmpty() || !QFileInfo(absPath).exists() ||
        !QFileInfo(absPath).isFile()) {
      return;
    }

    if (!seen.contains(absPath)) {
      seen.insert(absPath);
      csvFiles.push_back(absPath);
    }
  };

  const auto baseDir = outputDirPath.trimmed();
  if (!baseDir.isEmpty()) {
    const QDir dir(baseDir);
    const auto csvNames = dir.entryList({"*.csv"}, QDir::Files | QDir::Readable,
                                        QDir::Name);
    for (const auto &csvName : csvNames) {
      addIfValid(dir.absoluteFilePath(csvName));
    }

    const auto nestedDirPath = dir.filePath("output");
    const QDir nestedDir(nestedDirPath);
    if (nestedDir.exists()) {
      const auto nestedCsvNames = nestedDir.entryList({"*.csv"},
                                                    QDir::Files | QDir::Readable,
                                                    QDir::Name);
      for (const auto &csvName : nestedCsvNames) {
        addIfValid(nestedDir.absoluteFilePath(csvName));
      }
    }
  }

  for (const auto &fallbackPath : fallbackPaths) {
    addIfValid(fallbackPath);
  }

  std::sort(csvFiles.begin(), csvFiles.end(), [](const QString &lhs,
                                               const QString &rhs) {
    const auto lhsBase = QFileInfo(lhs).baseName();
    const auto rhsBase = QFileInfo(rhs).baseName();
    if (lhsBase == rhsBase) {
      return lhs < rhs;
    }
    static const QSet<QString> prioritized = {"displacement_histogram",
                                              "electronic_loss_histogram",
                                              "nuclear_loss_histogram",
                                              "total_loss_histogram"};
    if (prioritized.contains(lhsBase) && !prioritized.contains(rhsBase)) {
      return true;
    }
    if (prioritized.contains(rhsBase) && !prioritized.contains(lhsBase)) {
      return false;
    }
    return lhsBase < rhsBase;
  });

  return csvFiles;
}

HistogramDataKind detectHistogramKind(const QString &filePath) {
  const auto name = QFileInfo(filePath).baseName().toLower();
  if (name.contains("displacement")) return HistogramDataKind::Displacement;
  if (name.contains("electronic")) return HistogramDataKind::ElectronicLoss;
  if (name.contains("nuclear")) return HistogramDataKind::NuclearLoss;
  if (name.contains("total")) return HistogramDataKind::TotalLoss;
  return HistogramDataKind::Generic;
}

QVector<QPointF> loadHistogramPoints(const QString &filePath,
                                    QString *errorMessage = nullptr) {
  QVector<QPointF> points;
  QFile csvFile(filePath);
  if (!csvFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
    if (errorMessage) {
      *errorMessage = QString("Unable to open %1").arg(filePath);
    }
    return points;
  }

  QTextStream stream(&csvFile);
  if (!stream.atEnd()) {
    stream.readLine();
  }

  while (!stream.atEnd()) {
    const auto line = stream.readLine();
    if (line.trimmed().isEmpty()) continue;

    const auto fields = line.split(',');
    if (fields.size() < 3) continue;

    const double left = fields[0].toDouble();
    const double right = fields[1].toDouble();
    const double value = fields[2].toDouble();
    if (std::isfinite(left) && std::isfinite(right) && std::isfinite(value)) {
      points.append(QPointF((left + right) / 2.0, value));
    }
  }
  csvFile.close();

  if (points.isEmpty() && errorMessage) {
    *errorMessage = QString("No valid numeric rows found in %1").arg(filePath);
  }
  return points;
}

QString histogramAxisLabel(HistogramDataKind kind, bool xAxis) {
  switch (kind) {
    case HistogramDataKind::Displacement:
      return xAxis ? QStringLiteral("Depth (Å)") : QStringLiteral("Displacements");
    case HistogramDataKind::ElectronicLoss:
      return xAxis ? QStringLiteral("Depth (Å)") : QStringLiteral("Electronic loss (eV)");
    case HistogramDataKind::NuclearLoss:
      return xAxis ? QStringLiteral("Depth (Å)") : QStringLiteral("Nuclear loss (eV)");
    case HistogramDataKind::TotalLoss:
      return xAxis ? QStringLiteral("Depth (Å)") : QStringLiteral("Total loss (eV)");
    case HistogramDataKind::Generic:
    default:
      return xAxis ? QStringLiteral("Depth") : QStringLiteral("Value");
  }
}

QString buildSimulationSummaryHtml(const QString &outputDir,
                                  const QStringList &filePaths,
                                  const QString &logFile,
                                  bool loggingEnabled,
                                  const QString &incidentType,
                                  bool damageCascade,
                                  double incidentEnergy,
                                  int incidentMass,
                                  int incidentCharge,
                                  double incidentAngle,
                                  double range,
                                  int bombardmentCount,
                                  int threadCount,
                                  int layerCount) {
  QStringList rows;
  rows << QStringLiteral("<h2>Simulation Summary</h2>")
       << QStringLiteral("<h3>Input settings</h3>")
       << QStringLiteral("<table border='0' cellspacing='4' cellpadding='4'>")
       << QStringLiteral("<tr><td><b>Output directory</b></td><td>%1</td></tr>").arg(outputDir)
       << QStringLiteral("<tr><td><b>Incident type</b></td><td>%1</td></tr>").arg(incidentType)
       << QStringLiteral("<tr><td><b>Damage cascade</b></td><td>%1</td></tr>").arg(damageCascade ? "Enabled" : "Disabled")
       << QStringLiteral("<tr><td><b>Incident energy</b></td><td>%1 keV</td></tr>").arg(QString::number(incidentEnergy, 'f', 1))
       << QStringLiteral("<tr><td><b>Incident mass</b></td><td>%1 amu</td></tr>").arg(incidentMass)
       << QStringLiteral("<tr><td><b>Incident charge</b></td><td>%1</td></tr>").arg(incidentCharge)
       << QStringLiteral("<tr><td><b>Incident angle</b></td><td>%1°</td></tr>").arg(QString::number(incidentAngle, 'f', 2))
       << QStringLiteral("<tr><td><b>Range</b></td><td>%1 Å</td></tr>").arg(QString::number(range, 'f', 1))
       << QStringLiteral("<tr><td><b>Bombardments</b></td><td>%1</td></tr>").arg(bombardmentCount)
       << QStringLiteral("<tr><td><b>Threads</b></td><td>%1</td></tr>").arg(threadCount)
       << QStringLiteral("<tr><td><b>Layers</b></td><td>%1</td></tr>").arg(layerCount)
       << QStringLiteral("<tr><td><b>Displacement histogram</b></td><td>%1</td></tr>").arg(filePaths.value(0, QStringLiteral("(not set)")))
       << QStringLiteral("<tr><td><b>Electronic loss histogram</b></td><td>%1</td></tr>").arg(filePaths.value(1, QStringLiteral("(not set)")))
       << QStringLiteral("<tr><td><b>Nuclear loss histogram</b></td><td>%1</td></tr>").arg(filePaths.value(2, QStringLiteral("(not set)")))
       << QStringLiteral("<tr><td><b>Total loss histogram</b></td><td>%1</td></tr>").arg(filePaths.value(3, QStringLiteral("(not set)")))
       << QStringLiteral("<tr><td><b>Log file</b></td><td>%1</td></tr>").arg(logFile)
       << QStringLiteral("<tr><td><b>Logging enabled</b></td><td>%1</td></tr>").arg(loggingEnabled ? "Yes" : "No")
       << QStringLiteral("</table>");

  const auto summaryPath = QDir(outputDir).filePath("tallies.json");
  QFile talliesFile(summaryPath);
  if (talliesFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
    const auto json = talliesFile.readAll();
    talliesFile.close();
    const auto doc = QJsonDocument::fromJson(json);
    if (doc.isObject()) {
      const auto obj = doc.object();
      rows << QStringLiteral("<h3>Statistics</h3>")
           << QStringLiteral("<table border='0' cellspacing='4' cellpadding='4'>")
           << QStringLiteral("<tr><td><b>Incidents</b></td><td>%1</td></tr>").arg(obj.value("incidents").toVariant().toString())
           << QStringLiteral("<tr><td><b>Total displacements</b></td><td>%1</td></tr>").arg(obj.value("total_displacements").toVariant().toString())
           << QStringLiteral("<tr><td><b>Sputtered</b></td><td>%1</td></tr>").arg(obj.value("sputtered").toVariant().toString())
           << QStringLiteral("<tr><td><b>Mean range</b></td><td>%1</td></tr>").arg(QString::number(obj.value("range_mean").toDouble(), 'g', 6))
           << QStringLiteral("<tr><td><b>Mean energy loss</b></td><td>%1</td></tr>").arg(QString::number(obj.value("energy_loss_mean").toDouble(), 'g', 6))
           << QStringLiteral("</table>");
    }
  }

  return rows.join(QStringLiteral("\n"));
}
} // namespace

class HistogramPlotWidget : public QWidget {
  Q_OBJECT
public:
  explicit HistogramPlotWidget(const QString &filePath, const QString &title,
                               HistogramDataKind kind = HistogramDataKind::Generic,
                               QWidget *parent = nullptr)
      : QWidget(parent), title_(title), kind_(kind), series_(new QLineSeries(this)) {
    const auto points = loadHistogramPoints(filePath);
    for (const auto &point : points) {
      series_->append(point);
    }
  }

protected:
  void paintEvent(QPaintEvent *event) override {
    QWidget::paintEvent(event);

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    painter.fillRect(rect(), QColor(245, 245, 245));

    const auto points = series_->points();
    if (points.size() < 2) {
      painter.setPen(Qt::gray);
      painter.drawText(rect(), Qt::AlignCenter, tr("No data available"));
      return;
    }

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    for (const auto &point : points) {
      minX = std::min(minX, point.x());
      maxX = std::max(maxX, point.x());
      minY = std::min(minY, point.y());
      maxY = std::max(maxY, point.y());
    }

    if (qFuzzyCompare(maxX, minX)) {
      minX -= 1.0;
      maxX += 1.0;
    }
    if (qFuzzyCompare(maxY, minY)) {
      minY = std::min(0.0, minY) - 1.0;
      maxY = std::max(1.0, maxY) + 1.0;
    }

    const auto margin = 24.0;
    auto transform = [&](const QPointF &pt) {
      const double sx = (pt.x() - minX) / std::max(maxX - minX, 1e-12);
      const double sy = (pt.y() - minY) / std::max(maxY - minY, 1e-12);
      const double x = margin + sx * (width() - 2.0 * margin);
      const double y = height() - margin - sy * (height() - 2.0 * margin);
      return QPointF(x, y);
    };


    QPen pen(Qt::blue, 2.0);
    painter.setPen(pen);
    QPolygonF polygon;
    for (const auto &point : points) {
      polygon << transform(point);
    }
    painter.drawPolyline(polygon);

    painter.setPen(Qt::black);
    painter.drawLine(QPointF(margin, height() - margin),
                     QPointF(width() - margin, height() - margin));
    painter.drawLine(QPointF(margin, margin),
                     QPointF(margin, height() - margin));

    const auto xLabel = histogramAxisLabel(kind_, true);
    const auto yLabel = histogramAxisLabel(kind_, false);
    painter.drawText(QRectF(0, height() - 20, width(), 20), Qt::AlignCenter, xLabel);
    painter.save();
    painter.rotate(-90.0);
    painter.drawText(QRectF(-height(), 0, height(), 20), Qt::AlignCenter, yLabel);
    painter.restore();
    painter.drawText(QRectF(0, 0, width(), 18), Qt::AlignCenter, title_);
  }

private:
  QString title_;
  HistogramDataKind kind_;
  QLineSeries *series_;
};

// TODO: Implement proper error handling
void MainWindow::handleSimulationResult(int result) {
  activeSimulation = nullptr;
  if (result != 0) {
    std::cerr << "Something went wrong" << std::endl;
    std::terminate();
  } else {
    std::cerr << "Success" << std::endl;
    refreshResultsPane();
  }
};

void MainWindow::launchSimulationThread(
    std::unique_ptr<common::Simulation> simulation) {

  simulationThread = new SimulationThread(std::move(simulation));
  connect(simulationThread, &SimulationThread::resultReady, this,
          &MainWindow::handleSimulationResult);
  simulationThread->start();
}

QString MainWindow::buildOutputPath(const QString &filename) const {
  const auto outputDir = outputDirectoryEdit ? outputDirectoryEdit->text().trimmed() : QString();
  const auto normalized = filename.trimmed();

  if (normalized.isEmpty()) {
    return outputDir;
  }

  const QFileInfo fileInfo(normalized);
  if (fileInfo.isAbsolute()) {
    return normalized;
  }

  QString cleanPath = normalized;
  const QString outputDirPrefix = QStringLiteral("output/");
  if (cleanPath.startsWith(outputDirPrefix)) {
    cleanPath = cleanPath.mid(outputDirPrefix.size());
  } else if (cleanPath.startsWith("./")) {
    cleanPath = cleanPath.mid(2);
  }

  if (outputDir.isEmpty()) {
    return cleanPath;
  }

  return QDir(outputDir).absoluteFilePath(cleanPath);
}

void MainWindow::updateResultPane(const QString &message) {
  if (!resultPane) {
    return;
  }
  const auto text = resultPane->toPlainText();
  if (!text.isEmpty()) {
    resultPane->appendPlainText("");
  }
  resultPane->appendPlainText(message);
}

void MainWindow::browseOutputDirectory(QLineEdit *lineEdit,
                                      const QString &title) {
  if (!lineEdit) return;
  const auto directory = QFileDialog::getExistingDirectory(
      this, title, lineEdit->text().trimmed(), QFileDialog::ShowDirsOnly);
  if (!directory.isEmpty()) {
    lineEdit->setText(directory);
  }
}

void MainWindow::browseOutputFile(QLineEdit *lineEdit, const QString &title,
                                 const QString &defaultSuffix) {
  if (!lineEdit) return;
  const auto selected = QFileDialog::getSaveFileName(
      this, title, lineEdit->text().trimmed(),
      defaultSuffix.isEmpty() ? QStringLiteral("All Files (*)")
                              : QStringLiteral("%1 (*.%2)").arg(defaultSuffix.toUpper(), defaultSuffix),
      nullptr, QFileDialog::DontUseNativeDialog);
  if (!selected.isEmpty()) {
    lineEdit->setText(selected);
  }
}

void MainWindow::populateResultCharts() {
  auto populateTabWidget = [this](QTabWidget *target) {
    if (!target) {
      return;
    }

    // Avoid calling QTabWidget::clear() directly — some Qt builds can crash
    // when clearing tabs that own complex child layouts. Remove tabs one by
    // one and delete widgets safely instead. Guard against an uninitialized
    // internal QTabBar which makes `count()` crash in certain Qt builds.
    if (target->tabBar()) {
      while (target->tabBar()->count() > 0) {
        QWidget *w = target->widget(0);
        target->removeTab(0);
        if (w) {
          w->deleteLater();
        }
      }
    }

    auto addFileChart = [&](const QString &title, const QString &filePath) {
      const QFileInfo fileInfo(filePath);
      if (!fileInfo.exists() || !fileInfo.isFile()) {
        return;
      }
      const auto kind = detectHistogramKind(filePath);
      // Parent the plot widget to the tab widget to avoid reparenting/layout
      // issues when `QTabWidget::clear()` is called.
      auto *plotWidget = new HistogramPlotWidget(filePath, title, kind, target);
      target->addTab(plotWidget, title);
    };

    const auto outputDirPath = outputDirectoryEdit ? outputDirectoryEdit->text().trimmed() : QStringLiteral("output");
    const auto fallbackFilePaths = QStringList{
        displacementOutputEdit ? displacementOutputEdit->text().trimmed() : QString(),
        electronicOutputEdit ? electronicOutputEdit->text().trimmed() : QString(),
        nuclearOutputEdit ? nuclearOutputEdit->text().trimmed() : QString(),
        totalOutputEdit ? totalOutputEdit->text().trimmed() : QString(),
        QStringLiteral("output/collisions.csv"),
    };

    auto csvFiles = collectCsvFilesForDisplay(outputDirPath, fallbackFilePaths);
    if (!csvFiles.isEmpty()) {
      std::sort(csvFiles.begin(), csvFiles.end(), [](const QString &lhs, const QString &rhs) {
        const auto lhsBase = QFileInfo(lhs).baseName();
        const auto rhsBase = QFileInfo(rhs).baseName();
        if (lhsBase == rhsBase) {
          return lhs < rhs;
        }
        if (lhsBase == "collisions") {
          return true;
        }
        if (rhsBase == "collisions") {
          return false;
        }
        return lhsBase < rhsBase;
      });

      for (const auto &csvPath : csvFiles) {
        const auto title = QFileInfo(csvPath).baseName();
        addFileChart(title, csvPath);
      }
      return;
    }

    const int existingTabs = target->tabBar() ? target->tabBar()->count() : 0;
    if (existingTabs == 0) {
      auto *placeholder = new QLabel("No histogram CSV files found.", target);
      placeholder->setAlignment(Qt::AlignCenter);
      placeholder->setWordWrap(true);
      target->addTab(placeholder, "No data");
    }
  };

  if (resultsTabs) {
    populateTabWidget(resultsTabs);
  }
  if (histogramTabs) {
    populateTabWidget(histogramTabs);
  }
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
    google::SetLogDestination(google::INFO, stdPath.c_str());
    google::SetLogDestination(google::WARNING, stdPath.c_str());
    google::SetLogDestination(google::ERROR, stdPath.c_str());
    google::SetLogDestination(google::FATAL, stdPath.c_str());
    updateResultPane(QString("Logging enabled to %1").arg(logFile));
  } else {
    google::SetLogDestination(google::INFO, "");
    google::SetLogDestination(google::WARNING, "");
    google::SetLogDestination(google::ERROR, "");
    google::SetLogDestination(google::FATAL, "");
  }

  const auto outputDir = outputDirectoryEdit->text().trimmed();
  if (!outputDir.isEmpty()) {
    QDir().mkpath(outputDir);
  }

  counter->reset();

  // Avoid edge case where clearing the screen takes a long time.
  if (use2DView->isChecked()) {
    canvas->clear();
  }
  if (use3DView->isChecked()) {
    view3D->clearScene();
  }

  // Create AMCSET settings from the GUI configuration.
  const auto displacementOutput = buildOutputPath(displacementOutputEdit->text().trimmed());
  const auto electronicOutput = buildOutputPath(electronicOutputEdit->text().trimmed());
  const auto nuclearOutput = buildOutputPath(nuclearOutputEdit->text().trimmed());
  const auto totalOutput = buildOutputPath(totalOutputEdit->text().trimmed());
  const auto particleTreeOutput = buildOutputPath(particleTreeEdit->text().trimmed());

  common::energy_quantity electron_stopping_energy =
      40.0 * common::electron_volt; // Energy at which electrons stop.
  common::energy_quantity ion_stopping_energy =
      40.0 * common::electron_volt; // Typical for iron.
  common::energy_quantity ion_displacement_energy =
      40.0 * common::electron_volt; // Typical for iron
  const bool log_single_displacement =
      logSingleDisplacementCheck ? logSingleDisplacementCheck->isChecked() : false;
  const size_t divisor_angle_number =
      angleDivisorSpinBox ? static_cast<size_t>(angleDivisorSpinBox->value()) : 1000;
  const size_t flying_distance_number =
      flyingDistanceSpinBox ? static_cast<size_t>(flyingDistanceSpinBox->value()) : 1000;

  const bool log_end_of_flying_distance_only =
      logEndOfFlyingDistanceOnlyCheck ? logEndOfFlyingDistanceOnlyCheck->isChecked() : false;
  const bool log_stopping_point_only =
      logStoppingPointOnlyCheck ? logStoppingPointOnlyCheck->isChecked() : false;
  const std::string electron_screening_path =
      electronScreeningParametersEdit ? electronScreeningParametersEdit->text().trimmed().toStdString()
                                      : "common/assets/electron_screening_parameters.csv";
  const std::string mott_path =
      mottScatteringParametersEdit ? mottScatteringParametersEdit->text().trimmed().toStdString()
                                  : "common/assets/mott_scattering_parameters.csv";

  bool is_electron =
      (incidentElectron->checkState() ==
       Qt::Checked); // Start with ions. TEST: Test with electrons later.
  bool enable_damage_cascade = (damageCascade->checkState() == Qt::Checked);
  size_t mass_number = incidentMass->value();
  size_t z_number = incidentCharge->value();
  common::angle_quantity incident_angle =
      common::degree * incidentAngle->value();
  common::length_quantity range = rangeSpinBox->value() * common::angstrom;
  size_t bombardment_count =
      static_cast<size_t>(bombardmentCountSpinBox->value());
  common::energy_quantity incident_energy =
      incidentEnergy->value() * common::kilo_electron_volt;
  size_t thread_count = threadCountSpinBox->value();

  // PERF: Add an option to use GPU (Use Kokkos for cross-platform
  // compatability. Note that it is not available in conan. Either make a
  // conan recipe for it or use another installation method.

    // Map GUI selection to backend ScreeningModel
    common::Simulation::Settings::ScreeningModel selected_model =
      (screeningModelCombo && screeningModelCombo->currentIndex() == 1)
        ? common::Simulation::Settings::ScreeningModel::QuickBasic
        : common::Simulation::Settings::ScreeningModel::ThomasFermi;

    auto settings = common::Simulation::Settings(
      electron_stopping_energy, z_number, mass_number, enable_damage_cascade,
      ion_stopping_energy, ion_displacement_energy, log_single_displacement,
      divisor_angle_number, flying_distance_number, range, bombardment_count,
      is_electron, incident_energy, thread_count, incident_angle, 100, 100, 100,
      100, displacementOutput.toStdString(), electronicOutput.toStdString(),
      nuclearOutput.toStdString(), totalOutput.toStdString(),
      enableLoggingCheck->isChecked(), particleTreeOutput.toStdString(),
      electron_screening_path, mott_path, selected_model,
      log_end_of_flying_distance_only, log_stopping_point_only);

  // Create AMCSET volume (TODO: Read this from a file and allow a user to
  // configure)

  const size_t layerCount = layerModel->rowCount();

  if (layerCount == 0) {
    throw std::logic_error("Layer count cannot be zero.");
  }

  std::vector<common::Layer> layers;
  for (size_t layerRow = 0; layerRow < layerCount; layerRow++) {
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
    // TODO: Check for duplicate elements
    size_t elementCount = materialModel->rowCount();
    common::Layer::material_vector material;
    auto proportionSum = 0.0;
    for (size_t elementRow = 0; elementRow < elementCount; elementRow++) {
      const auto atomicSymbolColumn = 0;
      const auto proportionColumn = 1;

      auto atomicSymbol = materialModel->index(elementRow, atomicSymbolColumn)
                              .data(Qt::DisplayRole)
                              .toString();
      auto atomicProportion = materialModel->index(elementRow, proportionColumn)
                                  .data(Qt::DisplayRole)
                                  .toDouble();
      // TODO: Use displacement energy
      // auto displacementEnergy =
      //     materialModel->index(elementRow, displacementEnergyColumn)
      //         .data(Qt::DisplayRole)
      //         .toDouble();
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
    }
    layers.push_back(common::Layer(std::move(material), depth, density));
  }

  auto volume = common::Volume(std::move(layers));

  // Create the AMCSET simulation object. PERF: Replace the hacky lambda drawing
  // function with proper queue feed to the frontend so drawing doesn't slow
  // down the simulation.
  // Use an empty draw function; frontend will poll the simulation queue instead.
  auto simulation = std::make_unique<common::Simulation>(settings, std::move(volume),
                                                         [](const std::vector<common::Coordinate> &, size_t) {});

  activeSimulation = simulation.get();

  try {
    launchSimulationThread(std::move(simulation));
    // Start a short-interval timer to poll for completed bombardment
    // collision vectors and draw them in the UI. This yields for 1ms
    // after each processed bombardment to give the backend time to publish.
    if (!bombardmentPollTimer) {
      bombardmentPollTimer = new QTimer(this);
      connect(bombardmentPollTimer, &QTimer::timeout, this, &MainWindow::pollBombardmentQueue);
      bombardmentPollTimer->start(1);
    }
  } catch (const std::exception &e) {
    common::print_exception(e);
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

void MainWindow::pollBombardmentQueue() {
  if (!activeSimulation) return;
  std::vector<common::CollisionPoint> vec;
  if (!activeSimulation->try_pop_bombardment_results(vec)) {
    // no data
    return;
  }

  const auto canvasPtr = QPointer<Canvas>(canvas);
  const auto viewPtr = QPointer<Widget3DView>(view3D);
  const auto counterPtr = QPointer<Counter>(counter);

  const auto range = activeSimulation->get_settings(&common::Simulation::Settings::range_);

  if (canvasPtr && use2DView && use2DView->checkState() == Qt::Checked) {
    const auto width = canvasPtr->width();
    const auto half_height = canvasPtr->height() / 2;
    for (size_t i = 0; i + 1 < vec.size(); ++i) {
      const auto &p1 = vec[i].coordinate;
      const auto &p2 = vec[i + 1].coordinate;
      const auto start_x = (p1.x_ / range).value() * width + half_height;
      const auto start_z = (p1.z_ / range).value() * width;
      const auto finish_x = (p2.x_ / range).value() * width + half_height;
      const auto finish_z = (p2.z_ / range).value() * width;
      canvasPtr->addLine({(int)start_z, (int)start_x, (int)finish_z, (int)finish_x});
    }
  }

  if (viewPtr && use3DView && use3DView->checkState() == Qt::Checked) {
    for (size_t i = 0; i + 1 < vec.size(); ++i) {
      const auto &p1 = vec[i].coordinate;
      const auto &p2 = vec[i + 1].coordinate;
      const auto start_x = p1.x_ / common::angstrom;
      const auto start_y = p1.y_ / common::angstrom;
      const auto start_z = p1.z_ / common::angstrom;
      const auto finish_x = p2.x_ / common::angstrom;
      const auto finish_y = p2.y_ / common::angstrom;
      const auto finish_z = p2.z_ / common::angstrom;
      QVector3D point1(start_x, start_y, start_z);
      QVector3D point2(finish_x, finish_y, finish_z);
      viewPtr->requestDrawLine(point1, point2, Qt::green);
    }
  }

  if (counterPtr && activeSimulation != nullptr) {
    counterPtr->setValue(static_cast<int>(activeSimulation->get_completed_bombardment_count()));
  }

  // Yield for 1ms to allow backend to publish
  QThread::msleep(1);
}

MainWindow::~MainWindow() {
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

  // TODO: Uncomment all
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
  auto *cancelButton = new QPushButton("Cancel Run ");
  cancelButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  // TODO: Figure out how to implement cancellation.

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
      "QCheckBox::indicator:checked { background: #3d7de8; border-color: #2d66c8; }"
  );
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
  // Defer population until the event loop runs to ensure internal tab
  // structures are fully initialized and avoid crashes in some Qt builds.
  QTimer::singleShot(0, this, &MainWindow::populateResultCharts);
  mainTabs->addTab(histogramTabs, "Histograms");

  simulationSummaryView = new QTextBrowser(this);
  simulationSummaryView->setOpenExternalLinks(false);
  simulationSummaryView->setReadOnly(true);
  simulationSummaryView->setStyleSheet("QTextBrowser { background: #f8f8f8; }");
  // refreshResultsPane() is deferred until after settings widgets are created.
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

  createLayerRow();
  show();
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

  connect(incidentElectron, &QCheckBox::checkStateChanged, this,
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
  const auto maxRange = std::numeric_limits<decltype(minRange)>::max();
  const auto defaultDepth = 1000;
  rangeSpinBox->setRange(minRange, maxRange);
  rangeSpinBox->setValue(defaultDepth);
  connect(rangeSpinBox, &QDoubleSpinBox::valueChanged, view3D,
          &Widget3DView::setRange);
  // view3D->setRange(range / common::angstrom);

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

  // TODO: Put this somewhere else, around the view window maybe.
  counterLabel = new QLabel(this);
  counterLabel->setText(
      QString("Bombardment %1 out of %2").arg(0).arg(defaultBombardments));
  counter = new Counter(this);
  QTimer *timer = new QTimer(this);
  connect(timer, &QTimer::timeout, this, [this]() {
    if (activeSimulation != nullptr) {
      const auto completed = static_cast<int>(
          activeSimulation->get_completed_bombardment_count());
      if (completed != counter->getValue()) {
        counter->setValue(completed);
      }
    }

    if (counter->hasChanged()) {
      int value = counter->getValue();
      auto totalBombardments = bombardmentCountSpinBox->value();
      counterLabel->setText(QString("Bombardment %1 out of %2")
                                .arg(value)
                                .arg(totalBombardments));
    }
  });
  timer->start(10); // check every 100ms
  incidentFormLayout->addRow("Bombardment Counter", counterLabel);

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
  const auto minDepth = 0.01; // TODO: More elegant solution (optional)
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
  const auto defaultDensity =
      1000.0; // TODO: Pick some proper density for a default material
  auto *densitySpinBox = new QDoubleSpinBox();
  densitySpinBox->setRange(minDensity, maxDensity);
  densitySpinBox->setValue(defaultDensity);
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
  layerModel->setItem(row, layerLineEditColumn, layerLineEdit);
  layerTable->setIndexWidget(startingDepthLabelIndex, startingDepthLabel);
  layerTable->setIndexWidget(depthSpinBoxIndex, depthSpinBox);
  layerTable->setIndexWidget(densitySpinBoxIndex, densitySpinBox);
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
void MainWindow::refreshResultsPane() {
  QStringList summaryPaths;
  summaryPaths << (displacementOutputEdit ? displacementOutputEdit->text().trimmed() : QString())
               << (electronicOutputEdit ? electronicOutputEdit->text().trimmed() : QString())
               << (nuclearOutputEdit ? nuclearOutputEdit->text().trimmed() : QString())
               << (totalOutputEdit ? totalOutputEdit->text().trimmed() : QString());

  const auto outputDir = outputDirectoryEdit ? outputDirectoryEdit->text().trimmed() : QStringLiteral("output");
  if (simulationSummaryView) {
    simulationSummaryView->setHtml(buildSimulationSummaryHtml(
        outputDir, summaryPaths,
        logFileEdit ? logFileEdit->text().trimmed() : QStringLiteral("amcset.log"),
        enableLoggingCheck && enableLoggingCheck->isChecked(),
        incidentElectron && incidentElectron->isChecked() ? QStringLiteral("Electron") : QStringLiteral("Ion"),
        damageCascade && damageCascade->isChecked(),
        incidentEnergy ? incidentEnergy->value() : 0.0,
        incidentMass ? incidentMass->value() : 1,
        incidentCharge ? incidentCharge->value() : 1,
        incidentAngle ? incidentAngle->value() : 0.0,
        rangeSpinBox ? rangeSpinBox->value() : 0.0,
        bombardmentCountSpinBox ? bombardmentCountSpinBox->value() : 1,
        threadCountSpinBox ? threadCountSpinBox->value() : 1,
        layerModel ? static_cast<int>(layerModel->rowCount()) : 0));
  }

  if (!resultPane) return;
  resultPane->clear();
  QStringList summary;
  summary << "Output directory: " + outputDir;
  summary << "Displacement histogram: " + summaryPaths.value(0, QString());
  summary << "Electronic loss histogram: " + summaryPaths.value(1, QString());
  summary << "Nuclear loss histogram: " + summaryPaths.value(2, QString());
  summary << "Total loss histogram: " + summaryPaths.value(3, QString());
  summary << "Particle tree: " + (particleTreeEdit ? particleTreeEdit->text().trimmed() : QString());
  summary << "Log file: " + (logFileEdit ? logFileEdit->text().trimmed() : QString());
  summary << "Logging enabled: " + QString(enableLoggingCheck && enableLoggingCheck->isChecked() ? "Yes" : "No");
  summary << "Log single displacement: " + QString(logSingleDisplacementCheck && logSingleDisplacementCheck->isChecked() ? "Yes" : "No");
  summary << "Angle divisors: " + QString::number(angleDivisorSpinBox ? angleDivisorSpinBox->value() : 1000);
  summary << "Flying distances: " + QString::number(flyingDistanceSpinBox ? flyingDistanceSpinBox->value() : 1000);
  summary << "Log end-of-flight only: " + QString(logEndOfFlyingDistanceOnlyCheck && logEndOfFlyingDistanceOnlyCheck->isChecked() ? "Yes" : "No");
  summary << "Log stopping point only: " + QString(logStoppingPointOnlyCheck && logStoppingPointOnlyCheck->isChecked() ? "Yes" : "No");
  resultPane->setPlainText(summary.join("\n"));
  populateResultCharts();
}

void MainWindow::createSettingsDialog() {
  settingsDialog = new QDialog(this);
  settingsDialog->setWindowTitle("Simulation Settings");
  settingsDialog->resize(960, 680);
  settingsDialog->setMinimumWidth(900);
  auto *dialogLayout = new QVBoxLayout(settingsDialog);

  auto *settingsTabs = new QTabWidget(settingsDialog);
  auto *settingsTab = new QWidget(settingsTabs);
  auto *settingsLayout = new QVBoxLayout(settingsTab);
  auto *formLayout = new QFormLayout();
  formLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  formLayout->setLabelAlignment(Qt::AlignLeft);

  auto *electronSettingsGroup = new QGroupBox("Electron transport", settingsTab);
  auto *electronSettingsLayout = new QFormLayout(electronSettingsGroup);
  electronSettingsLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  electronSettingsLayout->setLabelAlignment(Qt::AlignLeft);

  auto *outputDirRow = new QWidget(settingsTab);
  auto *outputDirLayout = new QHBoxLayout(outputDirRow);
  outputDirLayout->setContentsMargins(0, 0, 0, 0);
  outputDirectoryEdit = new QLineEdit("output", outputDirRow);
  auto *outputDirBrowse = new QPushButton("Browse", outputDirRow);
  connect(outputDirBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputDirectory(outputDirectoryEdit, "Select output directory");
  });
  outputDirLayout->addWidget(outputDirectoryEdit);
  outputDirLayout->addWidget(outputDirBrowse);
  formLayout->addRow("Output directory", outputDirRow);

  const auto fileFields = std::vector<std::pair<QLineEdit **, const char *>>{
      {&displacementOutputEdit, "Displacement output"},
      {&electronicOutputEdit, "Electronic loss output"},
      {&nuclearOutputEdit, "Nuclear loss output"},
      {&totalOutputEdit, "Total loss output"},
      {&particleTreeEdit, "Particle tree output"},
  };

  for (const auto &entry : fileFields) {
    auto **field = entry.first;
    const char *label = entry.second;
    auto *row = new QWidget(settingsTab);
    auto *layout = new QHBoxLayout(row);
    layout->setContentsMargins(0, 0, 0, 0);
    *field = new QLineEdit(QString(QString(label).toLower().replace(" ", "_")) + ".csv", row);
    auto *browseButton = new QPushButton("File", row);
    connect(browseButton, &QPushButton::clicked, this, [this, field, label]() {
      browseOutputFile(*field, QString("Select %1").arg(label), "csv");
    });
    layout->addWidget(*field);
    layout->addWidget(browseButton);
    formLayout->addRow(label, row);
  }

  displacementOutputEdit->setText("displacement_histogram.csv");
  electronicOutputEdit->setText("electronic_loss_histogram.csv");
  nuclearOutputEdit->setText("nuclear_loss_histogram.csv");
  totalOutputEdit->setText("total_loss_histogram.csv");
  particleTreeEdit->setText("particle_tree.json");

  auto *logFileRow = new QWidget(settingsTab);
  auto *logFileLayout = new QHBoxLayout(logFileRow);
  logFileLayout->setContentsMargins(0, 0, 0, 0);
  logFileEdit = new QLineEdit("amcset.log", logFileRow);
  auto *logBrowse = new QPushButton("File", logFileRow);
  connect(logBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(logFileEdit, "Select log file", "log");
  });
  logFileLayout->addWidget(logFileEdit);
  logFileLayout->addWidget(logBrowse);
  formLayout->addRow("Log file", logFileRow);

  enableLoggingCheck = new QCheckBox(settingsTab);
  enableLoggingCheck->setChecked(false);
  formLayout->addRow("Enable logging", enableLoggingCheck);

  logSingleDisplacementCheck = new QCheckBox(settingsTab);
  logSingleDisplacementCheck->setChecked(false);
  formLayout->addRow("Log single displacement", logSingleDisplacementCheck);

  logEndOfFlyingDistanceOnlyCheck = new QCheckBox(settingsTab);
  logEndOfFlyingDistanceOnlyCheck->setChecked(false);
  electronSettingsLayout->addRow("Log end-of-flight only", logEndOfFlyingDistanceOnlyCheck);

  logStoppingPointOnlyCheck = new QCheckBox(settingsTab);
  logStoppingPointOnlyCheck->setChecked(false);
  electronSettingsLayout->addRow("Log stopping point only", logStoppingPointOnlyCheck);

  angleDivisorSpinBox = new QSpinBox(settingsTab);
  angleDivisorSpinBox->setRange(1, 1000000);
  angleDivisorSpinBox->setValue(1000);
  electronSettingsLayout->addRow("Angle divisors", angleDivisorSpinBox);

  flyingDistanceSpinBox = new QSpinBox(settingsTab);
  flyingDistanceSpinBox->setRange(1, 1000000);
  flyingDistanceSpinBox->setValue(1000);
  electronSettingsLayout->addRow("Flying distances", flyingDistanceSpinBox);

  auto *electronScreeningRow = new QWidget(settingsTab);
  auto *electronScreeningLayout = new QHBoxLayout(electronScreeningRow);
  electronScreeningLayout->setContentsMargins(0, 0, 0, 0);
  electronScreeningParametersEdit = new QLineEdit("common/assets/electron_screening_parameters.csv", settingsTab);
  auto *electronScreeningBrowse = new QPushButton("Browse", settingsTab);
  connect(electronScreeningBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(electronScreeningParametersEdit, "Select electron screening parameter file", "csv");
  });
  electronScreeningLayout->addWidget(electronScreeningParametersEdit);
  electronScreeningLayout->addWidget(electronScreeningBrowse);
  electronSettingsLayout->addRow("Electron screening parameters", electronScreeningRow);

  auto *mottRow = new QWidget(settingsTab);
  auto *mottLayout = new QHBoxLayout(mottRow);
  mottLayout->setContentsMargins(0, 0, 0, 0);
  mottScatteringParametersEdit = new QLineEdit("common/assets/mott_scattering_parameters.csv", settingsTab);
  auto *mottBrowse = new QPushButton("Browse", settingsTab);
  connect(mottBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(mottScatteringParametersEdit, "Select Mott scattering parameter file", "csv");
  });
  mottLayout->addWidget(mottScatteringParametersEdit);
  mottLayout->addWidget(mottBrowse);
  electronSettingsLayout->addRow("Mott scattering parameters", mottRow);

  // Screening model selection
  auto *screeningRow = new QWidget(settingsTab);
  auto *screeningLayout = new QHBoxLayout(screeningRow);
  screeningLayout->setContentsMargins(0, 0, 0, 0);
  screeningModelCombo = new QComboBox(settingsTab);
  screeningModelCombo->addItem("Thomas-Fermi (default)");
  screeningModelCombo->addItem("QuickBasic (ion.bas)");
  screeningLayout->addWidget(screeningModelCombo);
  electronSettingsLayout->addRow("Screening model", screeningRow);

  settingsLayout->addLayout(formLayout);
  settingsLayout->addWidget(electronSettingsGroup);
  settingsTabs->addTab(settingsTab, "Settings");

  auto *resultsPage = new QWidget(settingsTabs);
  auto *resultsLayout = new QVBoxLayout(resultsPage);
  resultsTabs = new QTabWidget(resultsPage);
  resultsLayout->addWidget(resultsTabs);
  settingsTabs->addTab(resultsPage, "Results");

  dialogLayout->addWidget(settingsTabs);

  auto *buttons = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel,
                                       Qt::Horizontal, settingsDialog);
  connect(buttons, &QDialogButtonBox::accepted, settingsDialog, &QDialog::accept);
  connect(buttons, &QDialogButtonBox::rejected, settingsDialog, &QDialog::reject);
  connect(buttons, &QDialogButtonBox::accepted, this, &MainWindow::refreshResultsPane);
  dialogLayout->addWidget(buttons);

  resultPane = new QPlainTextEdit(settingsDialog);
  resultPane->setReadOnly(true);
  resultPane->setPlaceholderText("Simulation results appear here...");
  resultPane->setMaximumHeight(180);
  dialogLayout->addWidget(resultPane);

  // refreshResultsPane() is deferred until after settings widgets are created.
}

void MainWindow::openSettingsDialog() {
  if (settingsDialog) {
    settingsDialog->show();
    settingsDialog->raise();
    settingsDialog->activateWindow();
  }
}

void MainWindow::createActions() { /*{{{*/
  settingsAction = new QAction(tr("Settings"), this);
  settingsAction->setStatusTip(tr("Open simulation settings"));
  connect(settingsAction, &QAction::triggered, this,
          &MainWindow::openSettingsDialog);

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
  fileMenu->addAction(settingsAction);
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

  settingsMenu = menuBar()->addMenu(tr("&Settings"));
  settingsMenu->addAction(settingsAction);

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

#include "mainwindow.moc"

} // namespace gui
} // namespace amcset
