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

#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QJsonDocument>
#include <QJsonObject>
#include <QPlainTextEdit>
#include <QTabWidget>
#include <QTextBrowser>
#include <QVBoxLayout>

namespace amcset {
namespace gui {

enum class HistogramDataKind {
  Generic,
  Displacement,
  ElectronicLoss,
  NuclearLoss,
  TotalLoss
};

namespace {
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
}

void MainWindow::openSettingsDialog() {
  if (settingsDialog) {
    settingsDialog->show();
    settingsDialog->raise();
    settingsDialog->activateWindow();
  }
}

} // namespace gui
} // namespace amcset

#include "mainwindow_results.moc"
