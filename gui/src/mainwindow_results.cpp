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
#include <QPushButton>
#include <QImage>
#include <QPainter>
#include <QPen>
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
  , Range
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
  if (name.contains("range")) return HistogramDataKind::Range;
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
    case HistogramDataKind::Range:
      return xAxis ? QStringLiteral("Range (Å)") : QStringLiteral("Count");
    case HistogramDataKind::Generic:
    default:
      return xAxis ? QStringLiteral("Depth") : QStringLiteral("Value");
  }
}

QString buildSimulationSummaryHtml(const QString &outputDir,
                                  const QStringList &filePaths,
                                  const QString &logFile,
                                  bool loggingEnabled,
                                  const QColor &backgroundColor,
                                  const QColor &textColor,
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
  rows << QStringLiteral("<body style='background-color:%1;color:%2;'>")
              .arg(backgroundColor.name(), textColor.name())
       << QStringLiteral("<h2>Simulation Summary</h2>")
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
      << QStringLiteral("<tr><td><b>Range histogram</b></td><td>%1</td></tr>").arg(filePaths.value(4, QStringLiteral("(not set)")))
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
               << QStringLiteral("<tr><td><b>Range straggling</b></td><td>%1</td></tr>").arg(QString::number(obj.value("range_straggling").toDouble(), 'g', 6))
           << QStringLiteral("<tr><td><b>Mean energy loss</b></td><td>%1</td></tr>").arg(QString::number(obj.value("energy_loss_mean").toDouble(), 'g', 6))
           << QStringLiteral("</table>");
    }
  }

  rows << QStringLiteral("</body>");
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
    setMinimumSize(320, 240);
    const auto points = loadHistogramPoints(filePath);
    for (const auto &point : points) {
      series_->append(point);
    }
  }

  void setHistogram(const common::Histogram &histogram) {
    QVector<QPointF> points;
    if (!histogram.bins.empty()) {
      const auto width = histogram.max - histogram.min;
      points.reserve(static_cast<int>(histogram.bins.size()));
      for (size_t i = 0; i < histogram.bins.size(); ++i) {
        const auto x = histogram.min + width * (static_cast<double>(i) + 0.5) /
                                      static_cast<double>(histogram.bins.size());
        points.append(QPointF(x, histogram.bins[i]));
      }
    }
    normalizedPerBombardment_ = false;
    series_->replace(points);
    update();
  }

  void setNormalizedHistogram(const common::Histogram &histogram,
                              std::size_t bombardmentCount) {
    QVector<QPointF> points;
    if (!histogram.bins.empty() && bombardmentCount > 0) {
      const auto width = histogram.max - histogram.min;
      points.reserve(static_cast<int>(histogram.bins.size()));
      for (size_t i = 0; i < histogram.bins.size(); ++i) {
        const auto x = histogram.min + width * (static_cast<double>(i) + 0.5) /
                                      static_cast<double>(histogram.bins.size());
        points.append(QPointF(
            x, histogram.bins[i] / static_cast<double>(bombardmentCount)));
      }
    }
    normalizedPerBombardment_ = true;
    series_->replace(points);
    update();
  }

  bool savePng(const QString &filePath) const {
    const QSize imageSize(1280, 800);
    QImage image(imageSize, QImage::Format_ARGB32_Premultiplied);
    image.fill(palette().color(QPalette::Base));
    QPainter painter(&image);
    drawPlot(painter, imageSize);
    painter.end();
    return image.save(filePath, "PNG");
  }

protected:
  void paintEvent(QPaintEvent *event) override {
    QWidget::paintEvent(event);
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, true);
    drawPlot(painter, size());
  }

private:
  static double niceStep(double range) {
    if (!std::isfinite(range) || range <= 0.0) {
      return 1.0;
    }
    const double rawStep = range / 5.0;
    const double magnitude = std::pow(10.0, std::floor(std::log10(rawStep)));
    const double fraction = rawStep / magnitude;
    const double niceFraction = fraction <= 1.0 ? 1.0
                                : fraction <= 2.0 ? 2.0
                                : fraction <= 5.0 ? 5.0
                                                  : 10.0;
    return niceFraction * magnitude;
  }

  static QString formatTick(double value, double step) {
    const int decimals = std::clamp(
        static_cast<int>(std::ceil(-std::log10(step))) + 1, 0, 8);
    return QString::number(value, 'f', decimals);
  }

  void drawPlot(QPainter &painter, const QSize &canvasSize) const {
    painter.setRenderHint(QPainter::Antialiasing, true);
    const auto colors = palette();
    painter.fillRect(QRect(QPoint(0, 0), canvasSize), colors.color(QPalette::Base));

    const auto points = series_->points();
    if (points.isEmpty()) {
      painter.setPen(colors.color(QPalette::Text));
      painter.drawText(QRect(QPoint(0, 0), canvasSize), Qt::AlignCenter,
                       tr("No data available"));
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

    minY = std::min(0.0, minY);
    const double xStep = niceStep(maxX - minX);
    const double yStep = niceStep(maxY - minY);
    minX = std::floor(minX / xStep) * xStep;
    maxX = std::ceil(maxX / xStep) * xStep;
    minY = std::floor(minY / yStep) * yStep;
    maxY = std::ceil(maxY / yStep) * yStep;
    if (maxX <= minX) maxX = minX + xStep;
    if (maxY <= minY) maxY = minY + yStep;

    const double leftMargin = 76.0;
    const double rightMargin = 18.0;
    const double topMargin = 30.0;
    const double bottomMargin = 58.0;
    const QRectF plotArea(leftMargin, topMargin,
                          std::max(1.0, canvasSize.width() - leftMargin - rightMargin),
                          std::max(1.0, canvasSize.height() - topMargin - bottomMargin));
    auto transform = [&](const QPointF &pt) {
      const double sx = (pt.x() - minX) / std::max(maxX - minX, 1e-12);
      const double sy = (pt.y() - minY) / std::max(maxY - minY, 1e-12);
      const double x = plotArea.left() + sx * plotArea.width();
      const double y = plotArea.bottom() - sy * plotArea.height();
      return QPointF(x, y);
    };

    const int xTickCount = std::min(10, static_cast<int>((maxX - minX) / xStep) + 1);
    const int yTickCount = std::min(10, static_cast<int>((maxY - minY) / yStep) + 1);
    painter.setFont(QFontDatabase::systemFont(QFontDatabase::GeneralFont));
    painter.setPen(colors.color(QPalette::Text));
    for (int i = 0; i < xTickCount; ++i) {
      const double value = minX + i * xStep;
      const double x = plotArea.left() + (value - minX) / (maxX - minX) * plotArea.width();
      painter.setPen(QPen(colors.color(QPalette::AlternateBase), 1.0));
      painter.drawLine(QPointF(x, plotArea.top()), QPointF(x, plotArea.bottom()));
      painter.setPen(colors.color(QPalette::Text));
      painter.drawText(QRectF(x - 42.0, plotArea.bottom() + 4.0, 84.0, 18.0),
                       Qt::AlignHCenter | Qt::AlignTop, formatTick(value, xStep));
    }
    for (int i = 0; i < yTickCount; ++i) {
      const double value = minY + i * yStep;
      const double y = plotArea.bottom() - (value - minY) / (maxY - minY) * plotArea.height();
      painter.setPen(QPen(colors.color(QPalette::AlternateBase), 1.0));
      painter.drawLine(QPointF(plotArea.left(), y), QPointF(plotArea.right(), y));
      painter.setPen(colors.color(QPalette::Text));
      painter.drawText(QRectF(0.0, y - 9.0, leftMargin - 8.0, 18.0),
                       Qt::AlignRight | Qt::AlignVCenter, formatTick(value, yStep));
    }

    painter.setPen(QPen(colors.color(QPalette::Mid), 1.2));
    painter.drawRect(plotArea);
    painter.setPen(QPen(colors.color(QPalette::Highlight), 2.0));
    QPolygonF polygon;
    for (const auto &point : points) {
      polygon << transform(point);
    }
    painter.drawPolyline(polygon);

    painter.setPen(colors.color(QPalette::Text));
    const auto xLabel = histogramAxisLabel(kind_, true);
    const auto yLabel = normalizedPerBombardment_
                            ? histogramAxisLabel(kind_, false) +
                                  tr(" per bombardment")
                            : histogramAxisLabel(kind_, false);
    painter.drawText(QRectF(leftMargin, canvasSize.height() - 24.0,
                            canvasSize.width() - leftMargin - rightMargin, 20.0),
                     Qt::AlignCenter, xLabel);
    painter.save();
    painter.rotate(-90.0);
    painter.drawText(QRectF(-canvasSize.height(), 0, canvasSize.height(), 20),
                     Qt::AlignCenter, yLabel);
    painter.restore();
    painter.drawText(QRectF(0, 4.0, canvasSize.width(), 20.0), Qt::AlignCenter,
                     title_);
  }

  QString title_;
  HistogramDataKind kind_;
  QLineSeries *series_;
  bool normalizedPerBombardment_ = false;
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

  QString cleanPath = QDir::fromNativeSeparators(normalized);
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
      const auto kind = detectHistogramKind(filePath);
      const auto csvBaseName = fileInfo.completeBaseName();
      auto *plotWidget = new HistogramPlotWidget(filePath, title, kind, target);
      plotWidget->setProperty("graphFileName", csvBaseName);
      if (!fileInfo.exists() || !fileInfo.isFile()) {
        plotWidget->setHistogram(common::Histogram());
      }
      plotWidget->setProperty("histogramKind", static_cast<int>(kind));
      auto *chartPage = new QWidget(target);
      auto *chartLayout = new QVBoxLayout(chartPage);
      auto *chartActions = new QHBoxLayout();
      auto *saveButton = new QPushButton(tr("Save PNG"), chartPage);
      chartActions->addStretch();
      chartActions->addWidget(saveButton);
      chartLayout->addWidget(plotWidget, 1);
      chartLayout->addLayout(chartActions);
      target->addTab(chartPage, title);

      connect(saveButton, &QPushButton::clicked, this,
                [this, plotWidget, csvBaseName]() {
                auto outputDir = outputDirectoryEdit
                                     ? outputDirectoryEdit->text().trimmed()
                                     : QStringLiteral("output");
                if (outputDir.isEmpty()) {
                  outputDir = QStringLiteral("output");
                }
                if (!QDir().mkpath(outputDir)) {
                  QMessageBox::warning(
                      this, tr("Save graph"),
                      tr("Unable to create output directory %1").arg(outputDir));
                  return;
                }
                const auto suggestedPath =
                    QDir(outputDir).absoluteFilePath(
                        csvBaseName + QStringLiteral(".png"));
                auto imagePath = QFileDialog::getSaveFileName(
                    this, tr("Save graph"), suggestedPath,
                    tr("PNG image (*.png)"));
                if (imagePath.isEmpty()) {
                  return;
                }
                const QFileInfo selectedFile(imagePath);
                if (selectedFile.suffix().compare(QStringLiteral("png"),
                                                  Qt::CaseInsensitive) != 0) {
                  imagePath = selectedFile.dir().filePath(
                      selectedFile.completeBaseName() + QStringLiteral(".png"));
                }
                if (!plotWidget->savePng(imagePath)) {
                  QMessageBox::warning(this, tr("Save graph"),
                                       tr("Unable to save graph to %1").arg(imagePath));
                }
              });
    };

    const auto outputDirPath = outputDirectoryEdit ? outputDirectoryEdit->text().trimmed() : QStringLiteral("output");
    const auto fallbackFilePaths = QStringList{
        displacementOutputEdit ? displacementOutputEdit->text().trimmed() : QString(),
        electronicOutputEdit ? electronicOutputEdit->text().trimmed() : QString(),
        nuclearOutputEdit ? nuclearOutputEdit->text().trimmed() : QString(),
        totalOutputEdit ? totalOutputEdit->text().trimmed() : QString(),
        buildOutputPath(rangeOutputEdit ? rangeOutputEdit->text()
                                        : QStringLiteral("range_histogram.csv")),
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
    }

    const auto livePaths = QStringList{
        displacementOutputEdit ? displacementOutputEdit->text().trimmed() : QString(),
        electronicOutputEdit ? electronicOutputEdit->text().trimmed() : QString(),
        nuclearOutputEdit ? nuclearOutputEdit->text().trimmed() : QString(),
        totalOutputEdit ? totalOutputEdit->text().trimmed() : QString(),
        buildOutputPath(rangeOutputEdit ? rangeOutputEdit->text()
                                        : QStringLiteral("range_histogram.csv"))};
    for (const auto &path : livePaths) {
      const auto title = QFileInfo(path).baseName();
      bool exists = false;
      for (int i = 0; i < target->count(); ++i) {
        if (target->tabText(i) == title) {
          exists = true;
          break;
        }
      }
      if (!exists && !title.isEmpty()) {
        addFileChart(title, path);
      }
    }
    if (target->count() > 0) {
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

  if (histogramTabs) {
    populateTabWidget(histogramTabs);
    if (latestLiveTalliesBombardmentCount > 0) {
      refreshLiveHistograms(latestLiveTallies,
                            latestLiveTalliesBombardmentCount);
    }
  }
}

void MainWindow::saveProjection(Canvas *projection, const QString &fileName) {
  if (!projection) {
    return;
  }

  auto outputDir = outputDirectoryEdit
                       ? outputDirectoryEdit->text().trimmed()
                       : QStringLiteral("output");
  if (outputDir.isEmpty()) {
    outputDir = QStringLiteral("output");
  }
  if (!QDir().mkpath(outputDir)) {
    QMessageBox::warning(
        this, tr("Save graph"),
        tr("Unable to create output directory %1").arg(outputDir));
    return;
  }

  auto imagePath = QFileDialog::getSaveFileName(
      this, tr("Save graph"),
      QDir(outputDir).absoluteFilePath(fileName), tr("PNG image (*.png)"));
  if (imagePath.isEmpty()) {
    return;
  }
  const QFileInfo selectedFile(imagePath);
  if (selectedFile.suffix().compare(QStringLiteral("png"),
                                    Qt::CaseInsensitive) != 0) {
    imagePath = selectedFile.dir().filePath(
        selectedFile.completeBaseName() + QStringLiteral(".png"));
  }
  if (!projection->savePng(imagePath)) {
    QMessageBox::warning(this, tr("Save graph"),
                         tr("Unable to save graph to %1").arg(imagePath));
  }
}

void MainWindow::saveGraphsToOutputDirectory() {
  auto outputDir = outputDirectoryEdit
                       ? outputDirectoryEdit->text().trimmed()
                       : QStringLiteral("output");
  if (outputDir.isEmpty()) {
    outputDir = QStringLiteral("output");
  }
  if (!QDir().mkpath(outputDir)) {
    const auto message =
        tr("Unable to create graph output directory %1").arg(outputDir);
    updateResultPane(message);
    QMessageBox::warning(this, tr("Save graphs"), message);
    return;
  }

  QStringList failedImages;
  int savedCount = 0;
  auto saveProjectionIfEnabled =
      [&](Canvas *projection, QCheckBox *check, const QString &fileName) {
    if (!projection || !check || !check->isChecked()) {
      return;
    }
    const auto imagePath = QDir(outputDir).filePath(fileName);
    if (projection->savePng(imagePath)) {
      ++savedCount;
    } else {
      failedImages.append(imagePath);
    }
  };
  saveProjectionIfEnabled(xyCanvas, xyProjectionCheck,
                          QStringLiteral("xy_projection.png"));
  saveProjectionIfEnabled(yzCanvas, yzProjectionCheck,
                          QStringLiteral("yz_projection.png"));
  saveProjectionIfEnabled(xzCanvas, xzProjectionCheck,
                          QStringLiteral("xz_projection.png"));

  if (histogramTabs) {
    for (int tabIndex = 0; tabIndex < histogramTabs->count(); ++tabIndex) {
      auto *page = histogramTabs->widget(tabIndex);
      auto *plotWidget =
          page ? page->findChild<HistogramPlotWidget *>() : nullptr;
      if (!plotWidget || !plotWidget->isEnabled()) {
        continue;
      }

      auto graphFileName =
          plotWidget->property("graphFileName").toString().trimmed();
      if (graphFileName.isEmpty()) {
        graphFileName = histogramTabs->tabText(tabIndex);
      }
      const auto imagePath =
          QDir(outputDir).filePath(graphFileName + QStringLiteral(".png"));
      if (plotWidget->savePng(imagePath)) {
        ++savedCount;
      } else {
        failedImages.append(imagePath);
      }
    }
  }

  if (savedCount > 0) {
    updateResultPane(
        tr("Saved %1 graph PNG(s) to %2").arg(savedCount).arg(outputDir));
  }
  if (!failedImages.isEmpty()) {
    const auto message =
        tr("Unable to save graph images:\n%1").arg(failedImages.join('\n'));
    updateResultPane(message);
    QMessageBox::warning(this, tr("Save graphs"), message);
  }
}

void MainWindow::refreshLiveHistograms(const common::Tallies &tallies,
                                      std::size_t bombardmentCount) {
  if (!histogramTabs) {
    return;
  }
  for (auto *plot : histogramTabs->findChildren<HistogramPlotWidget *>()) {
    const auto kind = static_cast<HistogramDataKind>(
        plot->property("histogramKind").toInt());
    switch (kind) {
      case HistogramDataKind::Displacement:
        plot->setNormalizedHistogram(tallies.displacement_histogram,
                                     bombardmentCount);
        break;
      case HistogramDataKind::ElectronicLoss:
        plot->setNormalizedHistogram(
            tallies.electronic_energy_loss_histogram, bombardmentCount);
        break;
      case HistogramDataKind::NuclearLoss:
        plot->setNormalizedHistogram(tallies.nuclear_energy_loss_histogram,
                                     bombardmentCount);
        break;
      case HistogramDataKind::TotalLoss:
        plot->setNormalizedHistogram(tallies.total_energy_loss_histogram,
                                     bombardmentCount);
        break;
      case HistogramDataKind::Range:
        plot->setNormalizedHistogram(tallies.range_histogram,
                                     bombardmentCount);
        break;
      case HistogramDataKind::Generic:
        break;
    }
  }
}

void MainWindow::refreshSummaryView() {
  if (!simulationSummaryView) {
    return;
  }

  QStringList summaryPaths;
  summaryPaths << (displacementOutputEdit ? displacementOutputEdit->text().trimmed() : QString())
               << (electronicOutputEdit ? electronicOutputEdit->text().trimmed() : QString())
               << (nuclearOutputEdit ? nuclearOutputEdit->text().trimmed() : QString())
               << (totalOutputEdit ? totalOutputEdit->text().trimmed() : QString());
  summaryPaths << buildOutputPath(rangeOutputEdit ? rangeOutputEdit->text()
                                                  : QStringLiteral("range_histogram.csv"));
  const auto outputDir = outputDirectoryEdit ? outputDirectoryEdit->text().trimmed() : QStringLiteral("output");
  const auto summaryPalette = simulationSummaryView->palette();
  simulationSummaryView->setHtml(buildSimulationSummaryHtml(
      outputDir, summaryPaths,
      logFileEdit ? logFileEdit->text().trimmed() : QStringLiteral("amcset.log"),
      enableLoggingCheck && enableLoggingCheck->isChecked(),
      summaryPalette.color(QPalette::Base),
      summaryPalette.color(QPalette::Text),
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

void MainWindow::refreshResultsPane() {
  QStringList summaryPaths;
  summaryPaths << (displacementOutputEdit ? displacementOutputEdit->text().trimmed() : QString())
               << (electronicOutputEdit ? electronicOutputEdit->text().trimmed() : QString())
               << (nuclearOutputEdit ? nuclearOutputEdit->text().trimmed() : QString())
               << (totalOutputEdit ? totalOutputEdit->text().trimmed() : QString());

  const auto outputDir = outputDirectoryEdit ? outputDirectoryEdit->text().trimmed() : QStringLiteral("output");
  refreshSummaryView();

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
  summary << "Angle divisors: " + QString::number(angleDivisorSpinBox ? angleDivisorSpinBox->value() : 1000);
  summary << "Flying distances: " + QString::number(flyingDistanceSpinBox ? flyingDistanceSpinBox->value() : 1000);
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
  auto *generalTab = new QWidget(settingsTabs);
  auto *histogramsTab = new QWidget(settingsTabs);
  auto *electronTab = new QWidget(settingsTabs);
  auto *preferencesTab = new QWidget(settingsTabs);
  auto *settingsLayout = new QVBoxLayout(generalTab);
  auto *formLayout = new QFormLayout();
  formLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  formLayout->setLabelAlignment(Qt::AlignLeft);
  auto *histogramLayout = new QFormLayout(histogramsTab);
  histogramLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  histogramLayout->setLabelAlignment(Qt::AlignLeft);

  auto *electronSettingsGroup = new QGroupBox("Electron transport", electronTab);
  auto *electronSettingsLayout = new QFormLayout(electronSettingsGroup);
  electronSettingsLayout->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
  electronSettingsLayout->setLabelAlignment(Qt::AlignLeft);

  auto *outputDirRow = new QWidget(generalTab);
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
      {&displacementOutputEdit, "Displacement output file"},
      {&electronicOutputEdit, "Electronic loss output file"},
      {&nuclearOutputEdit, "Nuclear loss output file"},
      {&totalOutputEdit, "Total loss output file"},
      {&rangeOutputEdit, "Range output file"},
  };

  // Per-bombardment CSV output
  auto *perBombardmentRow = new QWidget(generalTab);
  auto *perBombardmentLayout = new QHBoxLayout(perBombardmentRow);
  perBombardmentLayout->setContentsMargins(0,0,0,0);
  perBombardmentOutputDirEdit = new QLineEdit("bombardments/", perBombardmentRow);
  auto *perBombardmentBrowse = new QPushButton("Browse", perBombardmentRow);
  connect(perBombardmentBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputDirectory(perBombardmentOutputDirEdit, "Select per-bombardment output directory");
  });
  perBombardmentLayout->addWidget(perBombardmentOutputDirEdit);
  perBombardmentLayout->addWidget(perBombardmentBrowse);
  formLayout->addRow("Per-bombardment output directory", perBombardmentRow);

  savePerBombardmentCheck = new QCheckBox(generalTab);
  savePerBombardmentCheck->setChecked(false);
  formLayout->addRow("Save per-bombardment CSVs", savePerBombardmentCheck);

  // Show a warning when the user enables the option from the settings UI.
  // Use "clicked" so programmatic changes (loading settings at startup)
  // do not trigger the warning — only user clicks do.
  connect(savePerBombardmentCheck, &QCheckBox::clicked, this, [this](bool checked) {
    if (!checked || (app_ && !app_->showLargeCsvWarning())) return;
    QMessageBox warning(QMessageBox::Warning, "Large Output Warning",
      "Saving per-bombardment CSVs can generate very large files. It is not recommended for more than 100 bombardments.",
      QMessageBox::Ok, this);
    auto *dontTellAgain = new QCheckBox(tr("Do not tell me again"), &warning);
    warning.setCheckBox(dontTellAgain);
    warning.exec();
    if (app_ && dontTellAgain->isChecked()) {
      app_->setShowLargeCsvWarning(false);
      if (showLargeCsvWarningCheck) showLargeCsvWarningCheck->setChecked(false);
    }
  });

  for (const auto &entry : fileFields) {
    auto **field = entry.first;
    const char *label = entry.second;
    auto *row = new QWidget(histogramsTab);
    auto *layout = new QHBoxLayout(row);
    layout->setContentsMargins(0, 0, 0, 0);
    *field = new QLineEdit(QString(QString(label).toLower().replace(" ", "_")) + ".csv", row);
    auto *browseButton = new QPushButton("Browse", row);
    connect(browseButton, &QPushButton::clicked, this, [this, field, label]() {
      browseOutputFile(*field, QString("Select %1").arg(label), "csv");
    });
    layout->addWidget(*field);
    layout->addWidget(browseButton);
    histogramLayout->addRow(label, row);
  }

  displacementOutputEdit->setText("displacement_histogram.csv");
  electronicOutputEdit->setText("electronic_loss_histogram.csv");
  nuclearOutputEdit->setText("nuclear_loss_histogram.csv");
  totalOutputEdit->setText("total_loss_histogram.csv");
  rangeOutputEdit->setText("range_histogram.csv");
  if (perBombardmentOutputDirEdit) perBombardmentOutputDirEdit->setText("bombardments/");
  if (savePerBombardmentCheck) savePerBombardmentCheck->setChecked(false);

  auto *particleTreeRow = new QWidget(generalTab);
  auto *particleTreeLayout = new QHBoxLayout(particleTreeRow);
  particleTreeLayout->setContentsMargins(0, 0, 0, 0);
  particleTreeEdit = new QLineEdit("particle_tree.json", particleTreeRow);
  auto *particleTreeBrowse = new QPushButton("Browse", particleTreeRow);
  connect(particleTreeBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(particleTreeEdit, "Select particle tree output file", "json");
  });
  particleTreeLayout->addWidget(particleTreeEdit);
  particleTreeLayout->addWidget(particleTreeBrowse);
  formLayout->addRow("Particle tree output file", particleTreeRow);

  auto addBinCount = [histogramsTab, histogramLayout](QSpinBox **field,
                                                       const QString &label) {
    *field = new QSpinBox(histogramsTab);
    (*field)->setRange(1, 1000000);
    (*field)->setValue(100);
    histogramLayout->addRow(label, *field);
  };
  addBinCount(&displacementBinCountSpinBox, tr("Displacement bins"));
  addBinCount(&electronicLossBinCountSpinBox, tr("Electronic loss bins"));
  addBinCount(&nuclearLossBinCountSpinBox, tr("Nuclear loss bins"));
  addBinCount(&totalLossBinCountSpinBox, tr("Total loss bins"));
  addBinCount(&rangeBinCountSpinBox, tr("Range bins"));

  auto *logFileRow = new QWidget(generalTab);
  auto *logFileLayout = new QHBoxLayout(logFileRow);
  logFileLayout->setContentsMargins(0, 0, 0, 0);
  logFileEdit = new QLineEdit("amcset.log", logFileRow);
  auto *logBrowse = new QPushButton("Browse", logFileRow);
  connect(logBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(logFileEdit, "Select log file", "log");
  });
  logFileLayout->addWidget(logFileEdit);
  logFileLayout->addWidget(logBrowse);
  formLayout->addRow("Log file", logFileRow);

  enableLoggingCheck = new QCheckBox(generalTab);
  enableLoggingCheck->setChecked(false);
  formLayout->addRow("Enable logging", enableLoggingCheck);

  angleDivisorSpinBox = new QSpinBox(electronTab);
  angleDivisorSpinBox->setRange(1, 1000000);
  angleDivisorSpinBox->setValue(1000);
  electronSettingsLayout->addRow("Angle divisors", angleDivisorSpinBox);

  flyingDistanceSpinBox = new QSpinBox(electronTab);
  flyingDistanceSpinBox->setRange(1, 1000000);
  flyingDistanceSpinBox->setValue(1000);
  electronSettingsLayout->addRow("Flying distances", flyingDistanceSpinBox);

  auto *electronScreeningRow = new QWidget(electronTab);
  auto *electronScreeningLayout = new QHBoxLayout(electronScreeningRow);
  electronScreeningLayout->setContentsMargins(0, 0, 0, 0);
  electronScreeningParametersEdit = new QLineEdit("common/assets/electron_screening_parameters.csv", electronTab);
  auto *electronScreeningBrowse = new QPushButton("Browse", electronTab);
  connect(electronScreeningBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(electronScreeningParametersEdit, "Select electron screening parameter file", "csv");
  });
  electronScreeningLayout->addWidget(electronScreeningParametersEdit);
  electronScreeningLayout->addWidget(electronScreeningBrowse);
  electronSettingsLayout->addRow("Electron screening parameters", electronScreeningRow);

  auto *mottRow = new QWidget(electronTab);
  auto *mottLayout = new QHBoxLayout(mottRow);
  mottLayout->setContentsMargins(0, 0, 0, 0);
  mottScatteringParametersEdit = new QLineEdit("common/assets/mott_scattering_parameters.csv", electronTab);
  auto *mottBrowse = new QPushButton("Browse", electronTab);
  connect(mottBrowse, &QPushButton::clicked, this, [this]() {
    browseOutputFile(mottScatteringParametersEdit, "Select Mott scattering parameter file", "csv");
  });
  mottLayout->addWidget(mottScatteringParametersEdit);
  mottLayout->addWidget(mottBrowse);
  electronSettingsLayout->addRow("Mott scattering parameters", mottRow);

  auto *screeningRow = new QWidget(electronTab);
  auto *screeningLayout = new QHBoxLayout(screeningRow);
  screeningLayout->setContentsMargins(0, 0, 0, 0);
  screeningModelCombo = new QComboBox(electronTab);
  screeningModelCombo->addItem("Thomas-Fermi (default)");
  screeningModelCombo->addItem("QuickBasic (ion.bas)");
  screeningLayout->addWidget(screeningModelCombo);
  electronSettingsLayout->addRow("Screening model", screeningRow);

  auto *preferencesGroup = new QGroupBox("Preferences", preferencesTab);
  auto *preferencesLayout = new QFormLayout(preferencesGroup);
  themeModeCombo = new QComboBox(preferencesGroup);
  themeModeCombo->addItem(tr("System"), static_cast<int>(Application::ThemeMode::System));
  themeModeCombo->addItem(tr("Light"), static_cast<int>(Application::ThemeMode::Light));
  themeModeCombo->addItem(tr("Dark"), static_cast<int>(Application::ThemeMode::Dark));
  if (app_) {
    const auto themeIndex =
        themeModeCombo->findData(static_cast<int>(app_->getThemeMode()));
    themeModeCombo->setCurrentIndex(themeIndex < 0 ? 0 : themeIndex);
  }
  connect(themeModeCombo, qOverload<int>(&QComboBox::currentIndexChanged),
          this, [this](int index) {
            if (app_ && themeModeCombo) {
              app_->setThemeMode(static_cast<Application::ThemeMode>(
                  themeModeCombo->itemData(index).toInt()));
              refreshSummaryView();
            }
          });
  preferencesLayout->addRow(tr("Color theme"), themeModeCombo);

  saveAllGraphsCheck = new QCheckBox(preferencesGroup);
  saveAllGraphsCheck->setToolTip(
      tr("Save a PNG for each enabled graph to the output directory after a "
         "simulation completes."));
  preferencesLayout->addRow(
      tr("Automatically save graphs after simulation"), saveAllGraphsCheck);

  showLargeCsvWarningCheck = new QCheckBox(preferencesGroup);
  showLargeCsvWarningCheck->setChecked(!app_ || app_->showLargeCsvWarning());
  connect(showLargeCsvWarningCheck, &QCheckBox::toggled, this, [this](bool show) {
    if (app_) app_->setShowLargeCsvWarning(show);
  });
  preferencesLayout->addRow(tr("Show large CSV warning"), showLargeCsvWarningCheck);

  showOutputOverwriteWarningCheck = new QCheckBox(preferencesGroup);
  showOutputOverwriteWarningCheck->setChecked(!app_ || app_->showOutputOverwriteWarning());
  connect(showOutputOverwriteWarningCheck, &QCheckBox::toggled, this, [this](bool show) {
    if (app_) app_->setShowOutputOverwriteWarning(show);
  });
  preferencesLayout->addRow(tr("Show output directory overwrite warning"),
                            showOutputOverwriteWarningCheck);

  settingsLayout->addLayout(formLayout);
  settingsLayout->addStretch();
  auto *electronTabLayout = new QVBoxLayout(electronTab);
  electronTabLayout->addWidget(electronSettingsGroup);
  electronTabLayout->addStretch();
  auto *preferencesTabLayout = new QVBoxLayout(preferencesTab);
  preferencesTabLayout->addWidget(preferencesGroup);
  preferencesTabLayout->addStretch();
  settingsTabs->addTab(generalTab, tr("General"));
  settingsTabs->addTab(histogramsTab, tr("Histograms"));
  settingsTabs->addTab(electronTab, tr("Electron transport"));
  settingsTabs->addTab(preferencesTab, tr("Preferences"));

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
