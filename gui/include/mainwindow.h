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

// Qt includes
#include <QCheckBox>
#include <QContextMenuEvent>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QStandardItemModel>
#include <QTabWidget>
#include <QTableView>
#include <QTextBrowser>
#include <QUndoStack>
#include <QVBoxLayout>
#include <QtGraphs/QLineSeries>
#include <deque>
#include <atomic>
#include <mutex>
#include <qnamespace.h>
#include <qpainter.h>
#include <qthread.h>

// Local includes
#include "3dview.h"
#include "amcset_common.h"
#include "application.h"
#include "canvas.h"
#include "counter.h"
#include "layerwindow.h"
#include "simulation.h"

namespace amcset {
namespace gui {

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MainWindow(Application *app);
  ~MainWindow() override;

protected:
#ifndef QT_NO_CONTEXTMENU
  void contextMenuEvent(QContextMenuEvent *event) override;
#endif // QT_NO_CONTEXTMENU

private slots:
  void newFile();
  void open();
  void save();
  void print();
  void undo();
  void redo();
  void cut();
  void copy();
  void paste();
  void bold();
  void italic();
  void leftAlign();
  void rightAlign();
  void justify();
  void center();
  void setLineSpacing();
  void setParagraphSpacing();
  void about();
  void aboutQt();
  void openLayerWindow(size_t row);
  void handleLayerButtonClicked(size_t row);
  void toggleIncidentIonFormVisibility(int state);
  void createLayerRow();
  void runSimulation();
  void cancelSimulation();

private:
  void createActions();
  void createMenus();
  void createSettingsDialog();
  void openSettingsDialog();
  void browseOutputDirectory(QLineEdit *lineEdit, const QString &title);
  void browseOutputFile(QLineEdit *lineEdit, const QString &title,
                       const QString &defaultSuffix = QString());
  void populateResultCharts();
  void drainTrackQueue();
  QString buildOutputPath(const QString &filename) const;
  void updateResultPane(const QString &message);
  void refreshResultsPane();
  QGroupBox *createIncidentParticleForms(QWidget *parent = nullptr);
  QGroupBox *createLayerForms(QWidget *parent = nullptr);
  double calculateLayerDensity();
  void update2DViewMode(Qt::CheckState state);
  void update3DViewMode(Qt::CheckState state);
  void launchSimulationThread(std::unique_ptr<common::Simulation> simulation);
  void handleSimulationResult(int result);

  Application *app_;
  QMenu *fileMenu;
  QMenu *editMenu;
  QMenu *formatMenu;
  QMenu *helpMenu;
  QMenu *settingsMenu;
  QActionGroup *alignmentGroup;
  QAction *settingsAction;
  QAction *newAct;
  QAction *openAct;
  QAction *saveAct;
  QAction *printAct;
  QAction *exitAct;
  QAction *undoAct;
  QAction *redoAct;
  QAction *cutAct;
  QAction *copyAct;
  QAction *pasteAct;
  QAction *boldAct;
  QAction *italicAct;
  QAction *leftAlignAct;
  QAction *rightAlignAct;
  QAction *justifyAct;
  QAction *centerAct;
  QAction *setLineSpacingAct;
  QAction *setParagraphSpacingAct;
  QAction *aboutAct;
  QAction *aboutQtAct;
  QLabel *infoLabel = nullptr;

  // Incident particle forms
  QCheckBox *incidentElectron = nullptr;
  QCheckBox *damageCascade = nullptr;
  QDoubleSpinBox *incidentEnergy = nullptr;
  QDoubleSpinBox *electronStoppingEnergySpinBox = nullptr;
  QSpinBox *electronAngleDivisorSpinBox = nullptr;
  QSpinBox *electronFlyingDistanceSpinBox = nullptr;
  QLineEdit *electronScreeningFileEdit = nullptr;
  QLineEdit *mottScatteringFileEdit = nullptr;
  QSpinBox *electronScreeningCountSpinBox = nullptr;
  QSpinBox *mottScatteringCountSpinBox = nullptr;
  QCheckBox *logSingleDisplacementCheck = nullptr;
  QCheckBox *logEndOfFlyingDistanceCheck = nullptr;
  QCheckBox *logStoppingPointCheck = nullptr;
  QSpinBox *incidentMass = nullptr;
  QSpinBox *incidentCharge = nullptr;
  QDoubleSpinBox *incidentAngle = nullptr;
  QFormLayout *incidentFormLayout = nullptr;
  size_t incidentChargeFormIndex;
  size_t incidentMassFormIndex;
  size_t electronScreeningFileFormIndex;
  size_t mottScatteringFileFormIndex;
  size_t electronScreeningCountFormIndex;
  size_t mottScatteringCountFormIndex;
  size_t electronStoppingEnergyFormIndex;
  size_t electronAngleDivisorFormIndex;
  size_t electronFlyingDistanceFormIndex;
  size_t logSingleDisplacementFormIndex;
  size_t logEndOfFlyingDistanceFormIndex;
  size_t logStoppingPointFormIndex;
  QSpinBox *bombardmentCountSpinBox = nullptr;
  QSpinBox *threadCountSpinBox = nullptr;
  QDoubleSpinBox *rangeSpinBox = nullptr;
  QTabWidget *mainTabs = nullptr;
  QDialog *settingsDialog = nullptr;
  QTabWidget *resultsTabs = nullptr;
  QTabWidget *histogramTabs = nullptr;
  QTextBrowser *simulationSummaryView = nullptr;
  QLineEdit *outputDirectoryEdit = nullptr;
  QLineEdit *logFileEdit = nullptr;
  QLineEdit *displacementOutputEdit = nullptr;
  QLineEdit *electronicOutputEdit = nullptr;
  QLineEdit *nuclearOutputEdit = nullptr;
  QLineEdit *totalOutputEdit = nullptr;
  QLineEdit *particleTreeEdit = nullptr;
  QCheckBox *enableLoggingCheck = nullptr;
  QPlainTextEdit *resultPane = nullptr;
  QGroupBox *electronSettingsGroup = nullptr;

  // View controls
  QCheckBox *use3DView = nullptr;
  QCheckBox *use2DView = nullptr;
  QCheckBox *viewModeToggle = nullptr;
  QGroupBox *canvasBox = nullptr;
  QGroupBox *view3DBox = nullptr;
  Counter *counter = nullptr;
  QLabel *counterLabel = nullptr;

  // Layer forms
  QStandardItemModel *layerModel = nullptr;
  QTableView *layerTable = nullptr;
  QStringList layerHeaders;

  // For particle track painting
  Canvas *canvas = nullptr;
  Widget3DView *view3D = nullptr;
  QPainter *painter = nullptr;
  // boost::lockfree::spsc_queue<common::Simulation::DataVariant *> &queue;

  common::Simulation *activeSimulation = nullptr;
  SimulationThread *simulationThread = nullptr;
  QPushButton *cancelButton = nullptr;
  std::atomic<bool> simulationFinished_{true};
  struct CompletedTrack {
    std::vector<common::Coordinate> points;
    size_t bombardmentId = 0;
  };
  std::mutex trackQueueMutex;
  std::deque<CompletedTrack> completedTrackQueue;
  // For Undo/Redo
  // QUndoStack *undoStack; // TODO: Implement undo/redo

public slots:
  void updateLayerDepths();

private slots:

  void onOkClicked() {
    QMessageBox::information(this, "Info",
                             "OK button was clicked in the subwindow!");
  }

  void onCancelClicked() {
    QMessageBox::information(this, "Info",
                             "Cancel button was clicked in the subwindow!");
  }
};

} // namespace gui
} // namespace amcset
