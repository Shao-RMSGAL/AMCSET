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
#include <QDoubleSpinBox>
#include <QFormLayout>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QMainWindow>
#include <QPushButton>
#include <QSpinBox>
#include <QStandardItemModel>
#include <QTableView>
#include <QUndoStack>
#include <QVBoxLayout>
#include <boost/lockfree/spsc_queue.hpp>
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

private:
  void createActions();
  void createMenus();
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
  QActionGroup *alignmentGroup;
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
  QLabel *infoLabel;

  // Incident particle forms
  QCheckBox *incidentElectron;
  QCheckBox *damageCascade;
  QDoubleSpinBox *incidentEnergy;
  QSpinBox *incidentMass;
  QSpinBox *incidentCharge;
  QDoubleSpinBox *incidentAngle;
  QFormLayout *incidentFormLayout;
  size_t incidentChargeFormIndex;
  size_t incidentMassFormIndex;
  QSpinBox *bombardmentCountSpinBox;
  QSpinBox *threadCountSpinBox;
  QDoubleSpinBox *rangeSpinBox;

  // View controls
  QCheckBox *use3DView;
  QCheckBox *use2DView;
  QGroupBox *canvasBox;
  QGroupBox *view3DBox;
  Counter *counter;
  QLabel *counterLabel;

  // Layer forms
  QStandardItemModel *layerModel;
  QTableView *layerTable;
  QStringList layerHeaders;

  // For particle track painting
  Canvas *canvas;
  Widget3DView *view3D;
  QPainter *painter;
  // boost::lockfree::spsc_queue<common::Simulation::DataVariant *> &queue;

  SimulationThread *simulationThread;
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
