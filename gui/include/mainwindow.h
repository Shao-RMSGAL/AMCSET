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
#include <QVBoxLayout>

// Local includes
#include "application.h"
#include "layerwindow.h"

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

private:
  void createActions();
  void createMenus();
  QGroupBox *createIncidentParticleForms(QWidget *parent = nullptr);
  QGroupBox *createLayerForms(QWidget *parent = nullptr);
  void openLayerWindow();
  void toggleIncidentIonFormVisibility(int state);

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

  // Layer forms
  QStandardItemModel *layerModel;
  QTableView *layerTable;

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
