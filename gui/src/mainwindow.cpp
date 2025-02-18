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

// Qt includes
#include <QFormLayout>

#include <QMenu>
#include <QMenuBar>
#include <QSplitter>
#include <QStatusBar>
#include <QToolBar>
#include <QVBoxLayout>
#include <QtWidgets>
// #include <QBoxLayout>
#include <QCheckBox>
#include <QFormLayout>
#include <qcontainerfwd.h>
#include <qdebug.h>
#include <qtableview.h>
#include <qtablewidget.h>

// Local includes
#include "application.h"
#include "canvas.h"
#include "layerwindow.h"
#include "mainwindow.h"

namespace amcset {
namespace gui {

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

  auto *runButton = new QPushButton("Run Simulation");
  runButton->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
  connect(runButton, &QPushButton::clicked, app_, &Application::runSimulation);

  auto *canvas = new Canvas;
  canvas->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Expanding);

  // Layouts
  auto *topLayout = new QVBoxLayout; // The outermost layout.
  topLayout->setContentsMargins(5, 5, 5, 5);

  auto *formBox = new QGroupBox(""); // Contains forms for user input.
  auto *formLayout = new QVBoxLayout(formBox);

  formLayout->addWidget(infoLabel);
  formLayout->addWidget(
      createIncidentParticleForms(formBox)); // Layout within the widget.
  formLayout->addWidget(createLayerForms());
  formLayout->addWidget(runButton);

  auto *hPane = new QSplitter(Qt::Horizontal); // Split the forms and canvas
  hPane->addWidget(formBox);
  hPane->addWidget(canvas);
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

  incidentElectron = new QCheckBox();
  incidentElectron->setChecked(false); // Default to ion simulation
  incidentFormLayout->addRow("Incident electron?", incidentElectron);

  damageCascade = new QCheckBox();
  damageCascade->setChecked(true); // Default to ion simulation
  incidentFormLayout->addRow("Simulate damage cascade?", damageCascade);

  incidentEnergy = new QDoubleSpinBox();
  incidentEnergy->setRange(0.0,
                           10000.0); // TODO: Ensure that range is reasonable;
  incidentEnergy->setDecimals(1);
  incidentFormLayout->addRow("Incident Energy (keV)", incidentEnergy);

  incidentMass = new QSpinBox();
  incidentMass->setRange(0, 295);
  incidentMassFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Mass (amu)", incidentMass);

  incidentCharge = new QSpinBox();
  incidentCharge->setRange(0, 118);
  incidentChargeFormIndex = incidentFormLayout->rowCount();
  incidentFormLayout->addRow("Incident Z-number", incidentCharge);

  connect(incidentElectron, &QCheckBox::stateChanged, this,
          &MainWindow::toggleIncidentIonFormVisibility);

  incidentAngle = new QDoubleSpinBox();
  incidentAngle->setRange(0.0, 90.0);
  incidentAngle->setDecimals(2);
  incidentFormLayout->addRow("Incident Angle (degree)", incidentAngle);

  return formGroupBox;
}

void MainWindow::toggleIncidentIonFormVisibility(int state) {
  incidentFormLayout->setRowVisible(incidentMassFormIndex,
                                    state != Qt::Checked);
  incidentFormLayout->setRowVisible(incidentChargeFormIndex,
                                    state != Qt::Checked);
  incidentMass->setVisible(state != Qt::Checked);
}

void handleButtonClick(size_t row) {
  qDebug() << "Button in row " << row << " clicked!";
}

QGroupBox *MainWindow::createLayerForms(QWidget *parent) {
  auto layerGroupBox = new QGroupBox("Layers", parent);
  auto layerVBox = new QVBoxLayout(layerGroupBox);
  auto addLayerButton = new QPushButton("Add Layer");
  connect(addLayerButton, &QPushButton::clicked, this,
          &MainWindow::openLayerWindow);

  layerModel = new QStandardItemModel(3, 6);
  layerModel->setHorizontalHeaderLabels(
      QStringList() << "Index" << "Layer" << "Density" << "Composition"
                    << "Phase" << "Edit");
  layerTable = new QTableView();
  layerTable->setModel(layerModel);
  layerTable->horizontalHeader()->setSectionResizeMode(
      layerModel->columnCount() - 1, QHeaderView::ResizeToContents);

  for (size_t row = 0; row < (size_t)layerModel->rowCount(); row++) {
    auto *button = new QPushButton();
    button->setIcon(QApplication::style()->standardIcon(QStyle::SP_FileIcon));
    button->setToolTip("Edit");
    connect(button, &QPushButton::clicked, this,
            [row]() { handleButtonClick(row); });
    layerTable->setIndexWidget(layerModel->index(row, 5), button);
  }

  layerVBox->addWidget(addLayerButton);
  layerVBox->addWidget(layerTable);

  return layerGroupBox;
}

void MainWindow::openLayerWindow() {
  // Create and show the layerWindow
  LayerWindow *layerWindow = new LayerWindow(this); // Pass 'this' as the parent
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
