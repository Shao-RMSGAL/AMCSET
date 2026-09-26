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

#include <QAction>
#include <QActionGroup>
#include <QApplication>
#include <QContextMenuEvent>
#include <QMenu>
#include <QMessageBox>

namespace amcset {
namespace gui {

#ifndef QT_NO_CONTEXTMENU
void MainWindow::contextMenuEvent(QContextMenuEvent *event) {
  QMenu menu(this);
  menu.addAction(cutAct);
  menu.addAction(copyAct);
  menu.addAction(pasteAct);
  menu.exec(event->globalPos());
}
#endif // QT_NO_CONTEXTMENU

void MainWindow::newFile() {
  statusBar()->showMessage(tr("Invoked File|New"), 3000);
}

void MainWindow::open() { statusBar()->showMessage(tr("Invoked File|Open"), 3000); }

void MainWindow::save() { statusBar()->showMessage(tr("Invoked File|Save"), 3000); }

void MainWindow::print() {
  statusBar()->showMessage(tr("Invoked File|Print"), 3000);
}

void MainWindow::undo() { statusBar()->showMessage(tr("Invoked Edit|Undo"), 3000); }

void MainWindow::redo() { statusBar()->showMessage(tr("Invoked Edit|Redo"), 3000); }

void MainWindow::cut() { statusBar()->showMessage(tr("Invoked Edit|Cut"), 3000); }

void MainWindow::copy() { statusBar()->showMessage(tr("Invoked Edit|Copy"), 3000); }

void MainWindow::paste() {
  statusBar()->showMessage(tr("Invoked Edit|Paste"), 3000);
}

void MainWindow::bold() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Bold"), 3000);
}

void MainWindow::italic() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Italic"), 3000);
}

void MainWindow::leftAlign() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Left Align"), 3000);
}

void MainWindow::rightAlign() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Right Align"), 3000);
}

void MainWindow::justify() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Justify"), 3000);
}

void MainWindow::center() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Center"), 3000);
}

void MainWindow::setLineSpacing() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Set Line Spacing"), 3000);
}

void MainWindow::setParagraphSpacing() {
  statusBar()->showMessage(tr("Invoked Edit|Format|Set Paragraph Spacing"), 3000);
}

void MainWindow::about() {
  statusBar()->showMessage(tr("Invoked Help|About"), 3000);
  QMessageBox::about(this, tr("About Menu"),
                     tr("The <b>Menu</b> example shows how to create "
                        "menu-bar menus and context menus."));
}

void MainWindow::aboutQt() {
  statusBar()->showMessage(tr("Invoked Help|About Qt"), 3000);
}

void MainWindow::createActions() {
  settingsAction = new QAction(tr("Settings"), this);
  settingsAction->setStatusTip(tr("Open simulation settings"));
  connect(settingsAction, &QAction::triggered, this,
          &MainWindow::openSettingsDialog);

  newAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::DocumentNew),
                       tr("&New"), this);
  newAct->setShortcuts(QKeySequence::New);
  newAct->setStatusTip(tr("Create a new file"));
  connect(newAct, &QAction::triggered, this, &MainWindow::newFile);

  openAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::DocumentOpen),
                        tr("&Open..."), this);
  openAct->setShortcuts(QKeySequence::Open);
  openAct->setStatusTip(tr("Open an existing file"));
  connect(openAct, &QAction::triggered, this, &MainWindow::open);

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
  cutAct->setStatusTip(tr("Cut the current selection's contents to the clipboard"));
  connect(cutAct, &QAction::triggered, this, &MainWindow::cut);

  copyAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditCopy),
                        tr("&Copy"), this);
  copyAct->setShortcuts(QKeySequence::Copy);
  copyAct->setStatusTip(tr("Copy the current selection's contents to the clipboard"));
  connect(copyAct, &QAction::triggered, this, &MainWindow::copy);

  pasteAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::EditPaste),
                         tr("&Paste"), this);
  pasteAct->setShortcuts(QKeySequence::Paste);
  pasteAct->setStatusTip(tr("Paste the clipboard's contents into the current selection"));
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
  setLineSpacingAct->setStatusTip(tr("Change the gap between the lines of a paragraph"));
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

  leftAlignAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyLeft),
                  tr("&Left Align"), this);
  leftAlignAct->setCheckable(true);
  leftAlignAct->setShortcut(tr("Ctrl+L"));
  leftAlignAct->setStatusTip(tr("Left align the selected text"));
  connect(leftAlignAct, &QAction::triggered, this, &MainWindow::leftAlign);

  rightAlignAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyRight),
                  tr("&Right Align"), this);
  rightAlignAct->setCheckable(true);
  rightAlignAct->setShortcut(tr("Ctrl+R"));
  rightAlignAct->setStatusTip(tr("Right align the selected text"));
  connect(rightAlignAct, &QAction::triggered, this, &MainWindow::rightAlign);

  justifyAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyFill),
                  tr("&Justify"), this);
  justifyAct->setCheckable(true);
  justifyAct->setShortcut(tr("Ctrl+J"));
  justifyAct->setStatusTip(tr("Justify the selected text"));
  connect(justifyAct, &QAction::triggered, this, &MainWindow::justify);

  centerAct = new QAction(QIcon::fromTheme(QIcon::ThemeIcon::FormatJustifyCenter),
                  tr("&Center"), this);
  centerAct->setCheckable(true);
  centerAct->setShortcut(tr("Ctrl+E"));
  centerAct->setStatusTip(tr("Center the selected text"));
  connect(centerAct, &QAction::triggered, this, &MainWindow::center);

  alignmentGroup = new QActionGroup(this);
  alignmentGroup->addAction(leftAlignAct);
  alignmentGroup->addAction(rightAlignAct);
  alignmentGroup->addAction(justifyAct);
  alignmentGroup->addAction(centerAct);
  leftAlignAct->setChecked(true);
}

void MainWindow::createMenus() {
  fileMenu = menuBar()->addMenu(tr("&File"));
  fileMenu->addAction(newAct);
  fileMenu->addAction(settingsAction);
  fileMenu->addAction(openAct);
  fileMenu->addAction(saveAct);
  fileMenu->addAction(printAct);
  fileMenu->addSeparator();
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
}

} // namespace gui
} // namespace amcset
