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

// Subwindow class
#include <QApplication>
#include <QMessageBox>
#include <QPushButton>
#include <QVBoxLayout>
#include <QWidget>

class LayerWindow : public QWidget {
  Q_OBJECT

public:
  explicit LayerWindow(QWidget *parent = nullptr);
signals:
  void okClicked();     // Signal emitted when "OK" is clicked
  void cancelClicked(); // Signal emitted when "Cancel" is clicked

private slots:
  void onOkClicked() {
    emit okClicked(); // Emit the okClicked signal
    close();          // Close the window
  }

  void onCancelClicked() {
    emit cancelClicked(); // Emit the cancelClicked signal
    close();              // Close the window
  }
};
