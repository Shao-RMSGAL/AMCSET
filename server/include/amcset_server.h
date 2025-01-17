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

#include <QApplication>
#include <QFrame>
#include <QLabel>
#include <QMainWindow>
#include <QScreen>

// #include "amcset_common.h"

namespace amcset {
namespace server {

class Window1 : public QMainWindow {
public:
  Window1() {
    label1.setText("\U0001F44B, \U0001F30E\U00002757");
    label1.setFont({label1.font().family(), 72});
    label1.resize(label1.sizeHint());

    setCentralWidget(&frame);
    setWindowTitle("Hello world (emoticons)");
    resize(label1.sizeHint());
  }

private:
  QFrame frame;
  QLabel label1{&frame};
};

// common::Simulation start_simulation();

} // namespace server
} // namespace amcset
