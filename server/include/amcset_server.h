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

#include <stdexcept>

#include <QApplication>
#include <QCheckBox>
#include <QFrame>
#include <QLabel>
#include <QListWidget>
#include <QMainWindow>
#include <QMessageBox>
#include <QPushButton>
#include <QRadioButton>
#include <QScreen>

#include <glog/logging.h>

// #include "amcset_common.h"

namespace amcset {
namespace server {
class Application : public QApplication {

public:
  Application(int &argc, char **argv) : QApplication(argc, argv) {
    enableDarkMode();
  }

  void enableLightMode() {
#ifndef Q_OS_MACOS
    LOG(INFO) << "Enabling Light Mode";
    dark_mode = 0;
    qApp->setStyle("Fusion");
    auto lightPalette = QPalette();
    auto lightColor = QColor(255, 255, 255);
    auto disabledColor = QColor(127, 127, 127);
    lightPalette.setColor(QPalette::Window, lightColor);
    lightPalette.setColor(QPalette::WindowText, Qt::black);
    lightPalette.setColor(QPalette::Base, QColor(18, 18, 18));
    lightPalette.setColor(QPalette::AlternateBase, lightColor);
    lightPalette.setColor(QPalette::ToolTipBase, Qt::white);
    lightPalette.setColor(QPalette::ToolTipText, Qt::white);
    lightPalette.setColor(QPalette::Text, Qt::white);
    lightPalette.setColor(QPalette::Disabled, QPalette::Text, disabledColor);
    lightPalette.setColor(QPalette::Button, lightColor);
    lightPalette.setColor(QPalette::ButtonText, Qt::white);
    lightPalette.setColor(QPalette::Disabled, QPalette::ButtonText,
                          disabledColor);
    lightPalette.setColor(QPalette::BrightText, Qt::red);
    lightPalette.setColor(QPalette::Link, QColor(42, 130, 218));
    lightPalette.setColor(QPalette::Highlight, QColor(42, 130, 218));
    lightPalette.setColor(QPalette::HighlightedText, Qt::black);
    lightPalette.setColor(QPalette::Disabled, QPalette::HighlightedText,
                          disabledColor);

    lightPalette.setColor(QPalette::ButtonText, Qt::black);

    qApp->setPalette(lightPalette);
    qApp->setStyleSheet("QToolTip { color: #000000; background-color: #ffffff; "
                        "border: 1px solid black; }");
#else
    LOG(INFO) << "Light mode cannot be set. You are on MacOS."
#endif
  }

  void enableDarkMode() {
#ifndef Q_OS_MACOS
    LOG(INFO) << "Enabling Dark Mode";
    dark_mode = 1;
    qApp->setStyle("Fusion");
    auto darkPalette = QPalette();
    auto darkColor = QColor(38, 38, 38);
    auto disabledColor = QColor(127, 127, 127);
    darkPalette.setColor(QPalette::Window, darkColor);
    darkPalette.setColor(QPalette::WindowText, Qt::white);
    darkPalette.setColor(QPalette::Base, QColor(18, 18, 18));
    darkPalette.setColor(QPalette::AlternateBase, darkColor);
    darkPalette.setColor(QPalette::ToolTipBase, Qt::white);
    darkPalette.setColor(QPalette::ToolTipText, Qt::white);
    darkPalette.setColor(QPalette::Text, Qt::white);
    darkPalette.setColor(QPalette::Disabled, QPalette::Text, disabledColor);
    darkPalette.setColor(QPalette::Button, darkColor);
    darkPalette.setColor(QPalette::ButtonText, Qt::white);
    darkPalette.setColor(QPalette::Disabled, QPalette::ButtonText,
                         disabledColor);
    darkPalette.setColor(QPalette::BrightText, Qt::red);
    darkPalette.setColor(QPalette::Link, QColor(42, 130, 218));
    darkPalette.setColor(QPalette::Highlight, QColor(42, 130, 218));
    darkPalette.setColor(QPalette::HighlightedText, Qt::black);
    darkPalette.setColor(QPalette::Disabled, QPalette::HighlightedText,
                         disabledColor);

    darkPalette.setColor(QPalette::ButtonText, Qt::white);

    qApp->setPalette(darkPalette);
    qApp->setStyleSheet("QToolTip { color: #ffffff; background-color: #2a82da; "
                        "border: 1px solid white; }");
#else
    LOG(INFO) << "Dark mode cannot be set. You are on MacOS."
#endif
  }

  bool dark_mode;
};

class Window1 : public QMainWindow {
  Q_OBJECT
public:
  Application *app;
  Window1(Application *app_in) {

    app = app_in;

    generateHandledExceptionButton.setText("Generated handled exception");
    generateHandledExceptionButton.move(10, 10);
    connect(&generateHandledExceptionButton, &QPushButton::clicked, [&] {
      try {
        throw std::invalid_argument("Exception handled generated");
      } catch (const std::exception &e) {
        QMessageBox::information(this, "Exception handled", e.what());
      }
    });

    generateExceptionButton.setText("Generate exception");
    generateExceptionButton.move(10, 40);
    connect(&generateExceptionButton, &QPushButton::clicked,
            [&] { throw std::invalid_argument("Exception generated"); });

    generateUnknownExceptionButton.setText("Generate unknown exception");
    generateUnknownExceptionButton.move(10, 70);
    connect(&generateUnknownExceptionButton, &QPushButton::clicked, [&] {
      try {
        throw "Unkown exception generated";
      } catch (char const *e) {
        QMessageBox::information(this, "Exception with char", e);
      }
    });

    label1.setText("u da best!");
    label1.setFont({label1.font().family(), 72});
    label1.resize(label1.sizeHint());
    label1.move(0, 400);

    label.move(200, 0);
    label.setText("Dark Mode");

    pushButton.move(210, 10);

    if (app->dark_mode) {
      pushButton.setText("Switch to Light Mode");
    } else {
      pushButton.setText("Switch to Dark Mode");
    }
    connect(&pushButton, &QPushButton::clicked, [&] {
      try {
        if (app->dark_mode) {
          LOG(INFO) << "Changing to Light Mode";
          app->enableLightMode();
          pushButton.setText("Switch to Dark Mode");
        } else {
          LOG(INFO) << "Changing to Dark Mode";
          app->enableDarkMode();
          pushButton.setText("Switch to Light Mode");
        }
      } catch (char const *e) {
        QMessageBox::information(this, "Exception with char", e);
      }
    });

    listWidget.move(200, 50);
    listWidget.resize(120, 100);
    listWidget.addItems({"Item 1", "Item 2", "Item 3", "Item 4", "Item 5",
                         "Item 6", "Item 7", "Item 8", "Item 9", "Item 10"});
    listWidget.setCurrentRow(0);

    radioButton1.move(200, 170);
    radioButton1.setText("Radio 1");
    radioButton1.setChecked(true);

    radioButton2.move(310, 170);
    radioButton2.setText("Radio 2");

    checkBox1.move(200, 200);
    checkBox1.setText("Check 1");
    checkBox1.setChecked(true);

    checkBox2.move(310, 200);
    checkBox2.setText("Check 2");

    setCentralWidget(&frame);
    setWindowTitle("AMCSET");
    // resize(label1.sizeHint());
  }

  QFrame frame;
  QPushButton generateHandledExceptionButton{&frame};
  QPushButton generateExceptionButton{&frame};
  QPushButton generateUnknownExceptionButton{&frame};
  int button1Clicked = 0;
  int button2Clicked = 0;
  QLabel label1{&frame};
  QLabel label{&frame};
  QPushButton pushButton{&frame};
  QListWidget listWidget{&frame};
  QRadioButton radioButton1{&frame};
  QRadioButton radioButton2{&frame};
  QCheckBox checkBox1{&frame};
  QCheckBox checkBox2{&frame};
};
// common::Simulation start_simulation();

} // namespace server
} // namespace amcset
