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

/*! \file ion_panel.h
 *  \brief The incident particle and the particle source (beam or TRIM.DAT).
 */

#pragma once

#include <QWidget>

class QButtonGroup;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;

namespace amcset {
namespace config {
struct RunConfig;
}
namespace gui {

class FormBinder;
class Scenario;

class IonPanel : public QWidget {
  Q_OBJECT

public:
  IonPanel(Scenario *scenario, QWidget *parent = nullptr);

signals:
  void message(const QString &text);

private:
  void load_particles();

  Scenario *scenario_;
  FormBinder *binder_;
  QComboBox *element_;
  QLineEdit *mass_;
  QLineEdit *energy_;
  QLineEdit *angle_;
  QButtonGroup *source_;
  QLineEdit *particles_file_;
  QLabel *particles_info_;
  QPushButton *load_;
  // One recoil atom (a PKA) at a point: source "recoils" with one particle
  // named "PKA", repeated (TrimConfig::particle_repeat).
  QGroupBox *pka_box_;
  QComboBox *pka_element_;
  QLineEdit *pka_energy_;
  QLineEdit *pka_depth_;
  QLineEdit *pka_theta_;
  QLineEdit *pka_phi_;
  QSpinBox *pka_count_;
};

//! True if \p c runs one recoil atom set up on the Ion tab (a PKA).
bool is_single_recoil(const config::RunConfig &c);

} // namespace gui
} // namespace amcset
