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

#include "app/guard.h"
#include "setup/ion_panel.h"

#include "element_data.h"
#include "model/quantity_text.h"
#include "setup/form_binder.h"
#include "setup/periodic_table.h"
#include "trim/trim_dat.h"

#include <QButtonGroup>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QRadioButton>
#include <QToolButton>
#include <QVBoxLayout>

namespace amcset {
namespace gui {

using common::Element;

IonPanel::IonPanel(Scenario *scenario, QWidget *parent)
    : QWidget(parent), scenario_(scenario), binder_(new FormBinder(scenario, this)) {
  auto *ion = new QGroupBox(tr("Incident ion"), this);
  element_ = new QComboBox(ion);
  for (unsigned z = 1; z <= Element::max_z; ++z) {
    const auto &el = Element::by_z(z);
    element_->addItem(QStringLiteral("%1  %2 (%3)")
                          .arg(QString::fromStdString(std::string(el.symbol())))
                          .arg(QString::fromStdString(std::string(el.name())))
                          .arg(z),
                      z);
  }
  element_->setMaxVisibleItems(20);
  element_->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
  element_->setMinimumContentsLength(14);
  auto *pt = new QToolButton(ion);
  pt->setText(tr("PT"));
  pt->setToolTip(tr("Choose from the periodic table"));
  auto *el_row = new QHBoxLayout();
  el_row->addWidget(element_, 1);
  el_row->addWidget(pt);
  mass_ = new QLineEdit(ion);
  mass_->setToolTip(tr("Ion mass; the default is the most abundant isotope (TRIM's default)"));
  energy_ = new QLineEdit(ion);
  energy_->setToolTip(tr("Kinetic energy, e.g. \"1 MeV\", \"250 keV\" (a bare number is keV)"));
  angle_ = new QLineEdit(ion);
  angle_->setToolTip(tr("Angle of incidence to the surface normal, tilted in the depth-Y plane (0-89.9 deg)"));
  auto *f = new QFormLayout(ion);
  f->addRow(tr("Element"), el_row);
  f->addRow(tr("Mass"), mass_);
  f->addRow(tr("Energy"), energy_);
  f->addRow(tr("Angle"), angle_);

  auto *src = new QGroupBox(tr("Particle source"), this);
  source_ = new QButtonGroup(this);
  auto *beam = new QRadioButton(tr("Ion beam"), src);
  beam->setToolTip(tr("Every ion is the ion above, starting at the surface"));
  auto *ions = new QRadioButton(tr("Ions from a list (TRIM.DAT)"), src);
  ions->setToolTip(tr("Each line of the list starts one ion with its own energy, position and direction"));
  auto *recoils = new QRadioButton(tr("Recoil atoms from a list (TRIM.DAT)"), src);
  recoils->setToolTip(tr("Each line sets one target atom in motion, without an ion: damage from neutrons, "
                         "electrons or photons computed by another program"));
  source_->addButton(beam, 0);
  source_->addButton(ions, 1);
  source_->addButton(recoils, 2);
  particles_file_ = new QLineEdit(src);
  particles_file_->setReadOnly(true);
  particles_file_->setPlaceholderText(tr("no list loaded"));
  load_ = new QPushButton(tr("Load..."), src);
  load_->setToolTip(tr("TRIM.DAT, or BACKSCAT.txt / TRANSMIT.txt / SPUTTER.txt / TRIMOUT.txt"));
  particles_info_ = new QLabel(src);
  particles_info_->setWordWrap(true);
  auto *file_row = new QHBoxLayout();
  file_row->addWidget(particles_file_, 1);
  file_row->addWidget(load_);
  auto *sl = new QVBoxLayout(src);
  sl->addWidget(beam);
  sl->addWidget(ions);
  sl->addWidget(recoils);
  sl->addLayout(file_row);
  sl->addWidget(particles_info_);

  auto *lay = new QVBoxLayout(this);
  lay->addWidget(ion);
  lay->addWidget(src);
  lay->addStretch(1);

  // --- bindings ---------------------------------------------------------------
  connect(binder_, &FormBinder::error, this, guarded(this, tr("Ion settings (error)"), [this](auto &&...a) { message(std::forward<decltype(a)>(a)...); }));
  binder_->combo(element_, tr("Ion element"),
                 [](const config::RunConfig &c) { return static_cast<int>(c.trim.ion.z) - 1; },
                 [](config::RunConfig &c, int i) {
                   const unsigned z = static_cast<unsigned>(i + 1);
                   if (z == c.trim.ion.z) return;
                   c.trim.ion.z = z;
                   c.trim.ion.mass = Element::by_z(z).most_abundant_isotope_mass();
                 });
  connect(pt, &QToolButton::clicked, this, guarded(this, tr("Ion settings (clicked)"), [this] {
    const unsigned z = PeriodicTableDialog::pick(this, scenario_->config().trim.ion.z);
    if (z != 0) element_->setCurrentIndex(static_cast<int>(z) - 1);
  }));
  binder_->text(mass_, tr("Ion mass"),
                [](const config::RunConfig &c) {
                  return QStringLiteral("%1 amu").arg((c.trim.ion.mass / common::atomic_mass_unit).value(), 0, 'g', 8);
                },
                [](config::RunConfig &c, const QString &t) {
                  const auto m = parse_mass(t);
                  if (!(m > 0.0 * common::kilogram)) throw std::invalid_argument("must be positive");
                  c.trim.ion.mass = m;
                });
  binder_->text(energy_, tr("Ion energy"), [](const config::RunConfig &c) { return format_energy(c.trim.ion.energy); },
                [](config::RunConfig &c, const QString &t) {
                  const auto e = parse_energy(t, "keV");
                  if (!(e > 0.0 * common::electron_volt)) throw std::invalid_argument("must be positive");
                  c.trim.ion.energy = e;
                });
  binder_->text(angle_, tr("Angle"),
                [](const config::RunConfig &c) {
                  return QStringLiteral("%1 deg").arg(common::in_units(c.trim.ion.angle, common::degree), 0, 'g', 6);
                },
                [](config::RunConfig &c, const QString &t) {
                  const auto a = parse_angle(t);
                  if (!(a >= 0.0 * common::degree) || !(a < 90.0 * common::degree)) {
                    throw std::invalid_argument("must be in [0, 90) degrees");
                  }
                  c.trim.ion.angle = a;
                });
  binder_->custom([this](const config::RunConfig &c) {
    const int id = static_cast<int>(c.trim.source);
    if (auto *b = source_->button(id)) b->setChecked(true);
    const bool list = c.trim.source != trim::TrimConfig::Source::beam;
    load_->setEnabled(list);
    particles_file_->setEnabled(list);
    if (!list) {
      particles_info_->setText(tr("Every ion starts at the surface point with the energy and angle above."));
    } else if (c.trim.particles.empty()) {
      particles_info_->setText(tr("Load a particle list. The ion energy above must be at least the highest "
                                  "energy in the list (TRIM skips the others)."));
    } else {
      double emax = 0;
      for (const auto &p : c.trim.particles) emax = std::max(emax, common::in_units(p.energy, common::electron_volt));
      particles_info_->setText(
          tr("%L1 particles, up to %2. %3")
              .arg(c.trim.particles.size())
              .arg(format_energy(emax * common::electron_volt))
              .arg(c.trim.source == trim::TrimConfig::Source::particle_ions
                       ? tr("Each line starts one ion (species from above).")
                       : (c.trim.damage == trim::DamageCalculation::monolayer
                              ? tr("Each line moves one target atom like an ion, in monolayer steps (TRIM's code 7).")
                              : tr("Each line starts a cascade at one target atom; there is no ion (TRIM's code 6)."))));
    }
    if (!list) particles_file_->clear();
    else if (particles_file_->text().isEmpty() && !c.trim.particles_title.empty())
      particles_file_->setText(QString::fromStdString(c.trim.particles_title));
    energy_->setToolTip(list ? tr("Highest particle energy allowed (TRIM's maximum energy)")
                             : tr("Kinetic energy, e.g. \"1 MeV\", \"250 keV\" (a bare number is keV)"));
  });
  connect(source_, &QButtonGroup::idClicked, this, guarded(this, tr("Ion settings (idClicked)"), [this](int id) {
    binder_->apply(tr("Particle source"),
                   [&](config::RunConfig &c) {
                     c.trim.source = static_cast<trim::TrimConfig::Source>(id);
                     if (c.trim.source == trim::TrimConfig::Source::particle_ions) {
                       if (c.trim.damage == trim::DamageCalculation::monolayer) {
                         c.trim.damage = trim::DamageCalculation::full_cascade;
                       }
                     }
                   },
                   -1);
    if (id != 0 && scenario_->config().trim.particles.empty()) load_particles();
  }));
  connect(load_, &QPushButton::clicked, this, guarded(this, tr("Ion settings (clicked)"), [this] { load_particles(); }));
  binder_->refresh();
}

void IonPanel::load_particles() {
  const QString path = QFileDialog::getOpenFileName(
      this, tr("Load particle list"), {},
      tr("Particle lists (TRIM.DAT *.dat *.DAT *.txt);;All files (*)"));
  if (path.isEmpty()) return;
  trim::ParticleList list;
  try {
    const QString name = QFileInfo(path).fileName().toUpper();
    const bool exits = name == "BACKSCAT.TXT" || name == "TRANSMIT.TXT" || name == "SPUTTER.TXT" ||
                       name == "TRIMOUT.TXT";
    list = exits ? trim::read_exit_file_as_particles(std::filesystem::path(path.toStdWString()))
                 : trim::read_trim_dat(std::filesystem::path(path.toStdWString()));
  } catch (const std::exception &e) {
    emit message(tr("Cannot read %1: %2").arg(path, QString::fromUtf8(e.what())));
    return;
  }
  if (!list.warnings.empty()) {
    emit message(tr("%1: %2 lines skipped (first: %3)")
                     .arg(QFileInfo(path).fileName())
                     .arg(list.warnings.size())
                     .arg(QString::fromStdString(list.warnings.front())));
  }
  particles_file_->setText(path);
  binder_->apply(tr("Load particle list"),
                 [&](config::RunConfig &c) {
                   c.trim.particles = list.particles;
                   c.trim.particles_title = list.title.empty() ? QFileInfo(path).fileName().toStdString() : list.title;
                   if (c.trim.source == trim::TrimConfig::Source::beam) {
                     c.trim.source = trim::TrimConfig::Source::particle_ions;
                   }
                   c.trim.ion_count = 0;
                   // The ion energy is TRIM's upper limit of the list.
                   common::energy_quantity emax = c.trim.ion.energy;
                   for (const auto &p : list.particles) emax = std::max(emax, p.energy);
                   c.trim.ion.energy = emax;
                 },
                 -1);
}

} // namespace gui
} // namespace amcset
