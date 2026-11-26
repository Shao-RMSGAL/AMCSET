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
#include <QSpinBox>
#include <QToolButton>
#include <QVBoxLayout>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

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
  auto *pka = new QRadioButton(tr("One recoil atom in the target (PKA)"), src);
  pka->setToolTip(tr("A target atom set in motion at a point inside the target, without an ion (a primary "
                     "knock-on atom), run as many times as you choose: its cascades are shown in the path "
                     "views and tallied like TRIM's. In a crystal layer it starts from the nearest lattice "
                     "site of its element. (Scenario: source type recoils, one particle named PKA, "
                     "source.repeat; CLI --recoil.)"));
  source_->addButton(beam, 0);
  source_->addButton(ions, 1);
  source_->addButton(recoils, 2);
  source_->addButton(pka, 3);

  pka_box_ = new QGroupBox(tr("Recoil atom (PKA)"), this);
  pka_element_ = new QComboBox(pka_box_);
  pka_element_->setToolTip(tr("Element of the recoil: an element of the target"));
  pka_energy_ = new QLineEdit(pka_box_);
  pka_energy_->setToolTip(tr("Its kinetic energy, e.g. \"10 keV\" (a bare number is keV). With full cascades it "
                             "is displaced if the energy exceeds its Ed and leaves with the energy less Eb."));
  pka_depth_ = new QLineEdit(pka_box_);
  pka_depth_->setToolTip(tr("Depth below the surface where it starts, e.g. \"50 nm\" (lateral position 0)"));
  pka_theta_ = new QLineEdit(pka_box_);
  pka_theta_->setToolTip(tr("Angle of its direction from the depth axis (0 deg: into the target, 90 deg: "
                            "parallel to the surface, 180 deg: towards the surface)"));
  pka_phi_ = new QLineEdit(pka_box_);
  pka_phi_->setToolTip(tr("Azimuth of its direction about the depth axis, from the y axis"));
  pka_count_ = new QSpinBox(pka_box_);
  pka_count_->setRange(1, 1'000'000'000);
  pka_count_->setGroupSeparatorShown(true);
  pka_count_->setToolTip(tr("How many times the recoil is run, each with its own random numbers "
                            "(source.repeat). The path views show the first ones (Output, paths to record)."));
  auto *pf = new QFormLayout(pka_box_);
  pf->addRow(tr("Element"), pka_element_);
  pf->addRow(tr("Energy"), pka_energy_);
  pf->addRow(tr("Depth"), pka_depth_);
  pf->addRow(tr("Polar angle"), pka_theta_);
  pf->addRow(tr("Azimuth"), pka_phi_);
  pf->addRow(tr("Recoils"), pka_count_);
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
  sl->addWidget(pka);
  sl->addLayout(file_row);
  sl->addWidget(particles_info_);

  auto *lay = new QVBoxLayout(this);
  lay->addWidget(ion);
  lay->addWidget(src);
  lay->addWidget(pka_box_);
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
    const bool single = is_single_recoil(c);
    const int id = single ? 3 : static_cast<int>(c.trim.source);
    if (auto *b = source_->button(id)) b->setChecked(true);
    const bool list = c.trim.source != trim::TrimConfig::Source::beam && !single;
    load_->setEnabled(list);
    particles_file_->setEnabled(list);
    pka_box_->setVisible(single);
    if (single) {
      // The target's elements (the recoil must be one of them).
      std::vector<unsigned> zs;
      for (const auto &l : c.trim.layers) {
        for (const auto &e : l.elements) {
          if (std::find(zs.begin(), zs.end(), e.z) == zs.end()) zs.push_back(e.z);
        }
      }
      const unsigned z = c.trim.particles.front().z;
      if (std::find(zs.begin(), zs.end(), z) == zs.end()) zs.push_back(z);
      pka_element_->clear();
      for (const unsigned zi : zs) {
        const auto &el = Element::by_z(zi);
        pka_element_->addItem(QStringLiteral("%1  %2").arg(QString::fromStdString(std::string(el.symbol())),
                                                           QString::fromStdString(std::string(el.name()))),
                              zi);
      }
      pka_element_->setCurrentIndex(pka_element_->findData(z));
    }
    if (single) {
      particles_info_->setText(
          tr("%1 recoil atom(s) of %2, without an ion. %3")
              .arg(c.trim.particle_repeat)
              .arg(format_energy(c.trim.particles.front().energy))
              .arg(c.trim.damage == trim::DamageCalculation::kinchin_pease
                       ? tr("With the quick damage calculation only the Kinchin-Pease estimate is made: choose full "
                            "cascades (Calculation) to follow and see its cascade.")
                       : tr("Its cascades appear in the path views.")));
    } else if (!list) {
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
    if (id == 3) {
      if (is_single_recoil(scenario_->config())) return;
      binder_->apply(tr("Recoil atom (PKA)"),
                     [&](config::RunConfig &c) {
                       // A first atom of the target, 10 keV, halfway through the
                       // first layer (at most 50 nm deep), into the target.
                       trim::SourceParticle p;
                       p.name = "PKA";
                       p.z = c.trim.layers.empty() || c.trim.layers.front().elements.empty()
                                 ? c.trim.ion.z
                                 : c.trim.layers.front().elements.front().z;
                       p.energy = 10.0 * common::kilo_electron_volt;
                       p.x = c.trim.layers.empty() ? 500.0 * common::angstrom
                                                   : std::min(0.5 * c.trim.layers.front().width, 500.0 * common::angstrom);
                       p.y = 0.0 * common::angstrom;
                       p.z_position = 0.0 * common::angstrom;
                       c.trim.source = trim::TrimConfig::Source::particle_recoils;
                       c.trim.particles = {p};
                       c.trim.particles_title = "PKA";
                       c.trim.particle_repeat = 100;
                       c.trim.ion_count = 0;
                       if (c.trim.ion.energy < p.energy) c.trim.ion.energy = p.energy;
                       if (c.trim.damage == trim::DamageCalculation::kinchin_pease) {
                         c.trim.damage = trim::DamageCalculation::full_cascade; // (to follow the cascade)
                       }
                       if (c.trim.transport == trim::TrimConfig::Transport::refined) {
                         c.trim.transport = trim::TrimConfig::Transport::trim2013; // (beam only)
                       }
                     },
                     -1);
      return;
    }
    binder_->apply(tr("Particle source"),
                   [&](config::RunConfig &c) {
                     if (is_single_recoil(c)) {
                       // Leaving the PKA: its particle is not a list.
                       c.trim.particles.clear();
                       c.trim.particles_title.clear();
                       c.trim.particle_repeat = 1;
                       if (id == 0) c.trim.ion_count = 1000;
                     }
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

  // --- The recoil atom (PKA) ------------------------------------------------------
  // (Its fields edit the one particle; they are shown only for a PKA.)
  const auto pka_of = [](config::RunConfig &c) -> trim::SourceParticle & {
    if (!is_single_recoil(c)) throw std::invalid_argument("not a single recoil atom");
    return c.trim.particles.front();
  };
  const auto pka_get = [](const config::RunConfig &c) -> const trim::SourceParticle * {
    return is_single_recoil(c) ? &c.trim.particles.front() : nullptr;
  };
  connect(pka_element_, &QComboBox::activated, this, guarded(this, tr("Recoil element"), [this, pka_of](int i) {
    const unsigned z = pka_element_->itemData(i).toUInt();
    binder_->apply(tr("Recoil element"), [&](config::RunConfig &c) { pka_of(c).z = z; }, -1);
  }));
  binder_->text(pka_energy_, tr("Recoil energy"),
                [pka_get](const config::RunConfig &c) {
                  const auto *p = pka_get(c);
                  return p ? format_energy(p->energy) : QString();
                },
                [pka_of](config::RunConfig &c, const QString &t) {
                  const auto e = parse_energy(t, "keV");
                  if (!(e > 0.0 * common::electron_volt)) throw std::invalid_argument("must be positive");
                  pka_of(c).energy = e;
                  if (c.trim.ion.energy < e) c.trim.ion.energy = e; // (the tables' upper end)
                });
  binder_->text(pka_depth_, tr("Recoil depth"),
                [pka_get](const config::RunConfig &c) {
                  const auto *p = pka_get(c);
                  return p ? format_length(p->x) : QString();
                },
                [pka_of](config::RunConfig &c, const QString &t) {
                  const auto x = parse_length(t, "A");
                  if (!(x >= 0.0 * common::angstrom) || x > c.trim.target_thickness()) {
                    throw std::invalid_argument("must be inside the target (0 to its thickness)");
                  }
                  pka_of(c).x = x;
                });
  // The direction as a polar angle from the depth axis and an azimuth from y.
  const auto angles = [](const trim::SourceParticle &p, double &theta_deg, double &phi_deg) {
    const double n = std::sqrt(p.cos_x * p.cos_x + p.cos_y * p.cos_y + p.cos_z * p.cos_z);
    theta_deg = n > 0.0 ? std::acos(std::clamp(p.cos_x / n, -1.0, 1.0)) * 180.0 / common::pi_v : 0.0;
    phi_deg = std::atan2(p.cos_z, p.cos_y) * 180.0 / common::pi_v;
    if (std::abs(p.cos_y) < 1e-12 && std::abs(p.cos_z) < 1e-12) phi_deg = 0.0;
  };
  const auto set_angles = [](trim::SourceParticle &p, double theta_deg, double phi_deg) {
    const double t = theta_deg * common::pi_v / 180.0;
    const double f = phi_deg * common::pi_v / 180.0;
    p.cos_x = std::cos(t);
    p.cos_y = std::sin(t) * std::cos(f);
    p.cos_z = std::sin(t) * std::sin(f);
  };
  binder_->text(pka_theta_, tr("Recoil polar angle"),
                [pka_get, angles](const config::RunConfig &c) {
                  const auto *p = pka_get(c);
                  if (!p) return QString();
                  double t = 0.0, f = 0.0;
                  angles(*p, t, f);
                  return QStringLiteral("%1 deg").arg(t, 0, 'g', 6);
                },
                [pka_of, angles, set_angles](config::RunConfig &c, const QString &text) {
                  const double t = common::in_units(parse_angle(text), common::degree);
                  if (!(t >= 0.0 && t <= 180.0)) throw std::invalid_argument("must be 0 to 180 degrees");
                  auto &p = pka_of(c);
                  double t0 = 0.0, f = 0.0;
                  angles(p, t0, f);
                  set_angles(p, t, f);
                });
  binder_->text(pka_phi_, tr("Recoil azimuth"),
                [pka_get, angles](const config::RunConfig &c) {
                  const auto *p = pka_get(c);
                  if (!p) return QString();
                  double t = 0.0, f = 0.0;
                  angles(*p, t, f);
                  return QStringLiteral("%1 deg").arg(f, 0, 'g', 6);
                },
                [pka_of, angles, set_angles](config::RunConfig &c, const QString &text) {
                  const double f = common::in_units(parse_angle(text), common::degree);
                  if (!std::isfinite(f)) throw std::invalid_argument("must be a number of degrees");
                  auto &p = pka_of(c);
                  double t = 0.0, f0 = 0.0;
                  angles(p, t, f0);
                  set_angles(p, t, f);
                });
  binder_->spin(pka_count_, tr("Recoils"),
                [](const config::RunConfig &c) {
                  return static_cast<int>(std::min<std::size_t>(c.trim.particle_repeat, 1'000'000'000));
                },
                [](config::RunConfig &c, int v) {
                  if (!is_single_recoil(c)) return;
                  c.trim.particle_repeat = static_cast<std::size_t>(std::max(v, 1));
                });
  binder_->refresh();
}

bool is_single_recoil(const config::RunConfig &c) {
  return c.trim.source == trim::TrimConfig::Source::particle_recoils && c.trim.particles.size() == 1 &&
         c.trim.particles.front().name == "PKA";
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
