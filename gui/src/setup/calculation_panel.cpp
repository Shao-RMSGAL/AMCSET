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
#include "setup/calculation_panel.h"

#include "model/quantity_text.h"
#include "setup/form_binder.h"

#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>

#include <cmath>

namespace amcset {
namespace gui {

namespace {

QLineEdit *path_row(QWidget *parent, QFormLayout *form, const QString &label, bool directory,
                    const QString &filter = {}) {
  auto *edit = new QLineEdit(parent);
  auto *browse = new QPushButton(QObject::tr("Browse..."), parent);
  auto *row = new QHBoxLayout();
  row->addWidget(edit, 1);
  row->addWidget(browse);
  form->addRow(label, row);
  QObject::connect(browse, &QPushButton::clicked, parent, guarded(nullptr, QObject::tr("Choosing a file"), [=] {
    const QString p = directory ? QFileDialog::getExistingDirectory(parent, label, edit->text())
                                : QFileDialog::getOpenFileName(parent, label, edit->text(), filter);
    if (p.isEmpty()) return;
    edit->setText(p);
    emit edit->editingFinished();
  }));
  return edit;
}

QSpinBox *count_spin(QWidget *parent, int max, int min = 0) {
  auto *s = new QSpinBox(parent);
  s->setRange(min, max);
  s->setGroupSeparatorShown(true);
  return s;
}

} // namespace

CalculationPanel::CalculationPanel(Scenario *scenario, QWidget *parent)
    : QWidget(parent), scenario_(scenario), binder_(new FormBinder(scenario, this)) {
  connect(binder_, &FormBinder::error, this, guarded(this, tr("Calculation settings (error)"), [this](auto &&...a) { message(std::forward<decltype(a)>(a)...); }));
  auto *lay = new QVBoxLayout(this);

  // --- Engine and calculation ----------------------------------------------
  auto *calc = new QGroupBox(tr("Calculation"), this);
  auto *cf = new QFormLayout(calc);
  auto *engine = new QComboBox(calc);
  engine->addItems({tr("SRIM mode (TRIM-compatible)"), tr("Native AMCSET engine (ions and electrons)")});
  auto *damage = new QComboBox(calc);
  damage->addItems({tr("Ion distribution and quick damage (Kinchin-Pease)"),
                    tr("Detailed calculation with full damage cascades"),
                    tr("Monolayer collision steps / surface sputtering")});
  auto *ions = count_spin(calc, 2'000'000'000, 1);
  ions->setToolTip(tr("Number of incident ions (for particle lists: 0 = every line)"));
  ions->setMinimum(0);
  auto *seed = count_spin(calc, 2'147'483'647);
  seed->setGroupSeparatorShown(false);
  seed->setToolTip(tr("Random-number seed (TRIM's default is 716381). Every ion has its own stream, so the "
                      "results do not depend on the number of threads."));
  auto *threads = count_spin(calc, 1024);
  threads->setSpecialValueText(tr("all cores"));
  for (auto *c : {engine, damage}) {
    c->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    c->setMinimumContentsLength(18);
  }
  cf->addRow(tr("Engine"), engine);
  cf->addRow(tr("Damage"), damage);
  cf->addRow(tr("Ions"), ions);
  cf->addRow(tr("Random seed"), seed);
  cf->addRow(tr("Threads"), threads);
  lay->addWidget(calc);

  binder_->combo(engine, tr("Engine"), [](const config::RunConfig &c) { return c.mode == config::Mode::srim ? 0 : 1; },
                 [](config::RunConfig &c, int i) { c.mode = i == 0 ? config::Mode::srim : config::Mode::native; });
  binder_->combo(damage, tr("Damage calculation"),
                 [](const config::RunConfig &c) { return static_cast<int>(c.trim.damage) - 1; },
                 [](config::RunConfig &c, int i) {
                   const auto d = static_cast<trim::DamageCalculation>(i + 1);
                   if (d == trim::DamageCalculation::monolayer &&
                       c.trim.source == trim::TrimConfig::Source::particle_ions) {
                     throw std::invalid_argument("TRIM.DAT ions run with the quick or full-cascade calculation");
                   }
                   c.trim.damage = d;
                 });
  binder_->spin(ions, tr("Number of ions"), [](const config::RunConfig &c) { return static_cast<int>(std::min<std::size_t>(c.trim.ion_count, 2'000'000'000)); },
                [](config::RunConfig &c, int v) {
                  if (v == 0 && c.trim.source == trim::TrimConfig::Source::beam) throw std::invalid_argument("at least one ion");
                  c.trim.ion_count = static_cast<std::size_t>(v);
                });
  binder_->spin(seed, tr("Seed"), [](const config::RunConfig &c) { return static_cast<int>(c.trim.seed & 0x7fffffff); },
                [](config::RunConfig &c, int v) { c.trim.seed = static_cast<std::uint32_t>(v); });
  binder_->spin(threads, tr("Threads"), [](const config::RunConfig &c) { return static_cast<int>(c.trim.threads); },
                [](config::RunConfig &c, int v) { c.trim.threads = static_cast<unsigned>(v); });

  // --- SRIM-mode physics -----------------------------------------------------
  srim_ = new QGroupBox(tr("Physics (SRIM mode)"), this);
  auto *sf = new QFormLayout(srim_);
  auto *transport = new QComboBox(srim_);
  transport->addItems({tr("TRIM-2013 (as SRIM-2013's TRIM.exe)"), tr("Refined (smaller steps)")});
  transport->setToolTip(tr("How ions and recoils are followed through amorphous layers. TRIM-2013 follows TRIM.exe "
                           "step by step; Refined limits flight paths to a 5 % energy loss and uses free flight "
                           "paths for recoils. (Scenario key physics.transport.) With the crystal transport "
                           "switched on below, this choice waits: it applies again when that is switched off."));
  auto *crystal = new QCheckBox(tr("Crystal transport (follows the lattice)"), srim_);
  crystal->setToolTip(tr("On: layers that have a crystal (Crystal tab; scenario key layers.N.crystal) are followed "
                         "atom by atom, with collision partners on lattice sites and thermal vibrations, so that "
                         "ions can channel; amorphous layers take random partners. Off: every layer is amorphous, "
                         "as in TRIM, and the transport chosen above is used.\n\n"
                         "Disclaimer: this is not TRIM. It takes none of TRIM.exe's bugs (the list below is then "
                         "all corrected), and the division of its electronic stopping between the path and the "
                         "collisions, which decides the depth of channeling tails, is not calibrated: measured "
                         "boron profiles in silicon allow a wide range (docs/CLI.md \"Crystals\"). Switching it on "
                         "also switches the automatic target depth off. (Scenario key physics.transport: "
                         "\"crystal\".)"));
  auto *stopping = new QComboBox(srim_);
  stopping->addItems({tr("SRIM-2013 stopping tables"), tr("Built-in (no data files)"), tr("Bethe-Bloch"),
                      tr("Lindhard-Scharff")});
  auto *table = path_row(srim_, sf, tr("Stopping tables"), false,
                         tr("AMCSET stopping tables (*.amcset *.amcset.gz);;All files (*)"));
  table->setPlaceholderText(tr("SRIM-2013 tables that come with AMCSET"));
  table->setToolTip(tr("Empty: the SRIM-2013 tables of every ion in every element that come with AMCSET "
                       "(data/srim2013), which TRIM's own stopping reproduces. Or a file made by tools/srim "
                       "(from SRIM's SR Module), or a folder of SR Module outputs."));
  for (auto *c : {transport, stopping}) {
    c->setSizeAdjustPolicy(QComboBox::AdjustToMinimumContentsLengthWithIcon);
    c->setMinimumContentsLength(18);
  }
  sf->insertRow(0, tr("Transport"), transport);
  sf->insertRow(1, QString(), crystal);
  sf->insertRow(2, tr("Electronic stopping"), stopping);
  auto *nickel = new QCheckBox(tr("Correct the stopping of Ni ions"), srim_);
  nickel->setToolTip(
      tr("SRIM-2013's electronic stopping of Ni projectiles is up to 40 % below that of Co and Cu at "
         "0.02-0.3 MeV/u, in every target, and measured Ni ranges in Fe and Mo are 1.3-1.5 times shorter than "
         "SRIM-2013's. On: Ni ions (and Ni recoils) take the geometric mean of Co's and Cu's stopping at the same "
         "velocity, which shortens their ranges by about 1.27 times at 5-12 MeV.\n\n"
         "Not TRIM: off by default, so that results agree with TRIM.exe. Evidence: docs/VALIDATION.md §3, §4.1, "
         "§9. (Scenario key physics.stopping_corrections.nickel_ions; SRIM tables only.)"));
  sf->insertRow(3, QString(), nickel);
  auto *srim_data = path_row(srim_, sf, tr("SRIM Data folder"), true);
  srim_data->setToolTip(tr("Optional: SRIM's Data folder (pair potentials, if requested)"));
  auto *straggling = new QCheckBox(tr("Electronic energy-loss straggling"), srim_);
  sf->addRow(QString(), straggling);
  auto *final_e = new QLineEdit(srim_);
  final_e->setToolTip(tr("Atoms below this energy stop (TRIM: 1 eV)"));
  sf->addRow(tr("Final energy"), final_e);
  auto *max_recoils = count_spin(srim_, 2'000'000'000);
  max_recoils->setSpecialValueText(tr("no limit"));
  max_recoils->setToolTip(tr("Largest number of recoils followed per ion (TRIM.cfg's cascade size)"));
  sf->addRow(tr("Recoils per ion"), max_recoils);
  auto *fast = new QCheckBox(tr("Fast approximations (about 2x faster)"), srim_);
  fast->setToolTip(
      tr("Takes the magic-formula scattering angle, the stopping powers and a few power functions from "
         "precomputed tables (relative error below 4e-5) and uses cheaper sines, cosines and "
         "rotations.\n\n"
         "Disclaimer: this may not faithfully replicate TRIM. Results agree with the exact calculation "
         "within statistics, but not number for number: the same seed gives different individual "
         "histories, so a run with this option does not reproduce a run without it, nor TRIM.exe step "
         "by step. Leave it off for TRIM comparisons and reference results. (No TRIM equivalent; "
         "scenario key physics.fast_math.)"));
  sf->addRow(QString(), fast);
  lay->addWidget(srim_);

  // The list holds the two transports of amorphous targets; the crystal
  // transport is the switch below it. While that is on, the list shows the
  // transport that comes back when it is switched off.
  binder_->combo(transport, tr("Transport"),
                 [this](const config::RunConfig &c) {
                   const auto t = c.trim.transport == trim::TrimConfig::Transport::crystal ? scenario_->amorphous_transport()
                                                                                            : c.trim.transport;
                   return t == trim::TrimConfig::Transport::refined ? 1 : 0;
                 },
                 [](config::RunConfig &c, int i) {
                   if (c.trim.transport == trim::TrimConfig::Transport::crystal) return; // the switch decides
                   c.trim.transport = i == 1 ? trim::TrimConfig::Transport::refined : trim::TrimConfig::Transport::trim2013;
                 });
  binder_->custom([transport, crystal](const config::RunConfig &c) {
    const bool on = c.trim.transport == trim::TrimConfig::Transport::crystal;
    crystal->setChecked(on); // (clicked() is not emitted for this)
    transport->setEnabled(!on);
  });
  connect(crystal, &QCheckBox::clicked, this, guarded(this, tr("Crystal transport"), [this](bool on) {
    const QString err = scenario_->set_crystal_transport(on);
    if (!err.isEmpty()) {
      emit message(tr("Crystal transport: %1").arg(err));
      binder_->refresh(); // the box goes back to what the scenario says
    }
  }));
  binder_->combo(stopping, tr("Electronic stopping"), [](const config::RunConfig &c) { return static_cast<int>(c.trim.stopping_model); },
                 [](config::RunConfig &c, int i) { c.trim.stopping_model = static_cast<physics::ElectronicStoppingModel>(i); });
  binder_->check(nickel, tr("Ni stopping correction"),
                 [](const config::RunConfig &c) { return c.trim.stopping_corrections.nickel_ions; },
                 [](config::RunConfig &c, bool on) { c.trim.stopping_corrections.nickel_ions = on; });
  binder_->custom([nickel](const config::RunConfig &c) {
    nickel->setEnabled(c.trim.stopping_model == physics::ElectronicStoppingModel::srim_tables);
  });
  binder_->text(table, tr("Stopping tables"), [](const config::RunConfig &c) { return QString::fromStdString(c.trim.stopping_table); },
                [](config::RunConfig &c, const QString &t) {
                  c.trim.stopping_table = t.trimmed().toStdString();
                  if (!c.trim.stopping_table.empty()) c.trim.stopping_model = physics::ElectronicStoppingModel::srim_tables;
                });
  binder_->text(srim_data, tr("SRIM Data folder"), [](const config::RunConfig &c) { return QString::fromStdString(c.trim.srim_data_directory); },
                [](config::RunConfig &c, const QString &t) { c.trim.srim_data_directory = t.trimmed().toStdString(); });
  binder_->check(straggling, tr("Straggling"), [](const config::RunConfig &c) { return c.trim.electronic_straggling; },
                 [](config::RunConfig &c, bool on) { c.trim.electronic_straggling = on; });
  binder_->check(fast, tr("Fast approximations"), [](const config::RunConfig &c) { return c.trim.fast_math; },
                 [](config::RunConfig &c, bool on) { c.trim.fast_math = on; });
  binder_->text(final_e, tr("Final energy"), [](const config::RunConfig &c) { return format_energy(c.trim.final_energy); },
                [](config::RunConfig &c, const QString &t) {
                  const auto e = parse_energy(t);
                  if (!(e > 0.0 * common::electron_volt)) throw std::invalid_argument("must be positive");
                  c.trim.final_energy = e;
                });
  binder_->spin(max_recoils, tr("Recoils per ion"),
                [](const config::RunConfig &c) { return static_cast<int>(std::min<std::size_t>(c.trim.max_recoils_per_ion, 2'000'000'000)); },
                [](config::RunConfig &c, int v) { c.trim.max_recoils_per_ion = static_cast<std::size_t>(v); });

  // --- Target depth (SRIM mode) --------------------------------------------------
  auto *depth = new QGroupBox(tr("Target depth (SRIM mode)"), this);
  auto *df = new QFormLayout(depth);
  auto *auto_depth = new QCheckBox(tr("Automatic target depth"), depth);
  auto_depth->setToolTip(
      tr("Before the calculation, quick test runs (Kinchin-Pease, fast approximations) size the last layer so "
         "that the deepest test ion stops inside the target, in its last part (Margin): the layer grows while "
         "ions are transmitted and shrinks while the deepest stops shallower. The first guess comes from the "
         "stopping powers (the ions' path length). The sized thickness replaces the last layer's and can be "
         "undone.\n\nWith the crystal transport, channeled ions stop far beyond that guess: the first test run "
         "is made in a box 32 times as deep, which is then trimmed to the deepest ion (two test runs as a "
         "rule; they take as long per ion as the calculation itself).\n\nLeave this off to study "
         "transmission through a given target, and for recoil sources "
         "(TRIM.DAT atoms start at given depths). (No TRIM equivalent; scenario key auto_depth.enabled.)"));
  auto *test_ions = count_spin(depth, static_cast<int>(10'000'000), 1);
  test_ions->setToolTip(tr("Ions per test run. The deepest of more ions lies deeper in the tail of the range "
                           "distribution; a long calculation can still have a few ions deeper than every test "
                           "ion. (auto_depth.test_ions)"));
  auto *margin = new QSpinBox(depth);
  margin->setRange(1, 50);
  margin->setSuffix(tr(" %"));
  margin->setToolTip(tr("The deepest test ion must stop within this fraction of the target at its back; a "
                        "larger margin leaves more room for the deepest ions of the full calculation. "
                        "(auto_depth.margin)"));
  auto *warn = new QCheckBox(tr("Warn when ions leave through the back"), depth);
  warn->setToolTip(tr("During a calculation, report ions transmitted through the target (the calculation goes "
                      "on) and offer to run again with the automatic target depth. Turn it off for "
                      "transmission studies. (auto_depth.warn_transmitted)"));
  df->addRow(QString(), auto_depth);
  df->addRow(tr("Test ions"), test_ions);
  df->addRow(tr("Margin"), margin);
  df->addRow(QString(), warn);
  lay->addWidget(depth);

  binder_->check(auto_depth, tr("Automatic target depth"), [](const config::RunConfig &c) { return c.trim.auto_depth.enabled; },
                 [](config::RunConfig &c, bool on) { c.trim.auto_depth.enabled = on; });
  binder_->spin(test_ions, tr("Test ions"),
                [](const config::RunConfig &c) {
                  return static_cast<int>(std::min<std::size_t>(c.trim.auto_depth.test_ions, 10'000'000));
                },
                [](config::RunConfig &c, int v) { c.trim.auto_depth.test_ions = static_cast<std::size_t>(v); });
  binder_->spin(margin, tr("Margin"),
                [](const config::RunConfig &c) { return static_cast<int>(std::lround(c.trim.auto_depth.margin * 100.0)); },
                [](config::RunConfig &c, int v) { c.trim.auto_depth.margin = v / 100.0; });
  binder_->check(warn, tr("Warn when ions escape"), [](const config::RunConfig &c) { return c.trim.auto_depth.warn_transmitted; },
                 [](config::RunConfig &c, bool on) { c.trim.auto_depth.warn_transmitted = on; });
  binder_->custom([depth, test_ions, margin](const config::RunConfig &c) {
    depth->setVisible(c.mode == config::Mode::srim);
    const bool on = c.trim.auto_depth.enabled && c.trim.source == trim::TrimConfig::Source::beam;
    test_ions->setEnabled(on);
    margin->setEnabled(on);
  });

  // --- Bugs of TRIM.exe ----------------------------------------------------------
  bugs_ = new QGroupBox(tr("TRIM.exe bugs (trim2013 transport)"), this);
  bugs_->setToolTip(tr("SRIM-2013's TRIM.exe has bugs that change its results. Checked = reproduce the bug, so that "
                       "AMCSET agrees with TRIM (the default). Unchecked = the corrected behaviour: results then "
                       "depart from TRIM's and should not be compared with TRIM runs or published TRIM numbers."));
  auto *bl = new QVBoxLayout(bugs_);
  bugs_note_ = new QLabel(tr("<b>All corrected.</b> The crystal transport is not TRIM and takes none of TRIM.exe's "
                             "bugs. The switches set for the other transports are kept and apply again when the "
                             "crystal transport is switched off."),
                          bugs_);
  bugs_note_->setWordWrap(true);
  bl->addWidget(bugs_note_);
  const QString caveat = tr("\n\nUnchecked results depart from TRIM's; keep it checked to compare with TRIM.");
  struct BugRow {
    bool trim::TrimBugs::*flag;
    QString label;
    QString tip;
  };
  const std::vector<BugRow> rows = {
      {&trim::TrimBugs::stale_primary_energy, tr("Stale energy of undisplaced primaries"),
       tr("TRIM: when an ion gives a target atom T <= Ed, the atom is not displaced, but TRIM still runs its "
          "cascade with the energy the recoil slot had at the end of the previous cascade "
          "(STARTRECOIL).\nFixed: no cascade; T stays in the ion's phonons.\nEffect of the fix: recoil phonons -8 to "
          "-11 %, recoil ionization -6 to -9 %, sputtering -3 to -5 % at 1 keV (Fe -> Fe, Ar -> Cu); about -1 % "
          "at 1 MeV.")},
      {&trim::TrimBugs::recoil_energy_created, tr("Energy created in slow recoil collisions"),
       tr("TRIM: a recoil's electronic loss over one step and the energy it gives in the collision are both "
          "computed at its energy before the step, can together exceed it, and are both booked in full; the "
          "mover keeps max(0, E - T - DEE) (0x436de5). Fe -> Fe at 1 keV ends with about 105 % of the incident "
          "energy.\nFixed: the struck atom gets at most E - DEE.\nEffect of the fix: replacements -3 to -7 %, "
          "sputtering -3 to -7 % at 1-10 keV, recoil ionization up to -4 %, vacancies up to -1 %.")},
      {&trim::TrimBugs::struck_atom_direction, tr("Direction of atoms struck by recoils"),
       tr("TRIM: the atom struck by a recoil is turned by DIRCOS with sin x = 0 (0x4377d5), so its direction "
          "depends only on the mover's x cosine and its transverse direction is random.\nFixed: turned from the "
          "mover's direction, as the primary of an ion collision is.\nEffect of the fix: sputtering -17 % (Ar -> "
          "Cu 1 keV), -4 % (Fe -> Fe 10 keV); vacancies up to -3 %; sputtered atoms slower and less oblique.")},
      {&trim::TrimBugs::replacement_interstitial, tr("Replacing atoms counted as interstitials"),
       tr("TRIM: an atom that makes a replacement collision fills the vacancy, but is still counted as an "
          "interstitial (RECSTOPDATA, 0x457970), so interstitials exceed vacancies by the replacements.\nFixed: "
          "it is not an interstitial.\nEffect of the fix: interstitials -2 to -5 % (only the interstitial "
          "count and profile change).")},
      {&trim::TrimBugs::subthreshold_sputter_counted, tr("Sputtered sub-threshold atoms counted as vacancies"),
       tr("TRIM: a sputtered atom that was set in motion below Ed (near the surface) adds a displacement and a "
          "vacancy to the totals of VACANCY.txt and COLLISON.txt, not to the depth tables (SPUTTER, "
          "0x441ff0).\nFixed: it is not counted; it was never displaced.\nEffect of the fix: total vacancies and "
          "displacements -9 % (Fe -> Fe 1 keV), -17 % (Ar -> Cu 1 keV), -2.5 % at 10 keV, none at MeV energies; "
          "the vacancy profiles do not change.")},
      {&trim::TrimBugs::ion_phonons_counted_twice, tr("Small ion transfers counted twice"),
       tr("TRIM: an ion collision that gives less than the stop energy (1 eV) is booked as ion phonons and again "
          "as recoil phonons (0x44d4e0, 0x4575e0).\nFixed: once, as ion phonons.\nEffect of the fix: recoil "
          "phonons about -0.02 %.")},
      {&trim::TrimBugs::cascade_random_replay, tr("Cascades replay the ion's random numbers"),
       tr("TRIM: the random-number state is saved before each cascade and restored after it (0x434e89, "
          "0x438140), so a cascade uses the numbers the ion draws next and consecutive cascades are "
          "correlated.\nFixed: each cascade has its own random-number stream.\nEffect of the fix: within "
          "statistics for the totals (below 1 %), but correlated cascades make TRIM's run-to-run scatter "
          "larger.")},
      {&trim::TrimBugs::zero_energy_defaults_skip_first_atom, tr("Zero energies kept for the first atom"),
       tr("TRIM replaces a zero displacement, surface or lattice binding energy by 20, 4 and 2 eV for target "
          "atoms 2 to n, but its loop (0x44efc0) skips the first atom, which keeps zeros.\nFixed: every atom "
          "with a zero energy gets TRIM's default.\nEffect of the fix: only for inputs with zero energies; for "
          "the first atom a zero Ed displaces every struck atom.")},
      {&trim::TrimBugs::integer_layer_z, tr("Layer atomic number rounded"),
       tr("TRIM rounds the layer's mean atomic number at every element (z = CInt(c Z + z), 0x44f8c0) for the "
          "ion's flight path, straggling and the quick calculation.\nFixed: the stoichiometric mean.\nEffect of "
          "the fix: compounds and alloys only, below 0.5 % for the cases tried.")},
      {&trim::TrimBugs::exyz_skips_first_step, tr("EXYZ.txt skips its first energy step"),
       tr("TRIM sets the next EXYZ record a step too low after the record at the start (0x434745), so the "
          "first crossing of an energy step is not written.\nFixed: every crossing is written.\nEffect: EXYZ.txt "
          "only.")},
      {&trim::TrimBugs::recoil_event_counted_twice, tr("TRIM.DAT recoil events counted twice"),
       tr("TRIM's code 6 counts every event twice in the ion column of RANGE.txt (where the atom starts and "
          "where it stops).\nFixed: once, where it starts.\nEffect: the ion range table of TRIM.DAT recoil runs "
          "only (halved).")},
  };
  for (const auto &row : rows) {
    auto *cb = new QCheckBox(row.label, bugs_);
    cb->setToolTip(row.tip + caveat);
    bl->addWidget(cb);
    bug_widgets_.push_back(cb);
    const auto flag = row.flag;
    // (Shown as the run takes them: all corrected with the crystal transport.)
    binder_->check(cb, row.label, [flag](const config::RunConfig &c) { return c.trim.effective_bugs().*flag; },
                   [flag](config::RunConfig &c, bool on) { c.trim.bugs.*flag = on; });
  }
  {
    auto *sp = new QCheckBox(tr("Single-precision sums in the TRIM text tables"), bugs_);
    sp->setToolTip(tr("TRIM adds 1/(bin width) per event to its VACANCY, NOVAC and RANGE tables in single "
                      "precision, which shifts the table totals of long runs by up to +-1.5 %. Only the "
                      "TRIM-format .txt tables are affected; the totals, the CSV/JSON/Parquet tables and the "
                      "physics are exact either way. (Also on the Output tab.)") + caveat);
    bl->addWidget(sp);
    bug_widgets_.push_back(sp);
    binder_->check(sp, tr("TRIM single precision"),
                   [](const config::RunConfig &c) { return c.trim.effective_single_precision_tables(); },
                   [](config::RunConfig &c, bool on) { c.trim.single_precision_tables = on; });
  }
  auto *buttons = new QHBoxLayout;
  auto *as_trim = new QPushButton(tr("Reproduce all"), bugs_);
  as_trim->setToolTip(tr("Behave exactly as TRIM.exe (the default)."));
  auto *fix_all = new QPushButton(tr("Fix all"), bugs_);
  fix_all->setToolTip(tr("Correct every bug. Results then depart from TRIM's (see each entry).") + caveat);
  buttons->addWidget(as_trim);
  buttons->addWidget(fix_all);
  bug_widgets_.push_back(as_trim);
  bug_widgets_.push_back(fix_all);
  buttons->addStretch(1);
  bl->addLayout(buttons);
  connect(as_trim, &QPushButton::clicked, this, guarded(this, tr("Calculation settings (clicked)"), [this] {
    binder_->apply(tr("TRIM.exe bugs"), [](config::RunConfig &c) {
      c.trim.bugs.set_all(true);
      c.trim.single_precision_tables = true;
    }, -1);
    binder_->refresh();
  }));
  connect(fix_all, &QPushButton::clicked, this, guarded(this, tr("Calculation settings (clicked)"), [this] {
    binder_->apply(tr("TRIM.exe bugs"), [](config::RunConfig &c) {
      c.trim.bugs.set_all(false);
      c.trim.single_precision_tables = false;
    }, -1);
    binder_->refresh();
  }));
  lay->addWidget(bugs_);

  // --- Native engine -----------------------------------------------------------
  native_ = new QGroupBox(tr("Native engine"), this);
  auto *nf = new QFormLayout(native_);
  auto *electron = new QCheckBox(tr("Electrons instead of ions"), native_);
  auto *cascade = new QCheckBox(tr("Follow recoil cascades"), native_);
  auto *screening = new QComboBox(native_);
  screening->addItems({"thomas-fermi", "quickbasic", "zbl", "moliere", "custom"});
  auto *electronic = new QComboBox(native_);
  electronic->addItems({"zbl", "lss", "custom", "bethe-bloch", "builtin", "srim"});
  auto *ion_stop = new QLineEdit(native_);
  auto *ed = new QLineEdit(native_);
  auto *e_stop = new QLineEdit(native_);
  auto *bins = count_spin(native_, 100000, 1);
  auto *angles = count_spin(native_, 1'000'000, 1);
  auto *flights = count_spin(native_, 1'000'000, 1);
  auto *per_history = new QCheckBox(tr("Write one file per incident particle"), native_);
  nf->addRow(QString(), electron);
  nf->addRow(QString(), cascade);
  nf->addRow(tr("Screening"), screening);
  nf->addRow(tr("Electronic stopping"), electronic);
  nf->addRow(tr("Ion stop energy"), ion_stop);
  nf->addRow(tr("Displacement energy"), ed);
  nf->addRow(tr("Electron stop energy"), e_stop);
  nf->addRow(tr("Histogram bins"), bins);
  nf->addRow(tr("Angle divisions"), angles);
  nf->addRow(tr("Flight-length divisions"), flights);
  nf->addRow(QString(), per_history);
  lay->addWidget(native_);
  lay->addStretch(1);

  const auto energy_field = [this](QLineEdit *w, const QString &what, common::energy_quantity config::NativeOptions::*m) {
    binder_->text(w, what, [m](const config::RunConfig &c) { return format_energy(c.native.*m); },
                  [m](config::RunConfig &c, const QString &t) {
                    const auto e = parse_energy(t);
                    if (!(e >= 0.0 * common::electron_volt)) throw std::invalid_argument("must be >= 0");
                    c.native.*m = e;
                  });
  };
  binder_->check(electron, tr("Electrons"), [](const config::RunConfig &c) { return c.native.electron; },
                 [](config::RunConfig &c, bool on) { c.native.electron = on; });
  binder_->check(cascade, tr("Cascades"), [](const config::RunConfig &c) { return c.native.damage_cascade; },
                 [](config::RunConfig &c, bool on) { c.native.damage_cascade = on; });
  binder_->combo(screening, tr("Screening"),
                 [screening](const config::RunConfig &c) { return screening->findText(QString::fromStdString(c.native.screening_model)); },
                 [screening](config::RunConfig &c, int i) { c.native.screening_model = screening->itemText(i).toStdString(); });
  binder_->combo(electronic, tr("Electronic stopping"),
                 [electronic](const config::RunConfig &c) {
                   return electronic->findText(QString::fromStdString(c.native.electronic_stopping_model));
                 },
                 [electronic](config::RunConfig &c, int i) { c.native.electronic_stopping_model = electronic->itemText(i).toStdString(); });
  energy_field(ion_stop, tr("Ion stop energy"), &config::NativeOptions::ion_stop_energy);
  energy_field(ed, tr("Displacement energy"), &config::NativeOptions::displacement_energy);
  energy_field(e_stop, tr("Electron stop energy"), &config::NativeOptions::electron_stop_energy);
  binder_->spin(bins, tr("Bins"), [](const config::RunConfig &c) { return static_cast<int>(c.native.bins); },
                [](config::RunConfig &c, int v) { c.native.bins = static_cast<std::size_t>(v); });
  binder_->spin(angles, tr("Angle divisions"), [](const config::RunConfig &c) { return static_cast<int>(c.native.angle_divisors); },
                [](config::RunConfig &c, int v) { c.native.angle_divisors = static_cast<std::size_t>(v); });
  binder_->spin(flights, tr("Flight divisions"), [](const config::RunConfig &c) { return static_cast<int>(c.native.flying_distances); },
                [](config::RunConfig &c, int v) { c.native.flying_distances = static_cast<std::size_t>(v); });
  binder_->check(per_history, tr("Per-particle files"), [](const config::RunConfig &c) { return c.native.per_bombardment_output; },
                 [](config::RunConfig &c, bool on) { c.native.per_bombardment_output = on; });

  binder_->custom([this, damage](const config::RunConfig &c) {
    const bool srim = c.mode == config::Mode::srim;
    srim_->setVisible(srim);
    bugs_->setVisible(srim);
    // The switches belong to the trim2013 transport. (The group itself stays
    // enabled, so that its note for the crystal transport can be read.)
    for (QWidget *w : bug_widgets_) w->setEnabled(c.trim.transport == trim::TrimConfig::Transport::trim2013);
    bugs_note_->setVisible(c.trim.transport == trim::TrimConfig::Transport::crystal);
    native_->setVisible(!srim);
    damage->setEnabled(srim);
  });
  binder_->refresh();
}

} // namespace gui
} // namespace amcset
