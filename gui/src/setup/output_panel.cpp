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
#include "setup/output_panel.h"

#include "model/quantity_text.h"
#include "setup/form_binder.h"
#include "setup/output_spec_dialog.h"

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QFormLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QPushButton>
#include <QStringList>
#include <QVBoxLayout>

#include <algorithm>

namespace amcset {
namespace gui {

OutputPanel::OutputPanel(Scenario *scenario, std::function<std::shared_ptr<const output::Dataset>()> sample,
                         QWidget *parent)
    : QWidget(parent), scenario_(scenario), binder_(new FormBinder(scenario, this)), sample_(std::move(sample)) {
  connect(binder_, &FormBinder::error, this, guarded(this, tr("Output settings (error)"), [this](auto &&...a) { message(std::forward<decltype(a)>(a)...); }));
  auto *lay = new QVBoxLayout(this);

  auto *files = new QGroupBox(tr("Output files"), this);
  auto *ff = new QFormLayout(files);
  auto *dir = new QLineEdit(files);
  dir->setPlaceholderText(tr("(no files are written)"));
  auto *browse = new QPushButton(tr("Browse..."), files);
  auto *drow = new QHBoxLayout();
  drow->addWidget(dir, 1);
  drow->addWidget(browse);
  ff->addRow(tr("Folder"), drow);
  auto *fmt_row = new QHBoxLayout();
  auto *f_trim = new QCheckBox(tr("TRIM text"), files);
  f_trim->setToolTip(tr("RANGE.txt, VACANCY.txt, IONIZ.txt, ... in TRIM's own layout"));
  auto *f_csv = new QCheckBox(tr("CSV"), files);
  auto *f_json = new QCheckBox(tr("JSON"), files);
  auto *f_pq = new QCheckBox(tr("Parquet"), files);
  for (auto *c : {f_trim, f_csv, f_json, f_pq}) fmt_row->addWidget(c);
  fmt_row->addStretch(1);
  ff->addRow(tr("Formats"), fmt_row);
  auto *single = new QCheckBox(tr("TRIM's rounding in the TRIM text tables"), files);
  single->setToolTip(tr("TRIM sums its vacancy, replacement and range tables in single precision, which shifts "
                        "them by up to about 1.5 %; reproduce that in the TRIM text files (CSV, JSON and "
                        "Parquet always hold the exact values). Not with the crystal transport, which is not "
                        "TRIM: its tables are always exact."));
  ff->addRow(QString(), single);
  auto *customize = new QPushButton(tr("Tables, units, numbers..."), files);
  customize->setToolTip(tr("Choose tables, columns, rows, units and number formats of the output files"));
  summary_ = new QLabel(files);
  summary_->setWordWrap(true);
  ff->addRow(QString(), customize);
  ff->addRow(QString(), summary_);
  lay->addWidget(files);

  auto *tally = new QGroupBox(tr("Depth distributions"), this);
  auto *tf = new QFormLayout(tally);
  auto *w0 = new QLineEdit(tally);
  auto *w1 = new QLineEdit(tally);
  w1->setPlaceholderText(tr("0 = end of the target"));
  auto *bins = new QSpinBox(tally);
  bins->setRange(1, 100000);
  bins->setToolTip(tr("TRIM uses 100"));
  auto *wrow = new QHBoxLayout();
  wrow->addWidget(w0);
  wrow->addWidget(new QLabel(tr("to"), tally));
  wrow->addWidget(w1);
  tf->addRow(tr("Depth window"), wrow);
  tf->addRow(tr("Bins"), bins);
  lay->addWidget(tally);

  auto *lists = new QGroupBox(tr("Particle lists"), this);
  auto *lf = new QFormLayout(lists);
  auto *range3d = new QCheckBox(tr("Ion end points (RANGE_3D)"), lists);
  range3d->setToolTip(tr("Final position of every ion that stops in the target"));
  auto *exits = new QCheckBox(tr("Atoms leaving the target"), lists);
  exits->setToolTip(tr("Backscattered and transmitted ions and sputtered atoms (BACKSCAT, TRANSMIT, SPUTTER)"));
  auto *coll = new QSpinBox(lists);
  coll->setRange(0, 2'000'000'000);
  coll->setSpecialValueText(tr("off"));
  coll->setToolTip(tr("Ions whose collisions and cascades are listed (COLLISON.txt)"));
  auto *coll_atoms = new QCheckBox(tr("List the atoms of each cascade"), lists);
  coll_atoms->setToolTip(tr("COLLISON.txt also lists each displacement and replacement of the cascades "
                            "(TRIM's \"Ion + Recoils\"); the file grows with the damage"));
  auto *exyz = new QLineEdit(lists);
  exyz->setToolTip(tr("Record the ion every time it has lost this much energy (EXYZ.txt); 0 = off"));
  auto *paths = new QSpinBox(lists);
  paths->setRange(0, 1'000'000);
  paths->setToolTip(tr("Number of ions whose paths (with their recoils) are drawn"));
  lf->addRow(QString(), range3d);
  lf->addRow(QString(), exits);
  lf->addRow(tr("Collision details for"), coll);
  lf->addRow(QString(), coll_atoms);
  lf->addRow(tr("Energy step (EXYZ)"), exyz);
  lf->addRow(tr("Paths drawn"), paths);
  auto *fluence = new QLineEdit(lists);
  fluence->setToolTip(tr("The fluence (ions per cm2) the dpa table and plot are given for: dpa = "
                         "displacements per ion and depth / atomic density x fluence (output.dpa_fluence; "
                         "default 1e15). Not a TRIM setting."));
  lf->addRow(tr("dpa at fluence (ions/cm2)"), fluence);
  auto *all_paths = new QCheckBox(tr("Record every ion's path"), lists);
  all_paths->setToolTip(tr("Record the trajectories of all ions (output.trajectories = \"all\"), e.g. for the "
                           "trajectory file. The path views keep at most their segment limit."));
  lf->addRow(QString(), all_paths);
  lay->addWidget(lists);

  // --- Trajectory and MD event files (not TRIM; trim2013 transport) ---------
  auto *md = new QGroupBox(tr("Trajectory and MD event files (not TRIM outputs)"), this);
  md->setToolTip(tr("Files for coupling the cascades to molecular dynamics (hybrid BCA-MD) and for "
                    "analysing trajectories. SRIM mode with the TRIM-2013 transport only. File names "
                    "ending in .csv or .parquet; relative names go into the output folder. See "
                    "docs/OUTPUT.md."));
  auto *mf = new QFormLayout(md);
  auto *traj_file = new QLineEdit(md);
  traj_file->setPlaceholderText(tr("off (e.g. trajectories.parquet)"));
  traj_file->setToolTip(tr("Every recorded point of the paths above (Paths drawn / every ion): ion, point, "
                           "track, parent track and point, kind (birth, collision, stop, exit), Z, mass, "
                           "position, direction, energy, energy before the collision, energy transferred, "
                           "electronic loss, time. Large: about 100 000 points per 5 MeV Fe ion; Parquet "
                           "is the compact choice. (output.trajectory_file)"));
  auto *births = new QCheckBox(tr("Only where atoms start to move (births)"), md);
  births->setToolTip(tr("Write only the first point of every track (output.trajectory_file.points = births)"));
  auto *traj_emin = new QLineEdit(md);
  traj_emin->setToolTip(tr("Write only points with at least this energy (0 = all)"));
  auto *md_file = new QLineEdit(md);
  md_file->setPlaceholderText(tr("off (e.g. events.csv)"));
  md_file->setToolTip(tr("Recoil events for MD (output.md_events): for each transition energy E_c, a recoil "
                         "set in motion with E_min <= T < E_c, and a moving atom at the collision where its "
                         "energy first falls below E_c. Needs full cascades (or monolayer steps). Columns: "
                         "ion_id, event_id, E_c_eV, kind, t_ps, x_A, y_A, z_A (depth), E_eV, ux, uy, uz, gen, "
                         "Z, mass_amu, track, parent_track, parent_event."));
  auto *md_ec = new QLineEdit(md);
  md_ec->setToolTip(tr("Transition energies E_c, ascending: \"1 keV\" or \"0.5, 1, 2, 5 keV\" "
                       "(a list only in record mode)"));
  auto *md_emin = new QLineEdit(md);
  md_emin->setToolTip(tr("Recoils set in motion below this are not events; they stay with the BCA "
                         "(\"Ed\" = each struck element's displacement energy)"));
  auto *md_mode = new QComboBox(md);
  md_mode->addItems({tr("Record (transport and tallies unchanged)"), tr("Hand off (MD takes the events)")});
  md_mode->setToolTip(tr("Record: only write the events; every tally is the same as without them.\n"
                         "Hand off: the BCA stops following an event's atom, so nothing is counted twice "
                         "and the run is faster; the tallies then cover the BCA part only, and summary.json "
                         "gives the energy handed to MD. One E_c. Not TRIM: for MD coupling, fix TRIM's "
                         "bugs as well (Calculation tab), or energy is created in slow collisions."));
  mf->addRow(tr("Trajectory file"), traj_file);
  mf->addRow(QString(), births);
  mf->addRow(tr("Minimum energy"), traj_emin);
  mf->addRow(tr("MD event file"), md_file);
  mf->addRow(tr("E_c"), md_ec);
  mf->addRow(tr("E_min"), md_emin);
  mf->addRow(tr("Mode"), md_mode);
  lay->addWidget(md);
  lay->addStretch(1);

  binder_->text(dir, tr("Output folder"), [](const config::RunConfig &c) { return QString::fromStdString(c.trim.output_directory); },
                [](config::RunConfig &c, const QString &t) { c.trim.output_directory = t.trimmed().toStdString(); });
  connect(browse, &QPushButton::clicked, this, guarded(this, tr("Output settings (clicked)"), [this, dir] {
    const QString p = QFileDialog::getExistingDirectory(this, tr("Output folder"), dir->text());
    if (p.isEmpty()) return;
    dir->setText(p);
    emit dir->editingFinished();
  }));
  const auto format_box = [this](QCheckBox *c, output::Format f, const QString &name) {
    binder_->check(c, name, [f](const config::RunConfig &r) { return r.trim.outputs.formats.count(f) > 0; },
                   [f](config::RunConfig &r, bool on) {
                     if (on) r.trim.outputs.formats.insert(f);
                     else r.trim.outputs.formats.erase(f);
                   });
  };
  format_box(f_trim, output::Format::trim, tr("TRIM files"));
  format_box(f_csv, output::Format::csv, tr("CSV files"));
  format_box(f_json, output::Format::json, tr("JSON file"));
  format_box(f_pq, output::Format::parquet, tr("Parquet files"));
  binder_->check(single, tr("TRIM rounding"),
                 [](const config::RunConfig &c) { return c.trim.effective_single_precision_tables(); },
                 [](config::RunConfig &c, bool on) { c.trim.single_precision_tables = on; });
  binder_->custom([single](const config::RunConfig &c) {
    single->setEnabled(c.trim.transport != trim::TrimConfig::Transport::crystal);
  });
  connect(customize, &QPushButton::clicked, this, guarded(this, tr("Output settings (clicked)"), [this] { this->customize(); }));
  binder_->text(w0, tr("Depth window start"), [](const config::RunConfig &c) { return format_length(c.trim.window_start); },
                [](config::RunConfig &c, const QString &t) { c.trim.window_start = parse_length(t); });
  binder_->text(w1, tr("Depth window end"), [](const config::RunConfig &c) { return format_length(c.trim.window_end); },
                [](config::RunConfig &c, const QString &t) { c.trim.window_end = parse_length(t); });
  binder_->spin(bins, tr("Bins"), [](const config::RunConfig &c) { return static_cast<int>(c.trim.depth_bins); },
                [](config::RunConfig &c, int v) { c.trim.depth_bins = static_cast<std::size_t>(v); });
  binder_->check(range3d, tr("RANGE_3D"), [](const config::RunConfig &c) { return c.trim.record_ion_positions; },
                 [](config::RunConfig &c, bool on) { c.trim.record_ion_positions = on; });
  binder_->check(exits, tr("Exit lists"), [](const config::RunConfig &c) { return c.trim.record_exiting_atoms; },
                 [](config::RunConfig &c, bool on) { c.trim.record_exiting_atoms = on; });
  binder_->spin(coll, tr("Collision details"),
                [](const config::RunConfig &c) { return static_cast<int>(std::min<std::size_t>(c.trim.collision_details_ions, 2'000'000'000)); },
                [](config::RunConfig &c, int v) { c.trim.collision_details_ions = static_cast<std::size_t>(v); });
  binder_->check(coll_atoms, tr("Cascade atoms"), [](const config::RunConfig &c) { return c.trim.collision_recoil_details; },
                 [](config::RunConfig &c, bool on) { c.trim.collision_recoil_details = on; });
  binder_->text(exyz, tr("EXYZ step"), [](const config::RunConfig &c) { return format_energy(c.trim.exyz_interval); },
                [](config::RunConfig &c, const QString &t) {
                  const auto e = parse_energy(t);
                  if (!(e >= 0.0 * common::electron_volt)) throw std::invalid_argument("must be >= 0");
                  c.trim.exyz_interval = e;
                });
  binder_->spin(paths, tr("Paths drawn"),
                [](const config::RunConfig &c) {
                  return static_cast<int>(std::min<std::size_t>(c.trim.trajectory_ions, 1'000'000));
                },
                [](config::RunConfig &c, int v) { c.trim.trajectory_ions = static_cast<std::size_t>(v); });
  binder_->text(fluence, tr("dpa fluence"),
                [](const config::RunConfig &c) { return QString::number(c.trim.dpa_fluence_per_cm2, 'g', 6); },
                [](config::RunConfig &c, const QString &t) {
                  bool ok = false;
                  const double v = t.trimmed().toDouble(&ok);
                  if (!ok) throw std::invalid_argument("ions per cm2, e.g. 1e15");
                  c.trim.dpa_fluence_per_cm2 = v;
                });
  binder_->check(all_paths, tr("Every ion's path"),
                 [](const config::RunConfig &c) { return c.trim.trajectory_ions == trim::TrimConfig::all_ions; },
                 [](config::RunConfig &c, bool on) { c.trim.trajectory_ions = on ? trim::TrimConfig::all_ions : 10; });
  binder_->text(traj_file, tr("Trajectory file"),
                [](const config::RunConfig &c) { return QString::fromStdString(c.trim.trajectory_export.file); },
                [](config::RunConfig &c, const QString &t) {
                  c.trim.trajectory_export.file = t.trimmed().toStdString();
                  if (!c.trim.trajectory_export.file.empty() && c.trim.trajectory_ions == 0) {
                    c.trim.trajectory_ions = trim::TrimConfig::all_ions; // the file needs recorded paths
                  }
                });
  binder_->check(births, tr("Births only"), [](const config::RunConfig &c) { return c.trim.trajectory_export.births_only; },
                 [](config::RunConfig &c, bool on) { c.trim.trajectory_export.births_only = on; });
  binder_->text(traj_emin, tr("Trajectory minimum energy"),
                [](const config::RunConfig &c) { return format_energy(c.trim.trajectory_export.min_energy); },
                [](config::RunConfig &c, const QString &t) { c.trim.trajectory_export.min_energy = parse_energy(t); });
  binder_->text(md_file, tr("MD event file"), [](const config::RunConfig &c) { return QString::fromStdString(c.trim.md_events.file); },
                [](config::RunConfig &c, const QString &t) {
                  c.trim.md_events.file = t.trimmed().toStdString();
                  if (c.trim.md_events.enabled() && c.trim.md_events.transition_energies.empty()) {
                    c.trim.md_events.transition_energies = {1000.0 * common::electron_volt};
                  }
                });
  binder_->text(md_ec, tr("E_c"),
                [](const config::RunConfig &c) {
                  QStringList parts;
                  for (const auto &e : c.trim.md_events.transition_energies) parts << format_energy(e);
                  return parts.join(QStringLiteral(", "));
                },
                [](config::RunConfig &c, const QString &t) {
                  c.trim.md_events.transition_energies = config::parse_energy_list(t.toStdString(), "E_c");
                });
  binder_->text(md_emin, tr("E_min"),
                [](const config::RunConfig &c) {
                  return c.trim.md_events.min_energy > 0.0 * common::electron_volt ? format_energy(c.trim.md_events.min_energy)
                                                                                    : QStringLiteral("Ed");
                },
                [](config::RunConfig &c, const QString &t) {
                  const QString v = t.trimmed();
                  c.trim.md_events.min_energy = v.compare(QStringLiteral("ed"), Qt::CaseInsensitive) == 0
                                                    ? 0.0 * common::electron_volt
                                                    : parse_energy(v);
                });
  binder_->combo(md_mode, tr("MD mode"), [](const config::RunConfig &c) { return c.trim.md_events.handoff ? 1 : 0; },
                 [](config::RunConfig &c, int i) { c.trim.md_events.handoff = i == 1; });
  binder_->custom([paths, md, births, traj_emin, md_ec, md_emin, md_mode](const config::RunConfig &c) {
    paths->setEnabled(c.trim.trajectory_ions != trim::TrimConfig::all_ions);
    md->setVisible(c.mode == config::Mode::srim);
    const bool traj = !c.trim.trajectory_export.file.empty();
    births->setEnabled(traj);
    traj_emin->setEnabled(traj);
    const bool ev = c.trim.md_events.enabled();
    md_ec->setEnabled(ev);
    md_emin->setEnabled(ev);
    md_mode->setEnabled(ev);
  });
  binder_->custom([this](const config::RunConfig &c) {
    const auto &s = c.trim.outputs;
    QStringList parts;
    if (!s.units.empty()) {
      QStringList u;
      for (const auto &kv : s.units) u << QString::fromStdString(kv.second);
      parts << tr("units %1").arg(u.join(QStringLiteral(", ")));
    }
    parts << tr("numbers: %1, %2 digits").arg(QString::fromStdString(output::to_string(s.format.notation))).arg(s.format.digits);
    int excluded = 0, shaped = 0;
    for (const auto &kv : s.tables) {
      if (!kv.second.include) ++excluded;
      else if (!kv.second.columns.empty() || !kv.second.filters.empty() || !kv.second.column_specs.empty() ||
               kv.second.drop_zero_rows || kv.second.every > 1 || kv.second.max_rows > 0)
        ++shaped;
    }
    if (excluded) parts << tr("%n table(s) left out", nullptr, excluded);
    if (shaped) parts << tr("%n table(s) customised", nullptr, shaped);
    summary_->setText(parts.join(QStringLiteral("; ")) + QLatin1Char('.'));
  });
  binder_->refresh();
}

void OutputPanel::customize() {
  OutputSpecDialog d(scenario_->config().trim.outputs, sample_ ? sample_() : nullptr, this);
  if (d.exec() != QDialog::Accepted) return;
  const auto spec = d.spec();
  binder_->apply(tr("Output settings"), [&](config::RunConfig &c) { c.trim.outputs = spec; }, -1);
}

} // namespace gui
} // namespace amcset
