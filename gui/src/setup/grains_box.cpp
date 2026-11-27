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

#include "setup/grains_box.h"

#include "app/guard.h"
#include "crystal/grain_io.h"
#include "crystal/grain_summary.h"
#include "model/quantity_text.h"
#include "model/scenario.h"
#include "setup/form_binder.h"
#include "views/axis_ticks.h"
#include "views/sample_scene.h"

#include <QCoreApplication>

#include <QCheckBox>
#include <QComboBox>
#include <QDir>
#include <QFileDialog>
#include <QFileInfo>
#include <QFormLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPixmap>
#include <QPushButton>
#include <QSignalBlocker>
#include <QSpinBox>
#include <QThread>
#include <QTimer>

#include <cmath>
#include <filesystem>
#include <stdexcept>

namespace amcset {
namespace gui {

namespace {

using crystal::Polycrystal;

QString angle_text(common::angle_quantity a) {
  return QStringLiteral("%1 deg").arg(format_number(common::in_units(a, common::degree), 6));
}

//! The grain file \p p names, relative to the scenario's folder.
std::filesystem::path grain_file_path(const Polycrystal &p, const QString &scenario_file) {
  const std::filesystem::path f = std::filesystem::path(p.file);
  if (f.is_absolute() || scenario_file.isEmpty()) return f;
  return std::filesystem::path(QFileInfo(scenario_file).absolutePath().toStdWString()) / f;
}

int source_index(Polycrystal::Source s) {
  switch (s) {
  case Polycrystal::Source::voronoi: return 1;
  case Polycrystal::Source::grain_file: return 2;
  case Polycrystal::Source::neper: return 3;
  case Polycrystal::Source::ebsd: return 4;
  }
  return 1;
}

} // namespace

GrainsBox::GrainsBox(Scenario *scenario, FormBinder *binder, std::function<int()> layer, QWidget *parent)
    : QGroupBox(tr("Grains (polycrystal)"), parent), scenario_(scenario), binder_(binder), layer_(std::move(layer)) {
  source_ = new QComboBox(this);
  source_->addItems({tr("Single crystal"), tr("Generated grains (Voronoi)"), tr("AMCSET grain file"),
                     tr("Neper tessellation (.tess)"), tr("EBSD map (.ang, .ctf)")});
  source_->setToolTip(tr("A single crystal, or grains of this unit cell in many orientations: generated (a "
                         "Poisson-Voronoi tessellation with a mean grain size and a texture), or read from a grain "
                         "file, a Neper tessellation or an EBSD map (whose grains become columns into the depth). "
                         "The crystal transport follows the lattice of the grain each atom is in. "
                         "(layers.N.crystal.polycrystal)"));
  grain_size_ = new QLineEdit(this);
  grain_size_->setToolTip(tr("Mean grain diameter (of the sphere of the mean grain volume; of the circle of the mean "
                             "section for columns), 10 A or more. (polycrystal.grain_size)"));
  columnar_ = new QCheckBox(tr("Columns through the layer"), this);
  columnar_->setToolTip(tr("Grains that are columns along the depth, as in many deposited films. (polycrystal.columnar)"));
  texture_ = new QComboBox(this);
  texture_->addItems({tr("Random"), tr("Fiber (the surface plane kept)"), tr("Mosaic (the orientation kept)")});
  texture_->setToolTip(tr("Orientations of the grains: uniformly random; a fiber texture, every grain with the "
                          "crystal's surface plane, turned at random about the normal; or a mosaic, every grain "
                          "near the crystal's orientation. The spread widens the last two. (polycrystal.texture)"));
  spread_ = new QLineEdit(this);
  spread_->setToolTip(tr("Fiber: rms tilt of the fiber axis; mosaic: rms rotation of the grains. (polycrystal.spread)"));
  seed_ = new QSpinBox(this);
  seed_->setRange(0, 2'000'000'000);
  seed_->setToolTip(tr("Seed of the grain generator: the same grains for every run; another seed, other grains of the "
                       "same kind. Not the run's random seed. (polycrystal.seed)"));
  period_ = new QLineEdit(this);
  period_->setPlaceholderText(tr("automatic"));
  period_->setToolTip(tr("The lateral period of the grains: a patch of the surface that repeats, one length or "
                         "\"y, z\". Empty: eight grain sizes for generated grains, the file's map otherwise. The "
                         "beam spot should cover it. (polycrystal.period)"));
  file_ = new QLabel(this);
  file_->setWordWrap(true);
  file_->setTextInteractionFlags(Qt::TextSelectableByMouse);
  choose_ = new QPushButton(tr("Choose..."), this);
  choose_->setToolTip(tr("The file of the grains (polycrystal.file)"));
  auto *file_row = new QWidget(this);
  auto *fr = new QHBoxLayout(file_row);
  fr->setContentsMargins(0, 0, 0, 0);
  fr->addWidget(file_, 1);
  fr->addWidget(choose_);
  scale_ = new QLineEdit(this);
  scale_->setToolTip(tr("Neper: the length of one unit of the tessellation. (polycrystal.scale)"));
  misorientation_ = new QLineEdit(this);
  misorientation_->setToolTip(tr("EBSD: neighbouring points closer in orientation than this (up to the lattice's own "
                                 "rotations) are one grain. (polycrystal.misorientation)"));
  phase_ = new QSpinBox(this);
  phase_->setRange(0, 1'000'000);
  phase_->setSpecialValueText(tr("every phase"));
  phase_->setToolTip(tr("EBSD: the phase of the map that is this crystal; the points of other phases take the "
                        "grain of the nearest point of it. (polycrystal.phase)"));
  boundary_ = new QLineEdit(this);
  boundary_->setToolTip(tr("Width of a disordered slab at every grain boundary, followed as amorphous (random "
                           "collision partners). 0: none; at most 100 A. (polycrystal.boundary_width)"));
  registry_ = new QLineEdit(this);
  registry_->setPlaceholderText(tr("off"));
  registry_->setToolTip(tr("Low-angle boundaries in registry: neighbouring grains misoriented by less than this "
                           "angle continue each other's lattice across their boundary, as at a subgrain boundary "
                           "(grains of one orientation are then one crystal). Off: each grain's lattice is its own, "
                           "and every boundary dechannels. At most 30 deg. (polycrystal.registry_below)"));
  map_ = new QLabel(this);
  map_->setAlignment(Qt::AlignCenter);
  map_->setMinimumHeight(40);
  map_->setToolTip(tr("The grains at the surface over one period, coloured by the beam direction in their crystal "
                      "(cubic: [001] red, [101] green, [111] blue)"));
  summary_ = new QLabel(this);
  summary_->setWordWrap(true);
  summary_->setTextInteractionFlags(Qt::TextSelectableByMouse);

  form_ = new QFormLayout(this);
  form_->addRow(tr("Grains"), source_);
  form_->addRow(tr("Grain size"), grain_size_);
  form_->addRow(QString(), columnar_);
  form_->addRow(tr("Texture"), texture_);
  form_->addRow(tr("Spread"), spread_);
  form_->addRow(tr("Generator seed"), seed_);
  form_->addRow(tr("File"), file_row);
  form_->addRow(tr("Unit length"), scale_);
  form_->addRow(tr("Same grain within"), misorientation_);
  form_->addRow(tr("Phase"), phase_);
  form_->addRow(tr("Period"), period_);
  form_->addRow(tr("Disordered boundaries"), boundary_);
  form_->addRow(tr("In registry below"), registry_);
  form_->addRow(map_);
  form_->addRow(summary_);

  timer_ = new QTimer(this);
  timer_->setSingleShot(true);
  timer_->setInterval(400);
  connect(timer_, &QTimer::timeout, this, guarded(this, tr("Grains (preview)"), [this] { start_preview(); }));
  connect(source_, &QComboBox::activated, this, guarded(this, tr("Grains (source)"), [this](int i) { choose_source(i); }));
  connect(choose_, &QPushButton::clicked, this, guarded(this, tr("Grains (file)"), [this] { choose_file(); }));

  // Fields of the polycrystal (empty while the layer is a single crystal).
  const auto text = [this](QLineEdit *w, const QString &what, std::function<QString(const Polycrystal &)> get,
                           std::function<void(Polycrystal &, const QString &)> set, bool rereads = false) {
    binder_->text(w, what,
                  [this, get](const config::RunConfig &c) {
                    const auto *p = current(c);
                    return p ? get(*p) : QString();
                  },
                  [this, set, rereads](config::RunConfig &c, const QString &t) {
                    auto &p = edited(c);
                    set(p, t);
                    if (rereads && p.source != Polycrystal::Source::voronoi) {
                      const auto &l = c.trim.layers.at(static_cast<std::size_t>(layer_()));
                      p.map = crystal::read_grain_map(p, grain_file_path(p, scenario_->file_path()), l.crystal->cell);
                    }
                  });
  };
  text(grain_size_, tr("Grain size"), [](const Polycrystal &p) { return format_length(p.grain_size); },
       [](Polycrystal &p, const QString &t) { p.grain_size = parse_length(t); });
  text(spread_, tr("Spread"), [](const Polycrystal &p) { return angle_text(p.spread); },
       [](Polycrystal &p, const QString &t) { p.spread = parse_angle(t); });
  text(period_, tr("Period"),
       [](const Polycrystal &p) {
         if (!(p.period_y > 0.0 * common::angstrom) && !(p.period_z > 0.0 * common::angstrom)) return QString();
         if (p.period_y == p.period_z) return format_length(p.period_y);
         return QStringLiteral("%1, %2").arg(format_length(p.period_y), format_length(p.period_z));
       },
       [](Polycrystal &p, const QString &t) {
         const QString s = t.trimmed();
         if (s.isEmpty() || s.compare(QStringLiteral("automatic"), Qt::CaseInsensitive) == 0) {
           p.period_y = p.period_z = 0.0 * common::angstrom;
           return;
         }
         const auto parts = s.split(QLatin1Char(','));
         if (parts.size() > 2) throw std::invalid_argument("one length, or \"y, z\"");
         p.period_y = parse_length(parts[0].trimmed());
         p.period_z = parts.size() == 2 ? parse_length(parts[1].trimmed()) : p.period_y;
       });
  text(scale_, tr("Unit length"), [](const Polycrystal &p) { return format_length(p.scale); },
       [](Polycrystal &p, const QString &t) { p.scale = parse_length(t, "um"); }, true);
  text(misorientation_, tr("Same grain within"), [](const Polycrystal &p) { return angle_text(p.misorientation); },
       [](Polycrystal &p, const QString &t) { p.misorientation = parse_angle(t); }, true);
  text(boundary_, tr("Disordered boundaries"), [](const Polycrystal &p) { return format_length(p.boundary_width); },
       [](Polycrystal &p, const QString &t) { p.boundary_width = t.trimmed().isEmpty() ? 0.0 * common::angstrom : parse_length(t); });
  text(registry_, tr("In registry below"),
       [](const Polycrystal &p) { return p.registry_below > 0.0 * common::degree ? angle_text(p.registry_below) : QString(); },
       [](Polycrystal &p, const QString &t) {
         const QString s = t.trimmed();
         p.registry_below = s.isEmpty() || s.compare(QStringLiteral("off"), Qt::CaseInsensitive) == 0
                                ? 0.0 * common::degree
                                : parse_angle(s);
       });
  binder_->check(columnar_, tr("Columns"), [this](const config::RunConfig &c) {
    const auto *p = current(c);
    return p != nullptr && p->columnar;
  }, [this](config::RunConfig &c, bool on) { edited(c).columnar = on; });
  binder_->combo(texture_, tr("Texture"), [this](const config::RunConfig &c) {
    const auto *p = current(c);
    return p ? static_cast<int>(p->texture) : 0;
  }, [this](config::RunConfig &c, int i) { edited(c).texture = static_cast<Polycrystal::Texture>(i); });
  binder_->spin(seed_, tr("Generator seed"), [this](const config::RunConfig &c) {
    const auto *p = current(c);
    return p ? static_cast<int>(std::min<std::uint64_t>(p->seed, 2'000'000'000)) : 1;
  }, [this](config::RunConfig &c, int v) { edited(c).seed = static_cast<std::uint64_t>(v); });
  binder_->spin(phase_, tr("Phase"), [this](const config::RunConfig &c) {
    const auto *p = current(c);
    return p ? p->phase : 0;
  }, [this](config::RunConfig &c, int v) {
    auto &p = edited(c);
    p.phase = v;
    if (p.source == Polycrystal::Source::ebsd) {
      const auto &l = c.trim.layers.at(static_cast<std::size_t>(layer_()));
      p.map = crystal::read_grain_map(p, grain_file_path(p, scenario_->file_path()), l.crystal->cell);
    }
  });
  binder_->custom([this](const config::RunConfig &c) { refresh(c); });
}

GrainsBox::~GrainsBox() {
  if (cancel_) cancel_->store(true);
  if (worker_) worker_->wait();
}

const crystal::Polycrystal *GrainsBox::current(const config::RunConfig &c) const {
  const int i = layer_();
  if (i < 0 || i >= static_cast<int>(c.trim.layers.size())) return nullptr;
  const auto &l = c.trim.layers[static_cast<std::size_t>(i)];
  return l.crystal && l.crystal->polycrystal ? &*l.crystal->polycrystal : nullptr;
}

crystal::Polycrystal &GrainsBox::edited(config::RunConfig &c) const {
  const int i = layer_();
  if (i < 0 || i >= static_cast<int>(c.trim.layers.size())) throw std::invalid_argument("no layer is selected");
  auto &l = c.trim.layers[static_cast<std::size_t>(i)];
  if (!l.crystal) throw std::invalid_argument("the layer has no crystal");
  if (!l.crystal->polycrystal) l.crystal->polycrystal = Polycrystal{};
  return *l.crystal->polycrystal;
}

void GrainsBox::refresh(const config::RunConfig &c) {
  const int i = layer_();
  const trim::TargetLayer *l = i >= 0 && i < static_cast<int>(c.trim.layers.size()) ? &c.trim.layers[static_cast<std::size_t>(i)] : nullptr;
  const bool crystal = l != nullptr && l->crystal.has_value();
  const auto *p = current(c);
  setEnabled(crystal);
  {
    const QSignalBlocker block(source_);
    source_->setCurrentIndex(p ? source_index(p->source) : 0);
  }
  const bool generated = p && p->source == Polycrystal::Source::voronoi;
  const bool from_file = p && !generated;
  for (QWidget *w : {static_cast<QWidget *>(grain_size_), static_cast<QWidget *>(columnar_),
                     static_cast<QWidget *>(texture_), static_cast<QWidget *>(spread_), static_cast<QWidget *>(seed_)}) {
    form_->setRowVisible(w, generated);
  }
  form_->setRowVisible(choose_->parentWidget(), from_file);
  form_->setRowVisible(scale_, p && p->source == Polycrystal::Source::neper);
  form_->setRowVisible(misorientation_, p && p->source == Polycrystal::Source::ebsd);
  form_->setRowVisible(phase_, p && p->source == Polycrystal::Source::ebsd);
  for (QWidget *w : {static_cast<QWidget *>(period_), static_cast<QWidget *>(boundary_), static_cast<QWidget *>(registry_),
                     static_cast<QWidget *>(map_), static_cast<QWidget *>(summary_)}) {
    form_->setRowVisible(w, p != nullptr);
  }
  spread_->setEnabled(generated && p->texture != Polycrystal::Texture::random);
  if (from_file) {
    file_->setText(p->file.empty() ? tr("(none)") : QDir::toNativeSeparators(QString::fromStdString(p->file)));
  }
  // The preview, when the grains (or what colours them) changed.
  if (p == nullptr) {
    previewed_.reset();
    map_->clear();
    summary_->clear();
    if (cancel_) cancel_->store(true);
    return;
  }
  const double width = common::in_units(l->width, common::angstrom);
  const double tilt = common::in_units(c.trim.ion.angle, common::degree);
  const bool same = previewed_ && *previewed_ == *l->crystal && previewed_->polycrystal->map == p->map &&
                    previewed_width_ == width && previewed_tilt_ == tilt;
  if (!same) {
    previewed_ = *l->crystal;
    previewed_width_ = width;
    previewed_tilt_ = tilt;
    summary_->setText(tr("Building the grains..."));
    timer_->start();
  }
}

void GrainsBox::choose_source(int index) {
  const int li = layer_();
  if (li < 0) return;
  if (index == 0) {
    binder_->apply(tr("Single crystal"), [this](config::RunConfig &c) {
      c.trim.layers.at(static_cast<std::size_t>(layer_())).crystal.value().polycrystal.reset();
    }, -1);
    return;
  }
  if (index == 1) {
    binder_->apply(tr("Generated grains"), [this](config::RunConfig &c) {
      auto &p = edited(c);
      p.source = Polycrystal::Source::voronoi;
      p.map.reset();
    }, -1);
    return;
  }
  // A file source: which file?
  static const Polycrystal::Source sources[] = {Polycrystal::Source::grain_file, Polycrystal::Source::neper,
                                                Polycrystal::Source::ebsd};
  const QString filters[] = {tr("Grain files (*.txt *.grains);;All files (*)"), tr("Neper tessellations (*.tess);;All files (*)"),
                             tr("EBSD maps (*.ang *.ctf);;All files (*)")};
  const QString f = QFileDialog::getOpenFileName(this, tr("File of the grains"), QFileInfo(scenario_->file_path()).absolutePath(),
                                                 filters[index - 2]);
  if (f.isEmpty()) {
    refresh(scenario_->config()); // (back to what the scenario says)
    return;
  }
  const auto source = sources[index - 2];
  binder_->apply(tr("Grains from %1").arg(QFileInfo(f).fileName()), [this, f, source](config::RunConfig &c) {
    auto &p = edited(c);
    p.source = source;
    p.file = QDir::toNativeSeparators(f).toStdString();
    const auto &l = c.trim.layers.at(static_cast<std::size_t>(layer_()));
    p.map = crystal::read_grain_map(p, std::filesystem::path(f.toStdWString()), l.crystal->cell);
  }, -1);
}

void GrainsBox::choose_file() {
  const auto *p = current(scenario_->config());
  if (p == nullptr || p->source == Polycrystal::Source::voronoi) return;
  choose_source(source_index(p->source));
}

void GrainsBox::start_preview() {
  const auto &c = scenario_->config();
  const int i = layer_();
  if (i < 0 || i >= static_cast<int>(c.trim.layers.size())) return;
  const auto &l = c.trim.layers[static_cast<std::size_t>(i)];
  if (!l.crystal || !l.crystal->polycrystal) return;
  if (cancel_) cancel_->store(true);
  if (worker_) {
    worker_->wait();
    worker_ = nullptr;
  }
  cancel_ = std::make_shared<std::atomic<bool>>(false);
  const quint64 generation = ++generation_;
  const crystal::Crystal cr = *l.crystal;
  const double width = common::in_units(l.width, common::angstrom);
  const double tilt = common::in_units(c.trim.ion.angle, common::radian);
  worker_ = QThread::create([this, cr, width, tilt, cancel = cancel_, generation] {
    Preview out;
    try {
      crystal::RegistryResult registry;
      const auto map = crystal::make_grain_map(cr, width, &registry);
      if (cancel->load()) return;
      // The surface over one period, coloured by the beam direction in each grain.
      const double py = map->period_y_A(), pz = map->period_z_A();
      const int w = 240;
      const int h = std::clamp(static_cast<int>(std::lround(w * pz / py)), 40, 360);
      QImage img(w, h, QImage::Format_RGB32);
      const crystal::Vec3 beam{std::cos(tilt), std::sin(tilt), 0.0};
      const bool cubic = crystal::lattice_rotations(cr.cell.lattice).size() == 24;
      std::vector<QRgb> colors(map->grain_count(), 0);
      std::vector<char> known(map->grain_count(), 0);
      for (int y = 0; y < h && !cancel->load(); ++y) {
        for (int x = 0; x < w; ++x) {
          const int g = map->at({0.0, (x + 0.5) / w * py, (h - y - 0.5) / h * pz}).grain;
          QRgb col = qRgb(150, 150, 150);
          if (g >= 0) {
            const auto gi = static_cast<std::size_t>(g);
            if (!known[gi]) {
              const auto &r = map->grain(g).rotation;
              colors[gi] = cubic ? ipf_color(r.transposed() * beam).rgb()
                                 : QColor::fromHsvF(static_cast<float>(std::fmod(0.11 + 0.618034 * g, 1.0)), 0.55f, 0.92f).rgb();
              known[gi] = 1;
            }
            col = colors[gi];
          }
          img.setPixel(x, y, col);
        }
      }
      if (cancel->load()) return;
      out.map = img;
      const auto s = crystal::summarize_grains(*map, cr.cell.lattice, cr.rotation(), width, 400'000);
      if (cancel->load()) return;
      QString text = QCoreApplication::translate("GrainsBox", "%L1 grains (%L2 in the layer); period %3 x %4; mean "
                                                              "diameter %5 (10-90 %: %6 to %7)%8")
                         .arg(s.grains)
                         .arg(s.grains_in_layer)
                         .arg(format_length_value(s.period_y_A), format_length_value(s.period_z_A),
                              format_length_value(s.diameter_mean_A), format_length_value(s.diameter_p10_A),
                              format_length_value(s.diameter_p90_A),
                              s.coarse ? QCoreApplication::translate("GrainsBox", " (coarse sample)") : QString());
      text += QCoreApplication::translate("GrainsBox", ". Misorientation from the crystal's orientation: mean %1 deg")
                  .arg(format_number(s.misorientation_mean_deg, 3));
      if (s.boundaries_counted) {
        text += QCoreApplication::translate("GrainsBox", "; %L1 boundaries, %2 % below 15 deg")
                    .arg(s.boundaries)
                    .arg(format_number(100.0 * s.low_angle_fraction, 3));
      }
      if (registry.groups > 0) {
        text += QCoreApplication::translate("GrainsBox", "; in registry: %L1 groups of %L2 grains").arg(registry.groups).arg(registry.grains);
      }
      out.summary = text + QStringLiteral(".");
    } catch (const std::exception &e) {
      out.summary = QCoreApplication::translate("GrainsBox", "The grains cannot be built: %1").arg(QString::fromUtf8(e.what()));
    }
    QMetaObject::invokeMethod(this, [this, generation, out] { preview_done(generation, out); }, Qt::QueuedConnection);
  });
  connect(worker_, &QThread::finished, worker_, &QObject::deleteLater);
  connect(worker_, &QThread::finished, this, [this, w = worker_] {
    if (worker_ == w) worker_ = nullptr;
  });
  worker_->start(QThread::LowPriority);
}

void GrainsBox::preview_done(quint64 generation, const Preview &p) {
  if (generation != generation_) return;
  if (!p.map.isNull()) map_->setPixmap(QPixmap::fromImage(p.map));
  summary_->setText(p.summary);
}

} // namespace gui
} // namespace amcset
