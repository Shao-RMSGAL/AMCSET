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
#include "views/track_panel.h"

#include "run/track_store.h"
#include "views/sample_scene.h"
#include "views/track_view.h"

#include <QCheckBox>
#include <QComboBox>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QPainter>
#include <QPushButton>
#include <QSettings>
#include <QSlider>
#include <QSpinBox>
#include <QThread>
#include <QTimer>
#include <QVBoxLayout>

namespace amcset {
namespace gui {

TrackPanel::TrackPanel(TrackStore *store, QWidget *parent) : QWidget(parent), store_(store) {
  view_ = new TrackView(store, this);
  arrangement_ = new QComboBox(this);
  arrangement_->setSizeAdjustPolicy(QComboBox::AdjustToContents);
  arrangement_->addItem(tr("Four views"));
  arrangement_->addItem(tr("Depth - Y"));
  arrangement_->addItem(tr("Depth - Z"));
  arrangement_->addItem(tr("Beam view (Y - Z)"));
  arrangement_->addItem(tr("3D"));
  color_by_ = new QComboBox(this);
  color_by_->addItem(tr("By element"));
  color_by_->addItem(tr("By generation"));
  color_by_->setToolTip(tr("Colour of the paths"));
  incident_ = new QCheckBox(tr("Incident"), this);
  incident_->setChecked(true);
  recoils_ = new QCheckBox(tr("Recoils to gen."), this);
  recoils_->setToolTip(tr("Show recoil paths up to this generation (1 = primary recoils)"));
  recoils_->setChecked(true);
  generation_ = new QSpinBox(this);
  generation_->setRange(1, 99);
  generation_->setValue(99);
  generation_->setToolTip(tr("1 = primary recoils only"));
  opacity_ = new QSlider(Qt::Horizontal, this);
  opacity_->setRange(2, 100);
  opacity_->setValue(80);
  opacity_->setMaximumWidth(120);
  opacity_->setToolTip(tr("Path opacity: lower values show where many paths overlap"));
  auto *reset = new QPushButton(tr("Fit"), this);
  reset->setToolTip(tr("Fit all views (double-click a view to fit it alone)"));
  auto *clear = new QPushButton(tr("Clear"), this);
  clear->setToolTip(tr("Remove the paths shown (the results are kept)"));
  auto *save = new QPushButton(tr("Save..."), this);
  save->setToolTip(tr("Save an image of the views"));
  spot_ = new QCheckBox(tr("Whole beam spot"), this);
  spot_->setToolTip(tr("Show the whole beam spot at true scale: every path at its true place on the sample (the "
                       "point where its ion entered the spot, marked on the surface), the spot's outline, and "
                       "the sample under it: the surface, sections through the spot's centre and a 3D block, "
                       "with the grains of a polycrystal coloured by the beam direction in their crystal. "
                       "Crystal transport only (the other transports have no entry points)."));
  spot_->setEnabled(false);
  scene_timer_ = new QTimer(this);
  scene_timer_->setSingleShot(true);
  scene_timer_->setInterval(300);

  auto *bar = new QHBoxLayout();
  bar->addWidget(arrangement_);
  bar->addWidget(color_by_);
  bar->addSpacing(8);
  bar->addWidget(spot_);
  bar->addSpacing(8);
  bar->addWidget(incident_);
  bar->addWidget(recoils_);
  bar->addWidget(generation_);
  bar->addSpacing(8);
  bar->addWidget(opacity_);
  bar->addStretch(1);
  bar->addWidget(reset);
  bar->addWidget(clear);
  bar->addWidget(save);
  auto *lay = new QVBoxLayout(this);
  lay->setContentsMargins(0, 0, 0, 0);
  lay->addLayout(bar);
  lay->addWidget(view_, 1);

  connect(arrangement_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Path views (currentIndexChanged)"), [this](int i) {
    using P = TrackView::Projection;
    if (i == 0) {
      view_->set_arrangement(TrackView::Arrangement::quad);
    } else {
      static const P ps[] = {P::depth_y, P::depth_z, P::beam, P::perspective};
      view_->set_single_projection(ps[i - 1]);
      view_->set_arrangement(TrackView::Arrangement::single);
    }
  }));
  connect(color_by_, &QComboBox::currentIndexChanged, this, guarded(this, tr("Path views (currentIndexChanged)"), [this](int i) {
    view_->set_color_by(i == 0 ? TrackView::ColorBy::element : TrackView::ColorBy::generation);
  }));
  connect(incident_, &QCheckBox::toggled, view_,
          guarded(this, tr("Path views (incident ions)"), [v = view_](bool on) { v->set_show_incident(on); }));
  connect(recoils_, &QCheckBox::toggled, view_,
          guarded(this, tr("Path views (recoils)"), [v = view_](bool on) { v->set_show_recoils(on); }));
  connect(recoils_, &QCheckBox::toggled, generation_, &QWidget::setEnabled);
  connect(generation_, &QSpinBox::valueChanged, view_,
          guarded(this, tr("Path views (generations)"), [v = view_](int g) { v->set_max_generation(g); }));
  connect(opacity_, &QSlider::valueChanged, this, guarded(this, tr("Path views (valueChanged)"), [this](int v) { view_->set_opacity(v / 100.0); }));
  connect(reset, &QPushButton::clicked, view_, guarded(this, tr("Path views (reset)"), [v = view_] { v->reset_views(); }));
  connect(clear, &QPushButton::clicked, this, guarded(this, tr("Path views (clicked)"), [this] { clear_requested(); }));
  connect(save, &QPushButton::clicked, this, guarded(this, tr("Path views (clicked)"), [this] {
    const QString f = QFileDialog::getSaveFileName(this, tr("Save paths image"), QStringLiteral("paths.png"),
                                                   tr("Images (*.png *.jpg *.bmp)"));
    if (f.isEmpty()) return;
    const QImage img = grab_image();
    if (!img.isNull()) img.save(f);
  }));
  connect(view_, &TrackView::hover_text, this, guarded(this, tr("Path views (hover text)"), [this](auto &&...a) { hover_text(std::forward<decltype(a)>(a)...); }));
  connect(spot_, &QCheckBox::toggled, this, guarded(this, tr("Path views (whole beam spot)"), [this](bool on) {
    view_->set_spot_mode(on);
    if (on && scene_stale_) build_scene();
  }));
  connect(scene_timer_, &QTimer::timeout, this, guarded(this, tr("Path views (sample)"), [this] { build_scene(); }));
}

TrackPanel::~TrackPanel() {
  // A sample still being drawn: stop it and wait (it checks often).
  if (cancel_) cancel_->store(true);
  if (worker_) worker_->wait();
}

void TrackPanel::set_target(const std::vector<double> &boundaries, unsigned incident_z) {
  view_->set_target(boundaries, incident_z);
}

void TrackPanel::set_scenario(const trim::TrimConfig &c, bool crystal_transport) {
  scenario_ = c;
  have_scenario_ = true;
  scene_stale_ = true;
  spot_->setEnabled(crystal_transport);
  if (!crystal_transport && spot_->isChecked()) spot_->setChecked(false);
  spot_->setToolTip(spot_->toolTip().section(QStringLiteral("\n\n"), 0, 0) +
                    (crystal_transport ? QString()
                                       : tr("\n\nNot available: this calculation does not use the crystal "
                                            "transport (Calculation tab), so its ions have no entry points.")));
  if (spot_->isChecked()) scene_timer_->start(); // (a burst of edits draws once)
}

void TrackPanel::build_scene() {
  if (!have_scenario_) return;
  scene_stale_ = false;
  if (cancel_) cancel_->store(true);
  if (worker_) {
    worker_->wait(); // (cancelled: returns within a row of cells)
    worker_ = nullptr;
  }
  cancel_ = std::make_shared<std::atomic<bool>>(false);
  const quint64 generation = ++scene_generation_;
  view_->set_scene(nullptr);
  const QColor line = palette().color(QPalette::WindowText);
  worker_ = QThread::create([this, c = scenario_, line, cancel = cancel_, generation] {
    auto scene = build_sample_scene(c, line, cancel.get());
    if (cancel->load()) return;
    // (Queued to the GUI thread; dropped if this panel is gone first: the
    // destructor waits for this thread.)
    QMetaObject::invokeMethod(this, [this, generation, scene] { scene_built(generation, scene); }, Qt::QueuedConnection);
  });
  connect(worker_, &QThread::finished, worker_, &QObject::deleteLater);
  connect(worker_, &QThread::finished, this, [this, w = worker_] {
    if (worker_ == w) worker_ = nullptr;
  });
  worker_->start(QThread::LowPriority);
}

QImage TrackPanel::grab_image() const {
  // The GPU image plus the axes painted over it.
  QImage img = view_->grabFramebuffer();
  if (img.isNull()) return img;
  img = img.convertToFormat(QImage::Format_ARGB32_Premultiplied);
  img.setDevicePixelRatio(view_->devicePixelRatioF());
  if (auto *overlay = view_->findChild<QWidget *>()) {
    QPainter p(&img);
    overlay->render(&p, QPoint(), QRegion(), QWidget::DrawChildren);
  }
  return img;
}

bool TrackPanel::whole_spot_available() const { return spot_->isEnabled(); }
bool TrackPanel::whole_spot() const { return spot_->isChecked(); }
void TrackPanel::set_whole_spot(bool on) { spot_->setChecked(on); }
int TrackPanel::arrangement_index() const { return arrangement_->currentIndex(); }
void TrackPanel::set_arrangement_index(int i) { arrangement_->setCurrentIndex(i); }

void TrackPanel::scene_built(quint64 generation, std::shared_ptr<const SampleScene> scene) {
  if (generation != scene_generation_) return; // a newer scenario is being drawn
  view_->set_scene(std::move(scene));
}

void TrackPanel::data_changed() { view_->data_changed(); }

void TrackPanel::save_settings() const {
  QSettings s;
  s.beginGroup(QStringLiteral("paths"));
  s.setValue(QStringLiteral("arrangement"), arrangement_->currentIndex());
  s.setValue(QStringLiteral("color_by"), color_by_->currentIndex());
  s.setValue(QStringLiteral("opacity"), opacity_->value());
  s.setValue(QStringLiteral("whole_spot"), spot_->isChecked());
  s.endGroup();
}

void TrackPanel::restore_settings() {
  QSettings s;
  s.beginGroup(QStringLiteral("paths"));
  arrangement_->setCurrentIndex(s.value(QStringLiteral("arrangement"), 0).toInt());
  color_by_->setCurrentIndex(s.value(QStringLiteral("color_by"), 0).toInt());
  opacity_->setValue(s.value(QStringLiteral("opacity"), 80).toInt());
  if (spot_->isEnabled()) spot_->setChecked(s.value(QStringLiteral("whole_spot"), false).toBool());
  s.endGroup();
}

} // namespace gui
} // namespace amcset
