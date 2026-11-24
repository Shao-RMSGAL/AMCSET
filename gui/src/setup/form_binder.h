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

/*! \file form_binder.h
 *  \brief Two-way binding of form widgets to the scenario.
 *
 * Each binding has a getter (scenario -> widget) and a setter (widget ->
 * scenario, through Scenario::modify so the change is undoable). The binder
 * refreshes every widget when the scenario changes (with its signals
 * blocked) and reports rejected input.
 */

#pragma once

#include "app/guard.h"
#include "model/scenario.h"

#include <QCheckBox>
#include <QComboBox>
#include <QDoubleSpinBox>
#include <QLineEdit>
#include <QObject>
#include <QSignalBlocker>
#include <QSpinBox>

#include <functional>
#include <vector>

namespace amcset {
namespace gui {

class FormBinder : public QObject {
  Q_OBJECT

public:
  using Get = std::function<QString(const config::RunConfig &)>;
  using Set = std::function<void(config::RunConfig &, const QString &)>;

  FormBinder(Scenario *scenario, QObject *parent) : QObject(parent), scenario_(scenario) {
    connect(scenario_, &Scenario::changed, this, &FormBinder::refresh);
  }

  //! A text field (any quantity with units); applied when editing finishes.
  void text(QLineEdit *w, const QString &what, Get get, Set set) {
    const int id = next_id_++;
    refreshers_.push_back([w, get](const config::RunConfig &c) {
      if (!w->hasFocus() || w->text().isEmpty()) w->setText(get(c));
    });
    connect(w, &QLineEdit::editingFinished, this, guarded(nullptr, what, [this, w, what, get, set, id] {
      const QString v = w->text();
      if (v == get(scenario_->config())) return;
      apply(what, [&](config::RunConfig &c) { set(c, v); }, id);
      w->setText(get(scenario_->config()));
    }));
  }
  void check(QCheckBox *w, const QString &what, std::function<bool(const config::RunConfig &)> get,
             std::function<void(config::RunConfig &, bool)> set) {
    refreshers_.push_back([w, get](const config::RunConfig &c) { w->setChecked(get(c)); });
    connect(w, &QCheckBox::toggled, this, guarded(nullptr, what, [this, what, set](bool on) {
      apply(what, [&](config::RunConfig &c) { set(c, on); }, -1);
    }));
  }
  //! A combo box; get/set work with the item index.
  void combo(QComboBox *w, const QString &what, std::function<int(const config::RunConfig &)> get,
             std::function<void(config::RunConfig &, int)> set) {
    refreshers_.push_back([w, get](const config::RunConfig &c) { w->setCurrentIndex(get(c)); });
    connect(w, &QComboBox::currentIndexChanged, this, guarded(nullptr, what, [this, what, set](int i) {
      if (i < 0) return;
      apply(what, [&](config::RunConfig &c) { set(c, i); }, -1);
    }));
  }
  void spin(QSpinBox *w, const QString &what, std::function<int(const config::RunConfig &)> get,
            std::function<void(config::RunConfig &, int)> set) {
    const int id = next_id_++;
    refreshers_.push_back([w, get](const config::RunConfig &c) { w->setValue(get(c)); });
    connect(w, &QSpinBox::valueChanged, this, guarded(nullptr, what, [this, what, set, id](int v) {
      apply(what, [&](config::RunConfig &c) { set(c, v); }, id);
    }));
  }
  //! Anything else that only needs refreshing.
  void custom(std::function<void(const config::RunConfig &)> refresh) { refreshers_.push_back(std::move(refresh)); }

  //! Show the scenario in every widget (a slot of Scenario::changed: it
  //! must not throw). A widget whose getter fails keeps its old text and the
  //! failure is logged; the others are still refreshed.
  void refresh() noexcept {
    // Widgets set here emit their change signals; apply() ignores them.
    // The flag is reset on every path (a stuck flag would ignore all edits).
    struct Flag {
      bool &f;
      explicit Flag(bool &x) : f(x) { f = true; }
      ~Flag() { f = false; }
    } flag(refreshing_);
    const auto &c = scenario_->config();
    for (auto &r : refreshers_) {
      try {
        r(c);
      } catch (...) {
        report_exception_quietly(QStringLiteral("Showing a setup field"), std::current_exception());
      }
    }
  }

  //! Change the scenario; reports a rejected value through error().
  bool apply(const QString &what, const std::function<void(config::RunConfig &)> &edit, int merge_id) {
    if (refreshing_) return false;
    const QString err = scenario_->modify(what, edit, merge_id < 0 ? -1 : 5000 + merge_id);
    if (!err.isEmpty()) {
      emit error(tr("%1: %2").arg(what, err));
      refresh();
      return false;
    }
    return true;
  }

signals:
  void error(const QString &message);

private:
  Scenario *scenario_;
  std::vector<std::function<void(const config::RunConfig &)>> refreshers_;
  bool refreshing_ = false;
  int next_id_ = 0;
};

} // namespace gui
} // namespace amcset
