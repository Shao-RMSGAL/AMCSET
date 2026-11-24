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

/*! \file guard.h
 *  \brief Exceptions must not reach Qt: guards for every Qt entry point.
 *
 * Qt does not support exceptions passing through its frames (event
 * dispatch, signal emission, model views, QRhi callbacks). Every place where
 * Qt calls AMCSET code is therefore guarded:
 *
 * - slots and lambdas given to connect() are wrapped with guarded() (or
 *   guarded_quiet() where a dialog would be wrong);
 * - overrides Qt calls directly (paint and mouse events, model data(),
 *   QRhiWidget render()) catch inside, with report_exception_quietly();
 * - Application::notify() catches whatever still escapes an event handler,
 *   as the last line of defence.
 *
 * A caught exception is written to the diagnostic log and shown to the user
 * in a message box, posted to the event loop (never opened inside a paint
 * or mouse handler), at most one at a time.
 */

#pragma once

#include <QApplication>
#include <QPointer>
#include <QString>
#include <QWidget>

#include <cstddef>
#include <exception>
#include <tuple>
#include <type_traits>
#include <utility>

namespace amcset {
namespace gui {

/*!
 * \brief Log the exception \p e caught at \p where and show it to the user
 * (a warning box over \p parent, posted to the event loop). Never throws.
 */
void report_exception(QWidget *parent, const QString &where, std::exception_ptr e) noexcept;

//! Log only (for places that run very often, e.g. paint events; the first
//! occurrence per place is also shown).
void report_exception_quietly(const QString &where, std::exception_ptr e) noexcept;

//! The message of an exception (for error strings).
QString exception_text(std::exception_ptr e) noexcept;

namespace detail {

template <class F, class Tuple, std::size_t... I>
constexpr bool invocable_with(std::index_sequence<I...>) {
  return std::is_invocable_v<F &, std::tuple_element_t<I, Tuple>...>;
}

/*!
 * Call \p f with the longest prefix of the signal's arguments it accepts
 * (Qt's own rule for slots that take fewer arguments than the signal, e.g.
 * a [this] { ... } lambda on QAction::triggered(bool)).
 */
template <std::size_t N, class F, class Tuple>
void call_prefix(F &f, Tuple &args) {
  if constexpr (invocable_with<F, Tuple>(std::make_index_sequence<N>{})) {
    [&]<std::size_t... I>(std::index_sequence<I...>) {
      f(std::forward<std::tuple_element_t<I, Tuple>>(std::get<I>(args))...);
    }(std::make_index_sequence<N>{});
  } else if constexpr (N > 0) {
    call_prefix<N - 1>(f, args);
  } else {
    static_assert(N > 0, "the slot cannot be called with the signal's arguments");
  }
}

} // namespace detail

/*!
 * \brief Wrap a slot so that no exception leaves it: connect(a, &A::sig,
 * this, guarded(this, tr("Open file"), [this] { ... })). The slot receives
 * as many of the signal's arguments as it takes.
 */
template <class F>
auto guarded(QWidget *parent, QString where, F f) {
  return [parent = QPointer<QWidget>(parent), where = std::move(where), f = std::move(f)](auto &&...args) mutable {
    try {
      auto t = std::forward_as_tuple(std::forward<decltype(args)>(args)...);
      detail::call_prefix<sizeof...(args)>(f, t);
    } catch (...) {
      report_exception(parent.data(), where, std::current_exception());
    }
  };
}

//! guarded() that only logs (slots that run on timers or very often).
template <class F>
auto guarded_quiet(QString where, F f) {
  return [where = std::move(where), f = std::move(f)](auto &&...args) mutable {
    try {
      auto t = std::forward_as_tuple(std::forward<decltype(args)>(args)...);
      detail::call_prefix<sizeof...(args)>(f, t);
    } catch (...) {
      report_exception_quietly(where, std::current_exception());
    }
  };
}

/*!
 * \brief The application object: QApplication whose notify() catches any
 * exception an event handler lets through (logs it, shows it, drops the
 * event) instead of letting it unwind through Qt.
 */
class Application : public QApplication {
public:
  Application(int &argc, char **argv) : QApplication(argc, argv) {}
  bool notify(QObject *receiver, QEvent *event) override;
};

} // namespace gui
} // namespace amcset
