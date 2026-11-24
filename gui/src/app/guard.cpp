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

#include "diagnostics.h"

#include <QEvent>
#include <QMessageBox>
#include <QMetaEnum>
#include <QSet>
#include <QTimer>

#include <new>
#include <stdexcept>

namespace amcset {
namespace gui {

QString exception_text(std::exception_ptr e) noexcept {
  try {
    try {
      if (e) std::rethrow_exception(e);
    } catch (const std::bad_alloc &) {
      return QStringLiteral("out of memory");
    } catch (const std::exception &x) {
      return QString::fromUtf8(x.what());
    } catch (...) {
      return QStringLiteral("unknown error (not a standard exception)");
    }
  } catch (...) {
    // Building the text failed (out of memory): nothing more to say.
  }
  return QStringLiteral("unknown error");
}

namespace {

bool g_box_open = false; // GUI thread only

void log_exception(const QString &where, const QString &text) noexcept {
  try {
    AMCSET_LOG(error) << "GUI: " << where.toStdString() << ": " << text.toStdString();
  } catch (...) {
  }
}

} // namespace

void report_exception(QWidget *parent, const QString &where, std::exception_ptr e) noexcept {
  try {
    const QString text = exception_text(e);
    log_exception(where, text);
    // Posted: the handler that caught it may be in the middle of painting
    // or of a mouse event, where a modal dialog must not open. One box at a
    // time, so that a repeating failure cannot flood the screen.
    QPointer<QWidget> p(parent);
    QTimer::singleShot(0, qApp, [p, where, text] {
      if (g_box_open) return;
      g_box_open = true;
      QMessageBox::warning(p ? p.data() : QApplication::activeWindow(), QObject::tr("AMCSET"),
                           QObject::tr("%1 failed:\n%2\n\nThe program continues; details are in the "
                                       "diagnostic log (Help > Open log folder).")
                               .arg(where, text));
      g_box_open = false;
    });
  } catch (...) {
    // Reporting must never throw out of a handler.
  }
}

void report_exception_quietly(const QString &where, std::exception_ptr e) noexcept {
  try {
    static QSet<QString> shown;
    const QString text = exception_text(e);
    log_exception(where, text);
    if (!shown.contains(where)) {
      shown.insert(where);
      report_exception(nullptr, where, e);
    }
  } catch (...) {
  }
}

bool Application::notify(QObject *receiver, QEvent *event) {
  try {
    return QApplication::notify(receiver, event);
  } catch (...) {
    QString where = QStringLiteral("An event handler");
    try {
      const char *type = QMetaEnum::fromType<QEvent::Type>().valueToKey(event->type());
      where = QStringLiteral("Handling a %1 event for %2")
                  .arg(QString::fromLatin1(type ? type : "Qt"),
                       QString::fromLatin1(receiver ? receiver->metaObject()->className() : "an object"));
    } catch (...) {
    }
    report_exception_quietly(where, std::current_exception());
    return false;
  }
}

} // namespace gui
} // namespace amcset
