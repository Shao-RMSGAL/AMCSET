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

#include "model/output_location.h"

#include <QDir>
#include <QFileInfo>
#include <QSettings>
#include <QStandardPaths>

#include <filesystem>

namespace amcset {
namespace gui {

namespace {
const QString key = QStringLiteral("paths/output_root");

std::string to_std(const QString &p) {
  return std::filesystem::path(QDir::toNativeSeparators(p).toStdWString()).string();
}
} // namespace

QString default_output_root() {
  // QStandardPaths::DocumentsLocation: Windows' known folder Documents, the
  // macOS Documents folder, XDG_DOCUMENTS_DIR on Linux; the home folder if
  // there is none.
  QString docs = QStandardPaths::writableLocation(QStandardPaths::DocumentsLocation);
  if (docs.isEmpty()) docs = QDir::homePath();
  return QDir(docs).filePath(QStringLiteral("AMCSET"));
}

QString output_root() {
  const QString chosen = QSettings().value(key).toString();
  if (!chosen.isEmpty() && QFileInfo(chosen).isAbsolute()) return QDir::cleanPath(chosen);
  return default_output_root();
}

void set_output_root(const QString &dir) {
  QSettings s;
  if (dir.isEmpty()) s.remove(key);
  else s.setValue(key, QDir::cleanPath(QFileInfo(dir).absoluteFilePath()));
}

std::string new_scenario_output_directory() { return to_std(QDir(output_root()).filePath(QStringLiteral("output"))); }

std::string absolute_output_directory(const std::string &dir, const QString &scenario_path) {
  if (dir.empty()) return dir; // (no files are written)
  const QString d = QString::fromStdWString(std::filesystem::path(dir).wstring());
  if (QFileInfo(d).isAbsolute()) return dir;
  const QString base = scenario_path.isEmpty() ? output_root() : QFileInfo(scenario_path).absolutePath();
  return to_std(QDir::cleanPath(QDir(base).filePath(d)));
}

} // namespace gui
} // namespace amcset
