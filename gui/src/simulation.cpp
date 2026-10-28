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

#include "simulation.h"

#include <QByteArray>
#include <QMessageBox>

namespace amcset {
namespace gui {

void SimulationThread::run() {
  try {
    const auto result = simulation_->run_simulation();
    emit resultReady(result);
  } catch (const std::exception &e) {
    QString error = QString::fromStdString(e.what());
    try {
      std::rethrow_if_nested(e);
    } catch (const std::exception &nested) {
      error = QString::fromStdString(std::string("Simulation failed:\n") + nested.what());
    } catch (...) {
    }
    emit errorReady(error);
  } catch (...) {
    emit errorReady(QStringLiteral("Simulation failed with an unknown error."));
  }
}

} // namespace gui
} // namespace amcset
