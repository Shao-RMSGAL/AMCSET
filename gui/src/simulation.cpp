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
