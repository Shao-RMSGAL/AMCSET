#include "simulation.h"

namespace amcset {
namespace gui {

void SimulationThread::run() {
  auto result = simulation_->run_simulation();
  emit resultReady(result);
}

} // namespace gui
} // namespace amcset
