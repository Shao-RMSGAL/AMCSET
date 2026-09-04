#pragma once

#include "amcset_common.h"
#include <QObject>
#include <QThread>

namespace amcset {
namespace gui {

class SimulationThread : public QThread {
  Q_OBJECT

public:
  SimulationThread(std::unique_ptr<amcset::common::Simulation> simulation)
      : simulation_(std::move(simulation)) {};

public slots:
  void run() override;

signals:
  void resultReady(int result);

private:
  std::unique_ptr<amcset::common::Simulation> simulation_;
};
} // namespace gui
} // namespace amcset
