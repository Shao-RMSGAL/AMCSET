#include <gtest/gtest.h>
#include <stdexcept>

#include "application.h"
#include "mainwindow.h"

TEST(MainWindowSimulationErrors, FormatsNestedExceptionsForUserDisplay) {
  try {
    try {
      throw std::invalid_argument("bad inputs");
    } catch (const std::exception &) {
      std::throw_with_nested(std::runtime_error("simulation failed"));
    }
  } catch (const std::exception &e) {
    const auto formatted = amcset::gui::MainWindow::formatSimulationError(e);
    EXPECT_TRUE(formatted.contains("simulation failed"));
    EXPECT_TRUE(formatted.contains("bad inputs"));
  }
}

TEST(ApplicationStatePersistence, SavesAndLoadsMainWindowState) {
  amcset::gui::Application app;

  QVariantMap state;
  state["incidentElectron"] = true;
  state["damageCascade"] = false;
  state["incidentEnergy"] = 42.5;
  state["incidentMass"] = 27;
  state["incidentCharge"] = 13;
  state["incidentAngle"] = 22.5;
  state["bombardmentCount"] = 9;
  state["threadCount"] = 8;
  state["range"] = 1500.0;

  QVariantList layers;
  QVariantMap layer;
  layer["startingDepth"] = 120.0;
  layer["depth"] = 80.0;
  layer["density"] = 850.0;
  layer["phase"] = "Solid";
  layer["name"] = "Layer 1";
  QVariantList composition;
  QVariantMap element;
  element["element"] = "H";
  element["atomicProportion"] = 0.75;
  element["atomicPercent"] = 75.0;
  element["displacementEnergy"] = 40.0;
  composition.append(element);
  layer["composition"] = composition;
  layers.append(layer);
  state["layers"] = layers;

  app.saveMainWindowState(state);

  const auto restored = app.loadMainWindowState();
  EXPECT_TRUE(restored["incidentElectron"].toBool());
  EXPECT_FALSE(restored["damageCascade"].toBool());
  EXPECT_DOUBLE_EQ(restored["incidentEnergy"].toDouble(), 42.5);
  EXPECT_EQ(restored["incidentMass"].toInt(), 27);
  EXPECT_EQ(restored["incidentCharge"].toInt(), 13);
  EXPECT_DOUBLE_EQ(restored["incidentAngle"].toDouble(), 22.5);
  EXPECT_EQ(restored["bombardmentCount"].toInt(), 9);
  EXPECT_EQ(restored["threadCount"].toInt(), 8);
  EXPECT_DOUBLE_EQ(restored["range"].toDouble(), 1500.0);
  EXPECT_EQ(restored["layers"].toList().size(), 1);
}
