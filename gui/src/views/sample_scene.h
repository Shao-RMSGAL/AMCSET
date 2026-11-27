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

/*! \file sample_scene.h
 *  \brief The sample under the beam spot as the path views draw it: the
 *  spot's outline, and the target's layers and grains on the surface, on two
 *  cross-sections and on the faces of a block.
 *
 * Built from the scenario (trim::TrimConfig) with the grain maps the crystal
 * transport uses (crystal::make_grain_map), so the grains drawn are the
 * grains the ions meet. Coordinates are TRIM's (x = depth, y and z lateral,
 * Angstrom) at their true places on the sample: the transport's positions
 * plus each ion's entry point on the spot (CollisionPoint::entry_x_A).
 */

#pragma once

#include "crystal/crystal.h"
#include "trim/trim.h"

#include <QColor>
#include <QString>

#include <array>
#include <atomic>
#include <cstdint>
#include <memory>
#include <vector>

namespace amcset {
namespace gui {

//! A vertex of the sample's drawing (TRIM axes, Angstrom; RGBA bytes, straight alpha).
struct SceneVertex {
  float x, y, z;
  std::uint32_t rgba;
};

//! What the views draw of the sample (index: TrackView::Projection).
struct SampleScene {
  //! Filled cells (triangles) and lines (pairs of vertices) per view:
  //! 0 depth-Y (the section through z = 0), 1 depth-Z (through y = 0),
  //! 2 beam view (the surface), 3 3D (the surface and two sides of a block).
  std::array<std::vector<SceneVertex>, 4> triangles;
  std::array<std::vector<SceneVertex>, 4> lines;
  //! The spot's outline on the surface (beam and 3D views), drawn over the paths.
  std::vector<SceneVertex> spot_lines;
  //! The colour of the spot's outline and of the entry points.
  static QColor spot_color() { return QColor(255, 196, 40); }
  //! The region drawn: depth 0 .. depth_A, y in [y0, y1], z in [z0, z1].
  double depth_A = 0.0, y0 = 0.0, y1 = 0.0, z0 = 0.0, z1 = 0.0;
  trim::BeamSpot spot;
  bool spot_drawn = false; //!< the spot's outline is in the lines (not for a point)
  //! Grains coloured by the beam direction in their crystal (inverse pole
  //! figure of cubic crystals: [001] red, [101] green, [111] blue); false:
  //! other crystals, coloured by grain.
  bool ipf = false;
  bool has_grains = false;
  QString message; //!< why nothing could be drawn (empty: fine)

  //! The layers, for the readout under the mouse.
  struct Layer {
    double top_A = 0.0, bottom_A = 0.0;
    QString name;
    bool crystal = false;
    crystal::Lattice lattice;
    crystal::Mat3 rotation;                      //!< single crystal: crystal -> TRIM axes
    std::shared_ptr<const crystal::GrainMap> grains; //!< polycrystal (null otherwise)
  };
  std::vector<Layer> layers;
  crystal::Vec3 beam{1.0, 0.0, 0.0}; //!< the incident direction (TRIM axes)

  //! What is at a point of the sample (TRIM axes, true place): layer, grain,
  //! the beam direction in its crystal and the angle to the nearest axis.
  QString describe(double x, double y, double z) const;
};

/*!
 * \brief The scene of the scenario \p c (its beam spot and target), with
 * lines in the colour \p line_color. Takes up to a second for large grain
 * maps: call it off the GUI thread. \p cancel (may be null) stops it early.
 * Never throws: a failure is in SampleScene::message.
 */
std::shared_ptr<const SampleScene> build_sample_scene(const trim::TrimConfig &c, QColor line_color,
                                                      const std::atomic<bool> *cancel = nullptr);

//! The colour of the beam direction \p d in a cubic crystal's frame (inverse pole figure).
QColor ipf_color(const crystal::Vec3 &d);

} // namespace gui
} // namespace amcset
