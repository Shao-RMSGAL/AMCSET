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


/*! \file target_model.h
 *  \brief SRIM mode's precomputed target (internal): per layer and element
 *  constants, TRIM's layer averages, stopping curves and collision pairs.
 *
 * Built by the TrimSimulation constructor (trim_engine.cpp) and read by the
 * transports and by the GPU backends (compute/kernel_input.cpp), which copy
 * it into the kernel's flat arrays.
 */

#pragma once

#include "trim/trim.h"

#include "crystal/grain_map.h"
#include "crystal/lattice_search.h"
#include "physics/scattering.h"
#include "physics/stopping.h"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <vector>

namespace amcset {
namespace trim {

using physics::CollisionPair;

namespace detail {

struct ElementModel {
  unsigned z = 0;
  mass_quantity mass;
  double cumulative = 1.0;  // cumulative stoichiometric fraction
  double fraction = 1.0;    // stoichiometric fraction
  energy_quantity ed, eb, es;
  std::size_t tally_index = 0;
  std::size_t species = 0;
  // Kinchin-Pease (Robinson) constants for this element as a recoil.
  double kp_kd = 0.0;
  double kp_eps_per_eV = 0.0;
  // arc-dpa / rpa constants (has_efficiency: known for this element).
  // arc_scale = (1 - c_arc)/(2 Ed/0.8)^b_arc, rpa_bc = b_rpa^c_rpa and
  // rpa_scale = rpa_bc/(2 Ed/0.8)^c_rpa + 1: the per-element factors of
  // Nordlund et al. 2018 Eqs. 7 and 9, computed once.
  bool has_efficiency = false;
  double arc_b = 0.0, arc_c = 1.0, rpa_c = 1.0;
  double arc_scale = 0.0, rpa_bc = 0.0, rpa_scale = 1.0;
};

/*!
 * Layer constants of SRIM-2013's TRIM.exe (TRIMSTART setup routine). TRIM
 * averages the layer into one "atom" for the ion's free flight path, the
 * Kinchin-Pease partition and the straggling: its mean atomic number is an
 * integer, accumulated as z = CInt(c_j Z_j + z) over the elements.
 */
struct Trim2013Layer {
  double z_int = 0.0;                //!< CInt-accumulated layer atomic number
  mass_quantity mass;                //!< Stoichiometric mean mass
  length_quantity pmax0;             //!< Monolayer impact parameter N^(-1/3)/sqrt(pi)
  length_quantity screening;         //!< 0.5292*0.8853 A / (Z1^0.23 + Z^0.23), ion on the layer
  common::inverse_energy_quantity fi; //!< Reduced energy per lab energy, ion on the layer
  double epsdg = 0.0;                //!< eps_min-type factor of TRIM's p_max (E_min = 5 eV)
  double kp_eps_per_eV = 0.0;        //!< 0.01 Z^(-7/3) per eV
  double kp_k = 0.0;                 //!< 0.1334 Z^(2/3) M^(-1/2)
  //! Straggling constants (TRIM's G754/G758, ion constants also used for
  //! recoils): Omega = min(k_low E, k_bohr f_rel) in eV/sqrt(Angstrom), with
  //! k_low in 1/sqrt(Angstrom) (times E in eV); includes TRIM's sqrt(12)
  //! for its uniform deviate.
  double k_low = 0.0;
  double k_bohr = 0.0;
};

struct LayerModel {
  length_quantity top, bottom;
  number_density_quantity density;
  length_quantity monolayer;  // N^(-1/3)
  length_quantity min_flight; // shortest free flight path
  bool gas = false;
  std::vector<ElementModel> elements;
  Trim2013Layer t13;
  //! Crystal transport: the layer's lattice (null: amorphous).
  std::shared_ptr<const crystal::LatticeSearch> lattice;
  /*!
   * Crystal transport, polycrystal: the grains and the lattice of each
   * (#lattice then stands for all of them: the same cell, impact parameter
   * and reach), and half the width of the slab at the grain boundaries
   * that is treated as amorphous.
   */
  std::shared_ptr<const crystal::GrainMap> grains;
  std::vector<crystal::LatticeSearch> grain_lattices;
  double boundary_half_A = 0.0;
  //! Number of the layer's lattice among all lattices of the target (a
  //! polycrystal's grains follow): tells their sites apart.
  std::uint32_t domain = 0;

  //! Index in #elements of element \p z (the first, if listed twice).
  std::size_t element_of(unsigned z) const {
    for (std::size_t j = 0; j < elements.size(); ++j) {
      if (elements[j].z == z) return j;
    }
    throw std::logic_error("crystal transport: an atom of the lattice is not among the layer's elements");
  }
};

struct SpeciesInfo {
  unsigned z = 0;
  mass_quantity mass;
};

struct SpeciesLayer {
  physics::StoppingCurve stopping;
  std::vector<CollisionPair> pairs; // per element of the layer
  std::vector<double> eps_min;      // reduced T_min per pair
  // Electronic straggling per path length is tabulated on the same grid
  // as the stopping: Omega^2/L in eV^2/Angstrom (log-spaced, linear interp).
  double strag_log_e0 = 0.0;
  double strag_inv_step = 1.0;
  std::vector<double> strag;
};

} // namespace detail

using detail::ElementModel;
using detail::LayerModel;
using detail::SpeciesInfo;
using detail::SpeciesLayer;

//! A TRIM.DAT particle checked against the target.
struct ResolvedParticle {
  SourceParticle particle;
  std::size_t layer = 0;   // layer at the starting depth
  std::size_t species = 0; // recoils: species of the atom
  std::size_t element_layer = 0; // recoils: layer whose element entry describes the atom
  std::size_t element = 0;       // recoils: element index in that layer
};

class TargetModel {
public:
  std::vector<ResolvedParticle> sources; // TRIM.DAT particles (valid ones)
  std::vector<SpeciesInfo> species; // index 0 = the ion
  std::vector<LayerModel> layers;
  std::vector<std::vector<SpeciesLayer>> data; // [species][layer]
  length_quantity thickness;
  energy_quantity lowest_ed;
  std::shared_ptr<const physics::SrimStoppingDatabase> stopping_db;
  std::shared_ptr<const physics::PairPotentialTable> potentials;

  std::size_t layer_index(length_quantity x) const {
    for (std::size_t i = 0; i < layers.size(); ++i) {
      if (x < layers[i].bottom) return i;
    }
    return layers.size() - 1;
  }
};

} // namespace trim
} // namespace amcset
