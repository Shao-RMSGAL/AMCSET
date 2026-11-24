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

/*! \file trim_tables.h
 *  \brief The results of a SRIM-mode run as output tables.
 *
 * Every quantity of TRIM's output files (and the data behind TRIM's plots)
 * as an output::Dataset, in TRIM's units, so it can be written as CSV, JSON
 * or Parquet with any units and number format, and plotted by the GUI:
 *
 * | Table key           | TRIM file / plot                                   |
 * |---------------------|----------------------------------------------------|
 * | range               | RANGE.txt, ion and recoil distribution plots       |
 * | lateral             | LATERAL.txt, lateral range plot                    |
 * | vacancies           | VACANCY.txt                                        |
 * | replacements        | NOVAC.txt                                          |
 * | damage              | collision events plot (displacements, vacancies,   |
 * |                     | replacements, interstitials)                       |
 * | ionization          | IONIZ.txt, ionization plot                         |
 * | phonons             | PHONON.txt, phonon plot                            |
 * | energy_to_recoils   | E2RECOIL.txt, energy to recoils / absorbed plots   |
 * | sputtering          | sputtering yields (TRIM calculation window)        |
 * | surface_energy      | integral and differential sputtering plots         |
 * | recoil_spectrum     | (AMCSET) energy spectrum of displaced atoms        |
 * | exits               | BACKSCAT.txt, TRANSMIT.txt, SPUTTER.txt, TRIMOUT.txt|
 * | ion_positions       | RANGE_3D.txt                                       |
 * | exyz                | EXYZ.txt                                           |
 * | collisions          | COLLISON.txt (cascade summaries)                   |
 */

#pragma once

#include "output/table.h"
#include "trim/trim.h"

namespace amcset {
namespace trim {

//! The finished run's results (TrimSimulation::tallies()).
output::Dataset make_dataset(const TrimSimulation &sim);
//! The results in \p tallies (e.g. TrimSimulation::snapshot() during a run).
output::Dataset make_dataset(const TrimSimulation &sim, const TrimTallies &tallies);

} // namespace trim
} // namespace amcset
