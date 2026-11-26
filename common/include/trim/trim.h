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

/*! \file trim.h
 *  \brief AMCSET's SRIM mode: a TRIM-compatible Monte Carlo transport engine.
 *
 * TrimSimulation reproduces the algorithms of TRIM (the transport program of
 * SRIM-2013) as described in Ziegler, Biersack & Ziegler, "SRIM - The
 * Stopping and Range of Ions in Matter" (2008/2015), chapters 7-9:
 *
 *  - amorphous multilayer targets (any number of layers and elements), with
 *    per-element displacement (Ed), lattice-binding (Eb) and surface-binding
 *    (Es) energies in every layer;
 *  - the free-flight-path method with TRIM's maximum impact parameter
 *    (energy-transfer and 5 % electronic-loss criteria), monolayer steps at
 *    low energy and near interfaces;
 *  - the magic formula (eps <= 10) / screened Rutherford (eps > 10)
 *    scattering with the ZBL universal potential, or SRIM's pair-specific
 *    potentials when a SRIM Data directory is given;
 *  - electronic stopping from SRIM-2013 tables (exact TRIM stopping) or the
 *    built-in models, Bragg's rule with a per-layer compound correction,
 *    gas/solid phase, TRIM's electronic straggling, and Bethe-Bloch
 *    continuation for high-energy particles;
 *  - damage in TRIM's three modes: quick Kinchin-Pease (Robinson's
 *    modification, as TRIM), full cascades (vacancies, replacement
 *    collisions, interstitials) and monolayer collision steps;
 *  - sputtering through a planar surface potential, backscattering and
 *    transmission;
 *  - separate random-number streams for ions and recoils (so the ion
 *    distributions of the quick and full-cascade modes are identical), a
 *    user seed, and results independent of the thread count;
 *  - every quantity TRIM tabulates, written in TRIM's file formats (RANGE.txt,
 *    VACANCY.txt, NOVAC.txt, IONIZ.txt, PHONON.txt, E2RECOIL.txt, LATERAL.txt,
 *    BACKSCAT.txt, TRANSMIT.txt, SPUTTER.txt, RANGE_3D.txt, TDATA.txt) plus a
 *    JSON summary.
 *
 * Coordinates follow TRIM: x is the depth below the front surface, y and z
 * are lateral; a tilted beam is tilted in the x-y plane. (The AMCSET-wide
 * SimulationEngine trajectories are converted to the +z-depth convention.)
 */

#pragma once

#include "amcset_units.h"
#include "crystal/crystal.h"
#include "output/table.h"
#include "physics/scattering.h"
#include "trim/trim_dat.h"
#include "physics/stopping.h"
#include "simulation_engine.h"
#include "tallies.h"

#include <array>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace amcset {
namespace common {
struct CollisionPoint;
}
namespace diag {
class Activity;
}

namespace trim {

using common::angle_quantity;
using common::energy_quantity;
using common::length_quantity;
using common::mass_density_quantity;
using common::mass_quantity;
using common::number_density_quantity;

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------

//! TRIM's "type of damage calculation".
enum class DamageCalculation {
  kinchin_pease = 1, //!< "Ion distribution and quick calculation of damage"
  full_cascade = 2,  //!< "Detailed calculation with full damage cascades"
  monolayer = 3      //!< "Monolayer collision steps / surface sputtering"
};

std::string to_string(DamageCalculation d);
//! Accepts "kp", "quick", "kinchin-pease", "full", "full-cascade",
//! "monolayer", "sputtering" or TRIM's numbers 1..3. \throws std::invalid_argument
DamageCalculation parse_damage_calculation(const std::string &s);

/*!
 * \brief Robinson's damage energy (the part of a recoil's energy that ends
 * up in nuclear motion) for a recoil of atomic number \p z and mass \p mass
 * slowing down in its own element:
 * E_d = T / (1 + k_d g(eps_d)), g = eps + 0.40244 eps^(3/4) + 3.4008 eps^(1/6),
 * eps_d = 0.01014 Z^(-7/3) T[eV], k_d = 0.1334 Z^(2/3) / sqrt(M[amu])
 * (Norgett, Robinson & Torrens, Nucl. Eng. Des. 33, 50 (1975)).
 */
energy_quantity robinson_damage_energy(unsigned z, mass_quantity mass,
                                       energy_quantity recoil_energy);

/*!
 * \brief Vacancies TRIM's quick (Kinchin-Pease) calculation assigns to a
 * primary recoil of energy \p recoil_energy: 0 if the damage energy is below
 * \p displacement_energy, 1 if it is below 2 Ed / 0.8, otherwise
 * 0.8 E_d / (2 Ed) (the NRT formula).
 */
double kinchin_pease_vacancies(unsigned z, mass_quantity mass,
                               energy_quantity displacement_energy,
                               energy_quantity recoil_energy);

/*!
 * \brief Material constants of the arc-dpa and rpa damage functions
 * (Nordlund et al., Nat. Commun. 9, 1084 (2018), Eqs. 7 and 9): the
 * defect production efficiency xi_arc(T_d) = (1 - c)/(2 Ed/0.8)^b T_d^b + c
 * and the mixing efficiency xi_rpa(T_d) = (b^c/(2 Ed/0.8)^c + 1)
 * T_d^c/(b^c + T_d^c), both multiplying NRT's 0.8 T_d/(2 Ed).
 */
struct DamageEfficiency {
  double arc_b = 0.0;               //!< b_arcdpa (< 0)
  double arc_c = 1.0;               //!< c_arcdpa (0 < c <= 1)
  energy_quantity rpa_b;            //!< b_rpa (an energy)
  double rpa_c = 1.0;               //!< c_rpa (> 0)
  //! The Ed the constants were fitted with (0: not known).
  energy_quantity fitted_ed;
};

/*!
 * \brief The constants Nordlund et al. (2018), Table 1, fitted to MD for Fe,
 * Cu, Ni, Pd, Pt and W (with Ed = 40, 33, 39, 41, 42 and 70 eV); nothing for
 * other elements.
 */
std::optional<DamageEfficiency> published_damage_efficiency(unsigned z);

//! One element of a target layer.
struct TargetElement {
  unsigned z = 0;                  //!< Atomic number
  mass_quantity mass;              //!< Atom mass (default: standard atomic weight)
  double stoichiometry = 1.0;      //!< Relative number of atoms (normalized per layer)
  energy_quantity displacement_energy;    //!< Ed
  energy_quantity lattice_binding_energy; //!< Eb
  energy_quantity surface_binding_energy; //!< Es
  std::string name;                //!< Display name (default: chemical symbol)
  //! Position (1-based) in a TRIM.IN's list of target atoms, or 0. TRIM
  //! numbers the atoms of the whole target; the number orders the output
  //! columns and decides TRIM's replacement of zero energies (see
  //! TrimSimulation). 0: numbered in order of first appearance.
  unsigned trim_atom = 0;
  //! arc-dpa/rpa constants for this element (not TRIM); empty: the
  //! published ones (published_damage_efficiency), if any.
  std::optional<DamageEfficiency> damage_efficiency;

  //! Element \p z with SRIM-like defaults (natural mass, Ed = 25 eV,
  //! Eb = 3 eV, Es = heat of sublimation).
  static TargetElement defaults(unsigned z, double stoichiometry = 1.0);
};

//! One target layer.
struct TargetLayer {
  std::string name;
  length_quantity width;           //!< Layer thickness
  mass_density_quantity density;   //!< Mass density
  bool gas = false;                //!< Gas phase (uses gas-phase stopping where available)
  double compound_correction = 1.0; //!< Multiplies Bragg's-rule stopping
  std::vector<TargetElement> elements;
  /*!
   * \brief The layer's crystal (not TRIM; empty: amorphous, as every TRIM
   * target). With a crystal the elements' stoichiometry and the density are
   * those of the unit cell (apply_crystal). The trim2013 and refined
   * transports ignore it: they treat every layer as amorphous.
   */
  std::optional<crystal::Crystal> crystal;
};

//! Limits of a layer's crystal (TrimConfig::validate).
inline constexpr double max_crystal_temperature_K = 10000.0;
inline constexpr double max_debye_temperature_K = 5000.0;
inline constexpr double max_vibration_rms_A = 1.0;

/*!
 * \brief Make \p layer's elements, stoichiometry and density those of its
 * crystal: elements of the unit cell are added (with TargetElement::defaults)
 * or keep their energies and masses, elements that are not in the cell are
 * removed, the stoichiometry is the number of atoms per cell and the density
 * the cell's. Does nothing without a crystal.
 * \throws std::invalid_argument if the cell has no atoms.
 */
void apply_crystal(TargetLayer &layer);
/*!
 * \brief Check a layer's crystal: a valid cell, the layer's composition and
 * density equal to the cell's (as apply_crystal leaves them), an orientation
 * that exists, the temperature and the thermal vibrations in range.
 * \throws std::invalid_argument starting with \p where.
 */
void check_layer_crystal(const TargetLayer &layer, const std::string &where);

struct TrimConfig;

/*!
 * \brief How a layer's crystal sits in the target and meets the beam, as the
 * outputs report it (summary.json "crystals", the "crystals" table, the
 * metadata of the trajectory and MD event files). Not TRIM.
 *
 * Two frames: TRIM's (x = depth, y, z lateral; the beam tilts in x-y) and
 * the output axes of the trajectory and event files (x = TRIM's y, y =
 * TRIM's z, z = depth).
 */
struct CrystalReport {
  std::size_t layer = 0; //!< 0-based
  std::string layer_name;
  std::string cell_name;
  bool followed = false;        //!< the crystal transport follows the lattice (else: amorphous)
  double top_A = 0.0;           //!< depth of the layer's top: the lattice origin's depth
  crystal::Mat3 lattice_A;      //!< rows: the lattice vectors a, b, c in the crystal frame
  crystal::Mat3 to_trim;        //!< crystal -> TRIM axes (v_trim = R v_crystal)
  crystal::Mat3 to_output;      //!< crystal -> output axes
  std::array<double, 3> euler_deg{0, 0, 0}; //!< Bunge (phi1, Phi, phi2) relative to the sample (z = outward normal)
  //! Crystal directions along the depth (inward normal), TRIM's y (the tilt
  //! plane) and TRIM's z, and of the beam: the nearest low-index [uvw] and
  //! the angle to it (degrees).
  std::array<int, 3> depth_uvw{0, 0, 1}, y_uvw{1, 0, 0}, z_uvw{0, 1, 0}, beam_uvw{0, 0, 1};
  double depth_off_deg = 0.0, y_off_deg = 0.0, z_off_deg = 0.0, beam_off_deg = 0.0;
  crystal::Vec3 beam_crystal;   //!< the beam's unit vector in the crystal frame
  double tilt_deg = 0.0;        //!< ion.angle
  double twist_deg = 0.0;
  double temperature_K = 0.0;
  //! rms thermal displacement along one axis per element as the transport
  //! uses it (0: static), and whether it came from Gemmell's table.
  struct Vibration {
    unsigned z = 0;
    double rms_A = 0.0;
    bool from_table = false;
  };
  std::vector<Vibration> vibrations;
  /*!
   * The low-index crystal axes nearest the beam (indices up to 2, nearest
   * first): the angle to the beam, the mean spacing d of the atoms along the
   * row, and Lindhard's critical angles for the ion at its energy E (Z2: the
   * cell's mean atomic number; Gemmell, Rev. Mod. Phys. 46, 129 (1974),
   * Eqs. 2.27-2.29): psi_1 = (2 Z1 Z2 e^2 / (E d))^1/2, valid while
   * psi_1 <= a/d, and below E' = 2 Z1 Z2 e^2 d / a^2 (psi_1 > a/d)
   * psi_2 = (C a psi_1 / (d sqrt 2))^1/2 with C = sqrt 3 and Lindhard's
   * screening length a = 0.8853 a0 (Z1^2/3 + Z2^2/3)^-1/2. An ion within
   * about the critical angle (psi_1 or psi_2, whichever applies) of an axis
   * channels along it.
   */
  struct Axis {
    std::array<int, 3> uvw{0, 0, 1};
    double angle_deg = 0.0;
    double spacing_A = 0.0;
    double psi1_deg = 0.0;
    double psi2_deg = 0.0;
    double critical_deg = 0.0; //!< psi_1 if psi_1 <= a/d, else psi_2
  };
  std::vector<Axis> axes;
};

//! The report of every layer that has a crystal (empty: none).
std::vector<CrystalReport> crystal_reports(const TrimConfig &config);
//! "[1 1 0]".
std::string uvw_text(const std::array<int, 3> &uvw);

/*!
 * \brief Where the ions of the beam meet the surface (not TRIM: TRIM's beam
 * is a point). Each ion enters at a random point of the spot. Only the
 * crystal transport uses it: an amorphous target looks the same everywhere.
 * Positions in the results stay relative to each ion's entry point.
 */
struct BeamSpot {
  enum class Shape {
    point,     //!< every ion at the same point of the surface (the lattice origin)
    rectangle, //!< uniform over width (y) x height (z)
    circle,    //!< uniform over a disc of diameter width
    gaussian   //!< normal, with full widths at half maximum width (y) and height (z)
  };
  Shape shape = Shape::rectangle;
  length_quantity width = 1.0 * common::micrometer;
  length_quantity height = 1.0 * common::micrometer;

  bool operator==(const BeamSpot &) const = default;
};

std::string to_string(BeamSpot::Shape s);
//! "point", "rectangle" ("uniform"), "circle", "gaussian". \throws std::invalid_argument
BeamSpot::Shape parse_beam_spot_shape(std::string_view s);

//! The incident beam.
struct IonBeam {
  unsigned z = 0;
  mass_quantity mass;     //!< Default: most-abundant-isotope mass
  energy_quantity energy; //!< Kinetic energy
  angle_quantity angle;   //!< Angle to the target normal (tilt in x-y plane)
  BeamSpot spot;          //!< Entry points on the surface (crystal transport only)

  static IonBeam defaults(unsigned z, energy_quantity energy,
                          angle_quantity angle = 0.0 * common::degree);
};

/*!
 * \brief Bugs of SRIM-2013's TRIM.exe that the trim2013 transport
 * reproduces (docs/TRIM2013_ALGORITHM.md, "Bugs of TRIM.exe").
 *
 * Each flag is true by default: AMCSET behaves as TRIM.exe does, so its
 * results agree with TRIM's within statistics. false selects the corrected
 * behaviour; results then depart from TRIM's (by the amounts given in the
 * documentation) and should not be compared with TRIM runs. The refined
 * transport has none of these bugs.
 */
struct TrimBugs {
  //! A primary recoil with T <= Ed is not displaced, but its cascade runs
  //! with the energy its slot had at the end of the previous cascade
  //! (STARTRECOIL in TRIMMain.Main). Fixed: no cascade; T stays in the ion's phonons.
  bool stale_primary_energy = true;
  //! A recoil's electronic loss over a step and the energy it gives in the
  //! collision can together exceed its energy; TRIM books both in full and
  //! sets the mover to max(0, E - T - DEE) (0x436de5), creating energy.
  //! Fixed: the struck atom gets at most E - DEE.
  bool recoil_energy_created = true;
  //! The struck atom of a recoil collision is turned by DIRCOS with sin x = 0
  //! (0x4377d5): its direction depends only on the mover's x cosine.
  //! Fixed: turned from the mover's direction like the primary of an ion
  //! collision. Changes sputtering by up to 20 %.
  bool struck_atom_direction = true;
  //! An atom that makes a replacement collision is also counted as an
  //! interstitial (RECSTOPDATA, 0x457970). Fixed: it fills the vacancy and
  //! is not an interstitial.
  bool replacement_interstitial = true;
  //! A sputtered atom that was set in motion below Ed counts as a
  //! displacement and a vacancy in the totals, not in the depth tables
  //! (0x441ff0). Fixed: it is not counted (it was never displaced).
  bool subthreshold_sputter_counted = true;
  //! An ion collision with T below the stop energy is booked as ion phonons
  //! and again as recoil phonons (0x44d4e0 and 0x4575e0). Fixed: once.
  bool ion_phonons_counted_twice = true;
  //! Cascades replay the ion's next random numbers (the generator state is
  //! saved before and restored after each cascade, 0x434e89/0x438140), which
  //! correlates consecutive cascades. Fixed: each cascade has its own stream.
  bool cascade_random_replay = true;
  //! Zero Ed, Es, Eb are replaced by 20, 4 and 2 eV for target atoms 2..n
  //! only; the loop of 0x44efc0 skips the first atom. Fixed: for every atom.
  bool zero_energy_defaults_skip_first_atom = true;
  //! The layer's atomic number in the ion's screening length, straggling and
  //! the quick calculation is rounded at every element (z = CInt(c Z + z),
  //! 0x44f8c0). Fixed: the stoichiometric mean sum c_j Z_j.
  bool integer_layer_z = true;
  //! EXYZ.txt: after the record at the start the next one is set a step too
  //! low (Int(E/I) - 1, 0x434745), so the first crossing is not recorded.
  //! Fixed: every crossing is recorded.
  bool exyz_skips_first_step = true;
  //! TRIM.DAT recoils (code 6) count each event twice in RANGE.txt's ion
  //! column (at the atom's start and where it stops). Fixed: once, at the start.
  bool recoil_event_counted_twice = true;

  //! Set every flag (true = reproduce TRIM, false = fix everything).
  void set_all(bool reproduce);
  //! True when every flag reproduces TRIM.
  bool all_reproduced() const;
  //! (name, flag) for every bug, in the order above; names are the
  //! configuration keys.
  static std::vector<std::pair<const char *, bool TrimBugs::*>> members();
};

//! Everything a SRIM-mode calculation needs.
struct TrimConfig {
  std::string title = "AMCSET SRIM-mode calculation";
  IonBeam ion;
  std::vector<TargetLayer> layers;
  DamageCalculation damage = DamageCalculation::kinchin_pease;
  /*!
   * \brief Where the particles come from.
   *
   * - beam: every ion is the IonBeam (TRIM's cascade codes 1-3);
   * - particle_ions: the ions of #particles, each with its own energy,
   *   starting point and direction (TRIM.DAT, codes 4 and 5); the ion species
   *   is always IonBeam::z, and IonBeam::energy is the highest energy allowed;
   * - particle_recoils: the target atoms of #particles are set in motion
   *   without an incident ion (damage from neutrons, electrons or photons;
   *   TRIM.DAT). Each atom must be an element of the layer at its starting
   *   depth. With full_cascade (TRIM's code 6) each atom is the primary
   *   recoil of a cascade; with monolayer (code 7) TRIM runs it through the
   *   ion loop, as an ion (of IonBeam's species) with monolayer steps;
   *   kinchin_pease (not in TRIM) books each atom's damage by the
   *   modified Kinchin-Pease formula.
   *
   * TRIM.DAT ions run with kinchin_pease or full_cascade. Every particle is
   * run (ion_count, if not zero, limits the number).
   */
  enum class Source { beam, particle_ions, particle_recoils };
  Source source = Source::beam;
  std::vector<SourceParticle> particles; //!< TRIM.DAT particles
  std::string particles_title;          //!< TRIM.DAT title (line 8)
  //! Number of ions (beam source); for particle sources the number of
  //! particles to run (0 = all).
  std::size_t ion_count = 1000;
  //! First ion number (0-based) of this run: a calculation continued from
  //! a previous one with first_ion = its ion count gives the same results as
  //! one longer run, because every ion has its own random-number stream.
  std::size_t first_ion = 0;
  std::uint32_t seed = 716381; //!< TRIM's default seed
  //! A new seed for every calculation (scenario "seed": "random", CLI
  //! --seed random): TrimSimulation replaces #seed by draw_random_seed()
  //! when it is built; its config() and every output (TDATA.txt,
  //! summary.json, results) give the seed it used, which reruns the same
  //! calculation. Continuing a calculation keeps its seed.
  bool random_seed = false;
  unsigned threads = 0;        //!< 0 = hardware concurrency

  // --- physics -----------------------------------------------------------
  physics::ElectronicStoppingModel stopping_model =
      physics::ElectronicStoppingModel::builtin;
  //! SRIM stopping tables (file or SR-Module output directory) for
  //! ElectronicStoppingModel::srim_tables; empty = the SRIM-2013 tables that
  //! come with AMCSET (data/srim2013).
  std::string stopping_table;
  /*!
   * \brief SRIM installation or its Data directory, for SRIM-2013's
   * pair-specific interatomic potentials (POTCOEF03.DAT). Optional.
   */
  std::string srim_data_directory;
  /*!
   * \brief Which collisions use the SRIM pair potentials (when
   * srim_data_directory is set). An AMCSET option: TRIM-2013 itself uses
   * the ZBL universal potential for every collision.
   */
  enum class PairPotentialUse { none, recoils, all };
  //! Default none: SRIM-2013's TRIM.exe never reads POTCOEF03.DAT and uses the
  //! ZBL universal potential for every collision.
  PairPotentialUse pair_potentials = PairPotentialUse::none;
  /*!
   * \brief Transport algorithm.
   *
   * trim2013 (default) follows SRIM-2013's TRIM.exe step by step, as found by
   * disassembling it (docs/TRIM2013_ALGORITHM.md): its free-flight-path rule
   * without an energy-loss limit, collisions at the energy after the
   * electronic loss of the flight, uniform (Bohr) straggling, monolayer
   * steps for recoils, sub-threshold surface atoms set in motion, recoils
   * reflected by the surface barrier, and its replacement/interstitial
   * bookkeeping. refined is AMCSET's earlier SRIM-mode transport (flight
   * paths limited to a 5 % electronic loss, collisions at the energy before
   * the loss, free flight paths for recoils), which converges better but
   * departs from TRIM where TRIM's steps are coarse (slow light ions).
   */
  /*!
   * trim2013: TRIM.exe step by step (default). refined: AMCSET's earlier
   * algorithm. crystal: follows the lattice of layers that have a crystal
   * (collision partners on lattice sites with thermal vibrations,
   * simultaneous collisions, impact-parameter-dependent electronic
   * stopping; after Robinson and Torrens, Phys. Rev. B 9, 5008 (1974));
   * amorphous layers get random partners. Not TRIM.
   */
  enum class Transport { trim2013, refined, crystal };
  Transport transport = Transport::trim2013;
  /*!
   * \brief Settings of the crystal transport (physics.crystal; not TRIM).
   */
  struct CrystalPhysics {
    /*!
     * Largest impact parameter of a collision. 0: for each crystal the
     * radius of the largest sphere that fits between its atoms
     * (crystal::largest_hole_radius_A), the smallest value that leaves no
     * part of the lattice free of collisions.
     */
    length_quantity max_impact_parameter = 0.0 * common::angstrom;
    /*!
     * Atoms whose closest approaches lie within this distance along the
     * path of the first are struck together: each collision is computed as
     * if alone and the momentum changes are added (as MARLOWE does).
     */
    length_quantity simultaneous_distance = 0.5 * common::angstrom;
    /*!
     * Fraction of the electronic stopping that is lost in the collisions,
     * as exp(-decay R/a) of the distance of closest approach R (a: the
     * pair's screening length), normalised to the stopping power for a
     * random path; the rest is lost along the path. Channeled particles
     * avoid close collisions and lose less. 0: no dependence on the path.
     */
    double local_stopping_fraction = 0.5;
    double local_stopping_decay = 0.3; //!< 0.3: the slowest exponential of the Moliere screening function
    bool thermal_vibrations = true;    //!< false: a static lattice whatever the crystals say

    bool operator==(const CrystalPhysics &) const = default;
  };
  CrystalPhysics crystal;
  physics::ScatteringMethod scattering = physics::ScatteringMethod::magic;
  physics::BetheBlochOptions bethe_bloch;
  //! Corrections of the SRIM-2013 stopping tables (not TRIM; off by
  //! default; physics.stopping_corrections): physics::StoppingCorrections.
  physics::StoppingCorrections stopping_corrections;
  bool electronic_straggling = true;   //!< TRIM's Eq. 7-47 straggling
  /*!
   * \brief Distribution of the straggled electronic loss of one flight.
   * truncated_gaussian (a normal deviate clipped at zero, the default)
   * reproduces TRIM, including its upward bias of the mean loss when the
   * straggling is comparable to the loss (slow recoils); gamma keeps the
   * mean loss exact. The refined transport only: the crystal transport
   * always uses gamma (its ~1 A steps would make the clipping bias large).
   */
  enum class StragglingDistribution { gamma, truncated_gaussian };
  StragglingDistribution straggling_distribution =
      StragglingDistribution::truncated_gaussian;
  bool free_flight_path = true;        //!< false: monolayer steps everywhere
  /*!
   * \brief Minimum energy transfer T_min defining TRIM's maximum impact
   * parameter (SRIM book Eq. 7-40). Collisions below it are skipped.
   */
  energy_quantity min_energy_transfer = 2.0 * common::electron_volt;
  /*!
   * \brief Upper limit of the maximum impact parameter, in Wigner-Seitz
   * radii (SRIM book p. 7-11: "up to two Seitz-Wigner radii"). It sets the
   * shortest free flight path, 1/(pi p_max^2 N) (about 0.5 Angstrom in
   * metals). Zero (default, which reproduces TRIM-2013's cascades) selects
   * the TRIM-85 rule of at least one monolayer, N^(-1/3), between
   * collisions.
   */
  double max_impact_parameter_ws = 0.0;
  /*!
   * \brief Final energy E_f: moving atoms (ion and recoils) are followed
   * until their energy falls below it, and a displacing collision is a
   * replacement collision when the moving atom is the same element as the
   * struck atom and keeps less than E_f. TRIM-2013 uses 1 eV.
   */
  energy_quantity final_energy = 1.0 * common::electron_volt;
  //! Override of final_energy for the ion (0 = use final_energy).
  energy_quantity ion_stop_energy = 0.0 * common::electron_volt;
  //! Override of final_energy for recoils (0 = use final_energy).
  energy_quantity recoil_stop_energy = 0.0 * common::electron_volt;
  /*!
   * \brief Subtract the lattice binding energy Eb from the recoil energy
   * before applying the Kinchin-Pease formula in quick mode.
   */
  bool kp_subtract_binding_energy = false;
  /*!
   * \brief Target atoms within this depth of the front surface are set in
   * motion (and may be sputtered) when they receive more than their surface
   * binding energy Es, even below Ed; such atoms leave no defect if they stay
   * in the target. Full-cascade and monolayer modes only. Zero disables it.
   */
  length_quantity surface_atom_depth = 0.0 * common::angstrom;
  //! Maximum recoils followed per incident ion (0 = unlimited).
  std::size_t max_recoils_per_ion = 0;
  /*!
   * \brief Fast approximations (off by default). The magic-formula
   * scattering angle comes from a precomputed table (physics::MagicTable,
   * relative error below 4e-5), and sines, cosines, a few powers and the
   * normalisation of directions use cheaper forms accurate to a few ulp.
   * Results agree with the exact calculation statistically but not number
   * for number: the same seed can give different individual histories (in
   * practice most histories stay identical), so a run with fast_math does
   * not reproduce a run without it, nor TRIM.exe step by step. Pair
   * potentials and ScatteringMethod::integral are always evaluated exactly.
   */
  bool fast_math = false;

  // --- tallies and output -------------------------------------------------
  //! Depth window of the tabulated distributions (end 0 = whole target).
  length_quantity window_start = 0.0 * common::angstrom;
  length_quantity window_end = 0.0 * common::angstrom;
  std::size_t depth_bins = 100;           //!< TRIM uses 100
  bool record_ion_positions = true;       //!< RANGE_3D.txt
  bool record_exiting_atoms = true;       //!< BACKSCAT / TRANSMIT / SPUTTER lists
  std::size_t max_exit_records = 1000000; //!< Cap for the exit lists
  //! Write COLLISON.txt cascade summaries for the first N ions (0 = none).
  std::size_t collision_details_ions = 0;
  //! COLLISON.txt also lists every displacement and replacement of the
  //! cascades (TRIM's Collisions = 2, "Ion + Recoils"). The file (and the
  //! memory it takes) grows with the number of displacements.
  bool collision_recoil_details = false;
  /*!
   * \brief TRIM's EXYZ.txt: record the ion's position, energy and stopping
   * every time its energy has dropped by this much (0 = off), at its start
   * and where it stops.
   */
  energy_quantity exyz_interval = 0.0 * common::electron_volt;
  //! Keep trajectories (for the GUI or #trajectory_export) of the first
  //! this many ions (0 = none, all_ions = every ion: "all" in a scenario).
  std::size_t trajectory_ions = 0;
  static constexpr std::size_t all_ions = static_cast<std::size_t>(-1);
  /*!
   * \brief Export of the recorded trajectories to a file (not a TRIM
   * output; every transport). One row per recorded point of the first
   * #trajectory_ions ions, in ion order for any number of threads: ion,
   * point, track, parent track and point, kind, generation, Z, mass,
   * position, direction, energy, energy before the collision, energy
   * transferred, electronic loss, time. docs/OUTPUT.md "Trajectory file".
   */
  struct TrajectoryExport {
    std::string file;          //!< "" = off; *.csv or *.parquet; relative paths: in output_directory
    bool births_only = false;  //!< Only the first point of every track (where atoms start to move)
    //! Only points whose energy is at least this (0 = every point).
    energy_quantity min_energy = 0.0 * common::electron_volt;
  };
  TrajectoryExport trajectory_export;
  /*!
   * \brief Recoil events for coupling to molecular dynamics (hybrid BCA-MD;
   * not TRIM; every transport, with full cascades or monolayer steps).
   * For each transition energy E_c: a recoil set in motion with
   * min_energy <= T < E_c, and a moving atom (the ion or a recoil) at the
   * collision where its energy first falls below E_c, are events; nothing
   * that descends from an event's atom is an event again for that E_c.
   * Record mode only writes them (the transport and every tally are
   * unchanged); handoff mode (one E_c) also stops following them, as MD
   * would take them over. docs/OUTPUT.md "MD event file".
   */
  struct MdEvents {
    std::string file;  //!< "" = off; *.csv or *.parquet; relative paths: in output_directory
    std::vector<energy_quantity> transition_energies; //!< E_c, ascending (1 to 32 values)
    //! Recoils set in motion below this are not events (they stay with the
    //! BCA); 0 = the struck element's displacement energy Ed.
    energy_quantity min_energy = 0.0 * common::electron_volt;
    bool handoff = false;  //!< Stop following events (exactly one E_c)
    bool enabled() const noexcept { return !file.empty(); }
  };
  MdEvents md_events;
  /*!
   * \brief The fluence the dpa table is given for, in ions per cm^2 (a
   * count per area kept as a number at this commented boundary; the dpa
   * scale linearly with it). Default 1e15 ions/cm^2.
   */
  double dpa_fluence_per_cm2 = 1e15;
  //! Write the output files here after run() (empty = do not write).
  std::string output_directory;
  /*!
   * \brief Which files to write and how (formats, tables, columns, rows,
   * units, number formats). Format::trim selects TRIM's text files (and
   * AMCSET's summary.json); CSV files go to `csv/`, Parquet files to
   * `parquet/` and the JSON document to `results.json` in the output
   * directory.
   */
  output::OutputSpec outputs;
  //! The scenario as JSON, copied into results.json (set by config::make_engine).
  std::string configuration_json;
  /*!
   * \brief Write the count tables of the TRIM-format files (VACANCY, NOVAC,
   * RANGE) as TRIM.exe does, summed in single precision (up to about 1.5 %
   * from the exact counts; see trim/single_precision.h). The totals, the
   * other files and the CSV/JSON/Parquet tables are exact either way.
   */
  bool single_precision_tables = true;
  //! Bugs of TRIM.exe reproduced by the trim2013 transport (true = as TRIM).
  TrimBugs bugs;
  /*!
   * \brief The bugs a run takes: #bugs, except with the crystal transport,
   * which is not TRIM and takes none of them whatever #bugs says (the
   * scenario keeps its switches for the other transports). TrimSimulation
   * runs with these, and results.json records them.
   */
  TrimBugs effective_bugs() const;
  //! #single_precision_tables; false with the crystal transport (as effective_bugs()).
  bool effective_single_precision_tables() const;

  /*!
   * \brief Automatic target depth (off by default; see size_target()). The
   * last layer (the substrate) is sized from quick-mode test runs so that
   * the deepest ion stops in its last #margin, before the calculation.
   * Leave it off to study transmission through a given thickness.
   */
  struct AutoDepth {
    bool enabled = false;
    std::size_t test_ions = 1000; //!< Ions per test run (quick mode)
    double margin = 0.1;          //!< Sizing puts the deepest test ion this fraction before the back
    //! Accepted while the deepest test ion reaches this fraction of the
    //! target; shallower, the last layer shrinks (used as at most 1 - margin).
    double shrink_below = 0.7;
    //! Front ends warn (without stopping) when ions leave through the back.
    bool warn_transmitted = true;
  };
  AutoDepth auto_depth;

  /*!
   * \brief Non-fatal problems: values that are allowed but probably not
   * meant (a density far from the material's, a gas flag that does not fit
   * the density, a target far thicker than the ions can reach, a
   * relativistic projectile). Never throws.
   */
  std::vector<std::string> plausibility_warnings() const;

  //! Upper limits checked by validate(): worker threads, and depth bins
  //! (every block of 16 ions holds ~30 histograms of depth_bins values).
  static constexpr unsigned max_threads = 4096;
  static constexpr std::size_t max_depth_bins = 10000;
  //! Highest ion energy: the stopping curves end at 100 GeV (SRIM tables,
  //! then relativistic Bethe-Bloch), beyond which they would be extrapolated.
  static constexpr double max_ion_energy_eV = 100e9;
  //! Upper limits of the target-atom energies (displacement energies are
  //! 5-100 eV, lattice and surface binding energies a few eV).
  static constexpr double max_displacement_energy_eV = 10000.0;
  static constexpr double max_binding_energy_eV = 100.0;
  //! Final (stop) energies and the refined transport's T_min: TRIM uses 1 eV
  //! and 2 eV; above 1 keV the low-energy physics is no longer followed.
  static constexpr double max_final_energy_eV = 1000.0;
  //! Densest material at normal conditions: osmium, 22.6 g/cm3.
  static constexpr double max_density_g_cm3 = 50.0;
  static constexpr double max_compound_correction = 10.0;

  /*!
   * \brief Masses accepted for an atom of atomic number \p z, in amu: from
   * 0.9 Z to 3 Z + 10 (every known isotope; an iron atom of 2 amu is a typo).
   */
  static std::pair<double, double> mass_range_amu(unsigned z) {
    return {0.9 * static_cast<double>(z), 3.0 * static_cast<double>(z) + 10.0};
  }

  //! TRIM's cascade code (TRIM.IN): 1-3 for the beam, 4/5 for TRIM.DAT
  //! ions, 6/7 for TRIM.DAT recoils (odd = full cascades except 1).
  int trim_code() const;
  //! Number of incident particles this configuration runs.
  std::size_t particle_count() const;
  /*!
   * \brief TRIM's list of target atoms: one element per distinct atom
   * (same name and Z), in TRIM.IN order (TargetElement::trim_atom), atoms
   * without a TRIM.IN number following in order of first appearance.
   */
  std::vector<const TargetElement *> trim_atoms() const;
  //! Total target thickness.
  length_quantity target_thickness() const;
  /*!
   * \brief Check consistency and the physical limits; throws
   * std::invalid_argument with a message naming the offending field, the
   * limit and what to use instead. With \p complete false (editing in the
   * GUI) a particle source may still lack its particle list.
   */
  void validate(bool complete = true) const;
  //! Human-readable summary.
  std::string describe() const;
};

std::string to_string(TrimConfig::Transport t);
//! A seed from the system's entropy source (std::random_device), mixed with
//! the clock and a call counter; 0 to 2^31 - 1.
std::uint32_t draw_random_seed();
//! "trim2013" ("trim"), "refined", "crystal". \throws std::invalid_argument
TrimConfig::Transport parse_transport(std::string_view s);

/*!
 * \brief The ions' path length (CSDA range) R = integral dE / (S_e(E) +
 * S_n(E)) from the ion energy to the stop energy, through the layers in
 * order (electronic stopping from the configured model, ZBL nuclear
 * stopping). It exceeds the mean projected range, but range straggling
 * takes the deepest ions beyond it (see estimate_deepest()).
 * \throws std::runtime_error if the stopping tables cannot be loaded
 */
length_quantity estimate_penetration(const TrimConfig &config);

/*!
 * \brief Rough depth of the deepest of ~1000 ions: R (1.2 + f_n), with R the
 * path length and f_n the path-averaged nuclear share of the stopping. The
 * first guess of size_target() and the suggestion of the plausibility
 * warnings.
 * \throws std::runtime_error if the stopping tables cannot be loaded
 */
length_quantity estimate_deepest(const TrimConfig &config);

//! What size_target() did.
struct AutoDepthResult {
  bool changed = false;              //!< The last layer's thickness was changed
  length_quantity last_layer_width;  //!< Its thickness now
  length_quantity deepest;           //!< Deepest stopped ion of the last test run
  std::size_t runs = 0;              //!< Test runs made
  std::size_t transmitted = 0;       //!< Ions transmitted in the last test run
  std::vector<std::string> notes;    //!< What was done, for the user
};

/*!
 * \brief Automatic target depth: size the last layer of \p config (in place)
 * so that the deepest ion of a quick-mode test run of
 * config.auto_depth.test_ions ions stops inside the target, in its last
 * config.auto_depth.margin. Starts from estimate_penetration(), grows the
 * layer while ions are transmitted and shrinks it while the deepest ion
 * stops shallower; at most 12 test runs. With the crystal transport the
 * first test run is made in a box 32 times as deep as that estimate
 * (channeled ions stop far beyond it), which is then trimmed to the deepest
 * ion. Not used for recoil sources (TRIM.DAT atoms start at given depths).
 * \p cancel (optional) stops early.
 * \throws std::invalid_argument, std::runtime_error (invalid configuration)
 */
AutoDepthResult size_target(TrimConfig &config, const std::atomic<bool> *cancel = nullptr);

// ---------------------------------------------------------------------------
// Results
// ---------------------------------------------------------------------------

//! A depth histogram with a fixed bin layout.
struct DepthHistogram {
  std::vector<double> values;
  void resize(std::size_t n) { values.assign(n, 0.0); }
  void add(const DepthHistogram &o);
};

//! Running power sums of a sample (for mean, variance, skewness, kurtosis).
struct Moments {
  double n = 0.0;
  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;
  double s4 = 0.0;
  void add(double x) {
    const double x2 = x * x;
    n += 1.0;
    s1 += x;
    s2 += x2;
    s3 += x2 * x;
    s4 += x2 * x2;
  }
  void add(const Moments &o) {
    n += o.n;
    s1 += o.s1;
    s2 += o.s2;
    s3 += o.s3;
    s4 += o.s4;
  }
  double mean() const { return n > 0 ? s1 / n : 0.0; }
  //! Population standard deviation (TRIM's "straggling").
  double stddev() const;
  //! Population skewness (TRIM's definition).
  double skewness() const;
  //! Population kurtosis (TRIM's definition, 3 for a Gaussian).
  double kurtosis() const;
  //! Root mean square.
  double rms() const { return n > 0 ? std::sqrt(s2 / n) : 0.0; }
};

/*!
 * \brief One hybrid BCA-MD recoil event (TrimConfig::md_events). Positions
 * and directions use AMCSET's trajectory axes: z is the depth, x and y are
 * lateral (TRIM's y and z).
 */
struct MdEvent {
  enum class Kind : std::uint8_t {
    birth,    //!< A recoil set in motion with E_min <= T < E_c (E = T)
    crossing  //!< A moving atom whose energy fell below E_c at a collision or on a crystal free step (E after it)
  };
  std::size_t ion = 0;              //!< 1-based incident particle
  std::size_t event = 0;            //!< 1-based within the ion and E_c, in time order
  double transition_energy_eV = 0.0; //!< The E_c it belongs to
  Kind kind = Kind::birth;
  double time_s = 0.0;              //!< Since the incident particle entered the target
  double x_A = 0.0, y_A = 0.0, z_A = 0.0; //!< z = depth
  double energy_eV = 0.0;
  double ux = 0.0, uy = 0.0, uz = 0.0;    //!< Unit direction (same axes)
  unsigned generation = 0;          //!< 0 = the incident ion, 1 = its primaries, ...
  unsigned z = 0;                   //!< Atomic number
  double mass_amu = 0.0;
  std::uint32_t track = 0;          //!< The atom's track (as in the trajectory file)
  std::uint32_t parent_track = 0xffffffffu; //!< Track of the atom that set it in motion (none: 0xffffffff)
};

//! Totals of the MD events of one transition energy (TrimTallies::md).
struct MdTransitionStats {
  double transition_energy_eV = 0.0;
  double events = 0.0;
  double birth_events = 0.0;        //!< Of which: recoils born below E_c
  double crossing_events = 0.0;     //!< Of which: movers falling below E_c
  double energy_eV = 0.0;           //!< Sum of the events' energies (handed to MD)
  //! Recoils set in motion below E_min (not events; left to the BCA).
  double below_min = 0.0;
  double below_min_energy_eV = 0.0;
  double ions_handed = 0.0;         //!< Incident ions that became events
  void add(const MdTransitionStats &o) noexcept {
    events += o.events;
    birth_events += o.birth_events;
    crossing_events += o.crossing_events;
    energy_eV += o.energy_eV;
    below_min += o.below_min;
    below_min_energy_eV += o.below_min_energy_eV;
    ions_handed += o.ions_handed;
  }
};

//! One atom leaving the target (TRIM's BACKSCAT/TRANSMIT/SPUTTER records).
struct ExitRecord {
  char kind = 'B';          //!< 'B' backscattered ion, 'T' transmitted, 'S' sputtered
  std::size_t ion = 0;      //!< 1-based ion number
  unsigned z = 0;           //!< Atomic number
  //! Kinetic energy leaving the target. With Transport::trim2013, sputtered
  //! atoms are listed as TRIM lists them: the energy inside the surface (the
  //! surface binding energy is not subtracted) and the cosines before the
  //! surface barrier.
  double energy_eV = 0.0;
  double x = 0.0, y = 0.0, z_pos = 0.0; //!< Last position (Angstrom)
  double cx = 0.0, cy = 0.0, cz = 0.0;  //!< Direction cosines
};

/*!
 * \brief One ion collision that displaced a target atom, with the summary of
 * the cascade it started (TRIM's COLLISON.txt).
 */
//! One displacement or replacement in a cascade (TRIM's COLLISON.txt with
//! "Ion + Recoils", in the order they happen).
struct CascadeAtom {
  unsigned number = 0;           //!< 1-based position in the cascade's list
  unsigned z = 0;                //!< Atomic number of the atom
  double energy_eV = 0.0;        //!< Energy of a displaced atom; 0 for a replacement
  double x = 0.0, y = 0.0, z_pos = 0.0; //!< Position (Angstrom)
  bool vacancy = false;          //!< The atom left a vacancy
  bool replacement = false;      //!< The atom replaced one of its own element
};

struct CascadeRecord {
  std::size_t ion = 0;          //!< 1-based ion number
  double ion_energy_keV = 0.0;  //!< Ion energy before the collision
  double x = 0.0, y = 0.0, z = 0.0; //!< Collision position (Angstrom)
  double se_eV_per_A = 0.0;     //!< Electronic stopping power of the ion there
  unsigned z_hit = 0;           //!< Atomic number of the struck atom
  double recoil_energy_eV = 0.0; //!< Energy transferred to the primary recoil
  double displacements = 0.0;   //!< Cascade displacements (incl. the primary)
  double vacancies = 0.0;
  double replacements = 0.0;
  double interstitials = 0.0;
  //! The cascade's displacements and replacements (only with
  //! TrimConfig::collision_recoil_details).
  std::vector<CascadeAtom> atoms;
};

//! Damage made by one ion (COLLISON.txt's summary after each ion).
struct IonCollisionSummary {
  std::size_t ion = 0;          //!< 1-based ion number
  double displacements = 0.0;
  double replacements = 0.0;
  double vacancies = 0.0;
  double interstitials = 0.0;
  double sputtered = 0.0;       //!< Target atoms sputtered
  double transmitted = 0.0;     //!< Target atoms leaving through the back
};

//! One line of TRIM's EXYZ.txt: the ion's state at an energy step.
struct EnergyPositionRecord {
  std::size_t ion = 0;          //!< 1-based ion number
  double energy_keV = 0.0;      //!< Ion energy
  double x = 0.0, y = 0.0, z = 0.0; //!< Position (Angstrom)
  double se_eV_per_A = 0.0;     //!< Electronic stopping power at that energy
  double last_recoil_eV = 0.0;  //!< Energy given to the last struck atom
};

//! Final position of one stopped ion (RANGE_3D.txt).
struct IonPosition {
  std::size_t ion = 0;
  double x = 0.0, y = 0.0, z = 0.0; //!< Angstrom
};

/*!
 * \brief Everything TRIM tabulates, in plain numbers (counts, eV, Angstrom).
 *
 * Depth histograms hold totals over all ions (not yet divided by the ion
 * count or the bin width); TrimResults' accessors and the writers normalize.
 * Per-element arrays are indexed by the "target element index" of
 * TrimSimulation::element_names().
 */
struct TrimTallies {
  std::size_t ions = 0;

  // Ion distributions
  DepthHistogram ion_stopped;          //!< Ions stopped per bin
  Moments range;                       //!< Depth of stopped ions (Angstrom)
  Moments lateral_y;                   //!< y of stopped ions
  Moments lateral_abs_y;               //!< |y|
  Moments radial;                      //!< sqrt(y^2 + z^2)
  //! Per-bin lateral statistics of stopped ions: n, sum|y|, sum y^2, sum r, sum r^2.
  std::vector<std::array<double, 5>> lateral_bins;

  // Recoil distributions (per element)
  std::vector<DepthHistogram> recoil_stopped; //!< Interstitial target atoms

  // Damage
  DepthHistogram knock_ons;                   //!< Primary knock-ons from the ion
  std::vector<DepthHistogram> vacancies;      //!< Per element (full cascade)
  //! Per element: vacancies removed again by replacements (counted in
  //! \c vacancies with -1; kept for TRIM's single-precision tables).
  std::vector<DepthHistogram> vacancies_removed;
  DepthHistogram vacancies_by_ion;            //!< Quick (KP) mode
  DepthHistogram vacancies_by_recoils;        //!< Quick (KP) mode
  DepthHistogram replacements;                //!< Replacement collisions
  DepthHistogram displacements;               //!< Displacements (all)
  DepthHistogram interstitials;               //!< Stopped recoils + ions that ended interstitial

  // Energy deposition (eV)
  DepthHistogram ionization_by_ion;
  DepthHistogram ionization_by_recoils;
  DepthHistogram phonons_by_ion;
  DepthHistogram phonons_by_recoils;
  DepthHistogram energy_to_recoils;               //!< Transferred by the ion (E2RECOIL)
  /*!
   * \brief Damage models of every primary knock-on atom (the atoms the
   * incident particle sets in motion, booked where they start; trim2013
   * transport): its damage energy T_d (TRIM's Robinson partition, as in the
   * quick calculation), and its displacements by NRT, by the original
   * Kinchin-Pease rule (T_d/(2 Ed)), by arc-dpa and its replacements by rpa
   * (the last two only for elements with DamageEfficiency constants).
   */
  DepthHistogram pka_damage_energy;
  DepthHistogram nrt_displacements;
  DepthHistogram kp_displacements;
  DepthHistogram arc_displacements;
  DepthHistogram rpa_replacements;
  double total_pka_damage_energy = 0.0;
  double total_nrt = 0.0;
  double total_kp = 0.0;
  double total_arc = 0.0;
  double total_rpa = 0.0;
  double pka_without_efficiency = 0.0; //!< PKAs of elements without arc-dpa/rpa constants
  std::vector<DepthHistogram> energy_absorbed;    //!< Per element (E2RECOIL)

  // Totals
  double total_displacements = 0.0;
  double total_vacancies = 0.0;
  double total_replacements = 0.0;
  double total_interstitials = 0.0;
  double total_ionization_ion = 0.0;
  double total_ionization_recoils = 0.0;
  double total_phonons_ion = 0.0;
  double total_phonons_recoils = 0.0;
  double total_energy_to_recoils = 0.0;
  double total_binding_energy_loss = 0.0;  //!< Eb lost by displaced atoms
  std::size_t backscattered = 0;
  double backscattered_energy = 0.0;
  std::size_t transmitted = 0;
  double transmitted_energy = 0.0;
  std::vector<double> sputtered;           //!< Per element (atoms)
  std::vector<double> sputtered_energy;    //!< Per element (eV)
  std::vector<double> transmitted_recoils; //!< Per element (atoms)
  double recoils_followed = 0.0;
  double recoil_limit_hits = 0.0;

  /*!
   * \brief Energy spectra of displaced atoms (recoil energy T - Eb), log
   * bins of spectrum_bins_per_decade from 1 eV: primary knock-ons made by
   * the ion, and all recoils (full cascades only).
   */
  std::vector<double> pka_spectrum;
  std::vector<double> recoil_spectrum;
  static constexpr int spectrum_bins_per_decade = 10;
  static constexpr int spectrum_decades = 10;

  /*!
   * \brief Atoms reaching the front surface from inside, by element and by
   * the energy of their motion normal to the surface (E cos^2 x), in bins of
   * surface_energy_bin_eV (TRIM's histogram behind its sputtering plots:
   * bin k counts round(10 E_n) = k, limited to 1..surface_energy_bins).
   * Atoms with E_n above the surface binding energy are sputtered, the
   * others reflected; summing the bins above any Es gives the yield the
   * target would have with that surface binding energy.
   */
  std::vector<std::vector<double>> surface_energy;
  static constexpr std::size_t surface_energy_bins = 200;
  static constexpr double surface_energy_bin_eV = 0.1;

  std::vector<ExitRecord> exits;
  std::vector<IonPosition> positions;
  std::vector<CascadeRecord> cascades; //!< For the first collision_details_ions ions
  std::vector<IonCollisionSummary> ion_summaries; //!< Same ions
  std::vector<EnergyPositionRecord> exyz; //!< EXYZ.txt (TrimConfig::exyz_interval)
  //! MD event totals, one per TrimConfig::md_events E_c (empty if off).
  std::vector<MdTransitionStats> md;

  void resize(std::size_t bins, std::size_t elements);
  void add(const TrimTallies &o, std::size_t max_exit_records);
  //! Copy without the per-particle lists (exits, positions, cascades, exyz).
  TrimTallies without_lists() const;
};

// ---------------------------------------------------------------------------
// The engine
// ---------------------------------------------------------------------------

class TargetModel; // internal precomputed tables (trim_engine.cpp)

/*!
 * \brief SRIM-mode simulation. Construct, run(), then read results() or
 * write_outputs().
 */
class TrimSimulation final : public common::SimulationEngine {
public:
  /*!
   * \brief Validate the configuration, load data files and precompute all
   * stopping and scattering tables. \throws std::invalid_argument,
   * std::runtime_error
   */
  explicit TrimSimulation(TrimConfig config);
  ~TrimSimulation() override;

  TrimSimulation(const TrimSimulation &) = delete;
  TrimSimulation &operator=(const TrimSimulation &) = delete;

  // SimulationEngine
  int run() override;
  void request_cancel() noexcept override;
  bool is_cancelled() const noexcept override { return cancel_.load(); }
  std::size_t completed_count() const noexcept override { return completed_.load(); }
  std::size_t total_count() const noexcept override { return total_; }
  bool supports_pause() const noexcept override { return true; }
  void set_paused(bool paused) override;
  bool is_paused() const noexcept override { return paused_.load(); }
  //! Statistics in the SimulationEngine form (derived from snapshot()).
  common::Tallies snapshot_tallies() const override;
  bool try_pop_trajectory(std::vector<common::CollisionPoint> &out) override;
  void clear_trajectories() override;
  length_quantity depth_extent() const override;
  std::string describe() const override { return config_.describe(); }

  //! The configuration.
  const TrimConfig &config() const noexcept { return config_; }
  //! Tallies of the finished run (valid after run() returns).
  const TrimTallies &tallies() const noexcept { return tallies_; }
  /*!
   * \brief Thread-safe copy of the tallies of the ions finished so far (in
   * ion order), for live plots while run() executes. The per-particle lists
   * are copied only if \p with_lists.
   */
  TrimTallies snapshot(bool with_lists = false) const;
  /*!
   * \brief Start from the tallies of an earlier run of the same target (call
   * before run()), to continue a calculation: set TrimConfig::first_ion to
   * the number of ions already run. \throws std::invalid_argument if the
   * tallies do not fit this target.
   */
  void continue_from(const TrimTallies &previous);
  //! Problems found in the input (particle list lines TRIM would skip,
  //! target energies TRIM replaces).
  const std::vector<std::string> &warnings() const noexcept { return warnings_; }
  /*!
   * \brief The target layers as the calculation uses them: with
   * Transport::trim2013, a displacement, surface or lattice binding energy
   * of zero (or less) of every target atom after the first (in TRIM's atom
   * order, TrimConfig::trim_atoms()) is replaced by 20, 4 or 2 eV, as
   * TRIM.exe does; the first atom keeps its values.
   */
  const std::vector<TargetLayer> &effective_layers() const noexcept { return layers_; }
  //! Names of the target elements, in the index order of the tallies.
  const std::vector<std::string> &element_names() const noexcept { return element_names_; }
  //! Atomic numbers of the target elements (tally index order).
  const std::vector<unsigned> &element_z() const noexcept { return element_z_; }
  //! Depth of the lower edge of tally bin \p i (Angstrom).
  double bin_start(std::size_t i) const;
  //! Width of a tally bin (Angstrom).
  double bin_width() const;
  //! Lowest displacement energy of the target (TRIM's recoil cut-off).
  energy_quantity lowest_displacement_energy() const;
  //! Wall-clock seconds spent in run().
  double elapsed_seconds() const noexcept { return elapsed_seconds_.load(); }
  //! Ions transmitted through the back so far (merged blocks; cheap to
  //! poll while run() executes, for the front ends' escape warning).
  std::size_t transmitted_so_far() const noexcept { return transmitted_.load(std::memory_order_relaxed); }

  /*!
   * \brief Write all TRIM-format output files and summary.json to
   * \p directory (created if needed).
   */
  void write_outputs(const std::filesystem::path &directory) const;

  //! Machine-readable summary of the results (JSON text).
  std::string summary_json() const;

  /*!
   * \brief Whether try_pop_trajectory() gets the recorded trajectories
   * (default true, for the GUI). A front end that never polls (the CLI)
   * turns it off, so that trajectories recorded for
   * TrimConfig::trajectory_export are not kept in memory as well.
   */
  void set_trajectory_polling(bool on) noexcept { poll_trajectories_ = on; }

  struct ExportBlock;   //!< A block's rows for the export files (trim_export.cpp)
  struct ExportWriters; //!< The open export files of a run

private:
  void run_block(std::size_t first_ion, std::size_t count, TrimTallies &out,
                 diag::Activity *activity, ExportBlock *exported = nullptr);

  TrimConfig config_;
  std::vector<TargetLayer> layers_; //!< effective_layers()
  std::unique_ptr<TargetModel> model_;
  std::size_t total_ = 0;
  std::vector<std::string> warnings_;
  std::atomic<bool> paused_{false};
  std::mutex pause_mutex_;
  std::condition_variable pause_cv_;
  mutable std::mutex tallies_mutex_;
  std::vector<std::string> element_names_;
  std::vector<unsigned> element_z_;
  TrimTallies tallies_;
  std::atomic<bool> cancel_{false};
  std::atomic<std::size_t> completed_{0};
  mutable std::mutex trajectory_mutex_;
  std::deque<std::vector<common::CollisionPoint>> trajectories_;
  bool poll_trajectories_ = true;
  std::atomic<double> elapsed_seconds_{0.0};
  std::atomic<std::size_t> transmitted_{0};
};

} // namespace trim
} // namespace amcset
