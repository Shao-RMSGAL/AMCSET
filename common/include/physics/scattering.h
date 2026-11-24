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

/*! \file scattering.h
 *  \brief Screened-Coulomb binary collisions: interatomic potentials, the
 *  Biersack-Haggmark "magic formula", the screened Rutherford limit, the exact
 *  scattering integral, kinematics and nuclear stopping.
 *
 * Reduced units (Lindhard): lengths in units of the screening length a,
 * energies as the reduced energy
 *   epsilon = a E_cm / (Z1 Z2 e^2),  E_cm = E M2 / (M1 + M2),
 * impact parameters as b = p / a. Reduced quantities are plain doubles.
 */

#pragma once

#include "amcset_units.h"

#include <array>
#include <cmath>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace amcset {
namespace physics {

using common::energy_quantity;
using common::length_quantity;
using common::mass_quantity;
using common::stopping_cross_section_quantity;

/*!
 * \brief Screening function phi(x) = sum_i c_i exp(-d_i x), x = r / a.
 *
 * The ZBL universal potential and the pair-specific SRIM-2013 potentials
 * both have this four-exponential form.
 */
struct ScreeningFunction {
  std::array<double, 4> c{}; //!< Amplitudes (sum to 1)
  std::array<double, 4> d{}; //!< Decay constants

  //! phi(x)
  double phi(double x) const noexcept {
    return c[0] * std::exp(-d[0] * x) + c[1] * std::exp(-d[1] * x) +
           c[2] * std::exp(-d[2] * x) + c[3] * std::exp(-d[3] * x);
  }
  //! dphi/dx
  double dphi(double x) const noexcept {
    return -(c[0] * d[0] * std::exp(-d[0] * x) + c[1] * d[1] * std::exp(-d[1] * x) +
             c[2] * d[2] * std::exp(-d[2] * x) + c[3] * d[3] * std::exp(-d[3] * x));
  }
  /*!
   * \brief phi(x) and dphi/dx together from the same four exponentials (half
   * the exp calls of phi() + dphi()). The sums are evaluated in the same
   * order as in phi() and dphi(), so the results are bit-identical to them.
   */
  void phi_and_derivative(double x, double &phi_x, double &dphi_x) const noexcept {
    const double e0 = std::exp(-d[0] * x);
    const double e1 = std::exp(-d[1] * x);
    const double e2 = std::exp(-d[2] * x);
    const double e3 = std::exp(-d[3] * x);
    phi_x = c[0] * e0 + c[1] * e1 + c[2] * e2 + c[3] * e3;
    dphi_x = -(c[0] * d[0] * e0 + c[1] * d[1] * e1 + c[2] * d[2] * e2 + c[3] * d[3] * e3);
  }

  //! The ZBL universal screening function (Ziegler, Biersack, Littmark 1985).
  static ScreeningFunction universal() noexcept;
};

//! The ZBL universal screening length a_U = 0.8854 a_0 / (Z1^0.23 + Z2^0.23).
length_quantity universal_screening_length(unsigned z1, unsigned z2);

/*!
 * \brief Pair-specific interatomic potentials and nuclear stopping fits.
 *
 * SRIM-2013 uses a four-exponential screening function fitted separately for
 * every pair (Z1, Z2) (file Data/POTCOEF03.DAT of a SRIM installation) and
 * the matching universal-form nuclear-stopping fits (Data/SNUC03.dat). They
 * are not distributed with AMCSET; load them from a local SRIM installation
 * to reproduce TRIM's nuclear scattering.
 *
 * Both files list one pair per line: "Z1 Z2 c1 d1 c2 d2 c3 d3 c4 d4" and
 * "Z1 Z2 A B C D" for Z1 <= Z2 (the potential is symmetric).
 */
class PairPotentialTable {
public:
  /*!
   * \brief Load from a SRIM "Data" directory (or the SRIM root directory, in
   * which case "Data" is appended). \throws std::runtime_error
   */
  static std::shared_ptr<const PairPotentialTable>
  load_srim(const std::filesystem::path &srim_data_dir);

  //! Screening function for the pair, if tabulated.
  std::optional<ScreeningFunction> screening(unsigned z1, unsigned z2) const;
  /*!
   * \brief Nuclear stopping fit coefficients (A, B, C, D) of
   * s_n(eps) = ln(1 + A eps) / (2 (eps + B eps^C + D eps^0.5)), if tabulated.
   */
  std::optional<std::array<double, 4>> nuclear_stopping_fit(unsigned z1,
                                                            unsigned z2) const;
  //! Number of tabulated potentials.
  std::size_t size() const noexcept { return potentials_.size(); }

  //! Add an entry (used by the loader and tests).
  void add_potential(unsigned z1, unsigned z2, const ScreeningFunction &f);
  void add_nuclear_fit(unsigned z1, unsigned z2, const std::array<double, 4> &f);

private:
  std::map<std::pair<unsigned, unsigned>, ScreeningFunction> potentials_;
  std::map<std::pair<unsigned, unsigned>, std::array<double, 4>> nuclear_;
};

/*!
 * \brief Reduced nuclear stopping s_n(epsilon).
 *
 * Universal ZBL form (or a pair fit when \p fit is given) for epsilon <= 30
 * and the unscreened limit ln(epsilon) / (2 epsilon) above, exactly as SRIM.
 */
double reduced_nuclear_stopping(double epsilon,
                                const std::array<double, 4> *fit = nullptr);

/*!
 * \brief A projectile/target-atom pair with all constants a collision needs.
 */
struct CollisionPair {
  unsigned z1 = 0;
  unsigned z2 = 0;
  mass_quantity m1;
  mass_quantity m2;
  length_quantity screening_length; //!< a (ZBL universal screening length)
  //! epsilon / E: reduced energy per laboratory energy (unit 1/energy).
  common::inverse_energy_quantity reduced_energy_per_energy;
  double mass_ratio = 1.0;  //!< M1 / M2
  double gamma = 1.0;       //!< 4 M1 M2 / (M1 + M2)^2
  ScreeningFunction screening = ScreeningFunction::universal();
  //! True when #screening is the ZBL universal function (no pair potential);
  //! only then can MagicTable stand in for the magic formula.
  bool universal_screening = true;
  std::optional<std::array<double, 4>> nuclear_fit; //!< SNUC fit, if any

  CollisionPair() = default;
  /*!
   * \brief Build the pair. \p potentials may be null (universal potential).
   */
  CollisionPair(unsigned z1, mass_quantity m1, unsigned z2, mass_quantity m2,
                const PairPotentialTable *potentials = nullptr);

  //! Reduced energy of a projectile with laboratory kinetic energy \p e.
  double reduced_energy(energy_quantity e) const {
    return (e * reduced_energy_per_energy).value();
  }

  //! Nuclear stopping cross section S_n(E) of this pair.
  stopping_cross_section_quantity nuclear_stopping(energy_quantity e) const;
};

//! How the centre-of-mass scattering angle is computed.
enum class ScatteringMethod {
  magic,     //!< Magic formula (eps <= 10) / screened Rutherford (eps > 10), as TRIM
  integral   //!< Numerical (Gauss-Mehler) evaluation of the classical integral
};

/*!
 * \brief sin^2(theta/2) from the magic formula (Biersack & Haggmark,
 * Nucl. Instr. Meth. 174, 257 (1980); SRIM book Eqs. 7-1 to 7-12) for the
 * screening function \p f, reduced energy \p eps and reduced impact
 * parameter \p b. The fit constants C1..C5 are those of the universal
 * potential (the ones TRIM uses).
 */
double magic_sin2_half_theta(const ScreeningFunction &f, double eps, double b);

/*!
 * \brief sin^2(theta/2) of Rutherford scattering with TRIM's screening
 * correction, 1 / (1 + (1 + b (1 + b)) (2 eps b)^2) (SRIM book Eq. 7-15).
 */
double rutherford_sin2_half_theta(double eps, double b) noexcept;

/*!
 * \brief sin^2(theta/2) from the classical scattering integral, evaluated by
 * Gauss-Mehler quadrature with \p nodes nodes (reference/validation method).
 */
double integral_sin2_half_theta(const ScreeningFunction &f, double eps,
                                double b, int nodes = 100);

/*!
 * \brief Reduced distance of closest approach R0 = r0/a for (eps, b).
 */
double closest_approach(const ScreeningFunction &f, double eps, double b);

/*!
 * \brief sin^2(theta/2) as TRIM computes it: Rutherford with screening
 * correction for eps > 10, otherwise the magic formula (or the exact
 * integral if \p method is ScatteringMethod::integral).
 */
double trim_sin2_half_theta(const ScreeningFunction &f, double eps, double b,
                            ScatteringMethod method = ScatteringMethod::magic);

/*!
 * \brief Laboratory-frame outcome of one elastic collision.
 */
struct CollisionKinematics {
  double sin2_half_theta = 0.0; //!< sin^2(theta_cm / 2)
  double cos_theta_cm = 1.0;    //!< cos(theta_cm)
  double cos_projectile = 1.0;  //!< cos of the projectile deflection (lab)
  double sin_projectile = 0.0;  //!< sin of the projectile deflection (lab)
  double cos_recoil = 0.0;      //!< cos of the recoil angle to the incident direction
  double sin_recoil = 1.0;      //!< sin of the recoil angle
  energy_quantity transferred;  //!< Energy given to the target atom
};

/*!
 * \brief Lab-frame kinematics from sin^2(theta_cm/2): energy transfer
 * T = gamma E sin^2(theta/2), projectile deflection
 * tan(psi) = sin(theta) / (cos(theta) + M1/M2) and recoil angle
 * (pi - theta) / 2 (SRIM book Eqs. 7-16, 7-17).
 */
CollisionKinematics collision_kinematics(const CollisionPair &pair,
                                         energy_quantity energy,
                                         double sin2_half_theta);

/*!
 * \brief TRIM's maximum reduced impact parameter for a minimum energy
 * transfer (SRIM book Eqs. 7-40, 7-41):
 * b_max = 1 / (xi + xi^0.5 + 0.125 xi^0.1), xi = (eps eps_min / gamma)^(1/2).
 */
double trim_max_reduced_impact_parameter(double eps, double eps_min,
                                         double gamma) noexcept;

/*!
 * \brief Kinematics from sin(theta/2) and cos(theta/2) given separately
 * (both with full relative precision, which 1 - sin^2 would lose near
 * head-on collisions). Same results as collision_kinematics up to rounding.
 */
CollisionKinematics collision_kinematics_half_angle(const CollisionPair &pair,
                                                    energy_quantity energy,
                                                    double sin_half, double cos_half);

/*!
 * \brief The magic formula for the ZBL universal potential, tabulated.
 *
 * The scattering angle depends only on (eps, b) for the universal potential,
 * so one table serves every pair. ln(theta) is tabulated on a grid of
 * ln(eps) (1e-9 <= eps <= 12, 800 points) and sqrt(b) (0 <= b <= 256, 800
 * points, which covers the monolayer impact parameters of solids and of
 * gases at normal density) and interpolated bicubically (Catmull-Rom); sin
 * and cos of theta/2 follow without cancellation. Against the converged
 * magic formula the largest relative error of sin^2(theta/2) and
 * cos^2(theta/2) is below 4e-5 (mean 5e-8); TRIM.exe's own Newton
 * iteration (|dR/R| < 3e-4) differs from the converged formula by more. A
 * lookup is about 5x faster than the formula. Used by SRIM mode's fast
 * approximations (TrimConfig::fast_math).
 */
class MagicTable {
public:
  //! The table of the universal potential (built on first use, in parallel,
  //! 5 MB; thread-safe).
  static const MagicTable &universal();

  /*!
   * \brief sin(theta/2) and cos(theta/2) for reduced energy \p eps <= 12
   * and impact parameter \p b. \returns false (outputs untouched) outside
   * the table; the caller then evaluates the formula.
   */
  bool half_angle(double eps, double b, double &sin_half, double &cos_half) const noexcept;

  static constexpr int energy_points = 800;
  static constexpr int impact_points = 800;
  static constexpr double min_eps = 1e-9;
  static constexpr double max_eps = 12.0;
  static constexpr double max_b = 256.0;

private:
  MagicTable();
  double u0_ = 0.0, inv_du_ = 1.0, inv_dt_ = 1.0;
  std::vector<double> ln_theta_; // [energy][impact]
};

/*!
 * \brief sin(theta/2) and cos(theta/2) as trim_sin2_half_theta computes the
 * angle, but through MagicTable where it applies (universal potential, magic
 * method, eps <= 10 inside the table); the screened Rutherford branch
 * (eps > 10) is evaluated directly.
 */
void fast_half_angle(const CollisionPair &pair, double eps, double b, ScatteringMethod method,
                     double &sin_half, double &cos_half);

} // namespace physics
} // namespace amcset
