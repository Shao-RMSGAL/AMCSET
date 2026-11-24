# Glossary

Terms used in the code and docs. "TRIM" means TRIM.exe of SRIM-2013.00.

| Term | Meaning | In the code |
|---|---|---|
| **BCA** | Binary-collision approximation: a moving atom interacts with one target atom at a time. Both engines are BCA codes. | |
| **SRIM / TRIM** | SRIM-2013 (J. F. Ziegler) is the package; TRIM is its Monte Carlo transport program, SR Module its stopping-table program. "SRIM mode" is AMCSET's TRIM-compatible engine. | `trim::` |
| **trim2013 / refined** | The two SRIM-mode transports: `trim2013` follows TRIM.exe step by step (default); `refined` is AMCSET's earlier algorithm. | `TrimConfig::Transport` |
| **Native engine** | AMCSET's original BCA engine (ions and electrons), not validated. | `common::Simulation` |
| **Ion / recoil / PKA** | The incident particle; a target atom set in motion; the primary knock-on atom (a recoil struck by the ion). | `Mover`, `generation` |
| **Cascade** | All recoils that follow from one PKA. TRIM follows them with a stack of slots. | `cascade_2013`, `slots_` |
| **Damage modes** | TRIM's "quick Kinchin–Pease" (KP, cascade code 1: recoils are not followed, damage from the KP formula), "full cascades" (FC, code 2: every recoil followed), "monolayer steps" (code 3: also the ion moves one monolayer per step; used for sputtering). | `DamageCalculation::{kinchin_pease, full_cascade, monolayer}` |
| **Cascade codes 4–7** | TRIM.DAT inputs: 4/5 ions from a list (KP/FC), 6/7 recoils from a list (neutron, electron or photon damage). | `TrimConfig::Source` |
| **Ed** | Displacement energy: the energy a target atom needs to leave its site and make a vacancy. SRIM proposes 25 eV for most elements; 40 eV is recommended for Fe (Stoller et al. 2013). A zero in TRIM.IN becomes 20 eV (except for atom 1, a TRIM bug). | `TargetElement::displacement_energy` |
| **Eb** | Lattice binding energy, subtracted from a displaced atom's energy (TRIM: only deeper than 2 monolayers). | `lattice_binding_energy` |
| **Es** | Surface binding energy: planar barrier an atom must overcome (normal energy) to sputter. | `surface_binding_energy` |
| **E_f / stop energy** | Atoms are followed until their energy falls below it; TRIM uses 1 eV for ions and recoils. | `final_energy` |
| **Vacancy / replacement / interstitial** | A displaced atom leaves a vacancy; if the mover stops right after displacing an atom of its own species, it fills the site (replacement collision, no vacancy); a stopped displaced atom is an interstitial. | `recoil_stop_2013` |
| **Displacements** | Vacancies + replacements. | |
| **Kinchin–Pease (modified)** | ν = 0.8 E_v / (2 Ed) vacancies for damage energy E_v > 2.5 Ed (1 for Ed ≤ E_v ≤ 2.5 Ed); E_v from Robinson's partition of T. | `kp_2013` |
| **Unit cell / site / occupancy** | A crystal's repeating cell: lattice vectors a, b, c and sites at fractional coordinates. A site holds each of its elements with a probability (occupancy); `lattice` sites belong to the host (substitutional alloying shares them), `interstitial` sites lie between. Followed by the crystal transport only. | `crystal::UnitCell`, `Site`, `Occupant` |
| **Miller indices, (h k l), [u v w]** | A plane by the reciprocals of its axis intercepts and a direction by its lattice-vector components; in cubic cells [h k l] is normal to (h k l), in others the normal comes from the reciprocal lattice. Hexagonal cells also use four indices (h k i l), [u v t w]. | `crystal::plane_normal`, `direction_vector` |
| **Tilt / twist** | The beam's angle to the surface normal (`ion.angle`) and the crystal's rotation about that normal (`orientation.twist`), as on an implanter (e.g. 7° tilt, 22° twist to avoid channeling in Si). | `crystal::tilt_twist_for_beam` |
| **Debye temperature, thermal amplitude** | Atoms vibrate about their sites; in the Debye model the rms displacement along one axis follows from the Debye temperature, the temperature and the mass (zero-point motion at 0 K). A few hundredths of an Å. | `crystal::debye_rms_displacement` |
| **Channeling** | Ions moving along open rows or planes of a crystal are steered by many glancing collisions and go deeper than in an amorphous target. Simulated by the crystal transport only (`physics.transport: "crystal"`; its stopping division is not calibrated; docs/CLI.md "Crystals"). The critical angle ψ₁ = (2 Z₁ Z₂ e²/(E d))^½ (Lindhard; d: atom spacing along the row) sets the width of the effect. | `crystal_follow`, `crystal::LatticeSearch` |
| **NRT-dpa / arc-dpa / rpa / dpa** | Damage functions of a PKA's damage energy T_d: NRT (Norgett-Robinson-Torrens, ASTM E521; equal to TRIM's modified Kinchin-Pease), arc-dpa (NRT times an athermal-recombination efficiency fitted to MD: surviving Frenkel pairs) and rpa (replaced atoms; atomic mixing), Nordlund et al. 2018. dpa: displacements per atom = displacements per volume / atomic density, at a fluence. | `book_pka`, `DamageEfficiency`, table `dpa` |
| **Robinson partition / damage energy** | E_v = T / (1 + k g(ε_d)): the part of a recoil's energy that goes to nuclear collisions (the rest to electrons). | `kp_eps_per_eV`, `kp_k` |
| **Ionization / phonons / energy to recoils** | TRIM's energy partition: energy to electrons, to lattice vibrations, and transferred to recoils, each by ions and by recoils. | `TrimTallies` |
| **Sputtering yield** | Target atoms leaving the front surface per ion. | `sputtered` |
| **Backscattered / transmitted** | Ions leaving through the front / back surface. | `leave_front`, `leave_back` |
| **Monolayer (l₀)** | N^(-1/3), the mean interatomic distance; TRIM's recoil step and minimum flight path. | `LayerModel::monolayer` |
| **Free flight path (L)** | Distance between ion collisions: L = 1/(π p_max² N). | `ion_flight_path_2013` |
| **Impact parameter (p), p_max** | Distance of closest approach of the straight-line paths; p = √(R/(πLN)) for ions, √R·p₀ for recoils. | |
| **Screening length (a)** | ZBL universal: 0.8853 a₀ / (Z₁^0.23 + Z₂^0.23). | `screening_length` |
| **Reduced energy ε** | ε = a M₂ E / (Z₁ Z₂ e² (M₁ + M₂)), dimensionless. | `reduced_energy` |
| **Magic formula** | Biersack–Haggmark closed-form approximation of the scattering angle (TRIM's Tmagic). | `magic_sin2_half_theta` |
| **ZBL universal potential** | Φ(r) = 0.18175e^(−3.1998r) + 0.50986e^(−0.94229r) + 0.28022e^(−0.4029r) + 0.02817e^(−0.20162r). TRIM uses it for every pair. | `ScreeningFunction` |
| **Electronic stopping S_e** | Energy loss to electrons per path length (eV/Å) or per areal density (eV/(10¹⁵ atoms/cm²)). | `StoppingCurve`, `stopping_cross_section_quantity` |
| **Straggling (Ω)** | Spread of the electronic loss: Bohr straggling at high energy; TRIM adds a uniform deviate (R − ½)·Ω·√L clipped to ±S·L. | `electronic_loss_2013` |
| **Bragg's rule, compound correction** | Stopping of a compound = sum of its elements' stopping, times a per-layer correction factor (TRIM.IN's Bragg column). | `compound_correction` |
| **Gas / solid phase** | SRIM has separate stopping for gas targets; TRIM uses the gas tables for every element of a gas layer. | `TargetLayer::gas` |
| **Depth window, bins** | The depth range and number of bins (TRIM: 100) of the tabulated profiles; totals include events outside the window. | `window_start/end`, `depth_bins` |
| **TRIM bugs** | Behaviours of TRIM.exe that are wrong (energy created, double counts, …); reproduced by default, switchable. | `TrimBugs` |
| **Unit (validation)** | All TRIM runs of one input (seed replicates `_s<k>`) compared with AMCSET's seed runs of that input. | `analyze.py` |
| **z** | Normal deviate of the two-sample t statistic of AMCSET − TRIM; |z| < 2 in ~95 % and < 3 in ~99.7 % of units if both are the same random process. | `analyze.py` |
| **Determinism harness** | Program printing every collision of fixed-seed native histories, to prove refactors bit-identical. | `tools/determinism_harness` |
