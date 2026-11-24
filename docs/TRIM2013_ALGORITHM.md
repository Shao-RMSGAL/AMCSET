# How SRIM-2013's TRIM.exe transports ions

This note records what TRIM.exe (SRIM-2013.00, 1,138,688 bytes) actually does,
as established by disassembling it, and how AMCSET's `trim2013` transport
(`physics.transport: "trim2013"`, the default) reproduces it. It describes the
algorithm and its constants only; no TRIM code is reproduced. Addresses refer
to that executable and are given so that each statement can be checked.

TRIM.exe is Visual Basic 5 compiled to native x86 code. Most numbers are
`Variant`s, so the arithmetic appears as calls into `MSVBVM50.DLL`
(`__vbaVarMul(res, rhs, lhs)` = lhs * rhs, `__vbaVarTstGt(rhs, lhs)` =
lhs > rhs, `__vbaVarMove` with ecx = destination and edx = source). The
routines below were read from disassembly turned into pseudo-code.

## Program structure

| Routine | Address | Role |
| --- | --- | --- |
| `TRIMMain.Main` | 0x433fd0 | Ion loop (IONCYCLE), primary recoil (STARTRECOIL), recoil loop (RECOILLOOP, NEXTRECOIL) |
| TRIMSTART setup | 0x44f8c0 | Layer averages, screening lengths, straggling constants |
| Ion free flight path | 0x44ab60 | Flight length of the ion |
| Layer interface / window edge | 0x44b330 / 0x44b7a0 | Shortening of flights near boundaries |
| `SELOSS` | 0x45d4e0 | Electronic energy loss and straggling of one flight |
| `Tmagic` | 0x467600 | Scattering angle (magic formula) and energy transfer |
| `Random` | 0x469050 | L'Ecuyer combined generator (m = 2147483563, 2147483399) |
| Target atom choice | 0x45f7c0 | One random number minus the stoichiometric fractions |
| `DIRCOS` | 0x45ce40 | New direction after a deflection (with a separate branch for motion along x) |
| Quick (Kinchin-Pease) damage | 0x44bae0 | Vacancies from one primary recoil |
| `SPUTTER` | 0x441ff0 | A recoil reaching the front surface |
| `RECSTOPDATA` | 0x457970 | A recoil that stops |
| Displacement / replacement / discard | 0x458b90 / 0x459140 / 0x459b70 | Recoil bookkeeping |
| Ion energy tallies | 0x44d4e0 | Ionization, E2RECOIL and ion phonons over a flight |

Energies are in eV, lengths in Angstrom. `L0 = N^(-1/3)` is the monolayer
spacing of a layer with atomic density N.

## Constants and layer averages (0x44f8c0)

* Layer mean mass M2 = sum c_j M_j (stoichiometric fractions c_j).
* Layer atomic number: an integer, accumulated as z = CInt(c_j Z_j + z) over
  the elements (VB rounding at every step).
* Screening lengths: a = 0.5292 x 0.8853 / (Z1^0.23 + Z2^0.23) for every pair
  and for the ion on the layer (with the integer layer Z). Reduced energy
  factors use e^2 = 14.41 eV A for pairs and 14.4 eV A for the layer.
* Only the ZBL universal potential is used. No executable of SRIM-2013 reads
  `POTCOEF03.DAT`.
* Robinson partition for the quick calculation, with the layer's integer Z
  and mean M: eps = 0.01 Z^(-7/3) T, k = 0.1334 Z^(2/3) M^(-1/2),
  g = eps + 0.4 eps^(3/4) + 3.4 eps^(1/6).
* Stop energy of ions and recoils: max(1 eV, min of Ed, Es, Eb), but the
  minimum starts from an uninitialised (zero) variant, so it is always 1 eV.

## Ion flight (0x44ab60)

1. First flight from the surface (X < L0 of layer 1 and E = E0): one
   monolayer, and the next two flights are randomised (below).
2. Surface-sputtering mode (cascade code 3): every flight is one monolayer.
3. Otherwise TRIM-85's rule with the layer averages:
   eps = E FI, EPSDG = 5 eV x FI (1 + M1/M2)^2 / (4 M1/M2),
   EEG = sqrt(eps EPSDG), p_max = a / (EEG + sqrt(EEG) + 0.125 EEG^0.1),
   L = 1 / (pi p_max^2 N). There is no limit on the energy lost over L and no
   lower limit other than those below.
4. If the flight would cross into another layer, leave the depth window, or
   the ion is within p_max of the front surface moving outwards, the flight
   is multiplied by random numbers until it fits (0x44b330, 0x44b7a0); the
   next two flights are then randomised.
5. Full-cascade mode: every flight is max(L0, R L).
6. Randomised flights (quick mode): max(L0, R L).

## Electronic loss of a flight (0x45d4e0)

DEE = S(E) L at the energy before the flight, plus a straggling deviate
(R - 0.5) Omega sqrt(L), limited to +-S(E) L, with DEE <= E. Omega is the
Bohr value times sqrt(12) (177 Z1 sqrt(Z2 N), so the uniform deviate has
Bohr's variance), with a relativistic factor (1/beta^2 - 1/2)/(1/beta^2 - 1),
or k_low E when that is smaller (k_low = sqrt(38.4 N (Z1 + Z2)^(8/3) /
(1000 M1))). Recoils use the ion's constants.

## Ion collision (IONCYCLE)

1. E -= DEE, then E is limited to [stop energy, E0].
2. The ion moves L; X < 0 is backscattering, X beyond the target transmission.
3. Target atom at the ion's position, chosen by stoichiometry; impact
   parameter p = sqrt(R / (pi L N)).
4. `Tmagic` at the energy after the electronic loss: magic formula for
   eps <= 10, 1/(1 + (1 + b(1 + b)) (2 eps b)^2) for eps > 10.
   The distance of closest approach starts from R = -2.7 ln(eps b) (taken
   again from that R, and kept only if larger than b; otherwise R = b) and
   takes at most 14 Newton steps, stopping when |dR/R| < 0.0003; R is
   limited to 200 screening lengths before each step. AMCSET's
   closest_approach converges fully (Newton with a bisection guard); for
   the impact parameters TRIM uses the two agree.
5. Tallies over the whole flight (0x44d4e0): ionization; E2RECOIL if T >= Ed,
   ion phonons if T < Ed.
6. T below the stop energy: also booked as recoil phonons (counted twice).
7. Quick mode: Kinchin-Pease with the layer constants; one vacancy if the
   damage energy is at least Ed, 0.8 E_d/(2 Ed) above 2.5 Ed.
8. Full cascades (STARTRECOIL): the primary takes recoil slot 1 at the ion's
   position. Its direction is turned from the ion's direction after the
   collision. If T > Ed it is displaced with energy T - Eb; otherwise the slot
   keeps the energy it had at the end of the previous cascade.

## Recoils (RECOILLOOP, NEXTRECOIL)

* The slot arrays (energy, position, direction cosines, sin x) are VB
  `Single`s, so recoil quantities are rounded to single precision (AMCSET
  keeps double precision; the difference is far below the statistics in
  solids; it is not yet excluded for the 10^6 A cascades of heavy gases,
  see SRIM_MODE.md §6.2).
* The stack holds at most the "Maximum number of secondary recoils" of
  Data/TRIM.cfg (200 000 in the installation used, 20 000 by default)
  minus one: a push beyond it leaves the stack pointer at the last slot,
  so the mover in the top slot is dropped, and a message box ("Minor TRIM
  Error: Cascade > N Atoms") is shown once. AMCSET's stack grows without
  limit.
* Recoils move one monolayer per step, with p = sqrt(R) L0/sqrt(pi). Light
  atoms (Z < 21) within 2 L0 of the surface get smaller impact parameters
  (p x R for Z < 14, an interpolation up to Z = 20).
* At the front surface (0x441ff0) an atom with E cos^2 >= Es is sputtered;
  otherwise its direction is reflected (cx -> -cx) and it continues at
  x = 1e-4 A. The sputtered energy tally is E - Es, but SPUTTER.txt lists E
  (the file says "SBE is NOT subtracted") and the direction cosines inside
  the surface, limited to +-0.9999999.
* Electronic loss over L0 at the energy before the collision; the collision
  also uses that energy, and the mover keeps max(0, E - T - DEE) (0x436de5).
  T + DEE can exceed E for slow atoms; TRIM still books all of DEE as
  ionization and all of T as the struck atom's energy, so energy is created
  there (about 3 % of the incident energy for Fe -> Fe at 1 keV).
* Directions after a recoil collision (0x43735e, 0x4377d5): the mover is
  turned by DIRCOS from its own direction. The struck atom is turned by a
  second DIRCOS call whose sin x argument is the constant 0, which selects
  DIRCOS's branch for a direction along the x axis: the struck atom leaves
  along (cx cos t, -sin t cos phi, -sin t sin phi), normalised, where cx is
  the mover's x cosine before the collision, t the recoil angle and phi the
  azimuth. Only the mover's x cosine enters: the struck atom's transverse
  direction is random, and its x cosine has the sign of the mover's. (The
  primary recoil of an ion collision is turned correctly, from the ion's
  new direction with its own sin x.) It matters most for sputtering: with
  the correct rotation, Ar -> Cu at 1 keV gives 2.95 sputtered atoms per
  ion instead of TRIM's 3.55, and the sputtered atoms are too fast and too
  oblique.
* A mover that strikes an atom moves up one slot and continues; the struck
  atom takes its old slot. When the mover stops (E < 1 eV), the next slot
  down is processed (one random number R is drawn for each, 0x437f62):
  * below Ed and deeper than 9 L0 R: not followed (its energy becomes
    phonons);
  * at or above Ed: displaced (a vacancy), losing Eb only if deeper than
    2 L0;
  * below Ed but near the surface: set in motion without being displaced;
    if it is sputtered TRIM adds a displacement and a vacancy to its totals
    only (VACANCY.txt's header and COLLISON.txt, not the depth tables;
    COLLISON.txt does not list the atom).
* A displaced mover that stops right after striking a displaced atom of its
  own element makes a replacement collision (the vacancy is removed). TRIM
  still counts the mover as an interstitial, so its interstitials exceed its
  vacancies by about the number of replacements.

## Ion energy file (EXYZ.txt, 0x43460b)

* At the top of each ion step: a record when the energy equals the TRIM.IN
  energy, after which the next record is due when Int(E / interval) falls
  below Int(E / interval) - 1 (0x434745), so the first crossing after the
  start is not recorded; otherwise a record whenever Int(E / interval) is
  below the last recorded value, which then becomes Int(E / interval).
  TRIM.DAT ions start with Int(E / interval) + 1 (0x449c10).

## Random numbers

* One L'Ecuyer generator for everything. Each ion first draws two numbers
  for its starting point, y = L0 cy R and z = L0 cz R (0x4344a2, 0x43453c;
  zero for normal incidence).
* Each flight, of the ion or of a recoil, draws in this order: the flight
  randomisation where it applies (ion only), the straggling deviate in
  SELOSS (always, also with no straggling to apply), the target atom
  (always, also for a pure element), the impact parameter, the light-atom
  numbers where they apply (recoils only), and the azimuth in Tmagic. Every
  atom taken up in NEXTRECOIL draws one more number.
* Cascades (full-cascade and surface-sputtering modes): the generator state
  is saved before STARTRECOIL and restored when the cascade is finished
  (0x434e89, 0x438140; not for codes 4-7). The cascade therefore draws the
  numbers the ion draws next, and each cascade replays the numbers of the
  previous one, shifted by the ion's draws in between. This does not change
  the ion, but it correlates consecutive cascades, and through the energy a
  primary inherits from the previous cascade (STARTRECOIL above) it shifts
  TRIM's recoil energy tallies and near-surface damage slightly. AMCSET
  reproduces it by running each cascade on a copy of the ion's generator
  and drawing the same numbers in the same order. In the
  surface-sputtering mode, where every ion collision starts a cascade, this
  removed a +1 % excess in recoil ionization and phonons and a +0.7 % excess
  in vacancies.

## Bugs of TRIM.exe

Several of the behaviours above are bugs rather than modelling choices: the
stale primary energy (STARTRECOIL), the energy created in slow recoil
collisions (0x436de5), the struck-atom direction (0x4377d5), interstitials
counted for replacements (0x457970), sub-threshold sputtered atoms counted
as vacancies (0x441ff0), ion transfers below the stop energy counted twice
(0x4575e0), the replayed random numbers of cascades (0x434e89/0x438140), the
zero-energy defaults that skip the first atom (0x44efc0), the rounded layer
Z (0x44f8c0), EXYZ.txt's skipped first step (0x434745) and the doubled
TRIM.DAT recoil events. AMCSET reproduces each by default and can correct
each (`physics.trim_bugs`); `docs/SRIM_MODE.md` §7 gives the corrections
and their effect.

## Effect on the comparison with TRIM

| Case | Refined transport | trim2013 transport |
| --- | --- | --- |
| H -> W 5 keV, range / vacancies / backscattering | -22.6 % / -27 % / 14.5 vs 20.5 % | -0.4 % / -0.9 % / 20.2 vs 20.5 % |
| H -> Fe 1 MeV, quick vacancies | -6.3 % | -1.4 % |
| Fe -> Fe 1 keV, quick range / vacancies | +3.8 % / +2.4 % | +0.1 % / -0.1 % |
| Fe -> Fe 1 keV, sputtering yield | 0.25 vs 1.98 | 1.95 vs 1.98 |
| Fe -> Fe 10 keV, sputtering yield | 1.27 vs 5.10 | 4.92 vs 5.10 |
| Ar -> Cu 1 keV, sputtering yield | 0.37 vs 3.55 | 3.59 vs 3.55 |
| Fe -> Fe 1 keV, recoil phonons (% of E0) | - | 72.0 vs 72.0 |

The sputtering yields came in two steps: the trim2013 flight and
bookkeeping rules (0.25 -> 2.00 at 1 keV, 1.27 -> 4.75 at 10 keV, Ar -> Cu
0.37 -> 2.95), then the struck-atom rotation above (Ar -> Cu 2.95 -> 3.54,
Fe 10 keV 4.75 -> 4.9). The energy spectrum and angular distribution of the
sputtered atoms then agree with TRIM's as well (Ar -> Cu: median energy
inside the surface 9.92 vs 9.94 eV, mean |cos x| 0.823 vs 0.821). The
max(0, ...) limit on the mover's energy brought the recoil phonons near the
surface to TRIM's values (Fe 1 keV: 68.9 -> 71.9 % vs 72.0 %).

Values for trim2013 are means of five runs with the number of ions of the
TRIM run. The full validation table, with statistical significance, is in
`docs/SRIM_MODE.md` §6.
