# 9. Conclusions

[Review contents](README.md) · previous: [8. Ionization by slow recoils](08-ionization.md)

## 9.1 What the review shows

**The code.** AMCSET's SRIM mode is a statistically faithful reproduction of
SRIM-2013's TRIM.exe: 1031 of 1032 integral quantities and the depth
profiles of 75 inputs agree within the statistics, as do 11 runs of Ni, Co
and Cu ions and published SRIM damage calculations (Chapter 1). Its bundled
stopping tables are SRIM-2013's (Chapter 2). Its damage models reproduce
published SRIM post-processing to 0.43 % (Chapter 6). The one open
difference from TRIM.exe is the replacements in Xe gas (+4.6 %).

**The physics it inherits.** Everything else is a test of SRIM-2013's
physics, since AMCSET reproduces it. Here the picture is the familiar one,
now with numbers:

| Quantity | Where it works | Where it does not | Chapter |
|---|---|---|---|
| Electronic stopping | median 0.996 over 38 268 measurements; H and He within 3 % of NIST above 0.1 MeV/u | slow ions below ~0.003 MeV/u (high); heavy ions at 0.025-0.25 MeV/u (often 10-15 % low); Ni ions (a dip of up to 39 %) | 2 |
| Ranges | light ions in Si within 1-7 %; heavy ions in Si and SiO2 on the median | Ni ions 1.2-1.9 × too deep; heavy ions in C, B, SiC 12-25 % too short; single heavy ions in Si ±10-22 % | 3, 4 |
| Sputtering | the shape of every yield curve; the median of all pairs (1.05) | each pair has its own factor (0.68-1.9); no yield near light-ion thresholds | 5 |
| Damage | NRT and Kinchin-Pease as defined; arc-dpa as fitted | full-cascade vacancies overestimate surviving defects 4-19 × | 6 |
| Ionization of keV recoils | — | Ge 25-32 % low, Si up to 3 × high | 8 |
| Channeling (AMCSET's own) | a channeling dip; B profiles at low dose | channeling tails of slow heavy ions near amorphization 9-20 × high | 7 |

## 9.2 The corrections that can be switched on

| Deviation | Switch | Default | Evidence |
|---|---|---|---|
| SRIM-2013's Ni-ion stopping | `physics.stopping_corrections.nickel_ions` (`--correct-nickel`, GUI Calculation tab) | off | measured Ni peaks at 3.5-9 MeV within 9 % (six of ten); Zhang 2002's stopping data |
| TRIM.exe's own bugs | 11 switches in `physics.trim_bugs` (`docs/SRIM_MODE.md` §7) | TRIM's behaviour | disassembly of TRIM.exe |
| Sputtering level of a pair | the surface binding energy `Es` per element | TRIM's table | a measured yield of that pair |
| Surviving defects | arc-dpa and rpa in the `dpa` table | reported alongside NRT | MD fits (Nordlund et al. 2018) |
| Native oxide on Si | an amorphous SiO2 top layer | not added | Chapter 7: it matters for slow ions |

The defaults reproduce TRIM.exe, so that results agree with the
literature's SRIM calculations. The corrections are opt-in and documented
with their evidence.

## 9.3 Recommendations for users

1. **Keep the SRIM-2013 tables** (the default). The built-in stopping is
   good only above about 2.5 MeV/u.
2. **For Ni-ion irradiations**, switch on `nickel_ions` and say so. Dose
   profiles from SRIM-2013 for Ni ions of 3.5-9 MeV put the implanted Ni
   20-42 % too deep (and the damage with it).
3. **For heavy ions in light targets** (C, B, SiC, polymers), expect ranges
   12-25 % short and profiles too narrow; check one measured range if the
   depth matters.
4. **For damage**, report NRT as the dose and arc-dpa for surviving defects
   in metals; do not use full-cascade vacancies as a dose.
5. **For sputtering**, calibrate Es against one measured yield of the pair.
6. **For channeling** at doses near amorphization, treat tails as upper
   estimates, and include the native oxide.
7. **State the inputs** (Ed, Eb, Es, density, quick or full cascades, the
   correction switches), as Weber and Zhang (2019) recommend for SRIM.

## 9.4 Open

- Replacements in Xe gas (+4.6 % against TRIM.exe; Chapter 1).
- The 2 MeV Ni peaks, and 5 MeV in Ni and 11.5 MeV in Fe, which stay deep
  with the correction; TRIM-98's Ni stopping (Song et al.'s supplementary
  tables, not obtained) or a fitted, energy-dependent correction would be
  the next test. More TRIM.exe seeds would settle the ~1 % range lean of
  5 MeV Ni in Fe.
- A model for heavy ions in light targets (Behar et al.'s correlated
  inelastic losses, or a revised electronic stopping).
- Dose-dependent dechanneling for the crystal transport.
- Measured reflection coefficients, sputtering at oblique incidence,
  compounds in the IAEA database, and the electron transport of the native
  engine, which is not validated.
- The full text of Gurbich et al. (2025).
