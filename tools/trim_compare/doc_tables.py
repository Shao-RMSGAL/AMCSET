# Copyright 2025, Texas A&M University
#
# This file is part of AMCSET.
#
# AMCSET is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# AMCSET. If not, see <https://www.gnu.org/licenses/>.

"""Markdown tables for docs/SRIM_MODE.md section 6.2 from analyze.py's JSON.

usage: doc_tables.py ANALYSIS_JSON [--units]

Prints the headline numbers (units, TRIM runs, AMCSET runs, comparisons
within 3 sigma, RMS of z), the summary table (one row per quantity: units,
|z| < 2, |z| < 3, chi^2/n of z, largest |z|) and, with --units, the table of
AMCSET minus TRIM in % (z) per unit. Paste them into docs/SRIM_MODE.md after
a new comparison run (tools/windows/compare_with_trim.bat writes
analysis.json).
"""

import argparse
import json
import math

QUANTITIES = [
    ('mean projected range (A)', 'Mean range'),
    ('range straggling (A)', 'Range straggling'),
    ('lateral range (A)', 'Lateral range'),
    ('backscattered (%)', 'Backscattered fraction'),
    ('transmitted (%)', 'Transmitted fraction'),
    ('vacancies / ion', 'Vacancies (tables)'),
    ('vacancies / ion (header)', 'Vacancies (header)'),
    ('displacements / ion (header)', 'Displacements (header)'),
    ('replacements / ion', 'Replacements'),
    ('sputtering yield (atoms/ion)', 'Sputtering yield'),
    ('ionization by ion (%)', 'Ionization by ion'),
    ('ionization by recoils (%)', 'Ionization by recoils'),
    ('phonons by ion (%)', 'Phonons by ion'),
    ('phonons by recoils (%)', 'Phonons by recoils'),
    ('energy to recoils (%)', 'Energy to recoils'),
    ('mean depth of vacancies (A)', 'Mean depth of vacancies'),
]
UNIT_COLUMNS = [
    ('mean projected range (A)', 'Range'),
    ('vacancies / ion', 'Vacancies'),
    ('replacements / ion', 'Replacements'),
    ('sputtering yield (atoms/ion)', 'Sputtered'),
    ('ionization by ion (%)', 'Ion ionization'),
    ('phonons by recoils (%)', 'Recoil phonons'),
]
PROFILES = [('vacancies', 'Vacancy profile (χ² per bin)'), ('ion ranges', 'Ion range profile (χ² per bin)')]


def z_values(data, quantity):
    out = []
    for name, unit in data.items():
        q = unit.get('quantities', {}).get(quantity)
        if q and q.get('z') is not None and math.isfinite(q['z']):
            out.append((q['z'], name))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('analysis_json')
    ap.add_argument('--units', action='store_true', help='also print the per-unit table')
    args = ap.parse_args()
    with open(args.analysis_json, encoding='utf-8') as f:
        data = json.load(f)

    trim_runs = sum(len(u['trim_runs']) for u in data.values())
    amcset_runs = sum(u.get('amcset_seeds', 0) for u in data.values())
    every = [z for q, _ in QUANTITIES for z, _ in z_values(data, q)]
    within = sum(abs(z) < 3 for z in every)
    rms = math.sqrt(sum(z * z for z in every) / len(every)) if every else float('nan')
    print(f'{len(data)} units, {trim_runs} TRIM runs, {amcset_runs} AMCSET runs; '
          f'{within} of {len(every)} comparisons within |z| < 3; RMS of z {rms:.2f}\n')

    print('| Quantity | Units | \\|z\\| < 2 | \\|z\\| < 3 | χ²/n of z | Largest \\|z\\| |')
    print('|---|---|---|---|---|---|')
    for q, label in QUANTITIES:
        zs = z_values(data, q)
        if not zs:
            continue
        worst = max(zs, key=lambda t: abs(t[0]))
        print(f'| {label} | {len(zs)} | {sum(abs(z) < 2 for z, _ in zs)} | {sum(abs(z) < 3 for z, _ in zs)} '
              f'| {sum(z * z for z, _ in zs) / len(zs):.2f} | {abs(worst[0]):.1f} (`{worst[1]}`) |')
    for key, label in PROFILES:
        rows = [(name, u['profiles'][key]) for name, u in data.items()
                if (u.get('profiles') or {}).get(key) and u['profiles'][key].get('ndf')]
        if not rows:
            continue
        ratios = sorted(p['chi2_per_ndf'] for _, p in rows)
        worst = max(rows, key=lambda t: t[1]['chi2_per_ndf'])
        # chi^2 is one-sided: a small chi^2 is not a disagreement.
        print(f'| {label} | {len(rows)} | {sum(p["z_equiv"] < 2 for _, p in rows)} '
              f'| {sum(p["z_equiv"] < 3 for _, p in rows)} | median χ²/ndf {ratios[len(ratios) // 2]:.2f} '
              f'| χ²/ndf {worst[1]["chi2_per_ndf"]:.2f} (`{worst[0]}`) |')

    if args.units:
        print()
        print('| Case | TRIM runs | AMCSET seeds | ' + ' | '.join(c + ' Δ % (z)' for _, c in UNIT_COLUMNS) + ' |')
        print('|---|---|---|' + '---|' * len(UNIT_COLUMNS))
        for name in sorted(data):
            u = data[name]
            cells = []
            for q, _ in UNIT_COLUMNS:
                v = u.get('quantities', {}).get(q)
                if not v or v.get('z') is None or not v.get('trim') or not math.isfinite(v['z']) or abs(v['z']) >= 39:
                    cells.append('–')
                else:
                    cells.append(f"{v['diff_percent']:+.1f} ({v['z']:+.1f})")
            print(f"| `{name}` | {len(u['trim_runs'])} | {u.get('amcset_seeds', 0)} | " + ' | '.join(cells) + ' |')


if __name__ == '__main__':
    main()
