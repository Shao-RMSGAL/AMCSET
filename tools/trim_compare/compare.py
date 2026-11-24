#!/usr/bin/env python3
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

"""Compare an AMCSET SRIM-mode run with a TRIM reference run.

Usage:
    compare.py TRIM_DIR AMCSET_DIR [--json] [--tolerance PCT]

TRIM_DIR is a directory with TRIM's output files (RANGE.txt, VACANCY.txt,
IONIZ.txt, PHONON.txt, E2RECOIL.txt, TDATA.txt, ...) and the TRIM.IN that
produced them. AMCSET_DIR is the output directory of
``amcset trim TRIM_DIR/TRIM.IN -o AMCSET_DIR`` (it must contain
summary.json).

The script prints one table of quantities (TRIM, AMCSET, relative
difference). Vacancies and replacements are integrated from the depth
tables, because TRIM rounds the totals in VACANCY.txt to integers.
With --tolerance, the exit status is 1 if any of the headline quantities
(mean range, straggling, vacancies) differ by more than PCT percent.
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trimout  # noqa: E402


def integrated_vacancies(case):
    """(vacancies per ion, replacements per ion or None) from the tables."""
    vac = case['vacancy']
    depth = vac['depth']
    cols = vac['columns']
    if not depth or not cols:
        return None, None
    if vac.get('sputtering_mode'):
        # Cascade code 3: "by IONS" are knock-ons, "by RECOILS" the vacancies.
        total = sum(trimout.integrate(depth, c) for c in cols[1:])
    elif vac['kp_mode']:
        # Kinchin-Pease layout: depth, vacancies by ions, vacancies by recoils.
        total = sum(trimout.integrate(depth, c) for c in cols)
    else:
        # Full-cascade layout: depth, knock-ons, vacancies per element.
        total = sum(trimout.integrate(depth, c) for c in cols[1:])
    repl = None
    novac = case.get('novac')
    if novac and len(novac) > 1:
        repl = sum(trimout.integrate(novac[0], c) for c in novac[1:])
        # TRIM writes NOVAC.txt with six decimals: in very wide bins (gas
        # targets) every value rounds to zero although the header counts
        # replacements; the table then says nothing.
        header = vac.get('total_replacements_per_ion')
        if repl == 0 and header and header >= 0.5:
            repl = None
    return total, repl


def energy_percent(case, column, e0_ev):
    table = case.get(column)  # list of columns; the first is the depth
    if not table or len(table) < 2:
        return None
    return [100.0 * trimout.integrate(table[0], c) / e0_ev for c in table[1:]]


def sputtered_count(trim_dir):
    """Number of sputtered atoms listed in TRIM's SPUTTER.txt (0 if none)."""
    path = os.path.join(trim_dir, 'SPUTTER.txt')
    if not os.path.exists(path):
        return 0
    with open(path, errors='replace') as f:
        return sum(1 for line in f if line.startswith('S ') and len(line.split()) >= 9)


def collect(trim_dir, amcset_dir):
    case = trimout.read_case(trim_dir)
    with open(os.path.join(amcset_dir, 'summary.json')) as f:
        a = json.load(f)
    e0 = a['ion']['energy_keV'] * 1000.0
    r = case['range']
    vac, repl = integrated_vacancies(case)
    rows = []

    def add(name, trim_value, amcset_value, headline=False):
        rows.append((name, trim_value, amcset_value, headline))

    add('mean projected range (A)', r['range'], a['range']['mean_A'], True)
    add('range straggling (A)', r['straggling'], a['range']['straggling_A'], True)
    add('lateral range (A)', r['lateral_range'], a['range']['lateral_A'])
    add('radial range (A)', r['radial_range'], a['range']['radial_A'])
    add('skewness', r['skewness'], a['range']['skewness'])
    add('kurtosis', r['kurtosis'], a['range']['kurtosis'])
    ions = r['ions'] or 1
    # RANGE.txt does not always carry the counts; TDATA.txt does.
    td = case.get('tdata', {})
    def count(key):
        v = td.get(key)
        return v if v is not None else (r[key] or 0)
    add('backscattered (%)', 100.0 * count('backscattered') / ions,
        100.0 * a['backscattered']['fraction'])
    add('transmitted (%)', 100.0 * count('transmitted') / ions,
        100.0 * a['transmitted']['fraction'])
    # Both sides integrate their VACANCY.txt (TRIM's header totals also count
    # sputtered atoms that had never been displaced).
    a_vac = a['per_ion']['vacancies']
    a_tab = os.path.join(amcset_dir, 'VACANCY.txt')
    if os.path.exists(a_tab):
        a_vac = integrated_vacancies({'vacancy': trimout.read_vacancy(a_tab)})[0]
    add('vacancies / ion', vac, a_vac, True)
    # The totals in VACANCY.txt's header (TRIM rounds them to integers; its
    # depth table does not always add up to them).
    if os.path.exists(a_tab):
        a_hdr = trimout.read_vacancy(a_tab)
        t_hdr = case['vacancy']
        if t_hdr.get('total_displacements_per_ion') and t_hdr['total_displacements_per_ion'] >= 100:
            add('vacancies / ion (header)', t_hdr['total_vacancies_per_ion'], a_hdr['total_vacancies_per_ion'])
            add('displacements / ion (header)', t_hdr['total_displacements_per_ion'],
                a_hdr['total_displacements_per_ion'])
    if repl is not None and not case['vacancy']['kp_mode']:
        add('replacements / ion', repl, a['per_ion']['replacements'])
    sput = sputtered_count(trim_dir)
    if sput:
        add('sputtering yield (atoms/ion)', sput / ions,
            sum(e.get('sputtering_yield', 0.0) for e in a.get('elements', [])))
    ioniz = energy_percent(case, 'ioniz', e0)
    phon = energy_percent(case, 'phonon', e0)
    e2r = energy_percent(case, 'e2recoil', e0)
    ep = a['energy_percent']
    if ioniz:
        add('ionization by ion (%)', ioniz[0], ep['ionization_ion'])
        if len(ioniz) > 1:
            add('ionization by recoils (%)', ioniz[1], ep['ionization_recoils'])
    if phon:
        add('phonons by ion (%)', phon[0], ep['phonons_ion'])
        if len(phon) > 1:
            add('phonons by recoils (%)', phon[1], ep['phonons_recoils'])
    if e2r:
        add('energy to recoils (%)', e2r[0], ep['to_recoils_from_ion'])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('trim_dir')
    ap.add_argument('amcset_dir')
    ap.add_argument('--json', action='store_true', help='print JSON instead of a table')
    ap.add_argument('--tolerance', type=float, default=None,
                    help='fail if a headline quantity differs by more than this many percent')
    args = ap.parse_args()
    rows = collect(args.trim_dir, args.amcset_dir)
    failed = []
    out = []
    for name, t, a, headline in rows:
        rel = None
        if t is not None and a is not None and t != 0 and math.isfinite(t):
            rel = 100.0 * (a - t) / abs(t)
        out.append({'quantity': name, 'trim': t, 'amcset': a, 'difference_percent': rel})
        if args.tolerance is not None and headline and rel is not None and abs(rel) > args.tolerance:
            failed.append(name)
    if args.json:
        print(json.dumps({'case': os.path.basename(os.path.normpath(args.trim_dir)), 'rows': out,
                          'failed': failed}, indent=2))
    else:
        print('%-30s %14s %14s %9s' % ('quantity', 'TRIM', 'AMCSET', 'diff %'))
        for row in out:
            fmt = lambda v: '%14s' % '-' if v is None else '%14.5g' % v
            d = '%9s' % '-' if row['difference_percent'] is None else '%+9.2f' % row['difference_percent']
            print('%-30s %s %s %s' % (row['quantity'], fmt(row['trim']), fmt(row['amcset']), d))
        if failed:
            print('outside tolerance: ' + ', '.join(failed))
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
