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

"""Run AMCSET on a directory of finished TRIM cases, several seeds each.

usage: batch.py CASES_DIR OUT_DIR --amcset EXE [options]

CASES_DIR holds one sub-directory per case with the TRIM.IN (and TRIM.DAT)
that TRIM ran and TRIM's output files. For every case whose TRIM run has
finished (RANGE.txt present) and every seed, AMCSET runs the same TRIM.IN with
TRIM's ion count into OUT_DIR/<case>/<seed>, keeping only the small
TRIM-format tables (the large per-ion lists are deleted), and compare.py's
result goes to OUT_DIR/cmp/<case>.<seed>.json. Finished runs are skipped, so
the script can be re-run as more TRIM cases finish. Uses only the standard
library (runs with the Windows Python too).
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402

# Per-ion lists and other large files that the analysis does not need.
DROP = ('RANGE_3D.txt', 'COLLISON.txt', 'EXYZ.txt', 'TRIMOUT.txt', 'BACKSCAT.txt', 'TRANSMIT.txt')


def ion_count(trim_in):
    with open(trim_in, errors='replace') as f:
        lines = f.read().splitlines()
    return int(float(lines[2].split()[4]))


def finished(case_dir):
    return os.path.exists(os.path.join(case_dir, 'RANGE.txt')) and \
        os.path.exists(os.path.join(case_dir, 'VACANCY.txt'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cases_dir')
    ap.add_argument('out_dir')
    ap.add_argument('--amcset', required=True, help='amcset executable')
    ap.add_argument('--seeds', default='101,102,103,104', help='comma-separated seeds')
    ap.add_argument('--threads', type=int, default=0, help='worker threads per run (0 = all cores)')
    ap.add_argument('--stopping-table')
    ap.add_argument('--srim-data')
    ap.add_argument('--only', help='regular expression on case names')
    ap.add_argument('--skip', help='regular expression on case names to skip')
    ap.add_argument('--extra', default='', help='extra amcset arguments (one string)')
    ap.add_argument('--replicate-seeds', default='',
                    help='TRIM seed replicates (cases named <id>_<name>_s<k> that differ only in the '
                         'seed) are run once, for the first finished one, with these seeds '
                         '(default: --seeds)')
    args = ap.parse_args()

    seeds = [int(s) for s in args.seeds.split(',') if s]
    rep_seeds = [int(s) for s in args.replicate_seeds.split(',') if s] or seeds
    groups_done = set()
    cmp_dir = os.path.join(args.out_dir, 'cmp')
    os.makedirs(cmp_dir, exist_ok=True)
    progress = os.path.join(args.out_dir, 'progress.txt')
    cases = sorted(d for d in os.listdir(args.cases_dir)
                   if os.path.exists(os.path.join(args.cases_dir, d, 'TRIM.IN')))
    for case in cases:
        if args.only and not re.search(args.only, case):
            continue
        if args.skip and re.search(args.skip, case):
            continue
        cdir = os.path.join(args.cases_dir, case)
        if not finished(cdir):
            continue
        n = ion_count(os.path.join(cdir, 'TRIM.IN'))
        case_seeds = seeds
        m = re.match(r'^\d+_(.*)_s\d+$', case)
        if m:
            if m.group(1) in groups_done:
                continue
            groups_done.add(m.group(1))
            case_seeds = rep_seeds
        for seed in case_seeds:
            out_json = os.path.join(cmp_dir, '%s.%d.json' % (case, seed))
            if os.path.exists(out_json):
                continue
            odir = os.path.join(args.out_dir, case, str(seed))
            cmd = [args.amcset, 'trim', os.path.join(cdir, 'TRIM.IN'), '-n', str(n), '-s', str(seed),
                   '-q', '--formats', 'trim', '-o', odir]
            if args.threads:
                cmd += ['-j', str(args.threads)]
            if args.stopping_table:
                cmd += ['--stopping-table', args.stopping_table]
            if args.srim_data:
                cmd += ['--srim-data', args.srim_data]
            if args.extra:
                cmd += args.extra.split()
            t0 = time.time()
            r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace')
            dt = time.time() - t0
            if r.returncode != 0:
                with open(progress, 'a') as f:
                    f.write('%s %d FAILED (%d): %s\n' % (case, seed, r.returncode, r.stdout.strip()[-500:]))
                break
            for name in DROP:
                p = os.path.join(odir, name)
                if os.path.exists(p):
                    os.remove(p)
            try:
                rows = compare.collect(cdir, odir)
                data = {'trim_dir': cdir, 'amcset_dir': odir, 'seed': seed, 'seconds': dt,
                        'rows': [{'quantity': q, 'trim': t, 'amcset': a,
                                  'difference_percent': (100.0 * (a - t) / abs(t)) if (t and a is not None) else None}
                                 for q, t, a, _ in rows]}
            except Exception as e:  # keep going; record the problem
                data = {'error': repr(e), 'seed': seed}
            with open(out_json, 'w') as f:
                json.dump(data, f, indent=1)
            with open(progress, 'a') as f:
                f.write('%s %d %.1f s\n' % (case, seed, dt))
    with open(progress, 'a') as f:
        f.write('PASS DONE %s\n' % time.strftime('%Y-%m-%d %H:%M:%S'))


if __name__ == '__main__':
    main()
