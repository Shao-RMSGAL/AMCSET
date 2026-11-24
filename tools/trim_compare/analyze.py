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

"""Statistical comparison of AMCSET seed runs with TRIM runs.

usage: analyze.py AMCSET_OUT TRIM_ROOT [TRIM_ROOT ...] [--json FILE] [-v]

AMCSET_OUT is the output directory of batch.py (<case>/<seed>/ with AMCSET's
TRIM-format tables and summary.json). TRIM_ROOT directories hold TRIM case
directories. Cases are matched by name after dropping the leading number and
a "rep_" prefix; TRIM runs of the same input with different seeds (names
ending in _s<k>, or an older run of the identical TRIM.IN) form one group, so
TRIM's own run-to-run scatter enters the comparison.

For every quantity of compare.py, plus the mean depth of the vacancy profile,
the script reports the TRIM and AMCSET means, their difference and a
t-statistic, using the pooled run-to-run standard deviation (from both codes
when both have several runs, otherwise from AMCSET's seeds alone), converted
to the normal deviate z with the same tail probability. For the
vacancy and ion-range depth profiles it reports chi^2 per bin between the two
mean profiles, with a variance proportional to the bin value (scaled to the
scatter of AMCSET's seeds) as the error.
"""

import argparse
import json
import math
import os
import re
import sys
from statistics import NormalDist

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import compare  # noqa: E402
import trimout  # noqa: E402


def unit_name(case):
    """Case directory name -> comparison unit (group) name."""
    name = re.sub(r'^\d+_', '', case)
    name = re.sub(r'^rep_', '', name)
    return re.sub(r'_s\d+$', '', name)


def input_signature(trim_in):
    """TRIM.IN lines without the seed line and the target name line."""
    with open(trim_in, errors='replace') as f:
        lines = [l.strip() for l in f.read().splitlines()]
    sig = []
    for i, l in enumerate(lines):
        if i == 4:      # cascades, seed, reminders
            sig.append(' '.join(l.split()[:1]))
        elif i == 8:    # "target name"  elements layers
            sig.append(' '.join(l.split()[-2:]))
        elif re.match(r'^\s*\d+\s+".*"', l):  # layer line: drop the layer name
            sig.append(re.sub(r'".*?"', '""', l))
        else:
            sig.append(' '.join(l.split()))
    return sig


def _betacf(a, b, x):
    """Continued fraction of the incomplete beta function (Numerical Recipes)."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        d = 1.0 / (d if abs(d) > tiny else tiny)
        c = 1.0 + aa / c
        c = c if abs(c) > tiny else tiny
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-12:
            break
    return h


def _betai(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1.0 - x)
    if x < (a + 1.0) / (a + b + 2.0):
        return math.exp(lbt) * _betacf(a, b, x) / a
    return 1.0 - math.exp(lbt) * _betacf(b, a, 1.0 - x) / b


def t_to_z(t, dof):
    """Normal deviate with the same two-sided tail probability as Student's t."""
    if not math.isfinite(t) or dof <= 0:
        return t
    p = _betai(0.5 * dof, 0.5, dof / (dof + t * t))  # two-sided tail probability
    if p <= 0.0:
        return math.copysign(40.0, t)
    if p >= 1.0:
        return 0.0
    # The lower tail keeps its precision for small p (1 - p/2 rounds to 1).
    return math.copysign(-NormalDist().inv_cdf(0.5 * p), t)


def mean_sd(v):
    n = len(v)
    if n == 0:
        return None, None
    m = sum(v) / n
    if n < 2:
        return m, None
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (n - 1))


def vacancy_profile(case_dir):
    v = trimout.read_vacancy(os.path.join(case_dir, 'VACANCY.txt'))
    depth, cols = v['depth'], v['columns']
    if not depth or not cols:
        return None, None
    if v.get('sputtering_mode') or not v['kp_mode']:
        cols = cols[1:]
    total = [sum(c[i] for c in cols) for i in range(len(depth))]
    return depth, total


def range_profile(case_dir):
    r = trimout.read_range(os.path.join(case_dir, 'RANGE.txt'))
    return r['depth'], r['ion_dist']


def profile_mean_depth(depth, values):
    if not depth or not values:
        return None
    w = sum(values)
    if w <= 0:
        return None
    return sum(d * v for d, v in zip(depth, values)) / w


def profile_chi2(t_profiles, a_profiles):
    """chi^2 between the mean TRIM and mean AMCSET profiles (same bins).

    The variance of one run in a bin is taken as k * (mean value), with one
    factor k for the whole profile (the median over bins of the AMCSET seeds'
    variance / mean): per-bin variances from a few seeds alone would inflate
    chi^2 by (n - 1) / (n - 3).
    """
    nt, na = len(t_profiles), len(a_profiles)
    if nt == 0 or na < 3:
        return None
    nb = min(min(len(p) for p in t_profiles), min(len(p) for p in a_profiles))
    ratios = []
    for b in range(nb):
        ma, sa = mean_sd([p[b] for p in a_profiles])
        # (Seeds that agree to rounding, e.g. a deterministic quick
        # calculation of a fixed particle list, give no scale.)
        if ma > 0 and sa > 1e-6 * ma:
            ratios.append(sa * sa / ma)
    if not ratios:
        return None
    k = sorted(ratios)[len(ratios) // 2]
    chi2, ndf, pulls = 0.0, 0, []
    for b in range(nb):
        mt = sum(p[b] for p in t_profiles) / nt
        ma = sum(p[b] for p in a_profiles) / na
        m = 0.5 * (mt + ma)
        if m <= 0:
            continue
        se = math.sqrt(k * m * (1.0 / na + 1.0 / nt))
        z = (ma - mt) / se
        chi2 += z * z
        ndf += 1
        pulls.append((abs(z), b, z, mt, ma))
    if ndf == 0:
        return None
    pulls.sort(reverse=True)
    # Wilson-Hilferty normal approximation of the chi^2 tail.
    zwh = ((chi2 / ndf) ** (1.0 / 3) - (1 - 2.0 / (9 * ndf))) / math.sqrt(2.0 / (9 * ndf))
    return {'chi2': chi2, 'ndf': ndf, 'chi2_per_ndf': chi2 / ndf, 'z_equiv': zwh, 'k': k,
            'worst_bins': [{'bin': b, 'pull': z, 'trim': mt, 'amcset': ma} for _, b, z, mt, ma in pulls[:3]]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('amcset_out')
    ap.add_argument('trim_roots', nargs='+')
    ap.add_argument('--json', help='write the full analysis to this file')
    ap.add_argument('--only', help='regular expression on unit names')
    ap.add_argument('-v', '--verbose', action='store_true')
    args = ap.parse_args()

    # TRIM runs per unit.
    trim_runs = {}
    for root in args.trim_roots:
        for d in sorted(os.listdir(root)):
            p = os.path.join(root, d)
            if not os.path.exists(os.path.join(p, 'RANGE.txt')) or '_diag_' in d:
                continue
            trim_runs.setdefault(unit_name(d), []).append(p)
    # AMCSET runs per unit (the seed directories with a summary.json).
    amc_runs = {}
    for d in sorted(os.listdir(args.amcset_out)):
        p = os.path.join(args.amcset_out, d)
        if not os.path.isdir(p) or d == 'cmp':
            continue
        seeds = sorted(os.path.join(p, s) for s in os.listdir(p)
                       if os.path.exists(os.path.join(p, s, 'summary.json')))
        if seeds:
            u = unit_name(d)
            if u not in amc_runs:
                amc_runs[u] = (os.path.join(p), seeds)

    results = {}
    for unit, (amc_case, seeds) in sorted(amc_runs.items()):
        if args.only and not re.search(args.only, unit):
            continue
        tdirs = trim_runs.get(unit, [])
        if not tdirs:
            continue
        # Keep only TRIM runs of the same input as the one AMCSET ran.
        ref_sig = None
        for t in tdirs:
            if os.path.basename(t).split('_', 1)[1] == os.path.basename(amc_case).split('_', 1)[1]:
                ref_sig = input_signature(os.path.join(t, 'TRIM.IN'))
        if ref_sig is None:
            ref_sig = input_signature(os.path.join(tdirs[-1], 'TRIM.IN'))
        tdirs = [t for t in tdirs if input_signature(os.path.join(t, 'TRIM.IN')) == ref_sig]
        if not tdirs:
            continue
        t_vals, a_vals = {}, {}
        for t in tdirs:
            for name, tv, _, _ in compare.collect(t, seeds[0]):
                if tv is not None and math.isfinite(tv):
                    t_vals.setdefault(name, []).append(tv)
        for s in seeds:
            for name, _, av, _ in compare.collect(tdirs[0], s):
                if av is not None and math.isfinite(av):
                    a_vals.setdefault(name, []).append(av)
        t_vp = [vacancy_profile(t) for t in tdirs]
        a_vp = [vacancy_profile(s) for s in seeds]
        t_rp = [range_profile(t) for t in tdirs]
        a_rp = [range_profile(s) for s in seeds]
        name = 'mean depth of vacancies (A)'
        t_vals[name] = [x for x in (profile_mean_depth(*p) for p in t_vp) if x is not None]
        a_vals[name] = [x for x in (profile_mean_depth(*p) for p in a_vp) if x is not None]
        quantities = {}
        for q in t_vals:
            tv, av = t_vals[q], a_vals.get(q, [])
            if not tv or not av:
                continue
            mt, st = mean_sd(tv)
            ma, sa = mean_sd(av)
            nt, na = len(tv), len(av)
            if st is not None and sa is not None:
                sp = math.sqrt(((nt - 1) * st * st + (na - 1) * sa * sa) / (nt + na - 2))
                dof = nt + na - 2
            else:
                sp = sa
                dof = na - 1
            if sp is None:
                continue
            se = sp * math.sqrt(1.0 / nt + 1.0 / na)
            t = (ma - mt) / se if se > 0 else (0.0 if ma == mt else float('inf'))
            z = t_to_z(t, dof)
            quantities[q] = {'trim': mt, 'trim_sd': st, 'n_trim': nt, 'amcset': ma, 'amcset_sd': sa,
                             'n_amcset': na, 'diff_percent': 100.0 * (ma - mt) / abs(mt) if mt else None,
                             'se_percent': 100.0 * se / abs(mt) if mt else None, 't': t, 'dof': dof,
                             'z': z}
        profiles = {}
        if all(p[0] for p in t_vp + a_vp):
            profiles['vacancies'] = profile_chi2([p[1] for p in t_vp], [p[1] for p in a_vp])
            if profiles['vacancies']:
                for w in profiles['vacancies']['worst_bins']:
                    w['depth_A'] = t_vp[0][0][w['bin']]
        if all(p[0] for p in t_rp + a_rp):
            profiles['ion ranges'] = profile_chi2([p[1] for p in t_rp], [p[1] for p in a_rp])
            if profiles['ion ranges']:
                for w in profiles['ion ranges']['worst_bins']:
                    w['depth_A'] = t_rp[0][0][w['bin']]
        results[unit] = {'trim_runs': [os.path.basename(t) for t in tdirs], 'amcset_case': os.path.basename(amc_case),
                         'amcset_seeds': len(seeds), 'quantities': quantities, 'profiles': profiles}

    # Summary per quantity over all units.
    names = []
    for r in results.values():
        for q in r['quantities']:
            if q not in names:
                names.append(q)
    print('%-32s %4s %6s %6s %8s  %s' % ('quantity', 'n', '|z|<2', '|z|<3', 'chi2/n', 'largest |z|'))
    for q in names:
        zs = [(abs(r['quantities'][q]['z']), u) for u, r in results.items()
              if q in r['quantities'] and math.isfinite(r['quantities'][q]['z'])]
        if not zs:
            continue
        chi = sum(z * z for z, _ in zs) / len(zs)
        worst = max(zs)
        print('%-32s %4d %6d %6d %8.2f  %.1f (%s)' % (q, len(zs), sum(z < 2 for z, _ in zs), sum(z < 3 for z, _ in zs),
                                                    chi, worst[0], worst[1]))
    for kind in ('vacancies', 'ion ranges'):
        ps = [(r['profiles'][kind]['chi2_per_ndf'], u) for u, r in results.items()
              if r['profiles'].get(kind)]
        if ps:
            print('%-32s %4d  median chi2/ndf %.2f, largest %.2f (%s)' % (
                'profile: ' + kind, len(ps), sorted(p for p, _ in ps)[len(ps) // 2], max(ps)[0], max(ps)[1]))
            zs = [r['profiles'][kind]['z_equiv'] for r in results.values() if r['profiles'].get(kind)]
            print('%-32s %4d  z<2: %d, z<3: %d (chi2 tail as a normal deviate)' % (
                '', len(zs), sum(z < 2 for z in zs), sum(z < 3 for z in zs)))
    if args.verbose:
        for u, r in results.items():
            print('\n%s  (TRIM runs %d, AMCSET seeds %d)' % (u, len(r['trim_runs']), r['amcset_seeds']))
            for q, x in r['quantities'].items():
                flag = ' <--' if abs(x['z']) >= 3 else ''
                print('  %-32s %12.5g %12.5g %+7.2f%% +-%5.2f%%  z %+6.1f%s' % (
                    q, x['trim'], x['amcset'], x['diff_percent'] or 0, x['se_percent'] or 0, x['z'], flag))
            for kind, p in r['profiles'].items():
                if p:
                    print('  profile %-24s chi2/ndf %.2f (ndf %d, z %.1f)  worst: %s' % (
                        kind, p['chi2_per_ndf'], p['ndf'], p['z_equiv'],
                        ', '.join('%.0f A %+.1f' % (w['depth_A'], w['pull']) for w in p['worst_bins'])))
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(results, f, indent=1)


if __name__ == '__main__':
    main()
