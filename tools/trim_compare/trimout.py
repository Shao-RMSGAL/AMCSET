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

"""Readers for TRIM-format output files (written by SRIM-2013 TRIM or by AMCSET's SRIM mode).

Every reader returns plain Python data; numbers keep the units of the file
(depths in Angstrom, distributions per (Angstrom*ion) or (atoms/cm3)/(atoms/cm2)).
"""

import math
import os
import re

_NUM = re.compile(r'^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$')


def _is_num(tok):
    return bool(_NUM.match(tok))


def _read(path):
    with open(path, encoding='latin-1') as f:
        return f.read().replace('\r', '').split('\n')


def numeric_table(path, min_cols=2):
    """All lines that consist only of numbers (at least min_cols); returns list of rows."""
    rows = []
    for line in _read(path):
        p = line.split()
        if len(p) >= min_cols and all(_is_num(t) for t in p):
            rows.append([float(t) for t in p])
    return rows


def _find(lines, pattern, conv=float, group=1):
    rx = re.compile(pattern)
    for line in lines:
        m = rx.search(line)
        if m:
            return conv(m.group(group))
    return None


def read_range(path):
    """RANGE.txt: summary statistics and the ion / recoil depth distributions."""
    lines = _read(path)
    num = r'([-+]?[\d.]+E[-+]?\d+|[-+]?[\d.]+)'
    out = {
        'ions': _find(lines, r'Total Ions calculated\s*=\s*' + num),
        'range': _find(lines, r'Ion Average Range\s*=\s*' + num),
        'straggling': _find(lines, r'Ion Average Range.*Straggling\s*=\s*' + num),
        'lateral_range': _find(lines, r'Ion Lateral Range\s*=\s*' + num),
        'lateral_straggling': _find(lines, r'Ion Lateral Range.*Straggling\s*=\s*' + num),
        'radial_range': _find(lines, r'Ion Radial\s+Range\s*=\s*' + num),
        'radial_straggling': _find(lines, r'Ion Radial\s+Range.*Straggling\s*=\s*' + num),
        'skewness': _find(lines, r'Range Skewne\S*\s*=\s*' + num),
        'kurtosis': _find(lines, r'Range Kurtosis\s*=\s*' + num),
        'backscattered': _find(lines, r'Backscattered Ions\s*=\s*(\d+)', int),
        'transmitted': _find(lines, r'Transmitted Ions\s*=\s*(\d+)', int),
    }
    rows = numeric_table(path, 2)
    out['depth'] = [r[0] for r in rows]
    out['ion_dist'] = [r[1] for r in rows]
    out['recoil_dist'] = [r[2:] for r in rows]
    return out


def read_columns(path, min_cols=2):
    rows = numeric_table(path, min_cols)
    return [list(c) for c in zip(*rows)] if rows else []


def read_vacancy(path):
    lines = _read(path)
    cols = read_columns(path, 2)
    return {
        'total_vacancies_per_ion': _find(lines, r'Total Target Vacancies\s*=\s*([\d.E+-]+)'),
        'total_displacements_per_ion': _find(lines, r'Total Target Displacements\s*=\s*([\d.E+-]+)'),
        'total_replacements_per_ion': _find(lines, r'Total Target Replacement Collisions\s*=\s*([\d.E+-]+)'),
        'kp_mode': any('VACANCIES' in l and 'by' in l.lower() for l in lines) or any('RECOILS' in l for l in lines),
        # TRIM's surface-sputtering mode (cascade code 3) writes the quick-mode
        # layout, but its columns are knock-ons and vacancies.
        'sputtering_mode': any('Vacancies not accurately calculated in SPUTTERING' in l for l in lines)
        and any('IONS' in l and 'RECOILS' in l for l in lines),
        'depth': cols[0] if cols else [],
        'columns': cols[1:] if cols else [],
    }


def read_tdata(path):
    lines = _read(path)
    return {
        'range': _find(lines, r'Average Range\s*=\s*([\d.E+-]+)'),
        'straggling': _find(lines, r'Average Straggling\s*=\s*([\d.E+-]+)'),
        'vacancies_per_ion': _find(lines, r'Average Vacancy/Ion\s*=\s*([\d.E+-]+)'),
        'backscattered': _find(lines, r'Total Backscattered Ions\s*=\s*(\d+)', int),
        'transmitted': _find(lines, r'Total Transmitted Ions\s*=\s*(\d+)', int),
        'ions': _find(lines, r'Total Ions calculated\s*=\s*(\d+)', int),
    }


def read_sputter_like(path):
    """BACKSCAT.txt / TRANSMIT.txt / TRIMOUT.txt / SPUTTER.txt atom records (code, ion, Z, E(eV), x, y, z, cx, cy, cz)."""
    recs = []
    if not os.path.exists(path):
        return recs
    for line in _read(path):
        p = line.split()
        if len(p) >= 10 and p[0] in ('B', 'T', 'S') and p[1].isdigit():
            try:
                recs.append({'kind': p[0], 'ion': int(p[1]), 'z': int(p[2]), 'energy_eV': float(p[3]),
                             'x': float(p[4]), 'y': float(p[5]), 'z_pos': float(p[6]),
                             'cx': float(p[7]), 'cy': float(p[8]), 'cz': float(p[9])})
            except ValueError:
                pass
    return recs


def read_range3d(path):
    xs = []
    started = False
    for line in _read(path):
        if line.startswith('-------'):
            started = True
            continue
        if started:
            p = line.split()
            if len(p) == 4:
                try:
                    xs.append((float(p[1]), float(p[2]), float(p[3])))
                except ValueError:
                    pass
    return xs


def read_case(directory):
    """Read every recognised TRIM output file in a directory into one dict."""
    d = {'dir': directory}
    def p(name):
        return os.path.join(directory, name)
    if os.path.exists(p('RANGE.txt')):
        d['range'] = read_range(p('RANGE.txt'))
    if os.path.exists(p('VACANCY.txt')):
        d['vacancy'] = read_vacancy(p('VACANCY.txt'))
    for key, name in (('ioniz', 'IONIZ.txt'), ('phonon', 'PHONON.txt'), ('e2recoil', 'E2RECOIL.txt'),
                      ('novac', 'NOVAC.txt'), ('lateral', 'LATERAL.txt')):
        if os.path.exists(p(name)):
            d[key] = read_columns(p(name), 2)
    if os.path.exists(p('TDATA.txt')):
        d['tdata'] = read_tdata(p('TDATA.txt'))
    for key, name in (('backscat', 'BACKSCAT.txt'), ('transmit', 'TRANSMIT.txt'), ('trimout', 'TRIMOUT.txt')):
        d[key] = read_sputter_like(p(name))
    if os.path.exists(p('RANGE_3D.txt')):
        d['range3d'] = read_range3d(p('RANGE_3D.txt'))
    return d


def integrate(depth, values):
    """Sum of a per-Angstrom distribution times the bin width (uniform bins)."""
    if len(depth) < 2:
        return 0.0
    w = depth[1] - depth[0]
    return sum(values) * w
