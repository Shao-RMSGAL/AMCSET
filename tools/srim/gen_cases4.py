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

"""TRIM reference cases, batch 4: Ni ions (the Ni stopping of SRIM-2013, docs/VALIDATION.md §3, §9).

Usage: python gen_cases4.py <atomdata.json> <out_dir>
  atomdata.json is written by gen_cases.py (step 1 of tools/srim/README.md);
  out_dir (normally <SRIM>/_amcset/trim_cases4) is replaced. Run the cases with
  run_trim_cases3.ps1 -CaseSet trim_cases4 -OutSet trim4 (results in results/trim4).

Ni ions in Fe, Ni, Mo and Si at 1-11.5 MeV (where SRIM-2013's Ni stopping is up
to 40 % below Co's and Cu's), with Co and Cu ions as controls; quick
(Kinchin-Pease) calculations except one with full cascades. They show that
AMCSET reproduces TRIM.exe for Ni ions with physics.stopping_corrections off,
and measure what the correction changes.
"""
import json, os, shutil, sys
ATOM = {int(k): v for k, v in json.load(open(sys.argv[1])).items()}
SYM = {v['sym']: k for k, v in ATOM.items()}
OUT = sys.argv[2]
shutil.rmtree(OUT, ignore_errors=True); os.makedirs(OUT)

def trim_in(name, ion, e_kev, layers, mode=1, n=1000, angle=0.0, seed=1, ed=None, eb=None, es=None,
            files=(1, 1, 1, 1, 0, 0), plot=(5, 0, 0), elems=None, bragg=None):
    """layers: (layer_name, width_A, density, {elem_key: stoich}, phase); elem_key = symbol or (symbol, tag)."""
    z1 = SYM[ion]; m1 = ATOM[z1]['mai']
    if elems is None:
        elems = []
        for _, _, _, comp, _ in layers:
            for s in comp:
                if s not in elems: elems.append(s)
    sym = lambda k: k if isinstance(k, str) else k[0]
    ed = ed or {}; eb = eb or {}; es = es or {}
    L = ['==> SRIM-2013.00 This file controls TRIM Calculations.',
         'Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.',
         f'     {z1}   {m1:.3f}   {e_kev:g}   {angle:g}   {n}   1   {n}',
         'Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders',
         f'                      {mode}                                   {seed}       0',
         'Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file',
         '                          ' + '       '.join(str(x) for x in files),
         'Target material : Number of Elements & Layers',
         f'"{name}"       {len(elems)}               {len(layers)}',
         'PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]',
         f'       {plot[0]}                         {plot[1]:g}            {plot[2]:g}',
         'Target Elements:    Z   Mass(amu)']
    for i, k in enumerate(elems):
        s = sym(k); L.append(f'Atom {i+1} = {s} =      {SYM[s]}  {ATOM[SYM[s]]["nat"]:.3f}')
    L.append('Layer   Layer Name /               Width Density    ' + '   '.join(f'{sym(k)}({SYM[sym(k)]})' for k in elems))
    L.append('Numb.   Description                (Ang) (g/cm3)    ' + '   '.join('Stoich' for _ in elems))
    for i, (ln, w, rho, comp, ph) in enumerate(layers):
        tot = sum(comp.values())
        L.append(f' {i+1}      "{ln}"           {w:g}  {rho:.6g}  ' + '  '.join(f'{comp.get(k, 0)/tot:.6f}' for k in elems))
    L.append('0  Target layer phases (0=Solid, 1=Gas)')
    L.append(' '.join(str(l[4]) for l in layers) + ' ')
    L.append('Target Compound Corrections (Bragg)')
    L.append(' ' + ' '.join(f'{b:g}' for b in (bragg or [1] * len(layers))) + '  ')
    L.append('Individual target atom displacement energies (eV)')
    L.append('   ' + '   '.join(f'{ed.get(k, ed.get(sym(k), 40.0)):.2f}' for k in elems))
    L.append('Individual target atom lattice binding energies (eV)')
    L.append('   ' + '   '.join(f'{eb.get(k, eb.get(sym(k), 0.0)):.2f}' for k in elems))
    L.append('Individual target atom surface binding energies (eV)')
    L.append('   ' + '   '.join(f'{es.get(k, es.get(sym(k), ATOM[SYM[sym(k)]]["es"])):.2f}' for k in elems))
    L.append('Stopping Power Version (1=2011, 0=2011)')
    L.append(' 0 ')
    return '\r\n'.join(L) + '\r\n'

cases = []

def add(cid, est, ion, e_kev, layers, **kw):
    files = {'TRIM.IN': trim_in(cid, ion, e_kev, layers, **kw)}
    meta = dict(id=cid, ion=ion, energy_keV=e_kev, est_s=est,
                kw={k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()}) for k, v in kw.items()})
    cases.append((est, cid, files, meta))

def el(sym, w):
    a = ATOM[SYM[sym]]
    return [(sym, w, a['rho'], {sym: 1.0}, 0)]

# name, ion, E (keV), target, thickness (A), ions, seed
for cid, ion, e, tgt, w, n, seed in [
        ('ni_fe_1000keV_kp', 'Ni', 1000, 'Fe', 15000, 4000, 7001),
        ('ni_fe_3000keV_kp', 'Ni', 3000, 'Fe', 30000, 4000, 7002),
        ('ni_fe_5000keV_kp', 'Ni', 5000, 'Fe', 40000, 4000, 7003),
        ('ni_fe_11500keV_kp', 'Ni', 11500, 'Fe', 60000, 3000, 7004),
        ('ni_mo_5000keV_kp', 'Ni', 5000, 'Mo', 40000, 4000, 7005),
        ('ni_ni_2000keV_kp', 'Ni', 2000, 'Ni', 25000, 4000, 7006),
        ('ni_ni_5000keV_kp', 'Ni', 5000, 'Ni', 40000, 4000, 7007),
        ('ni_si_3000keV_kp', 'Ni', 3000, 'Si', 50000, 4000, 7008),
        ('co_fe_5000keV_kp', 'Co', 5000, 'Fe', 40000, 4000, 7009),
        ('cu_fe_5000keV_kp', 'Cu', 5000, 'Fe', 40000, 4000, 7010)]:
    add(cid, 900, ion, e, el(tgt, w), mode=1, n=n, seed=seed)
add('ni_fe_5000keV_fc', 1800, 'Ni', 5000, el('Fe', 40000), mode=2, n=500, seed=7011)

cases.sort(key=lambda c: -c[0])
index = []
for k, (est, cid, files, meta) in enumerate(cases):
    d = os.path.join(OUT, f'{400 + k:03d}_{cid}')
    os.makedirs(d)
    for fn, txt in files.items():
        open(os.path.join(d, fn), 'w', newline='', encoding='latin-1').write(txt)
    meta['dir'] = os.path.basename(d); index.append(meta)
json.dump(index, open(os.path.join(OUT, 'cases.json'), 'w'), indent=1)
print(len(cases), 'cases')
