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

"""TRIM reference cases, batch 3: seed replicates, layer/alloy diagnostics, Ni, gas, TRIM.DAT, options.

Usage: python gen_cases3.py <atomdata.json> <out_dir>
  atomdata.json is written by gen_cases.py (step 1 of tools/srim/README.md);
  out_dir (normally <SRIM>/_amcset/trim_cases3) is replaced. Run the cases with
  run_trim_cases3.bat (results in results/trim3) and the TRIM.DAT recoil cases
  that TRIM holds with a message box with run_held_cases.bat.
"""
import json, math, os, random, shutil, sys
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

def trim_dat(title, rows):
    """rows: (name, z, E_eV, x, y, z, cx, cy, cz)."""
    L = ['TRIM.DAT written by AMCSET gen3.py'] + ['comment line'] * 6 + [title, 'Event  Atom  Energy  Depth   Lateral-Position   ----- Atom Direction ----',
         'Name   Numb   (eV)    _X_(A)   _Y_(A)  _Z_(A)   Cos(X)   Cos(Y)   Cos(Z)']
    for r in rows:
        L.append(f'{r[0]:<5s}  {r[1]:3d}  {r[2]:.6g}  {r[3]:.6g}  {r[4]:.6g}  {r[5]:.6g}  {r[6]:.7f}  {r[7]:.7f}  {r[8]:.7f}')
    return '\r\n'.join(L) + '\r\n\x1a\r\n'

cases = []   # (est_seconds, name, files)
def add(cid, est, ion, e_kev, layers, dat=None, **kw):
    files = {'TRIM.IN': trim_in(cid, ion, e_kev, layers, **kw)}
    if dat: files['TRIM.DAT'] = dat
    meta = dict(id=cid, ion=ion, energy_keV=e_kev, est_s=est, kw={k: (v if not isinstance(v, dict) else {str(a): b for a, b in v.items()}) for k, v in kw.items()})
    cases.append((est, cid, files, meta))

def el(sym, w, phase=0, rho=None):
    a = ATOM[SYM[sym]]
    return [(sym, w, rho if rho else (a['gasrho'] if phase else a['rho']), {sym: 1.0}, phase)]
rho = lambda s: ATOM[SYM[s]]['rho']

# ---- S: seed replicates of existing reference cases (same TRIM.IN as the original, new seed)
fecr = [('Fe-9Cr', 20000, 7.78, {'Fe': 0.91, 'Cr': 0.09}, 0)]
ss316 = [('316SS', 20000, 7.99, {'Fe': 0.655, 'Cr': 0.18, 'Ni': 0.12, 'Mo': 0.015, 'Mn': 0.02, 'Si': 0.01}, 0)]
crfe = [('Cr', 3000, rho('Cr'), {'Cr': 1.0}, 0), ('Fe', 20000, rho('Fe'), {'Fe': 1.0}, 0)]
crfe_h = [('Cr', 3000, rho('Cr'), {'Cr': 1.0}, 0), ('Fe', 150000, rho('Fe'), {'Fe': 1.0}, 0)]
reps = [  # name, est_per_run, ion, E, layers, mode, n, kw, nseeds
    ('fe_crfe_2000keV_fc', 1215, 'Fe', 2000, crfe, 2, 1000, {}, 6),
    ('fe_316ss_2000keV_fc', 1524, 'Fe', 2000, ss316, 2, 1000, {}, 6),
    ('h_fe_3000keV_kp', 2016, 'H', 3000, el('Fe', 1000000), 1, 20000, {}, 6),
    ('h_fe_1000keV_kp', 1478, 'H', 1000, el('Fe', 150000), 1, 20000, {}, 6),
    ('h_crfe_1000keV_kp', 825, 'H', 1000, crfe_h, 1, 20000, {}, 4),
    ('fe_fe_10keV_fc', 163, 'Fe', 10, el('Fe', 600), 2, 5000, {}, 4),
    ('fe_fe_100keV_fc', 868, 'Fe', 100, el('Fe', 2500), 2, 5000, {}, 4),
    ('fe_fe_1000keV_kp_ed25eb3', 1896, 'Fe', 1000, el('Fe', 12000), 1, 10000, dict(ed={'Fe': 25.0}, eb={'Fe': 3.0}), 4),
    ('he_fe_10keV_fc', 638, 'He', 10, el('Fe', 1500), 2, 5000, {}, 4),
    ('h_nb_1000keV_kp', 352, 'H', 1000, el('Nb', 150000), 1, 10000, {}, 4),
    ('h_ni_1000keV_kp', 389, 'H', 1000, el('Ni', 150000), 1, 10000, {}, 4),
    ('fe_fe_2000keV_fc', 2248, 'Fe', 2000, el('Fe', 18000), 2, 2000, {}, 6),
    ('h_fe_1000keV_fc', 703, 'H', 1000, el('Fe', 150000), 2, 5000, {}, 4),
    ('fe_fe_500keV_kp', 2765, 'Fe', 500, el('Fe', 8000), 1, 10000, {}, 2),
    ('fe_fe_500keV_fc', 1518, 'Fe', 500, el('Fe', 8000), 2, 2000, {}, 2),
    ('fe_fe_1000keV_fc_ed25eb3', 2037, 'Fe', 1000, el('Fe', 12000), 2, 2000, dict(ed={'Fe': 25.0}, eb={'Fe': 3.0}), 2),
    ('fe_fe_1000keV_fc', 1768, 'Fe', 1000, el('Fe', 12000), 2, 2000, {}, 2),
    ('fe_fe_100keV_kp', 699, 'Fe', 100, el('Fe', 2500), 1, 10000, {}, 2),
    ('fe_fe_1keV_sputt', 60, 'Fe', 1, el('Fe', 200), 3, 5000, {}, 4),
    ('ar_cu_1keV_sputt', 118, 'Ar', 1, el('Cu', 300), 3, 5000, {}, 4),
    ('ar_cu_1keV_fc', 118, 'Ar', 1, el('Cu', 300), 2, 5000, {}, 4),
    ('fe_fe_1keV_fc', 75, 'Fe', 1, el('Fe', 200), 2, 5000, {}, 4),
]
for name, est, ion, e, layers, mode, n, kw, ns in reps:
    for k in range(ns):
        add(f'rep_{name}_s{k+1}', est, ion, e, layers, mode=mode, n=n, seed=5000 + 37 * k + len(name), **kw)

# ---- L: layer and alloy diagnostics (Fe 2 MeV, full cascades)
add('lay_fe_fefe_2000keV_fc_s1', 2300, 'Fe', 2000, [('Fe1', 3000, rho('Fe'), {'Fe': 1.0}, 0), ('Fe2', 15000, rho('Fe'), {'Fe': 1.0}, 0)], mode=2, n=2000, seed=6001)
add('lay_fe_fefe_2000keV_fc_s2', 2300, 'Fe', 2000, [('Fe1', 3000, rho('Fe'), {'Fe': 1.0}, 0), ('Fe2', 15000, rho('Fe'), {'Fe': 1.0}, 0)], mode=2, n=2000, seed=6002)
add('lay_fe_fe2el_2000keV_fc', 2300, 'Fe', 2000, [('FeA', 3000, rho('Fe'), {('Fe', 'a'): 1.0}, 0), ('FeB', 15000, rho('Fe'), {('Fe', 'b'): 1.0}, 0)],
    mode=2, n=2000, seed=6003, elems=[('Fe', 'a'), ('Fe', 'b')])
for k in range(2):
    add(f'lay_fe_crfe_order_2000keV_fc_s{k+1}', 1250, 'Fe', 2000, crfe, mode=2, n=1000, seed=6010 + k, elems=['Fe', 'Cr'])
    add(f'lay_fe_fecr_rev_2000keV_fc_s{k+1}', 1250, 'Fe', 2000, [('Fe', 3000, rho('Fe'), {'Fe': 1.0}, 0), ('Cr', 20000, rho('Cr'), {'Cr': 1.0}, 0)], mode=2, n=1000, seed=6020 + k)
    add(f'lay_fe_cr_2000keV_fc_s{k+1}', 1250, 'Fe', 2000, el('Cr', 20000), mode=2, n=1000, seed=6030 + k)
add('lay_fe_crfe_2000keV_kp', 2400, 'Fe', 2000, crfe, mode=1, n=5000, seed=6040)
add('aly_fe_feni50_1000keV_fc', 1300, 'Fe', 1000, [('FeNi', 15000, 8.4, {'Fe': 0.5, 'Ni': 0.5}, 0)], mode=2, n=1000, seed=6050)
add('aly_fe_fecr50_1000keV_fc', 1300, 'Fe', 1000, [('FeCr', 15000, 7.5, {'Fe': 0.5, 'Cr': 0.5}, 0)], mode=2, n=1000, seed=6051)
add('aly_fe_fecr50_1000keV_kp', 1900, 'Fe', 1000, [('FeCr', 15000, 7.5, {'Fe': 0.5, 'Cr': 0.5}, 0)], mode=1, n=10000, seed=6052)

# ---- N: nickel
add('ni_fe_ni_1000keV_fc', 1300, 'Fe', 1000, el('Ni', 15000), mode=2, n=1000, seed=6100)
add('ni_ni_ni_100keV_fc', 1000, 'Ni', 100, el('Ni', 2500), mode=2, n=5000, seed=6101)
add('ni_h_ni_1000keV_fc', 700, 'H', 1000, el('Ni', 150000), mode=2, n=5000, seed=6102)
add('ni_he_ni_100keV_kp', 1500, 'He', 100, el('Ni', 6000), mode=1, n=20000, seed=6103)

# ---- G: gas targets (TRIM's gas phase) and their solid-phase twins
g = lambda s: ATOM[SYM[s]]['gasrho']
for ph in (1, 0):
    t = 'gas' if ph else 'solidph'
    add(f'gas_h_n_1000keV_kp_{t}', 600, 'H', 1000, [('N2', 4e8, g('N'), {'N': 1.0}, ph)], mode=1, n=10000, seed=6200 + ph)
    add(f'gas_he_ar_1000keV_kp_{t}', 1400, 'He', 1000, [('Ar', 1.5e8, g('Ar'), {'Ar': 1.0}, ph)], mode=1, n=10000, seed=6210 + ph)
    add(f'gas_h_he_100keV_kp_{t}', 700, 'H', 100, [('He', 2e8, g('He'), {'He': 1.0}, ph)], mode=1, n=10000, seed=6220 + ph)
add('gas_he_ar_1000keV_fc_gas', 700, 'He', 1000, [('Ar', 1.5e8, g('Ar'), {'Ar': 1.0}, 1)], mode=2, n=2000, seed=6230)
add('gas_ar_ar_100keV_fc_gas', 800, 'Ar', 100, [('Ar', 3e7, g('Ar'), {'Ar': 1.0}, 1)], mode=2, n=1000, seed=6231)
add('gas_n_n_500keV_kp_gas', 1000, 'N', 500, [('N2', 1e8, g('N'), {'N': 1.0}, 1)], mode=1, n=5000, seed=6232)
add('gas_xe_xe_1000keV_fc_gas', 1200, 'Xe', 1000, [('Xe', 2e7, g('Xe'), {'Xe': 1.0}, 1)], mode=2, n=500, seed=6233)
add('gas_he_co2_2000keV_kp_gas', 1400, 'He', 2000, [('CO2', 3e8, 1.977e-3, {'C': 1.0, 'O': 2.0}, 1)], mode=1, n=10000, seed=6234)
add('gas_he_alaral_5000keV_kp_mixed', 1200, 'He', 5000,
    [('Al', 100000, rho('Al'), {'Al': 1.0}, 0), ('Ar', 3e8, g('Ar'), {'Ar': 1.0}, 1), ('Al2', 300000, rho('Al'), {'Al': 1.0}, 0)], mode=1, n=5000, seed=6235)

# ---- D: TRIM.DAT calculations (codes 4-7)
rng = random.Random(42)
def iso_dir(fwd=None):
    while True:
        c = 2 * rng.random() - 1
        if fwd is None or (c > 0) == fwd: break
    phi = 2 * math.pi * rng.random(); s = math.sqrt(1 - c * c)
    return c, s * math.cos(phi), s * math.sin(phi)
rows = []
for i in range(2000):
    e = 10e3 + 490e3 * rng.random(); x = 2000 * rng.random()
    cx, cy, cz = iso_dir(True) if i % 5 else iso_dir(False)
    if x < 1 and cx < 0: cx = -cx
    rows.append((f'I{i % 10000:04d}', 26, e, x, 0.0, 0.0, cx, cy, cz))
add('dat_fe_fe_ions_kp', 300, 'Fe', 500, el('Fe', 10000), mode=4, n=2000, seed=6300, dat=trim_dat('Fe ions 10-500 keV at 0-2000 A, various angles', rows))
add('dat_fe_fe_ions_fc', 800, 'Fe', 500, el('Fe', 10000), mode=5, n=1000, seed=6301, dat=trim_dat('Fe ions 10-500 keV at 0-2000 A, various angles', rows[:1000]))
pka = []
for i in range(3000):
    e = math.exp(math.log(100) + (math.log(5e4) - math.log(100)) * rng.random()); x = 1000 + 8000 * rng.random()
    cx, cy, cz = iso_dir()
    pka.append((f'P{i:04d}', 26, e, x, 0.0, 0.0, cx, cy, cz))
add('dat_fe_fe_pka_fc', 900, 'Fe', 50, el('Fe', 10000), mode=7, n=3000, seed=6302, dat=trim_dat('Fe PKAs 0.1-50 keV in Fe', pka))
add('dat_fe_fe_pka_kp', 300, 'Fe', 50, el('Fe', 10000), mode=6, n=3000, seed=6303, dat=trim_dat('Fe PKAs 0.1-50 keV in Fe', pka))
pka2 = [(r[0], 24 if i % 4 == 0 else 26) + r[2:] for i, r in enumerate(pka[:2000])]
add('dat_fecr_pka_fc', 700, 'Fe', 50, fecr, mode=7, n=2000, seed=6304, dat=trim_dat('Fe and Cr PKAs in Fe-9Cr', pka2))

# ---- O: other inputs and scenarios
add('opt_fe_316ss_2000keV_kp_bragg09', 2400, 'Fe', 2000, ss316, mode=1, n=5000, seed=6400, bragg=[0.9])
add('opt_ar_cu_1keV_fc_75deg', 200, 'Ar', 1, el('Cu', 300), mode=2, n=5000, seed=6401, angle=75.0)
add('opt_ar_cu_5keV_fc_45deg', 400, 'Ar', 5, el('Cu', 600), mode=2, n=5000, seed=6402, angle=45.0)
add('opt_fe_fe_1000keV_kp_window', 1900, 'Fe', 1000, el('Fe', 12000), mode=1, n=10000, seed=6403, plot=(5, 1000, 3000))
add('opt_fe_fe_100keV_fc_eb3', 900, 'Fe', 100, el('Fe', 2500), mode=2, n=5000, seed=6404, eb={'Fe': 3.0})
add('opt_fe_fe_10keV_fc_es0ed20', 200, 'Fe', 10, el('Fe', 600), mode=2, n=5000, seed=6405, es={'Fe': 2.0}, ed={'Fe': 20.0})
add('opt_he_fe_2000keV_thin5000_kp', 600, 'He', 2000, el('Fe', 5000), mode=1, n=20000, seed=6406)
add('opt_h_w_3000keV_thin10000_kp', 600, 'H', 3000, el('W', 10000), mode=1, n=20000, seed=6407)
add('opt_he_w_2000keV_kp', 2000, 'He', 2000, el('W', 60000), mode=1, n=20000, seed=6408)
add('opt_h_w_100keV_kp', 1000, 'H', 100, el('W', 8000), mode=1, n=20000, seed=6409)
add('opt_h_w_1keV_kp', 800, 'H', 1, el('W', 300), mode=1, n=20000, seed=6410)
add('opt_he_fe_2keV_fc', 300, 'He', 2, el('Fe', 400), mode=2, n=5000, seed=6411)
add('opt_h_si_100keV_kp', 1000, 'H', 100, el('Si', 20000), mode=1, n=20000, seed=6412)
add('opt_he_c_1000keV_kp', 1400, 'He', 1000, el('C', 50000), mode=1, n=10000, seed=6413)
add('opt_si_sic_1000keV_fc', 1300, 'Si', 1000, [('SiC', 20000, 3.21, {'Si': 1.0, 'C': 1.0}, 0)], mode=2, n=1000, seed=6414)
add('opt_o_al2o3_500keV_kp', 1000, 'O', 500, [('Al2O3', 20000, 3.97, {'Al': 2.0, 'O': 3.0}, 0)], mode=1, n=5000, seed=6415)
add('opt_xe_uo2_1000keV_fc', 1200, 'Xe', 1000, [('UO2', 6000, 10.97, {'U': 1.0, 'O': 2.0}, 0)], mode=2, n=500, seed=6416)
add('opt_au_c_1000keV_fc', 800, 'Au', 1000, el('C', 15000), mode=2, n=500, seed=6417)
add('opt_fe_fe_50000keV_kp', 1500, 'Fe', 50000, el('Fe', 120000), mode=1, n=2000, seed=6418)
add('opt_h_fe_10000keV_kp', 1500, 'H', 10000, el('Fe', 4000000), mode=1, n=5000, seed=6419)
add('opt_fe_fe_100keV_fc_coll2', 60, 'Fe', 100, el('Fe', 2500), mode=2, n=30, seed=6420, files=(1, 1, 1, 1, 2, 0))
add('opt_fe_fe_1keV_sputt_coll2', 30, 'Fe', 1, el('Fe', 200), mode=3, n=100, seed=6421, files=(1, 1, 1, 1, 2, 0))
add('opt_he_fe_1000keV_fc_exyz', 30, 'He', 1000, el('Fe', 50000), mode=2, n=50, seed=6422, files=(1, 1, 1, 1, 0, 10000))

cases.sort(key=lambda c: -c[0])
index = []
for k, (est, cid, files, meta) in enumerate(cases):
    d = os.path.join(OUT, f'{300 + k:03d}_{cid}')
    os.makedirs(d)
    for fn, txt in files.items():
        open(os.path.join(d, fn), 'w', newline='', encoding='latin-1').write(txt)
    meta['dir'] = os.path.basename(d); index.append(meta)
json.dump(index, open(os.path.join(OUT, 'cases.json'), 'w'), indent=1)
tot = sum(c[0] for c in cases)
print(len(cases), 'cases, est', tot, 's =', round(tot / 14 / 3600, 2), 'h on 14 workers; longest', cases[0][0])
