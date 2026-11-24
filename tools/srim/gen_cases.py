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

"""Generate the TRIM reference-case matrix (TRIM.IN files) and SR Module jobs."""
import csv, json, math, os, sys
OUT = sys.argv[1]
ATOM = {}
for line in open(os.environ['ATOMDATA'], encoding='latin-1').read().replace('\r','').split('\n')[2:]:
    f = line.split()
    if len(f) < 12 or not f[0].isdigit(): continue
    z = int(f[0]); sym = f[1]
    # columns: Z sym name MAImass MAIweight natural density atomdens fermi heatsubl gasdens gasatomdens
    ATOM[z] = dict(z=z, sym=sym, name=f[2], mai_a=int(f[3]), mai=float(f[4]), nat=float(f[5]),
                   rho=float(f[6]), ndens=float(f[7]), vf=float(f[8]), es=float(f[9]),
                   gasrho=float(f[10]), gasn=float(f[11]))
SYM = {v['sym']: k for k, v in ATOM.items()}

def trim_in(name, ion, e_kev, layers, mode=1, n=1000, angle=0.0, seed=1, ed=None, eb=None, es=None,
            files=(1, 1, 1, 1, 0, 0), plot=(5, 0, 0)):
    """layers: list of (layer_name, width_A, density, {sym: stoich}, phase)"""
    z1 = SYM[ion]; m1 = ATOM[z1]['mai']
    elems = []
    for _, _, _, comp, _ in layers:
        for s in comp:
            if s not in elems: elems.append(s)
    ed = ed or {}; eb = eb or {}; es = es or {}
    L = []
    L.append('==> SRIM-2013.00 This file controls TRIM Calculations.')
    L.append('Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.')
    L.append(f'     {z1}   {m1:.3f}   {e_kev:g}   {angle:g}   {n}   1   {n}')
    L.append('Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders')
    L.append(f'                      {mode}                                   {seed}       0')
    L.append('Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file')
    L.append('                          ' + '       '.join(str(x) for x in files))
    L.append('Target material : Number of Elements & Layers')
    L.append(f'"{name}"       {len(elems)}               {len(layers)}')
    L.append('PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]')
    L.append(f'       {plot[0]}                         {plot[1]:g}            {plot[2]:g}')
    L.append('Target Elements:    Z   Mass(amu)')
    for i, s in enumerate(elems):
        L.append(f'Atom {i+1} = {s} =      {SYM[s]}  {ATOM[SYM[s]]["nat"]:.3f}')
    hdr = 'Layer   Layer Name /               Width Density    ' + '   '.join(f'{s}({SYM[s]})' for s in elems)
    L.append(hdr)
    L.append('Numb.   Description                (Ang) (g/cm3)    ' + '   '.join('Stoich' for _ in elems))
    for i, (ln, w, rho, comp, ph) in enumerate(layers):
        tot = sum(comp.values())
        st = '  '.join(f'{comp.get(s, 0)/tot:.6f}' for s in elems)
        L.append(f' {i+1}      "{ln}"           {w:g}  {rho:.4f}  {st}')
    L.append('0  Target layer phases (0=Solid, 1=Gas)')
    L.append(' '.join(str(l[4]) for l in layers) + ' ')
    L.append('Target Compound Corrections (Bragg)')
    L.append(' ' + ' '.join('1' for _ in layers) + '  ')
    L.append('Individual target atom displacement energies (eV)')
    L.append('   ' + '   '.join(f'{ed.get(s, 40.0):.2f}' for s in elems))
    L.append('Individual target atom lattice binding energies (eV)')
    L.append('   ' + '   '.join(f'{eb.get(s, 0.0):.2f}' for s in elems))
    L.append('Individual target atom surface binding energies (eV)')
    L.append('   ' + '   '.join(f'{es.get(s, ATOM[SYM[s]]["es"]):.2f}' for s in elems))
    L.append('Stopping Power Version (1=2011, 0=2011)')
    L.append(' 0 ')
    return '\r\n'.join(L) + '\r\n'

cases = []
def add(cid, ion, e_kev, layers, **kw):
    meta = dict(id=cid, ion=ion, energy_keV=e_kev, layers=[dict(name=l[0], width_A=l[1], density_g_cm3=l[2], composition=l[3], phase=l[4]) for l in layers])
    meta.update({k: v for k, v in kw.items() if k != 'files'})
    if 'files' in kw: meta['files'] = kw['files']
    cases.append((cid, trim_in(cid, ion, e_kev, layers, **kw), meta))

def elem_layer(sym, width):
    a = ATOM[SYM[sym]]
    return [(sym, width, a['rho'], {sym: 1.0}, 0)]

modes = {1: 'kp', 2: 'fc'}
# --- Priority: Fe -> Fe (Ed=40, Eb=0, Stoller et al. 2013 recommendations)
fe_w = {1: 200, 10: 600, 100: 2500, 500: 8000, 1000: 12000, 2000: 18000, 5000: 35000}
seed = 1
for e, w in fe_w.items():
    for m in (1, 2):
        n = (20000 if e <= 10 else 10000) if m == 1 else (5000 if e <= 100 else (2000 if e <= 2000 else 1000))
        add(f'fe_fe_{e}keV_{modes[m]}', 'Fe', e, elem_layer('Fe', w), mode=m, n=n, seed=seed); seed += 1
# Fe -> Fe 1 MeV with SRIM-style damage energies (Ed=25, Eb=3) and oblique incidence
for m in (1, 2):
    add(f'fe_fe_1000keV_{modes[m]}_ed25eb3', 'Fe', 1000, elem_layer('Fe', 12000), mode=m, n=10000 if m == 1 else 2000,
        seed=100 + m, ed={'Fe': 25.0}, eb={'Fe': 3.0})
    add(f'fe_fe_1000keV_{modes[m]}_60deg', 'Fe', 1000, elem_layer('Fe', 12000), mode=m, n=10000 if m == 1 else 2000,
        seed=110 + m, angle=60.0)
# --- H and He into Fe
for ion, ws in (('H', {10: 3000, 100: 15000, 1000: 150000, 3000: 1000000}), ('He', {10: 1500, 100: 6000, 1000: 50000, 3000: 150000})):
    for e, w in ws.items():
        for m in (1, 2):
            add(f'{ion.lower()}_fe_{e}keV_{modes[m]}', ion, e, elem_layer('Fe', w), mode=m, n=20000 if m == 1 else 5000, seed=seed); seed += 1
# --- Other transition metals (self-ion 1 MeV, H and He 1 MeV)
for sym in ('Cr', 'Ni', 'Cu', 'Ti', 'V', 'Co', 'Mn', 'Zr', 'Nb', 'Mo', 'Ta', 'W'):
    for m in (1, 2):
        add(f'{sym.lower()}_{sym.lower()}_1000keV_{modes[m]}', sym, 1000, elem_layer(sym, 15000), mode=m, n=5000 if m == 1 else 1000, seed=seed); seed += 1
    add(f'h_{sym.lower()}_1000keV_kp', 'H', 1000, elem_layer(sym, 150000), mode=1, n=10000, seed=seed); seed += 1
    add(f'he_{sym.lower()}_1000keV_kp', 'He', 1000, elem_layer(sym, 50000), mode=1, n=10000, seed=seed); seed += 1
    add(f'he_{sym.lower()}_1000keV_fc', 'He', 1000, elem_layer(sym, 50000), mode=2, n=2000, seed=seed); seed += 1
# --- Alloys and multilayers
fecr = [('Fe-9Cr', 20000, 7.78, {'Fe': 0.91, 'Cr': 0.09}, 0)]
ss316 = [('316SS', 20000, 7.99, {'Fe': 0.655, 'Cr': 0.18, 'Ni': 0.12, 'Mo': 0.015, 'Mn': 0.02, 'Si': 0.01}, 0)]
for m in (1, 2):
    add(f'fe_fecr_2000keV_{modes[m]}', 'Fe', 2000, fecr, mode=m, n=10000 if m == 1 else 1000, seed=seed); seed += 1
    add(f'fe_316ss_2000keV_{modes[m]}', 'Fe', 2000, ss316, mode=m, n=10000 if m == 1 else 1000, seed=seed); seed += 1
    add(f'ni_ni_2000keV_{modes[m]}', 'Ni', 2000, elem_layer('Ni', 18000), mode=m, n=10000 if m == 1 else 1000, seed=seed); seed += 1
add('he_fecr_1000keV_kp', 'He', 1000, [('Fe-9Cr', 50000, 7.78, {'Fe': 0.91, 'Cr': 0.09}, 0)], mode=1, n=20000, seed=seed); seed += 1
crfe = [('Cr', 3000, ATOM[24]['rho'], {'Cr': 1.0}, 0), ('Fe', 20000, ATOM[26]['rho'], {'Fe': 1.0}, 0)]
for m in (1, 2):
    add(f'fe_crfe_2000keV_{modes[m]}', 'Fe', 2000, crfe, mode=m, n=10000 if m == 1 else 1000, seed=seed); seed += 1
add('h_crfe_1000keV_kp', 'H', 1000, [('Cr', 3000, ATOM[24]['rho'], {'Cr': 1.0}, 0), ('Fe', 150000, ATOM[26]['rho'], {'Fe': 1.0}, 0)], mode=1, n=20000, seed=seed); seed += 1
# --- Backscattering, transmission, sputtering
add('he_w_10keV_kp', 'He', 10, elem_layer('W', 1500), mode=1, n=20000, seed=seed); seed += 1
add('h_w_5keV_kp', 'H', 5, elem_layer('W', 1500), mode=1, n=20000, seed=seed); seed += 1
add('h_fe_1000keV_thin_kp', 'H', 1000, elem_layer('Fe', 30000), mode=1, n=20000, seed=seed); seed += 1
add('he_fe_2000keV_thin_kp', 'He', 2000, elem_layer('Fe', 20000), mode=1, n=20000, seed=seed); seed += 1
add('ar_cu_1keV_fc', 'Ar', 1, elem_layer('Cu', 300), mode=2, n=5000, seed=seed); seed += 1
add('ar_cu_1keV_sputt', 'Ar', 1, elem_layer('Cu', 300), mode=3, n=5000, seed=seed); seed += 1
add('fe_fe_1keV_sputt', 'Fe', 1, elem_layer('Fe', 200), mode=3, n=5000, seed=seed); seed += 1
# --- Diagnostics: full collision listing and EXYZ trajectories (few ions)
add('diag_fe_fe_100keV_fc_coll', 'Fe', 100, elem_layer('Fe', 2500), mode=2, n=20, seed=7, files=(1, 1, 1, 1, 2, 0))
add('diag_fe_fe_100keV_kp_coll', 'Fe', 100, elem_layer('Fe', 2500), mode=1, n=200, seed=7, files=(1, 1, 1, 1, 1, 0))
add('diag_h_fe_1000keV_kp_exyz', 'H', 1000, elem_layer('Fe', 150000), mode=1, n=50, seed=7, files=(1, 1, 1, 1, 0, 1000))
add('diag_fe_fe_1000keV_kp_exyz', 'Fe', 1000, elem_layer('Fe', 12000), mode=1, n=50, seed=7, files=(1, 1, 1, 1, 0, 1000))
add('diag_he_fe_1000keV_kp_exyz', 'He', 1000, elem_layer('Fe', 50000), mode=1, n=50, seed=7, files=(1, 1, 1, 1, 0, 1000))

# priority order: diagnostics first, then Fe->Fe, H/He->Fe, then the rest
def prio(c):
    cid = c[0]
    if cid.startswith('diag'): return 0
    if cid.startswith('fe_fe'): return 1
    if cid.endswith('_fe_') or '_fe_' in cid[:6]: return 2
    return 3
cases.sort(key=lambda c: (prio(c), c[0]))
os.makedirs(os.path.join(OUT, 'trim_cases'), exist_ok=True)
index = []
for k, (cid, txt, meta) in enumerate(cases):
    d = os.path.join(OUT, 'trim_cases', f'{k:03d}_{cid}')
    os.makedirs(d, exist_ok=True)
    open(os.path.join(d, 'TRIM.IN'), 'w', newline='').write(txt)
    meta['dir'] = f'{k:03d}_{cid}'
    index.append(meta)
json.dump(index, open(os.path.join(OUT, 'trim_cases', 'cases.json'), 'w'), indent=1)

# --- SR Module jobs: every Z1 (1..92) into every Z2 (1..92) as a solid and as a gas.
# TRIM applies SRIM's gas-phase stopping to every element of a gas layer (e.g.
# the C of CO2), not only to the elements that are gases. The natural gases
# come first.
natural_gases = [1, 2, 7, 8, 10, 18, 36, 54]
with open(os.path.join(OUT, 'sr_jobs.csv'), 'w', newline='') as f:
    w = csv.writer(f, lineterminator='\r\n')
    w.writerow(['id', 'z1', 'm1', 'z2', 'name', 'm2', 'phase', 'density'])
    # order: elements of interest first
    focus = [26, 1, 2, 24, 28, 25, 27, 29, 22, 23, 40, 41, 42, 73, 74, 14, 6, 8, 18, 79, 47, 46, 78, 72, 75, 21, 30, 39, 44, 45, 77, 5, 4, 3, 13, 7]
    order = focus + [z for z in range(1, 93) if z not in focus]
    for z2 in order:
        for z1 in order:
            a1, a2 = ATOM[z1], ATOM[z2]
            w.writerow([f'sr_{z1:03d}_{z2:03d}_s', z1, f"{a1['mai']:.3f}", z2, a2['name'], f"{a2['nat']:.3f}", 0, f"{a2['rho']:.4f}"])
    for z2 in natural_gases + [z for z in order if z not in natural_gases]:
        for z1 in order:
            a1, a2 = ATOM[z1], ATOM[z2]
            # The density only scales SR Module's range columns, not the stopping used here.
            rho = a2['gasrho'] if a2['gasrho'] > 0 else a2['rho']
            w.writerow([f'sr_{z1:03d}_{z2:03d}_g', z1, f"{a1['mai']:.3f}", z2, a2['name'], f"{a2['nat']:.3f}", 1, f"{rho:.6g}"])
# energies (keV): 10 eV ... 2 GeV, 25 per decade  (SR Module accepts up to 1000 energies)
es = [10 ** (-2 + i / 25) for i in range(0, 25 * 8 + 9)]
open(os.path.join(OUT, 'sr_energies_keV.txt'), 'w', newline='').write('\r\n'.join(f' {e:.6g}' for e in es) + '\r\n')
print(len(cases), 'TRIM cases;', len(es), 'energies')
json.dump({str(z): ATOM[z] for z in ATOM}, open(os.path.join(OUT, 'atomdata.json'), 'w'))
