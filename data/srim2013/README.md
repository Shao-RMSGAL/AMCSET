# SRIM-2013 stopping tables

`srim2013_stopping.amcset` holds the electronic and nuclear stopping cross
sections of SRIM-2013.00 for every ion (Z1 = 1-92, most abundant isotope
mass) in every element (Z2 = 1-92), solid and gas phase: 16928 tables at
209 energies from 10 eV to 2.09 GeV (listed in the file),
in eV/(1e15 atoms/cm2). They were computed with SRIM-2013's SR Module, which
uses the same stopping coefficients (Data/SCOEF03.DAT) and the same
algorithm as TRIM.exe, and packed by `tools/srim` (`gen_cases.py` makes the
job list, `extract_stopping_tables.ps1` runs SR Module, `pack_stopping.py`
writes this file).

AMCSET's SRIM mode uses this file by default (stopping model `srim` with no
table given); a run reports "SRIM tables (bundled)". It looks for it in, in
order: `$AMCSET_DATA_DIR/srim2013/`, `data/srim2013/` next to the
executable or up to two folders above it, `share/amcset/srim2013/` next to
the executable's folder, and this source folder.

SRIM is (c) James F. Ziegler (www.srim.org) and free for non-commercial
use; these tables are derived from it and carry the same terms.
