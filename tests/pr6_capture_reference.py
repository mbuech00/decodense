#!/usr/bin/env python
"""
Not a test -- captures per-atom/per-orbital reference numbers for PR 6's
density reformulation, before any code change. Run from inside tests/:

    python3 pr6_capture_reference.py
"""
import numpy as np
from pyscf import gto, scf
import decodense

mol = gto.M(atom="geom/c5h5n.xyz", basis="pcseg1", verbose=0)
mf = scf.RKS(mol)
mf.xc = "pbe0"
mf.kernel()
mo_coeff = (mf.mo_coeff, mf.mo_coeff)

reference = {}
for part in ("atoms", "orbitals"):
    decomp = decodense.DecompCls(part=part)
    res = decodense.main(mol, decomp, mf, mo_coeff)
    reference[part] = dict(res.res_dict)

np.savez(
    "pr6_reference.npz",
    **{
        f"{part}_{key}": val
        for part, comps in reference.items()
        for key, val in comps.items()
    },
)
print("Saved reference values to pr6_reference.npz")
