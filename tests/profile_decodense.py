#!/usr/bin/env python
"""
Isolate decodense's own computation time from import/SCF overhead,
across all three partitioning modes, on a bigger molecule (pyridine).

Not a test -- a throwaway diagnostic script for profiling PR 6's
target bottleneck (_make_rho). Run from inside tests/:

    python3 profile_decodense.py
"""
import cProfile

# --- everything below this line happens BEFORE profiling starts ---
from pyscf import gto, scf
import decodense

mol = gto.M(
    atom="geom/c5h5n.xyz",
    basis="pcseg1",
    verbose=0,
)

mf = scf.RKS(mol)
mf.xc = "pbe0"
mf.kernel()

mo_coeff = (mf.mo_coeff, mf.mo_coeff)
# --- import/SCF overhead ends here ---


if __name__ == "__main__":
    for part in ("atoms", "eda", "orbitals"):
        decomp = decodense.DecompCls(part=part)
        # warm-up: trigger any lazy imports/one-time costs before profiling
        decodense.main(mol, decomp, mf, mo_coeff)
        # now profile a clean, warmed-up run
        cProfile.run(
            "decodense.main(mol, decomp, mf, mo_coeff)",
            f"profile_{part}.prof",
        )
