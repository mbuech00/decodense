#!/usr/bin/env python
# -*- coding: utf-8 -*

import numpy as np
import pytest
from unittest.mock import patch
from decodense.decodense import main
from decodense.orbitals import _population_mul, assign_rdm1s
from decodense.properties import (
    _dip_nuc,
    _get_nuc,
    _h_core,
    _make_rho,
    _make_rho_interm2,
    _trace,
    _e_xc,
    _xc_ao_deriv,
    _e_nuc,
    _vk_dft,
    _point_charges,
)
from decodense.tools import dim, make_rdm1, mf_info
from decodense.decomp import CompKeys, DecompCls, sanity_check
from decodense.results import atoms, orbs
from pyscf import gto, scf, dft
from pyscf.dft import numint


# shared fixtures
@pytest.fixture
def mf_h2o():
    mol = gto.M(
        verbose=0, output=None, basis="sto-3g", symmetry=True, atom="geom/h2o.xyz"
    )
    mf = scf.RHF(mol).run()
    return mf


@pytest.fixture
def mf_h2o_permute():
    mol = gto.M(
        verbose=0,
        output=None,
        basis="sto-3g",
        symmetry=True,
        atom="geom/h2o_permute.xyz",
    )
    mf = scf.RHF(mol).run()
    return mf


@pytest.fixture
def mf_h2o_dft():
    mol = gto.M(
        verbose=0, output=None, basis="sto-3g", symmetry=True, atom="geom/h2o.xyz"
    )
    mf = dft.RKS(mol, xc="pbe0").run()
    return mf


@pytest.fixture
def mf_li():
    mol = gto.M(
        verbose=0, output=None, spin=1, basis="sto-3g", symmetry=True, atom="Li 0 0 0"
    )
    mf = scf.UHF(mol).run()
    return mf

@pytest.fixture
def mf_li_rohf():
    mol = gto.M(
        verbose=0, output=None, spin=1, basis="sto-3g", symmetry=True, atom="Li 0 0 0"
    )
    mf = scf.ROHF(mol).run()
    return mf


# testing functions in decodense.py


# permutation invariance
@pytest.mark.parametrize("part_method", ["mo", "ao"])
def test_permute(mf_h2o, mf_h2o_permute, part_method):
    mol1 = mf_h2o.mol
    mo_coeff1, mo_occ1 = mf_info(mf_h2o)
    decomp = DecompCls(pop_method="iao", part="atoms", part_method=part_method)
    res = main(mol1, decomp, mf_h2o, mo_coeff1, mo_occ1)
    mol2 = mf_h2o_permute.mol
    mo_coeff2, mo_occ2 = mf_info(mf_h2o_permute)
    decomp_perm = DecompCls(pop_method="iao", part="atoms", part_method=part_method)
    res_perm = main(mol2, decomp_perm, mf_h2o_permute, mo_coeff2, mo_occ2)
    # every component, not just el: mol1's atoms 0,1,2 are mol2's atoms 2,0,1
    for key, val in res.res_dict.items():
        val_perm = res_perm.res_dict[key]
        assert np.isclose(val[0], val_perm[2])
        assert np.isclose(val[1], val_perm[0])
        assert np.isclose(val[2], val_perm[1])


# raw pyscf input (2-D mo_coeff + 1-D mo_occ) must land where the mf_info route lands.
def test_main_accepts_raw_pyscf_input(mf_li_rohf):
    mol = mf_li_rohf.mol
    mo_coeff, mo_occ = mf_info(mf_li_rohf)
    decomp = DecompCls(pop_method="iao", part="atoms")
    res = main(mol, decomp, mf_li_rohf, mo_coeff, mo_occ)
    decomp_raw = DecompCls(pop_method="iao", part="atoms")
    res_raw = main(mol, decomp_raw, mf_li_rohf, mf_li_rohf.mo_coeff, mf_li_rohf.mo_occ)
    assert np.allclose(res.el, res_raw.el)


# orbitals mode with a non-aufbau occupation (a gap in the occupied indices):
def test_main_orbitals_non_aufbau(mf_h2o):
    mol = mf_h2o.mol
    mf_u = scf.UHF(mol).run()
    occ_a = mf_u.mo_occ[0].copy()
    occ_a[4] = 0.0  
    occ_a[5] = 1.0  
    occ_b = mf_u.mo_occ[1].copy()
    decomp = DecompCls(pop_method="iao", part="orbitals")
    res = main(mol, decomp, mf_u, mf_u.mo_coeff, (occ_a, occ_b))
    idx_a = np.where(occ_a > 0.0)[0]
    idx_b = np.where(occ_b > 0.0)[0]
    mo_coeff_sliced = (mf_u.mo_coeff[0][:, idx_a], mf_u.mo_coeff[1][:, idx_b])
    mo_occ_sliced = (np.ones(idx_a.size), np.ones(idx_b.size))
    decomp_sliced = DecompCls(pop_method="iao", part="orbitals")
    res_sliced = main(mol, decomp_sliced, mf_u, mo_coeff_sliced, mo_occ_sliced)
    assert res.el[0].size == 5
    assert np.allclose(res.el[0], res_sliced.el[0])
    assert np.allclose(res.el[1], res_sliced.el[1])


# E_ne(1) + E_ne(2) split
@pytest.mark.parametrize("part_method", ["mo", "ao"])
def test_nuc_sum(mf_h2o, part_method):
    mol = mf_h2o.mol
    mo_coeff, mo_occ = mf_info(mf_h2o)
    decomp = DecompCls(pop_method="iao", part="atoms", part_method=part_method)
    res = main(mol, decomp, mf_h2o, mo_coeff, mo_occ)
    rdm1 = mf_h2o.make_rdm1()
    total_nuc_att = np.sum(mol.intor("int1e_nuc") * rdm1)
    assert np.isclose((res.nuc_att_glob + res.nuc_att_loc).sum(), total_nuc_att)


# decomp.py


# sanity_check
@pytest.mark.parametrize(
    "attr,bad_value,exception,message",
    [
        ("minao", "bad", ValueError, "invalid minao basis"),
        ("mo_basis", "bad", ValueError, "invalid MO basis"),
        ("pop_method", "bad", ValueError, "invalid population scheme"),
        ("mo_init", "bad", ValueError, "invalid MO start guess"),
        ("loc_exp", 3, ValueError, "invalid localization exponent"),
        ("part", "bad", ValueError, "invalid partitioning"),
        ("part_method", "bad", ValueError, "invalid partitioning method"),
        ("prop", "bad", ValueError, "invalid property"),
        ("unit", "bad", ValueError, "invalid unit"),
        ("verbose", -1, ValueError, "invalid verbosity"),
        ("verbose", 6, ValueError, "invalid verbosity"),
        ("ndo", "not_a_bool", TypeError, "invalid NDO argument"),
        ("write", "bad", ValueError, "invalid write format"),
        ("write", 123, TypeError, "invalid write format argument"),
        ("writename", 123, TypeError, "invalid write name argument"),
        ("verbose", "not_an_int", TypeError, "invalid verbosity"),
        ("unit", 123, TypeError, "invalid unit"),
    ],
)
def test_sanity_check_rejects_invalid_attr(attr, bad_value, exception, message):
    decomp = DecompCls(**{attr: bad_value})
    with pytest.raises(exception, match=message):
        sanity_check(None, None, decomp, None, None)


@pytest.mark.parametrize(
    "bad_gauge_origin,exception",
    [
        ([0.0, 0.0], ValueError),  # too short
        ([0.0, 0.0, 0.0, 0.0], ValueError),  # too long
        (["a", "b", "c"], ValueError),  # right length, wrong element type
        ("not a list", TypeError),  # wrong container type entirely
    ],
)
def test_sanity_check_gauge_origin_invalid(bad_gauge_origin, exception):
    decomp = DecompCls(gauge_origin=bad_gauge_origin)
    with pytest.raises(exception, match="invalid gauge origin"):
        sanity_check(None, None, decomp, None, None)


@pytest.mark.parametrize(
    "mo_coeff,exception,message",
    [
        ([1, 2, 3], TypeError, "invalid mo coefficients"),
        (np.zeros(5), ValueError, "invalid mo coefficients"),
        (
            (np.zeros((2, 2)),),
            TypeError,
            "invalid mo coefficients",
        ), 
    ],
)
def test_sanity_check_rejects_invalid_mo_coeff(mo_coeff, exception, message):
    decomp = DecompCls()
    with pytest.raises(exception, match=message):
        sanity_check(None, None, decomp, mo_coeff, None)


@pytest.mark.parametrize(
    "mo_occ,exception,message",
    [
        ("not valid", TypeError, "invalid mo occupation"),  # not ndarray or tuple
        (
            np.zeros((2, 2, 2)),
            ValueError,
            "invalid mo occupation",
        ),
        ((np.zeros(2),), TypeError, "invalid mo occupation"),  # wrong tuple length
    ],
)
def test_sanity_check_rejects_invalid_mo_occ(mo_occ, exception, message):
    decomp = DecompCls()
    mo_coeff = np.zeros((2, 2))
    with pytest.raises(exception, match=message):
        sanity_check(None, None, decomp, mo_coeff, mo_occ)


@pytest.mark.parametrize(
    "kwargs",
    [
        {},  # fully default
        {"verbose": 0},
        {"verbose": 5},
        {"gauge_origin": [0.0, 0.0, 0.0]},
        {"part": "orbitals"},
        {"part": "atoms", "part_method": "ao"},
        {"part": "eda"},  # deprecated, rewritten to atoms/ao
    ],
)
def test_sanity_check_accepts_valid_input(kwargs):
    # no assert needed -- passes by not raising; confirms valid input isn't over-rejected
    mo_coeff = np.zeros((2, 2))
    decomp = DecompCls(**kwargs)
    sanity_check(None, None, decomp, mo_coeff, None)


# tools.py


# dim
def test_dim():
    mo_occ = (np.array([1.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0]))
    alpha, beta = dim(mo_occ)
    assert np.array_equal(alpha, [0, 1])
    assert np.array_equal(beta, [0])


# mf_info
def test_mf_info_h2o(mf_h2o):
    mo_coeff, mo_occ = mf_info(mf_h2o)
    assert np.array_equal(mo_occ[0], np.ones(5))
    assert np.array_equal(mo_occ[1], np.ones(5))
    assert mo_coeff[0].shape[1] == 5
    assert mo_coeff[1].shape[1] == 5


def test_mf_info_li(mf_li):
    mo_coeff, mo_occ = mf_info(mf_li)
    assert np.array_equal(mo_occ[0], np.ones(2))
    assert np.array_equal(mo_occ[1], np.ones(1))
    assert mo_coeff[0].shape[1] == 2
    assert mo_coeff[1].shape[1] == 1


def test_mf_info_li_rohf(mf_li_rohf):
    mo_coeff, mo_occ = mf_info(mf_li_rohf)
    assert np.array_equal(mo_occ[0], np.ones(2))
    assert np.array_equal(mo_occ[1], np.ones(1))
    assert mo_coeff[0].shape[1] == 2
    assert mo_coeff[1].shape[1] == 1


# make_rdm1
def test_make_rdm1():
    mo = np.array([[1.0], [0.0]])
    occup = np.array([2.0])
    rdm = make_rdm1(mo, occup)
    assert np.array_equal(rdm, [[2.0, 0.0], [0.0, 0.0]])


def test_make_rdm1_equal_electron_count(mf_h2o):
    mo = mf_h2o.mo_coeff[:, :5]
    occup = mf_h2o.mo_occ[:5]
    S = mf_h2o.get_ovlp()
    D = make_rdm1(mo, occup)
    assert np.isclose(np.trace(D @ S), 10.0)


# testing functions in properties.py


# _e_nuc
def test_e_nuc_h2():
    mol = gto.M(atom="H 0 0 0; H 0 0 1.0", basis="sto-3g", unit="bohr", verbose=0)
    h2nuc = _e_nuc(mol)
    assert np.array_equal(h2nuc, [0.5, 0.5])


def test_e_nuc(mf_h2o):
    mol = mf_h2o.mol
    nuc_energy = _e_nuc(mol).sum()
    expected_nuc_energy = mol.energy_nuc()
    assert np.isclose(nuc_energy, expected_nuc_energy, atol=1e-10)


# _dip_nuc
def test_dip_nuc():
    mol = gto.M(atom="H 0 0 0", basis="sto-3g", unit="bohr", spin=1, verbose=0)
    gauge_origin = np.array([0.0, 0.0, 1.0])
    hdip_nuc = _dip_nuc(mol, gauge_origin)
    assert np.array_equal(hdip_nuc, np.array([[0, 0, -1]]))


def test_dip_nuc_h2o(mf_h2o):
    mol = mf_h2o.mol
    gauge_origin1 = np.array([0.0, 0.0, 1.0])
    dip_nuc1 = _dip_nuc(mol, gauge_origin1)
    gauge_origin2 = np.array([0.0, 1.0, 2.0])
    dip_nuc2 = _dip_nuc(mol, gauge_origin2)
    dip_nuc1 = dip_nuc1.sum(axis=0)
    dip_nuc2 = dip_nuc2.sum(axis=0)
    total_charge = mol.atom_charges().sum()
    assert np.allclose(
        dip_nuc1 - dip_nuc2,
        -1 * total_charge * (gauge_origin1 - gauge_origin2),
        atol=1e-10,
    )


# _h_core
def test_h_core(mf_h2o):
    mol = mf_h2o.mol
    kin, nuc, sub_nuc = _h_core(mol, mf_h2o)
    hcore = mf_h2o.get_hcore()
    assert np.allclose(hcore, kin + nuc, atol=1e-10)
    assert np.all(np.diag(kin) > 0.0)
    assert np.all(np.diag(nuc) < 0.0)


# _get_nuc
def test_get_nuc_matches_pyscf():
    mol = gto.M(atom="Li 0 0 0; H 0 0 1.0", basis="sto-3g", unit="bohr", verbose=0)
    sub_nuc = _get_nuc(mol)
    total = sub_nuc.sum(axis=0)
    assert np.allclose(total, mol.intor("int1e_nuc"), atol=1e-10)
    assert sub_nuc.shape == (mol.natm, mol.nao_nr(), mol.nao_nr())


# _point_charges
def test_point_charges():
    mol = gto.M(atom="Li 0 0 0; H 0 0 1", basis="sto-3g", unit="bohr", verbose=0)
    mm_mol = gto.M(atom="O 0 0 3; H 0 0 5", basis="sto-3g", spin=1, unit="bohr", verbose=0)
    mm_pot, nuc_solv = _point_charges(mol, mm_mol)
    assert nuc_solv.shape == (mol.natm,)
    assert np.allclose(nuc_solv, [3 * (8 / 3 + 1 / 5), 1 * (8 / 2 + 1 / 4)])
    ref = np.zeros_like(mm_pot)
    for coord, charge in zip(mm_mol.atom_coords(), mm_mol.atom_charges()):
        with mol.with_rinv_origin(coord):
            ref += -1.0 * mol.intor("int1e_rinv") * charge
    assert np.allclose(mm_pot, ref, atol=1e-10)


# _e_xc
def test_e_xc_1d():
    eps_xc = np.array([1.0, 2.0, 3.0])
    grid_weights = np.array([2.0, 1.0, 0.5])
    rho = np.array([1.0, 2.0, 3.0])
    assert _e_xc(eps_xc, grid_weights, rho) == 10.5


def test_e_xc_2d():
    eps_xc = np.array([1.0, 1.0, 1.0])
    grid_weights = np.array([1.0, 1.0, 1.0])
    rho = np.array([[1.0, 2.0, 3.0], [99.0, 99.0, 99.0]])
    assert _e_xc(eps_xc, grid_weights, rho) == 6.0


# _xc_ao_deriv
@pytest.mark.parametrize(
    "xc_func,expected",
    [
        ("lda,vwn", ("LDA", 0)),
        ("hf", ("HF", 0)),
        ("pbe", ("GGA", 1)),
        ("tpss", ("MGGA", 2)),
    ],
)
def test_xc_ao_deriv(xc_func, expected):
    assert _xc_ao_deriv(xc_func) == expected


# this pins current *buggy* behaviour on purpose: _xc_ao_deriv has no else branch,
# so an unrecognised xc_type leaves ao_deriv unassigned and the return line raises
# UnboundLocalError instead of a clear message. if someone later adds a proper
# `else: raise ValueError(...)`, this test SHOULD go red -- update it to expect
# ValueError, don't revert the fix.
# the mock is needed because no real functional can produce this: pyscf's xc_type
# only ever returns HF/LDA/GGA/MGGA/UNKNOWN, and UNKNOWN is unreachable in practice
def test_xc_ao_deriv_unknown_type():
    with patch("pyscf.dft.libxc.xc_type", return_value="UNKNOWN"):
        with pytest.raises(UnboundLocalError):
            _xc_ao_deriv("anything")


# _make_rho
def test_make_rho(mf_h2o):
    mol = mf_h2o.mol
    grids = dft.Grids(mol)
    grids.build()
    ao_value = numint.eval_ao(mol, grids.coords, deriv=0)
    rdm1 = mf_h2o.make_rdm1()
    c0, c1, rho = _make_rho(ao_value, rdm1, "LDA")
    rho_ref = numint.eval_rho(mol, ao_value, rdm1, xctype="LDA")
    assert np.allclose(rho, rho_ref, atol=1e-10)


def test_make_rho_gga(mf_h2o_dft):
    mol = mf_h2o_dft.mol
    grids = dft.Grids(mol)
    grids.build()
    ao_value = numint.eval_ao(mol, grids.coords, deriv=1)
    rdm1 = mf_h2o_dft.make_rdm1()
    c0, c1, rho = _make_rho(ao_value, rdm1, "GGA")
    rho_ref = numint.eval_rho(mol, ao_value, rdm1, xctype="GGA")
    assert np.allclose(rho, rho_ref, atol=1e-10)


# _make_rho_interm2 (atom-slicing, the decodense-specific reuse trick)
def test_make_rho_atom_slicing(mf_h2o):
    mol = mf_h2o.mol
    grids = dft.Grids(mol)
    grids.build()
    ao_value = numint.eval_ao(mol, grids.coords, deriv=0)
    rdm1 = mf_h2o.make_rdm1()
    c0, c1, rho_total = _make_rho(ao_value, rdm1, "LDA")
    ao_labels = mol.ao_labels(fmt=None)
    rho_sum = np.zeros_like(rho_total)
    for atom_idx in range(mol.natm):
        select = np.where([label[0] == atom_idx for label in ao_labels])[0]
        rho_atom = _make_rho_interm2(c0[:, select], None, ao_value[:, select], "LDA")
        rho_sum += rho_atom
    assert np.allclose(rho_sum, rho_total, atol=1e-10)


# _vk_dft
def test_vk_dft_pbe0(mf_h2o_dft):
    mol = mf_h2o_dft.mol
    rdm1 = 0.5 * mf_h2o_dft.make_rdm1()
    vj, vk = mf_h2o_dft.get_jk(mol=mol, dm=rdm1)
    vk_result = _vk_dft(mol, mf_h2o_dft, "pbe0", rdm1, vk, vj)
    # PBE0 is defined as exactly 25% exact exchange -- not looked up via any
    # decodense or PySCF code, just the functional's own definition
    expected = 0.25 * mf_h2o_dft.get_k(mol=mol, dm=rdm1)
    assert np.allclose(vk_result, expected, atol=1e-10)


# _trace
def test_trace_identity():
    op = np.array([[1.0, 2.0], [3.0, 4.0]])
    rdm1 = np.eye(2)
    assert _trace(op, rdm1) == 5.0


def test_trace_symmetric():
    op = np.array([[1.0, 2.0], [3.0, 4.0]])
    rdm1 = np.array([[2.0, 1.0], [1.0, 3.0]])
    assert _trace(op, rdm1) == np.trace(op @ rdm1)


def test_trace_3d():
    op = np.array(
        [
            [[1.0, 0.0], [0.0, 0.0]],
            [[0.0, 0.0], [0.0, 1.0]],
            [[2.0, 0.0], [0.0, 3.0]],
        ]
    )
    rdm1 = np.eye(2)
    expected = np.array([1.0, 1.0, 5.0])
    assert np.array_equal(_trace(op, rdm1), expected)



# testing functions in orbitals.py


# _population_mul
def test_population_mul_h2o(mf_h2o):
    mol = mf_h2o.mol
    ovlp = mf_h2o.get_ovlp()
    ao_labels = mol.ao_labels(fmt=None)
    mo = mf_h2o.mo_coeff[:, :5]
    mocc = mf_h2o.mo_occ[:5]
    overlap_mo = np.einsum("ji,jp->ip", ovlp, mo)
    pop = mocc[None, :] * mo * overlap_mo
    populations = _population_mul(mol.natm, ao_labels, pop)
    assert np.isclose(populations.sum(), 10.0)
    # symmetry-equivalent
    per_atom = populations.sum(axis=0)
    assert np.isclose(per_atom[1], per_atom[2])
    assert per_atom[0] > per_atom[1]


# assign_rdm1s
def test_assign_rdm1s_h2o(mf_h2o):
    mo_coeff, mo_occ = mf_info(mf_h2o)
    weights = assign_rdm1s(
        mf_h2o.mol, mf_h2o, mo_coeff, mo_occ, "MINAO", "mulliken", False, 0
    )
    assert len(weights) == 2
    assert np.array_equal(weights[1], weights[0])
    assert np.allclose(weights[0].sum(axis=1), 1.0)


def test_assign_rdm1s_h2o_iao(mf_h2o):
    mo_coeff, mo_occ = mf_info(mf_h2o)
    weights = assign_rdm1s(
        mf_h2o.mol, mf_h2o, mo_coeff, mo_occ, "MINAO", "iao", False, 0
    )
    assert len(weights) == 2
    assert np.array_equal(weights[1], weights[0])
    assert np.allclose(weights[0].sum(axis=1), 1.0)


def test_assign_rdm1s_li(mf_li):
    mo_coeff, mo_occ = mf_info(mf_li)
    weights = assign_rdm1s(
        mf_li.mol, mf_li, mo_coeff, mo_occ, "MINAO", "mulliken", False, 0
    )
    assert len(weights) == 2
    assert np.allclose(weights[0].sum(axis=1), 1.0)
    assert np.allclose(weights[1].sum(axis=1), 1.0)
    assert weights[0].shape[0] == 2
    assert weights[1].shape[0] == 1


def test_assign_rdm1s_li_iao(mf_li):
    mo_coeff, mo_occ = mf_info(mf_li)
    weights = assign_rdm1s(
        mf_li.mol, mf_li, mo_coeff, mo_occ, "MINAO", "iao", False, 0
    )
    assert len(weights) == 2
    assert np.allclose(weights[0].sum(axis=1), 1.0)
    assert np.allclose(weights[1].sum(axis=1), 1.0)
    assert weights[0].shape[0] == 2
    assert weights[1].shape[0] == 1


# partial charges from the iao population weights.
# expected signs are chemistry, not code: oxygen is the more electronegative atom,
# so it draws electron density and comes out negative, both hydrogens positive, and
# the two hydrogens equal by symmetry. the total-charge line is implied by the
# partition of unity asserted above (10 spin-orbitals, each summing to 1, vs sum(Z)
# = 10) and is kept as documentation rather than as a load-bearing assert
def test_rdm1_charge_conservation(mf_h2o):
    mol = mf_h2o.mol
    mo_coeff, mo_occ = mf_info(mf_h2o)
    weights = assign_rdm1s(mol, mf_h2o, mo_coeff, mo_occ, "MINAO", "iao", False, 0)
    population = np.sum(weights[0], axis=0) + np.sum(weights[1], axis=0)
    charge_atom = mol.atom_charges() - population
    total_charge = charge_atom.sum()
    assert np.isclose(total_charge, 0.0)
    assert charge_atom[0] < 0.0
    assert charge_atom[1] > 0.0
    assert charge_atom[2] > 0.0
    assert np.isclose(charge_atom[2], charge_atom[1])


# testing functions in results.py


# atoms
def test_atoms():
    mol = gto.M(
        atom="H 0 0 0; O 0 0 1.0", basis="sto-3g", spin=1, unit="bohr", verbose=0
    )
    res = {
        CompKeys.el: np.array([1.0, 4.0]),  # made-up numbers, one per atom
    }
    df = atoms(mol, res, "au")
    assert len(df) == 2
    assert list(df.index) == ["H0", "O1"]
    assert np.allclose(df[CompKeys.el], [1.0, 4.0])


# to_dataframe
def test_to_dataframe(mf_h2o):
    mol = mf_h2o.mol
    mo_coeff, mo_occ = mf_info(mf_h2o)
    decomp1 = DecompCls(pop_method="iao", part="atoms")
    res1 = main(mol, decomp1, mf_h2o, mo_coeff, mo_occ)
    decomp2 = DecompCls(pop_method="iao", part="orbitals")
    res2 = main(mol, decomp2, mf_h2o, mo_coeff, mo_occ)
    assert len(res1.to_dataframe()) == 3
    assert len(res2.to_dataframe()) == 10
    # __str__ is the path every example uses via print(res), and is where the
    # leftover charge_atom reference crashed. nothing else in the suite runs it
    assert "O0" in str(res1)
    assert len(str(res2).splitlines()) > 10


# orbs --  NDO 
def test_orbs_ndo():
    mol = gto.M(atom="H 0 0 0", basis="sto-3g", spin=1, unit="bohr", verbose=0)
    res = {
        CompKeys.el: [np.array([1.0, 2.0, 3.0, 4.0, 5.0]), np.array([])],
        CompKeys.mo_occ: (np.array([-0.8, -0.3, 0.0, 0.3, 0.8]), np.array([])),
        CompKeys.orbsym: (np.array(["A"] * 5, dtype=object), np.array([], dtype=object)),
    }
    df = orbs(mol, res, "au", ndo=True)
    assert len(df) == 5
    assert df.iloc[-1][CompKeys.el] == 3.0
    df1 = orbs(mol, res, "ev", ndo=True)
    assert np.isclose(df1.iloc[-1][CompKeys.el], 3.0 * 27.211386245988)


# DecompCls 
def test_decomp_cls_part_method_defaults():
    decomp_eda = DecompCls(part="eda")
    assert decomp_eda.part == "atoms"
    assert decomp_eda.part_method == "ao"
    decomp_atoms = DecompCls(part="atoms")
    assert decomp_atoms.part_method == "mo"
    decomp_orb = DecompCls(part="orbitals")
    assert decomp_orb.part_method is None


# pbctools.py (no tests yet) 