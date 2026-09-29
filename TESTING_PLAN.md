# decodense testing plan

Test suite for **mainline decodense** (the AD branch is out of scope — mainline first; these
tests become the reference the AD path must later agree with). Detailed scoping of **PR 12**,
companion to `PLAN.md`.

**[verified]** = actually measured against this repo. **[hypothesis]** = worth a test, not a
confirmed bug.

## Status at a glance

| | |
|---|---|
| Tests written | 43 test functions = 72 test cases, in `tests/test_unit.py`, all passing on `adapt-tests-to-restructuring` (PR #28 + 2 local fixes, see Part VII) |
| Source functions with a verdict | all of them. Part II, re-audited against PR #28 in Part VII, **re-audited again with measured coverage + mutation testing in Part VIII** |
| Original 23-item plan | 15 done, 1 partial, 3 declined with reasons, 4 untouched |
| Open work | Part VII's two remaining new tests (4b, 5) are **done** (N4, N1 below). **Now at Part VIII.7**, going through the rest of the decision list: T1–T30 (per-test review), then N2–N13. |
| Mutation score | Of 38 planted bugs: unit tests catch 18, legacy tests catch 5 more, **15 get past everything** (Part VIII.3) |
| Weakest coverage | Anything that moves energy *between atoms* without changing the total. The legacy tests check only totals, and the unit tests have almost no per-atom value checks (VIII.1 #5) |

## How to use this document

**Working session → start at Part VIII.7** (decision list); Part VIII holds the reasoning,
Part VII the PR #28 adaptation, Part I the older worklist (partly stale, see VIII.8).
Everything else is reference you pull from as needed.

| Part | What it is | When you need it |
|---|---|---|
| **I** | The worklist: A1-A13 (fix existing tests), B1-B19 (write new ones) | Every session — this is the agenda |
| **II** | Every function in `decodense/`: tested?, should it be?, how? | Deciding whether/how to test something |
| **III** | Every one of the 38 existing tests, with a verdict | Understanding why a test is weak or strong |
| **IV** | Mechanics: pytest, fixtures, parametrize, tolerances | While actually writing a test |
| **V** | Why we're doing this; how PySCF/QuTiP/Qiskit/OTR do it | Background, read once |
| **VI** | Decisions, history, corrections, superseded plans | Checking why something was declined |
| **VII** | PR #28 restructuring: function re-audit, new tests, updates | Adapting the suite to PR #28 |
| **VIII** | Critical re-audit: measured coverage, mutation testing, per-test verdicts, missing tests | **Current agenda**, start at VIII.7 |

**One convention worth knowing**, because it decides most arguments in here: *a test is only
worth its line count if the expected answer comes from outside the code under test* — physics
derived by hand, a genuinely separate implementation, or an invariant that must hold however
the code is written. A test that recomputes the function's own steps and compares verifies
wiring, not correctness.

---

# Part I — The worklist

This is the agenda. Each row points back to the reasoning in Part II or III.

#### Combined checklist — everything Parts II and III turned up

Two lists, because they are different kinds of work. **A** = fixing tests that already exist
(from the Part III review). **B** = writing tests that do not exist yet (from the Part II audit).
Each item names where its reasoning lives.

#### A. Fixes to existing tests in `test_unit.py`

| # | Fix | Why | Ref |
|---|---|---|---|
| A1 | Re-enable `comp_key_dict` as `@pytest.mark.xfail` — **declined by user, final answer.** Was actually built and passing earlier in this session, then lost to an editor-autosave-reverting-a-stale-buffer issue; a redo was attempted twice more and rejected both times, and the user then explicitly said they don't want this test at all. Stays commented out. `comp_key_dict`'s bug therefore still has zero coverage anywhere in the repo -- this is a known, accepted gap, not an oversight. | Guards an already-shipped bug and currently has **zero coverage anywhere** since `test_comp_keys.py` was deleted | Part III's table #1, consolidation note |
| A2 | `test_e_xc_1d`: use non-trivial `eps_xc`/`grid_weights` instead of all-ones | Currently blind to the function silently dropping an argument | Part III, weakness #1 |
| A3 | Replace `test_make_rdm1_symmetric` with `Tr(D·S) == sum(occup)` | Symmetry is true by construction for wrong implementations too; also closes the "no multi-orbital value check" gap | Part III, weakness #2, table #12 |
| A4 | ~~`_trace`: non-symmetric `rdm1` case~~ — **decided: no test needed** | Discussed live: every density matrix decodense builds is `Σ n_p c_p c_p^T`, an outer product, which is symmetric *by construction* for any real-valued coefficients. A non-symmetric `rdm1` cannot arise from any real decodense computation, so a test for one wouldn't simulate a realistic bug — it would only probe `_trace` as an abstract, disconnected function. Left as-is. | Part III, weakness #3 |
| A5 | `test_mf_info_h2o`: check `mo_coeff[1]` too — **done**. The dtype assertion was discussed and **declined**: `mf_info` has zero internal callers, so the int64-vs-float64 quirk can't break anything inside decodense; not worth the line given that low blast radius, even though it would have closed an open question from the original plan. | Part III, weakness #4, #9 |
| A6 | `test_xc_ao_deriv_unknown_type`: add a comment saying it pins a known bug on purpose | Otherwise a future fixer sees red and "fixes the test" by reverting the improvement | Part III, weakness #5 |
| A7 | `test_h_core`: label the `nuc == sub_nuc.sum(axis=0)` line as a consistency check, not verification | It re-runs the function's own line and cannot fail | Part III, weakness #6 |
| A8 | `test_charge_conservation` — **done properly.** First added just a comment flagging the sum-only limitation, then came back and added the real per-atom checks: oxygen negative, both hydrogens positive, and the two hydrogens close to each other via `np.isclose` (not `==`, since these are computed floats). | A sum of zero survives any bug that moves charge *between* atoms | Part III, weakness #7 |
| A9 | `test_vk_dft_pbe0`: pass a deliberately perturbed `rdm1` — **done**, `rdm1 = 0.5 * mf_h2o_dft.make_rdm1()`, one-line change since every other line already reused the same variable. Confirmed passing. | Currently does not confirm the function uses the density matrix it was handed | Part III, weakness #8 |
| A10 | ~~Add the `"NLC"` case to `_xc_ao_deriv`~~ — **revised: impossible with a real functional.** Checked PySCF's own source (`libxc.py:_xc_type`): the underlying C call can only ever return one of `('HF', 'LDA', 'GGA', 'MGGA', 'UNKNOWN')` — `"NLC"` is never a real return value. VV10-type nonlocal correlation is tracked by PySCF through a separate `mf.nlc` attribute, not through `xc_type()` at all. So this branch is dead code exactly like `"UNKNOWN"` — only reachable by mocking `dft.libxc.xc_type` to return a value it can't really produce. Left undone; same shape as `test_xc_ao_deriv_unknown_type` if ever added. | Part III's table #6-9 |
| A11 | Decide whether `verbose=True` should be accepted, and pin it — **deferred, not decided**. Discussed live: recommendation was to pin current behavior (bool coerces to 0/1, both already valid levels, so it's benign) rather than change `decomp.py`'s source. User wanted more time before deciding; revisit later. | `bool` subclasses `int`, so it currently passes the type check silently | Part III's table #35 |
| A12 | Fold `test_sanity_check_invalid_pop_method` into the parametrized table as a row — **done**, confirmed on disk. | One attribute in a standalone function while nine live in a table | Part III's table #14 |
| A13 | One-line comments on the three assertion-free "should not raise" tests | They read as empty tests to anyone unfamiliar with the pattern | Part III, "should not raise" |

#### B. New tests worth writing

**High value, low cost — do these first:**

| # | Test | Why |
|---|---|---|
| B1 | `orbs()` NDO odd-orbital count (`results.py`) — **skipped, not declined.** Verified the fix already exists, unmerged, on `fix-comp-key-dict` (same branch as the `comp_key_dict` fix): after the pairing loop, `if sort_idx.size % 2 == 1: mo_idx = np.append(mo_idx, sort_idx[sort_idx.size // 2])`. User plans to merge that branch soon, so an `xfail` test against `main` now would be short-lived churn. Revisit if the merge is delayed. | Known **shipped** bug; no SCF needed; the integration tests are structurally blind to it (they assert on `res.tot`, the bug is in the DataFrame built afterwards) |
| B2 | `assign_rdm1s` on `mf_li` (unrestricted) | Fixture already built and unused; without it, `weights[1] == weights[0]` would pass even if the code *always* duplicated |
| B3 | `ResultsCls.__init__` attribute round-trip — **skipped for now**, same reasoning as B1: would currently fail against `main` for the same `comp_key_dict` mismatch, and the user is about to merge that fix. Unlike B1, this test's value is broader than just that one bug (it guards `ResultsCls`'s general attribute-wiring, not just today's specific mismatch), so worth revisiting after the merge. | Would have caught the original `comp_key_dict` bug directly |
| B4 | `_e_nuc(mol).sum() == mol.energy_nuc()` — **skipped after re-examination.** Checked PySCF's `mol.energy_nuc()` source: it uses the exact same `inter_distance` helper, the exact same diagonal-masking trick, and the exact same einsum formula as `_e_nuc`, just fully summed instead of keeping the per-atom axis. Not an independent implementation the way `_get_nuc` vs `int1e_nuc` is — summing `_e_nuc`'s output is close to guaranteed to reconstruct PySCF's total by simple arithmetic. Downgraded from 'strong differential' to 'weak consistency check', same category as `_h_core`'s labelled circular line. Not written. | One line, independent PySCF reference |

**Medium value, low cost:**

| # | Test | Why |
|---|---|---|
| B5 | `_dip_nuc` gauge-origin shift property | Tests that the origin is actually *used*, which one hand-computed case cannot |
| B6 | `atoms()` structural test from a hand-built dict | Row count, index labels, columns, no NaN, unit scaling applied |
| B7 | `_solvent` dispatch via dummy objects | Pure decodense branching; needs no real PCM/QM-MM object |
| B8 | `ewald_e_nuc(cell).sum() == cell.ewald()` — **downgraded, same issue as B4.** `ewald_e_nuc`'s own docstring says "adapted from: pbc/gto/cell.py:ewald() from PySCF v2.1" — not independent code, same weak-consistency-check category, not the strong differential this row originally claimed. Lower priority than listed; revisit the whole PBC row in Part VI's original scoping too, since it was based on the same now-corrected assumption. | Independent PySCF reference; needs only a `Cell` |
| B9 | Parametrise the invariant tests over `part` in `("atoms", "eda")` | `prop_eda` is currently reached only by the old integration tests |
| B10 | `_get_nuc` shape + per-slice symmetry | One-liners the plan already asked for |
| B11 | `dim` edge case at the strict `> 0.0` boundary — **skipped after checking real call sites.** `dim()` is called from `properties.py`/`orbitals.py`/`results.py`/`tools.py`, but for NDO specifically, the caller (e.g. `test_ch2_hf_ndo.py`) pre-filters occupations with `make_natorb`'s own `NATORB_THRES = 1e-12` *before* ever calling `main()` — so `dim()` itself basically never actually receives a genuinely tiny-but-nonzero value in real usage; that decision already happened one level up, via a separate mechanism. Same pattern as B4/B8: justification didn't survive checking. | NDO produces small fractional occupations |

**Medium value, medium cost:**

| # | Test | Why |
|---|---|---|
| B12 | MGGA for `_make_rho`; range-separated branch for `_vk_dft` — **declined by user.** Verified `wb97m_v` is genuinely MGGA + range-separated (`omega=0.3`, `hyb=0.15`, `alpha=1.0`), and both branches are technically decodense's own code (not pure PySCF delegation) -- but the user judged that the real numerical work in both cases is PySCF's (`get_k(omega=...)`, `eval_rho(xctype="MGGA")`), with decodense mostly just combining/scaling those outputs, and preferred to skip on that basis. | Both are branches real examples actually use (`wb97m_v`) |
| B13 | `_make_rho_interm2` per-atom slicing (`Σ_atoms rho_atom == rho_total`) — **done for LDA** (`test_make_rho_atom_slicing`). A GGA variant was verified to work too (confirmed `_make_rho_interm2`'s LDA and GGA branches are separate `if/elif` code, including the `*2.0` symmetry-doubled gradient terms — an LDA-only test can never run that branch at all), but the user chose to keep just the one LDA test rather than add a near-duplicate; the GGA branch of this specific invariant is a known, deliberately accepted gap, not an oversight. | The decodense-specific reason the function is split in two |
| B14 | `_make_rho` unrestricted branch (`rdm1.ndim == 3`) — **declined by user.** Checked PySCF's real `eval_rho` source: not a copy like `_e_nuc` (PySCF uses its own internal `_dot_ao_dm`/`_contract_rho`, decodense uses generic `contract`), so a differential test against it is valid. But the core LDA/GGA formula was already validated by the restricted-case tests -- this test's real incremental value narrows to just "does the alpha/beta branch correctly avoid swapping the two spins," confirmed real (would silently corrupt XC energy for every open-shell system) and confirmed the test does catch a swap (walked through mechanically). User judged that narrower scope not worth it. | Two further code paths behind `np.allclose(rdm1[0], rdm1[1])` |
| B15 | `orbsym` fallback branches — **skipped.** Verified the fallback genuinely triggers (molecules built without `symmetry=True` really do raise inside `symm.label_orb_symm`). Discussed the real consequence: not silently-wrong physics (the output is display-only, the `"Symm."` column, never used in calculation), but a full crash if the fallback itself were broken. Also clarified this is a different category from `sanity_check` (graceful degradation of an internal PySCF call, not upfront user-input validation) -- user decided to move on without writing it. | Three near-duplicate shape-handling branches, all untested |
| B16 | `write_rdm1`'s two guard raises — **declined by user**, judged not necessary. | Called by `main()`; the raises are pure decodense |
| B17 | `DecompCls.__init__` fresh `gauge_origin` per instance — **skipped by user.** Confirmed the fix is present on this branch (`np.zeros(3)` built fresh in `__init__`, not a shared default), so this would have been a regression guard, not documenting an active bug. | Regression for the PR 2 shared-mutable-default fix |
| B18 | `logger_config` vs `sanity_check` range agreement — **skipped, user judged not plausible.** The class of bug is real (proven by the `comp_key_dict` precedent), but this specific pair sits in small, stable, rarely-touched config code, unlike the actively-developed area `comp_key_dict` lived in -- lower real-world likelihood than most of what got tested today. | If they drift, `main()` dies with `KeyError` on its first line |
| B19 | `_point_charges` nuclei-MM Coulomb loop — **done, and found a real bug in the process.** `nuc_solv = np.zeros(len(mol.atom))` used the character-count of the raw input string instead of `mol.natm` -- array ended up oversized with unused trailing zeros (harmless in practice, since only the first `natm` entries are ever read downstream, but genuinely wrong sizing). Written first as `xfail` to confirm and document the bug, then **fixed on this branch** (`len(mol.atom)` -> `mol.natm`, decodense/properties.py:609) after discussing the one-PR-one-topic convention with the user -- user chose to fix inline rather than defer to a separate branch. `xfail` marker removed, test now genuinely passes (`test_point_charges_nuc_solv_shape`). | Hand-computable analytically for one QM atom + one point charge |

**Explicitly not recommended**, each for the reason in its Part II row: `contract`, `_ao_val`,
`fmt`, `to_dataframe`, `__str__`, `info`, `git_version`, the logger methods, `_pcm`,
`make_natorb`, `res_add`/`res_sub`, `_population_becke`, and the vendored `pbctools.py`
builders and internals.

---

# Part II — Reference: every function

### Every function in `decodense/`, by file

Every function in `decodense/`, with: current coverage, whether it *should* be tested, and if
so a concrete assertion strategy. Ordered by file.

**The filter applied throughout:** many of these functions are wrappers around PySCF. We are
not in the business of re-testing PySCF's correctness. For each one the question is
specifically *what logic is decodense's own* — partitioning, summation, bookkeeping,
dispatch, formatting, edge-case handling — and only that part is proposed for testing. Where
a function is a thin pass-through, it is called out as such and skipped.

Legend: ✅ tested · ⚠️ partially tested · ❌ untested · ⛔ deliberately not worth testing

#### `decomp.py`

| Function | State | Assessment |
|---|---|---|
| `CompKeys` | ⛔ | A namespace of string literals. Nothing to break *in itself*. The thing that *did* break is its relationship to `comp_key_dict`, which is item #1's test — **currently commented out in `test_unit.py`**. That should be re-enabled (as `@pytest.mark.xfail` with a reason, so it runs and reports rather than disappearing) or the bug fixed. A commented-out test is the worst of both worlds: no coverage, and no visible red. |
| `DecompCls.__init__` | ❌ | Plain attribute assignment with one exception that is real decodense logic: `gauge_origin` builds a **fresh** `np.zeros(3)` per instance. That is the PR 2 fix for a shared-mutable-default bug. Two-line regression test (item #2): two instances, mutate one's `gauge_origin`, assert the other is unchanged. Currently deprioritised; cheap enough that it's worth reconsidering. |
| `sanity_check` | ✅ | Now well covered: 9 `ValueError` rows, 5 `TypeError` rows, `gauge_origin` shape/type, `verbose` boundaries both sides, `mo_coeff`/`mo_occ` validation, and a valid-input case. Pure decodense logic, no PySCF involved. Remaining gap: the 4 PBC-specific checks (need a real `Cell` + `mf.kpt`), deliberately skipped. |

#### `tools.py`

| Function | State | Assessment |
|---|---|---|
| `DecodenseLogger.info2/info3` | ⛔ | Thin wrapper over stdlib `logging._log` at a custom level. Testing it tests stdlib logging. |
| `logger_config` | ❌ | Maps `verbose` 0–5 to logging levels via a 6-entry dict. Trivial on its own — **but there is a real cross-file invariant worth one test**: `sanity_check` independently validates `0 <= verbose <= 5`, and `logger_config`'s dict has exactly those 6 keys. If the two ever drift (someone widens the allowed range without adding a key), the result is a `KeyError` at runtime, in `main()`, on the first line. Assert the two agree by deriving one from the other rather than hand-copying — this is Qiskit's introspective-sweep idea (item #4) applied to something concrete. |
| `git_version` | ⛔ | Subprocess + environment plumbing, returns `"Unknown"` on any failure. Not where a physics bug hides. |
| `dim` | ✅ | Tested. One edge case worth adding (item #5): the comparison is `np.abs(mo_occ) > 0.0`, strictly — so an occupation of exactly `0.0` is excluded but `1e-300` is included. NDO produces small fractional occupations, so pinning this boundary is meaningful, not hypothetical. |
| `mf_info` | ⚠️ | Tested for restricted + unrestricted. **Zero internal callers** (verified by grep) — exported in `__init__.py` but never used by decodense itself. Keep the tests, don't extend, *except* adding the dtype assertion (see Part III finding #4) since that documents a real open question. |
| `orbsym` | ❌ | The `symm.label_orb_symm` call is PySCF's. **The fallback is decodense's own and is the interesting part**: on any exception it returns `["A"] * n` as a `dtype=object` array, and there are *three* near-duplicate branches (2-D ndarray, 3-D ndarray, tuple) each with their own shape handling. A shape bug in one branch is plausible and would be invisible today. Cheap test, no symmetry needed: build a mol with `symmetry=False`, assert the result is all `"A"`, right shape, `dtype=object`, for both a 2-D `mo_coeff` and a tuple of two. |
| `make_rdm1` | ⚠️ | Tested; replace the weak symmetry assertion with `Tr(D·S) == sum(occup)` (Part III finding #2). |
| `make_natorb` | ❌ | Zero internal callers (verified). Mostly PySCF (`intor_symmetric`, `np.linalg.eigh`). Decodense's own: the reshape branches (2-D `mo_coeff` → duplicated across spins; 2-D `rdm1` → `[rdm1, rdm1] * 0.5`) and the `thres` masking that drops insignificant NOs. **Recommend leaving untested** given no callers — but if ever needed, test the masking (occupations below `thres` are dropped) and that occupations sum to the electron count, not the eigendecomposition itself. |
| `write_rdm1` | ❌ | **Is** called by `main()` (when `decomp.write != ""`). The cube/npz writing is PySCF/numpy. Decodense's own: the two guard raises, the per-atom weighted rdm1 accumulation, and filename construction. **Worth testing the cheap half**: both raises (`part != "atoms"`, `fmt` not in `("cube","numpy")`) are pure decodense and cost two lines each. Optionally the numpy path into pytest's `tmp_path` fixture, asserting one array per atom with the expected key names and shapes. Don't test cube output. |
| `res_add` / `res_sub` | ❌ | Zero internal callers; previously declined. If revisited, test that mismatched keys raise `ValueError` and that `"Symm."` is handled specially — **and fix the hardcoded literal `"Symm."` (`tools.py:291,306`) to `CompKeys.orbsym` at the same time**, since that is the identical hand-maintained-string hazard that caused the original `comp_key_dict` bug. |
| `contract` | ⛔ | Dispatch between `opt_einsum` and `np.einsum`. Asserting the two agree tests opt_einsum, not decodense. |

#### `properties.py`

| Function | State | Assessment |
|---|---|---|
| `prop_tot` | ⛔ (covered indirectly) | Not unit-testable: it defines `prop_atom`/`prop_eda`/`prop_orb` as closures over its own locals and uses `global`. Covered by the 7 integration tests, `test_permute`, and the invariant tests. |
| `prop_atom` / `prop_eda` / `prop_orb` | ⛔ (covered indirectly) | Literally uncallable in isolation. **But a real gap exists in how they're reached**: every invariant test we wrote uses `part="atoms"`, so `prop_eda` and `prop_orb` are exercised only by the old integration tests. Parametrising `test_charge_conservation` / `test_nuc_sum` over `part` in `("atoms", "eda")` is nearly free and covers a second kernel. |
| `_e_nuc` | ⚠️ | Analytic H₂ case done. **Missing the differential half item #13 asked for**: `_e_nuc(mol).sum() == mol.energy_nuc()`. One line, independent PySCF reference, works on any molecule. Add it. |
| `_dip_nuc` | ⚠️ | Hand-computed case done. **Missing the gauge-origin property** (item #15): shifting the origin by `d` must shift the total by exactly `−Z_total·d`. Pure property, no reference values, and it tests the one thing the single hand-computed case can't — that the origin is actually *used* correctly rather than coincidentally. |
| `_h_core` | ✅ | Good, with the one circular line noted in Part III, weakness #6. |
| `_get_nuc` | ⚠️ | Sum differential done. Item #14 also asked for shape `(natm, nao, nao)` and per-slice symmetry — both one-liners, both still missing. |
| `_solvent` | ❌ | **Cheaper to test than previously assumed, and worth revisiting.** The PCM/QM-MM/OpenMM machinery is external, but `_solvent` itself is pure decodense *dispatch*: `hasattr(mf, "mm_mol")` → point charges, `hasattr(mf, "with_solvent")` → PCM, `hasattr(mf, "h1e_mmpol")` → OpenMM, else `(None, None, None)`. That dispatch is testable **without building any real solvent object** — hand it a dummy object carrying only the relevant attribute and assert the right branch fires. The `(None, None, None)` path for a plain `mf` is a one-liner. Recommended. |
| `_point_charges` | ❌ | Mostly PySCF (`df.incore.aux_e2`, `fakemol_for_charges`, `lib.unpack_tril`). Decodense's own: the blocked accumulation loop and the nuclei–MM Coulomb loop. **The latter has a genuine analytic reference**: `nuc_solv[j] = Z_j · Σ_k q_k/r_jk` is plain Coulomb, hand-computable for one QM atom and one point charge. Moderate setup (a fake `mm_mol`), real value — this is decodense's own arithmetic, not PySCF's. |
| `_pcm` | ⛔ | **Revised: this was wrong, and it's testable.** Originally judged as needing a real, converged `solvent.PCM` object since the substantive maths (`_get_vind`, K/R linear solves) is PySCF's. **Verified otherwise**: a fully mocked `solvent_model` (`unittest.mock.MagicMock`, same trick already used for `_xc_ao_deriv_unknown_type`) with small hand-picked numbers for `surface`, `_get_vind`, `_intermediates["K"]/["R"]`, `_get_v`, `v_grids_n` drives `_pcm` through its real code with no actual PCM calculation needed, and the result (`vmat_e = 0.5 * mocked_value`, `nuc_solv_pcm` from the `q_sym = (q+qt)/2` symmetrization) is fully hand-derivable — confirmed by running it (`vmat_e=5.0` from `0.5*10.0`, `nuc_solv_pcm=0.125`). **Left untouched on purpose**: the user doesn't have much footing in PCM/solvent theory yet and didn't want to add a test they couldn't independently follow, not because it's infeasible. Revisit if/when that changes. |
| `_xc_ao_deriv` | ✅ | All 5 cases. See Part III, weakness #5 about labelling the bug-pinning test. |
| `_make_rho_interm1` / `_make_rho_interm2` | ⚠️ | Covered indirectly through `_make_rho` for LDA and GGA. **Two real gaps: (a) MGGA is untested**, and `wb97m_v` — used by a real integration test — is a meta-GGA, so this is the same "the untested branch is the one actually used" pattern GGA was. **(b) The per-atom slicing reuse is now tested for LDA** (`test_make_rho_atom_slicing`) — `_make_rho_interm2(c0[:, select], ..., ao_value[..., select])` is the entire reason these are split into two functions, and it's decodense-specific (PySCF's `eval_rho` has no equivalent); the property checked is that summing `rho_atom` over all atoms reproduces the total `rho` on every grid point. The GGA branch (a genuinely different code path, not just different numbers) was verified to satisfy the same property but deliberately left untested — see Part I, B13. |
| `_make_rho` | ⚠️ | LDA + GGA done. Missing MGGA, and missing the **unrestricted** branch (`rdm1.ndim == 3`, where `np.allclose(rdm1[0], rdm1[1])` selects between two further code paths) — `mf_li` could drive it. |
| `_vk_dft` | ⚠️ | PBE0 (global hybrid) done. **The range-separated branch (`ks_omega != 0`) is untested** — and `wb97m_v` is exactly that. Same pattern a third time. Reference: `mf.get_k(mol, dm, omega=ks_omega)` scaled by `ks_alpha - ks_hyb`. |
| `_ao_val` | ⛔ | One-line `isinstance` dispatch to `numint.eval_ao` / `pbc_numint.eval_ao`. Nothing of decodense's own. |
| `_trace` | ⚠️ | Tested; transposition orientation untested by accident (Part III, weakness #3). |
| `_e_xc` | ⚠️ | Tested; the 1-D case is blind to a dropped `grid_weights` (Part III, weakness #1). |

#### `orbitals.py`

| Function | State | Assessment |
|---|---|---|
| `assign_rdm1s` | ⚠️ | Mulliken + IAO done, **restricted only**. Item #17 explicitly wanted restricted *and* unrestricted, and `mf_li` (UHF) is already built and sitting unused. **This is the single cheapest high-value addition in the audit**: run it on `mf_li`, assert partition-of-unity holds for alpha and beta *separately* (different sizes: 2 and 1), and assert `weights[1]` is **not** a copy of `weights[0]` — which is what finally gives the RHF-shortcut test its meaning, by contrast. Without an unrestricted case, `assert np.array_equal(weights[1], weights[0])` could pass even if the code *always* duplicated, which would be a serious bug. |
| `get_weights` | ⛔ | A closure inside `assign_rdm1s`, capturing `pop_method`/`natm`/`ao_labels`/`ovlp`. Never returned or exposed — cannot be imported or called. Exercised via its parent; that is the only possible route. |
| `_population_mul` | ✅ | Charge conservation + O/H symmetry. Solid. |
| `_population_becke` | ❌ | Reachable (via `pop_method="becke"`) but **zero examples use it**, and the existing integration files' `POP_METHOD` loop already covers it end-to-end. If tested: same partition-of-unity property at `1e-6`, not `1e-12` — **[verified]** earlier at 1.6e-7 deviation, because it is grid-integrated. Explicitly do **not** cross-compare against `_population_mul`: different definitions, not expected to agree. Low priority. |

#### `results.py` — the one file with zero coverage, and it contains a known bug site

| Function | State | Assessment |
|---|---|---|
| `ResultsCls.__init__` | ❌ | **High value, low cost.** Its one line of real logic is `setattr(self, comp_key_dict[key], value)` — which is *exactly* where the original `comp_key_dict` bug surfaced as `res.mo_occ` raising `AttributeError`. Build a `ResultsCls` from a real `res_dict` and assert every expected attribute exists and holds the right array. This would have caught the original shipped bug directly, and needs no SCF. |
| `__str__` | ⛔ | Delegates to `atoms()`/`orbs()`; formatting only. |
| `to_dataframe` | ⛔ | One-line delegation to `fmt`. Covered by testing `atoms`/`orbs`. |
| `info` | ⛔ | Pure string formatting, already verified byte-for-byte during its f-string rewrite. |
| `fmt` | ⛔ | One-line branch on whether `charge_atom` is present. Covered by testing the two it dispatches to. |
| `atoms` | ❌ | Real decodense logic: unit scaling, scalar-vs-dipole column construction (the `" (x)"/" (y)"/" (z)"` suffixes), and atom labelling. **Testable from a hand-built `res` dict and a small mol — no SCF needed.** Assert: one row per atom; index entries equal `f"{symbol}{i}"`; expected columns present; no `NaN`; and unit scaling actually applied (`kcal_mol` → ×627.5094740631). The unit-scaling fall-through (item #22 — asking for energy in `debye` silently returns atomic units) should be **pinned as documented current behaviour**, not reported as a bug; that was already decided in PR 4. |
| `orbs` | ❌ | **The single highest-value untested thing in the codebase**, because it is a known, previously-shipped bug: the NDO pairing loop `range(sort_idx.size // 2)` silently drops the middle orbital when the count is odd. Test from a hand-built `res` dict with `ndo=True` and an **odd** number of orbitals; assert the row count equals the orbital count and that the index set is a permutation of `range(n)` (item #20). **Note for whoever writes it:** there is no `_ndo_mo_idx` helper to import — the logic is inline inside `orbs()`, so it must be driven through `orbs()` itself with a hand-built dict. (An earlier draft assumed such a helper existed and broke the whole test module with an `ImportError`.) |

#### `pbctools.py` — previously scoped out wholesale; one correction

The blanket "skip `pbctools.py` entirely" from Part VI's original scoping is right for most of the file but **wrong for
`ewald_e_nuc`**, which should be separated out.

| Function | State | Assessment |
|---|---|---|
| `ewald_e_nuc` | ❌ | **Cheaply testable and worth doing.** It is the PBC analogue of `_e_nuc` — per-atom nuclear repulsion for a cell — and needs only a `Cell` object: no k-points, no density fitting, no builders. **[verified]** PySCF's own `cell.ewald()` exists and returns the total independently, so `ewald_e_nuc(cell).sum() == cell.ewald()` is exactly the `_e_nuc`/`energy_nuc()` differential pattern, against a genuinely separate implementation. Also `cell.natm == 0` → empty array is a decodense-specific edge case, one line. |
| `get_nuc_pbc` | ⚠️ | Mostly dispatch, but the two `NotImplementedError` raises are decodense's own and testable by passing a dummy non-`DF` object. Low value on its own — but see the formatting bug below. |
| `_get_all_e_atomic_df` / `_get_pp_atomic_df` | ⛔ | Require real GDF builders and a converged PBC density fit. |
| `_RSNucBuilder` / `_IntPPBuilder` | ⛔ | Subclasses of *private* PySCF classes (`_RSGDFBuilder`, `Int3cBuilder`), vendored wholesale. No clean public reference exists to test against — that is the core reason this file was scoped out. |
| `_get_pp_nl`, `_int_dd_block_at`, `_merge_dd_at` | ⛔ | Same: vendored PySCF internals. |

**Separate bug found while auditing this file (not a testing issue).** Eight `raise` sites in
`pbctools.py` use the pattern `raise NotImplementedError("... %s ...", value)` — passing the
value as a *second exception argument* instead of formatting it. The message then renders as a
tuple, e.g. `('No Ewald sum for dimension %s.', 3)`, rather than the intended text. **[verified]**
PR 11 (branch `fix-pbctools-docs`, unmerged) fixes exactly two of these — `_int_dd_block_at`'s
`"% intor type"` and `ewald_e_nuc`'s `"No Ewald sum for dimension %s."` — leaving **seven more
of the identical pattern unfixed, on that branch too**: `pbctools.py` lines 51, 58, 150, 260,
293, 420, 556. PR 11 fixed the two instances it happened to notice rather than the class of
bug. Worth folding into PR 11 before it merges, or into a follow-up.

---

# Part III — Reference: every existing test

### Every current test, reviewed

Every test currently in `test_unit.py`, reviewed against four questions: is it testing the
right thing; is it circular; are the assertions meaningful; is the setup realistic enough to
catch a real bug. Verdicts are deliberately harsh — the point is to find what to improve, not
to confirm the work.

The organising principle used throughout: **a test is only worth its line count if the
expected answer comes from outside the code under test** — hand-derived physics, a separate
implementation, or an invariant that must hold regardless of how the code is written.

#### Keep as-is — genuinely strong

| Test | Why it holds up |
|---|---|
| `test_e_nuc_h2` | Expected value hand-derived from `Z₁Z₂/2R` at R=1 bohr. No code consulted. |
| `test_get_nuc_matches_pyscf` | Differential against `mol.intor("int1e_nuc")` — a separate PySCF integral, computed by a different route (decodense builds it per-atom from `int1e_rinv` with shifted origins). Li+H deliberately chosen over H₂ so a symmetry-hidden bug can't pass. |
| `test_dip_nuc` | Hand-derived from eq. 12 of the theory notes. |
| `test_make_rho` / `test_make_rho_gga` | Differential against PySCF's own `numint.eval_rho`. Strongest pair in the file: same physics, genuinely independent implementations, neither copied from the other. |
| `test_permute` | Pure property test — no expected values exist at all, so it cannot be written circularly. Catches the atom-indexing bug class nothing else covers. |
| `test_nuc_sum` | The combined `Tr[V_nuc D]` total is never computed anywhere in decodense — only ever as two split halves — so recomputing it independently (plain `np.sum(nuc * rdm1)`, deliberately *not* decodense's `_trace`) is a real external reference. |
| `test_e_xc_2d` | Underrated: the second `rho` row is filled with garbage `99.0`s, so it genuinely proves only `rho[0]` is consumed. Mutation-proven earlier. |
| `test_population_mul_h2o` | Charge conservation (from eq. 6) *plus* the O/H symmetry check — the symmetry half catches atom-mislabelling, which the sum alone structurally cannot. |

#### Weaknesses found — concrete improvements

1. **`test_e_xc_1d` is blind to a dropped argument.** Both `eps_xc` and `grid_weights` are
   `[1,1,1]`. If `_e_xc` silently stopped using `grid_weights` altogether, multiplying by
   all-ones is invisible and the test would still pass. **Fix:** use non-trivial values, e.g.
   `grid_weights = [0.5, 1.0, 2.0]` and `eps_xc = [2.0, 1.0, 0.5]`, and hand-compute the
   expected sum. Cheap, and closes a real hole.

2. **`test_make_rdm1_symmetric` is close to trivially true.** `D = Σₚ nₚ cₚcₚᵀ` is symmetric
   by construction for essentially *any* implementation, including several plausibly wrong
   ones (a swapped contraction index still yields a symmetric matrix). **Fix:** replace with
   the invariant the plan already named in item #6 — `Tr(D·S) == sum(occup)`, using a real
   overlap matrix `S` from a molecule. That is a genuine, non-trivial constraint: it ties the
   density matrix to the electron count, and a mis-weighted `D` breaks it immediately.

3. **None of the three `_trace` tests can catch a transposition bug.** Every one passes either
   `np.eye(2)` or a symmetric `rdm1`, and `Tr(opᵀ·D) == Tr(op·D)` whenever `D` is symmetric —
   so the op/rdm1 orientation is untested. The plan itself flagged that this only holds for
   symmetric input. Two defensible resolutions, and the choice should be deliberate rather
   than accidental: (a) add one case with a **non-symmetric** `rdm1` and a hand-computed
   expected value, pinning the current orientation; or (b) write down explicitly that every
   real caller passes symmetric input, so the orientation genuinely cannot matter. Right now
   it's neither — it's untested by accident.

4. **`test_mf_info_h2o` silently papers over the plan's own open question.** Item #8 asks
   whether `mf_info` returning int64 occupations (vs PySCF's float64) is a real bug.
   **[verified]** it *is* int64: `dim()` returns `np.where(...)[0]` (int64 indices), and
   `np.ones_like(alpha)` inherits that dtype. But `np.array_equal(int_array, np.ones(5))` is
   `True` regardless of dtype, so the test passes without surfacing the question at all.
   **Fix:** add `assert mo_occ[0].dtype == np.int64` to *document current behaviour*, so that
   if someone later "fixes" it to float64 the test goes red and forces the conversation.

5. **`test_xc_ao_deriv_unknown_type` asserts on a known bug, and nothing says so in the test.**
   It expects `UnboundLocalError` — the current, confusing failure mode. If someone fixes
   `_xc_ao_deriv` properly (adding `else: raise ValueError(...)`, as was drafted and
   deliberately declined), this test goes red *for the right reason*, and the temptation will
   be to "fix the test" by reverting the improvement. **Fix:** a one-line comment in the test
   body saying it pins current buggy behaviour on purpose, and what to change it to when the
   bug is fixed.

6. **`test_h_core`'s second assertion is quietly circular.** `_h_core` computes `nuc` *as*
   `np.sum(sub_nuc, axis=0)`, so asserting `nuc == sub_nuc.sum(axis=0)` re-runs the
   function's own line. It cannot fail. It's harmless only because the *other* three
   assertions in that test (`kin + nuc == mf.get_hcore()`, and the two diagonal-sign checks)
   carry the real weight. Worth keeping only as a cheap consistency check on the returned
   tuple, and worth labelling as such so it isn't mistaken for verification.

7. **`test_charge_conservation` checks only the sum.** Exactly the weakness we already fixed
   in `test_population_mul_h2o`: a sum of zero survives any bug that moves charge *between*
   atoms. **Fix:** add the per-atom sanity checks — oxygen's partial charge negative, both
   hydrogens positive and equal to each other by symmetry.

8. **`test_vk_dft_pbe0` does not verify the `rdm1` pass-through.** `_vk_dft` internally calls
   `mf.get_k(mol=mol, dm=rdm1)` with the `rdm1` it was handed, and the test builds its
   reference with that same `rdm1` — so "did the function actually use the density matrix it
   was given" is unchecked. Judged low-risk (there is no alternate density in scope for it to
   confuse itself with), and recorded rather than silently assumed. A stronger version would
   pass a deliberately perturbed `rdm1` (e.g. `0.5 * mf.make_rdm1()`) and confirm the result
   tracks it.

9. **Fixture inconsistency.** `test_mf_info_h2o` checks `mo_coeff[0].shape[1]` but not
   `mo_coeff[1]`; `test_mf_info_li` checks both. Trivial, but the H₂O one is the weaker of the
   two for no reason.

10. **`mf_li` (UHF) is built but barely used.** It appears only in `test_mf_info_li`. This is
    the biggest *cheap* gap in the file — see item #17 in the coverage audit below.

#### Complete per-test verdict — all 38 tests

The sections above single out what needs changing. This table is the exhaustive pass, so no
test is left without a verdict. Tests already discussed above are cross-referenced rather than
repeated.

| # | Test | Verdict |
|---|---|---|
| 1-3 | `test_trace_identity`, `_symmetric`, `_3d` | Hand-computed/`np.trace` differential, mutation-proven. **But all three use identity or symmetric `rdm1`, so none can catch a transposition bug** — see weakness #3. |
| 4 | `test_e_xc_1d` | Weak — blind to a dropped `grid_weights`, see weakness #1. |
| 5 | `test_e_xc_2d` | Strong — garbage `99.0` second row proves only `rho[0]` is consumed. |
| 6-9 | `test_xc_ao_deriv_lda/hf/gga/mgga` | Cheap and fine, but note **half of each assertion is really PySCF's**: the type string (`"GGA"` etc.) comes from `dft.libxc.xc_type`, so only the derivative level (`0`/`1`/`2`) is decodense's own logic. If PySCF ever reclassified a functional these fail for a non-decodense reason — acceptable, arguably useful as an upstream-change alarm. **`"NLC"` is untested, but not a real gap** — revised, see Part I item A10: no real functional makes `dft.libxc.xc_type` return `"NLC"`; it's dead code, same as the `"UNKNOWN"` branch, only reachable by mocking. |
| 10 | `test_xc_ao_deriv_unknown_type` | Good use of mocking, but pins a known bug with nothing saying so — see weakness #5. |
| 11 | `test_dim` | Correct and hand-built. Edge case (strict `> 0.0` boundary, fractional NDO occupations) noted in Part II. |
| 12 | `test_make_rdm1` | Correct, but **very small**: one orbital, one nonzero entry. It cannot catch a bug in the summation *over* orbitals, because there is only one. Combined with #13's weak assertion, **no current test checks a multi-orbital density matrix against a known value** — which is exactly what the proposed `Tr(D·S) == sum(occup)` replacement fixes. |
| 13 | `test_make_rdm1_symmetric` | Near-trivially true — see weakness #2. |
| 14 | `test_sanity_check_invalid_pop_method` | Correct, and currently the *only* coverage of `pop_method` (deliberately left out of the table to avoid duplication). Stylistically inconsistent though: one attribute in a standalone function while nine others live in a parametrized table. Cleaner to fold it in as a row and delete the standalone. |
| 15 | `test_e_nuc_h2` | Strong (hand-derived physics). Missing the `.sum() == mol.energy_nuc()` differential half — Part II. |
| 16 | `test_get_nuc_matches_pyscf` | Strong differential. Missing shape + per-slice symmetry assertions — Part II. |
| 17 | `test_dip_nuc` | Strong (hand-derived). Missing the gauge-shift property — Part II. |
| 18 | `test_mf_info_h2o` | Correct but papers over the int64 dtype question, and checks only `mo_coeff[0]` — weaknesses #4 and #9. |
| 19 | `test_mf_info_li` | Correct, and the better of the pair (checks both spins, genuinely different sizes 2 vs 1). Low blast radius: `mf_info` has zero internal callers. |
| 20 | `test_population_mul_h2o` | Strong — conservation *and* symmetry. |
| 21-22 | `test_assign_rdm1s_h2o`, `_iao` | Good as far as they go, covering the two population schemes that matter in practice. **But `assert np.array_equal(weights[1], weights[0])` is currently unfalsifiable in a dangerous way**: with only restricted cases in the suite, it would still pass if the code *always* duplicated alpha into beta — a serious bug for any unrestricted calculation. Adding the `mf_li` case is what gives this assertion its meaning. |
| 23 | `test_h_core` | Three strong assertions plus one quietly circular one — weakness #6. |
| 24-25 | `test_make_rho`, `_gga` | Strongest pair in the file. MGGA and the unrestricted branch still uncovered — Part II. |
| 26 | `test_vk_dft_pbe0` | Good, with the `rdm1` pass-through caveat — weakness #8. Range-separated branch uncovered. |
| 27 | `test_permute` | Strong pure property test. |
| 28 | `test_charge_conservation` | Sum-only, survives charge moving between atoms — weakness #7. |
| 29 | `test_nuc_sum` | Strong — the combined total genuinely exists nowhere in decodense. |
| 30 | `test_sanity_check_rejects_invalid_options` | Good: 9 rows, each reported individually. Minor note: `match` fragments like `"invalid unit"` appear in *both* that attribute's `TypeError` and `ValueError` messages — still unambiguous because the expected exception class is also specified, but the fragments are doing less work than they appear to. |
| 31 | `test_sanity_check_verbose_valid_boundaries` | Valuable — the only thing confirming valid input isn't over-rejected. Implicit assertion, see below. |
| 32-33 | `test_sanity_check_gauge_origin_wrong_shape`, `_wrong_type` | Good. These cover the one check in `sanity_check` with real logic (length *and* per-element type), including the `TypeError`/`ValueError` split between the two separate guards. |
| 34 | `test_sanity_check_gauge_origin_valid` | Good counterpart to 32-33. |
| 35 | `test_sanity_check_rejects_wrong_types` | Good. **One untested subtlety worth deciding on:** Python's `bool` is a subclass of `int`, so `verbose=True` passes `isinstance(decomp.verbose, int)` and is accepted. Whether that is intended is an open question, not currently pinned either way. |
| 36-37 | `test_sanity_check_rejects_invalid_mo_coeff`, `_mo_occ` | Good, and among the more valuable of the group: these run on *every* real `main()` call, and cover the `None`-allowed asymmetry (`mo_occ` may be `None`, `mo_coeff` may not). |
| 38 | `test_sanity_check_accepts_valid_input` | Valuable as the baseline "a fully default `DecompCls` is actually accepted" check. Implicit assertion, see below. |

#### On the "should not raise" tests

`test_sanity_check_verbose_valid_boundaries`, `test_sanity_check_gauge_origin_valid` and
`test_sanity_check_accepts_valid_input` contain no `assert` at all — they pass by *not*
raising. That is legitimate, idiomatic pytest, and worth keeping (they're the only thing
confirming valid input isn't over-rejected), but it's an implicit assertion and reads as an
empty test to someone unfamiliar with the pattern. A one-line comment in each would help.

---

---

# Part IV — How to write the tests

### What a unit test actually is

A test is a function that runs your code and **crashes on purpose** if the answer is wrong.

```python
def test_two_plus_two():
    assert 2 + 2 == 4
```

A **unit test** tests **one function, in isolation**, with inputs you build by hand, where
you know the right answer independently of the code being tested.

Three properties, all about *diagnosis*:

- **Small** — one function. When it fails you know exactly what is broken.
- **Isolated** — no SCF, no files, no globals. Nothing else can be the cause.
- **Fast** — microseconds, so you run them after every edit.

**[verified]** the speed difference here is not subtle:

| | time |
|---|---|
| 1000 calls to `_trace` | 0.008 s (8.2 µs each) |
| one existing integration test | 2.9 s |
| whole current suite | ~16 s |

That is why a good suite is mostly unit tests — you can afford hundreds.

---

### The three kinds of test

The hard part of any test is *where the expected answer comes from*.

| Kind | Expected answer comes from | Example in decodense |
|---|---|---|
| **Unit** | You compute it by hand | `_trace` of a 2×2 matrix |
| **Property** | A rule the answer must always obey | weights sum to 1 |
| **Differential** | A second, independent implementation | `_e_nuc` vs `mol.energy_nuc()` |

Kinds 2 and 3 are what decodense is missing most, and they are the most valuable because
**they need no reference data at all**.

#### Property tests — the key idea

You do not need to know the right answer. You only need a rule it must obey.

The central example: **population weights are a partition of unity.** Each orbital's
weights must sum across atoms to its occupation — for every molecule, basis and method.

**[verified]** on water/STO-3G:

```
mulliken     max|deviation from 1| = 8.9e-16
lowdin                               1.1e-15
meta_lowdin                          8.9e-16
becke                                1.6e-07     <-- grid-integrated, a million times looser
iao                                  1.8e-15
```

Two lessons in one experiment:

- It is a tight invariant, and it constrains exactly what the sum rule cannot see: the
  *distribution*.
- **The right tolerance is a property of the method, not a global constant.** A single
  `TOL = 9` for everything is wrong in both directions.

---

### pytest in ten minutes

#### Install

**[verified]** pytest is **not installed** in your environment yet:

```bash
pip install pytest
```

That is the only new dependency, and it **runs your existing `unittest` files unchanged** —
which is exactly how PySCF works (439 `unittest` files, CI runs `pytest`).

#### Your first test

`tests/test_trace.py`:

```python
import numpy as np
from decodense.properties import _trace


def test_trace_of_identity_rdm1_returns_matrix_trace():
    op = np.array([[1.0, 2.0],
                   [2.0, 3.0]])
    rdm1 = np.eye(2)
    assert _trace(op, rdm1) == 4.0     # 1 + 3, by hand
```

- File name and function name **must start with `test_`**, or pytest will not find them.
- Plain `assert` — no class, no `self.assertEqual`.
- **The expected value is computed by hand.** If you call the code under test to build
  your expectation, you are testing nothing.

#### Running and reading output

```bash
cd tests                       # required: geometry files use relative paths
python3 -m pytest -v
```

```
tests/test_trace.py::test_trace_of_identity_rdm1_returns_matrix_trace PASSED   [100%]
                    ^ test id — copy this to re-run just this one
```

On failure you get both sides for free:

```
>       assert _trace(op, rdm1) == 4.0
E       assert 5.0 == 4.0
tests/test_trace.py:10: AssertionError
```

Commands you will actually use:

```bash
python3 -m pytest                  # everything
python3 -m pytest test_trace.py    # one file
python3 -m pytest -x               # stop at first failure
python3 -m pytest -k weights       # only tests with "weights" in the name
python3 -m pytest --durations=10   # the 10 slowest
```

#### Floats never compare equal

```python
assert 0.1 + 0.2 == 0.3      # False. Really.
```

So every numeric assertion needs a tolerance. Use PySCF's spelling — it is what your
existing tests use and it works for scalars and arrays alike:

```python
assert abs(result - expected) < 1e-10
assert abs(array - expected).max() < 1e-10
```

Absolute vs relative: absolute (`abs(a-b) < tol`) asks "within tol of the answer";
relative asks "within tol of the *magnitude*". For energies near −76 Hartree use absolute.
For values near zero, relative blows up — use absolute.

#### Fixtures — build expensive things once

A fixture is a function that builds what a test needs; tests request it by name.

```python
import pytest
from pyscf import gto, scf


@pytest.fixture(scope="module")
def h2o_rhf():
    mol = gto.M(atom="geom/h2o.xyz", basis="sto-3g", verbose=0)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.kernel()
    return mol, mf


def test_weights_sum_to_occupation(h2o_rhf):
    mol, mf = h2o_rhf
    ...
```

`scope=` controls rebuild frequency: `function` (default, per test), `module` (once per
file — right for SCF), `session` (once per run). Use `module` for anything expensive, but
only if no test mutates it.

Put a fixture in `conftest.py` and every test file in that directory gets it with no
import. That is the one thing `conftest.py` is for.

#### Parametrize — replaces your `subTest` loops

```python
@pytest.mark.parametrize("pop_method",
                         ["mulliken", "lowdin", "meta_lowdin", "becke", "iao"])
def test_weights_sum_to_one(h2o_rhf, pop_method):
    ...
```

This gives **five separately named tests**, not one with five hidden cases:

```
test_weights_sum_to_one[mulliken]  PASSED
test_weights_sum_to_one[becke]     FAILED
```

You see immediately *which* case broke, and can re-run only it.

---

### Tolerances

The tolerance should match the precision of the **weakest step** producing the number.

| Testing | Tolerance | Why |
|---|---|---|
| Pure algebra (`_trace`, `_e_xc`, index logic) | `1e-12` … `1e-14` | Only round-off |
| Linear-algebra invariants (mulliken/lowdin/iao weights) | `1e-12` | **[verified]** deviation ~1e-15 |
| Anything on a DFT grid (becke, XC energies) | `1e-6` … `1e-7` | **[verified]** becke deviates 1.6e-7 |
| SCF energy vs stored reference | `1e-8` … `1e-10` | needs `conv_tol` ≤ 1e-10 |

Two rules that matter more than the table:

1. **Tighten `conv_tol` inside the test.** If SCF is converged to 1e-6, an assertion at
   1e-9 is testing convergence noise. PySCF comments this explicitly in its own tests.
2. **Never loosen a tolerance to make a test pass.** Either the code is wrong or your
   justification was. Loosening to green is how a suite becomes decorative.

---

### Gradient / finite-difference tests — background only

Not needed for mainline (it computes no derivatives), but you asked how they differ.

A unit test checks a *value*; a finite-difference test checks a *derivative*, by computing
it numerically a second way:

```
f'(x)  ≈  ( f(x + h) − f(x − h) ) / (2h)
```

Choosing `h` is a genuine trade-off, which is why gradient tolerances are loose: too large
and the formula's own truncation error dominates (`h²`); too small and you subtract two
nearly-identical numbers and round-off explodes (`ε/h`). For float64 the sweet spot is
`h ≈ 1e-4`–`1e-3`, agreeing to about `1e-6`.

The generalisation you *will* use on mainline is the **differential test** — "compute it a
second, independent way" — which is Part IV kind 3 and appears all over Part VI's original scoping.

---

---

# Part V — Why, and how others do it

### Why we are doing this

Every real test decodense has today asserts exactly one thing: **the decomposed pieces add
back up to PySCF's total energy.** Call it the *sum rule*.

The problem: a partition adds up to its total **no matter how wrongly you split it**. If
atom A's share is too big by 0.3 and atom B's too small by 0.3, the sum is unchanged. So
the sum rule is blind to the whole class of bugs decodense exists to avoid — *energy on
the wrong atom or the wrong orbital*.

This is not theoretical. Two real bugs were found in this codebase in the last weeks, and
**the suite was green for both**:

1. **`comp_key_dict` mapped wrong attribute names** — `res.mo_occ` raised `AttributeError`.
   Every number was correct, so every sum rule passed.
2. **NDO dropped an orbital when the count was odd** — the pairing loop
   `range(sort_idx.size // 2)` silently discarded the middle orbital from the results
   table.

Both are "right total, wrong distribution". Both would have been caught in milliseconds by
a three-line test. That is the argument for this whole project.

**Important:** the existing tests are not badly written. They follow PySCF's house style
closely. What is missing is a *layer of small tests underneath them*, not a rewrite.

#### This is not just us — a quick outside check

Two real sources on testing scientific code specifically (not general software), checked
this week:

- ["Ten Simple Rules for Writing Clean and Reliable Open-Source Scientific
  Software"](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1009481)
  (PLOS Comp. Biology) — a peer-reviewed rules paper aimed exactly at researchers writing
  scientific code, not professional software engineers.
- ["Testing your code"](https://fabienmaussion.info/scientific_programming/week_06/02-Testing.html),
  a university course on scientific Python.
- ["Unit testing best practices"](https://www.ibm.com/think/insights/unit-testing-best-practices)
  (IBM) — general software advice, not scientific-specific, useful as a sanity check on
  the mechanics.
- ["Fundamentals of Unit Testing. With Pytest, A Python Testing
  Framework"](https://medium.com/@taeefnajib/writing-unit-tests-for-your-data-science-projects-c3d46e98beb2) —
  aimed specifically at data-science code, and pytest-based.

Things they say that directly matter for you, given how much code there is here:

1. **A rule of thumb of ~60% test coverage is considered good, not 100%.** "Writing tests
   for scientific code can be difficult... tests also require maintenance, so ensure
   tests are of high quality and adequate utility to merit inclusion." You do not need to
   test every function in every file — Part VI's original scoping below makes an explicit decision about which
   ones matter and which do not.
2. **Write the regression test the moment you find a bug**, not on some later cleanup
   pass — "writing tests when you encounter bugs will help to validate a fix and prevent
   such issues from reoccurring." This is exactly how `test_comp_keys.py` and
   `test_ndo_mo_idx.py` came about — **note, both have since been deleted**, their content
   consolidated into `test_unit.py`. See the consolidation note in Part II.
3. **For numeric results specifically**: never compare floats for exact equality, always
   with a tolerance (`np.testing.assert_allclose` or `abs(a-b) < tol`); check outputs for
   `NaN` and for being in a physically sane range; check that shapes and types are
   preserved. All of this already matches what Part VI's original scoping/Part IV's tolerance table below recommend.
4. **Every test has the same three-part shape**, IBM calls it Arrange-Act-Assert: build
   the inputs, call the one function, check the result. Every example in Part IV/Part VI's original scoping already
   follows this; it is just a name for the pattern, worth knowing so you recognise it.
5. **Coverage targets vary by source (60% vs. 70-80%) — the number matters less than the
   principle**: none of these three sources says "test everything". Testing the ~15
   functions decided on in Part VI's original scoping and stopping there is following the advice, not falling
   short of it.
6. **One IBM recommendation that does *not* transfer well here: mocking.** Standard
   advice for business software is to mock out dependencies so a unit test never touches
   a real database or network call. Decodense's version of that idea is different:
   instead of mocking PySCF, the differential tests in Part VI's original scoping Tier 1 call the *real* PySCF
   function and compare against it directly (`_e_nuc` vs `mol.energy_nuc()`), because
   PySCF's own result is the known-good answer, not something to fake. This mirrors how
   OTR only mocks its Fortran solver to test the *interface*, never the physics.

---

### How the projects you respect actually do it

Researched directly from the checkouts in `../pyscf`, `../pyscfad`, `../opentrustregion`,
plus two outside libraries fetched for comparison: QuTiP and Qiskit — both mature,
professionally-maintained numerical/physics packages, useful precisely because they are
*not* quantum chemistry, so what transfers is the methodology, not the domain.

| | PySCF | PySCFAD | OpenTrustRegion | QuTiP | Qiskit |
|---|---|---|---|---|---|
| Framework | `unittest`, **run by pytest** | **pure pytest** | `unittest` driving Fortran | **pure pytest** | `unittest` + stestr |
| Layout | per-module `test/` | one top-level `tests/` | flat `tests/`, split by role | `qutip/tests/`, ships with package | mirrors package, custom base class |
| Fixtures | rare; `setUpModule` | yes, the normal way | `setUpClass` | yes + a seeded-RNG fixture | n/a (unittest) |
| Float assert | `assertAlmostEqual` (9910×) | `abs(a-b).max() < tol` | `abs(a-b) > tol` | `assert_allclose`, tiered `atol` | `assertEqual` w/ tolerance baked into `__eq__` |
| Slow tests | name suffix `_high_cost` + `-k` | same convention | n/a | `@pytest.mark.slow` | env-var gated `@slow_test` |

#### QuTiP and Qiskit — what's genuinely new here, not already covered above

Five things neither PySCF/PySCFAD/OTR showed clearly, all directly usable:

**1. Tolerance belongs in the test case, not a module constant.** Both libraries tie the
tolerance to *why* a value might be imprecise, and keep it next to the case rather than as
one global. QuTiP: exact invariants (norm conservation) get `atol=1e-15`; two ODE
integrators compared against each other get `atol=5e-5`; Monte-Carlo trajectories get
`tol=0.25`. Qiskit does this with `self._rtol`/`self._precision` set once in `setUp()`.
For us: `assign_rdm1s`'s tolerance already varies by `pop_method` (Part VI's original scoping Tier 2) — this is the
same idea, generalized, and confirms it's the right call rather than an ad hoc exception.

**2. `atol` alone (no `rtol`) when a value can legitimately be ~0.** Qiskit's house style
for near-zero quantities is `atol=1e-12, rtol=0` — a relative tolerance is meaningless when
the true answer might be zero (a hydrogen's exchange contribution, an off-diagonal density
matrix element). Worth using explicitly wherever a component could genuinely vanish.

**3. Invariants over reference values, even more aggressively than OTR already argued.**
Qiskit's `Operator(circuit) == Operator(transpiled_circuit)` — checking two different code
paths compute the *same* thing — is structurally identical to `_e_nuc` vs
`mol.energy_nuc()` in Part VI's original scoping, just phrased as a house rule rather than an occasional trick.
QuTiP adds one we hadn't used: **permutation invariance** — nothing in either library
depends on element order, so relabeling and re-checking is a free, reference-free test.

**4. Catch API drift at test time, not at the next upgrade.** Qiskit's base test class
turns every `DeprecationWarning` from a dependency into a hard failure
(`warnings.filterwarnings("error", category=DeprecationWarning)`). This is exactly the
failure mode behind PR 17 (the examples broke silently against a newer PySCF) — cheap to
add, and it converts a future silent break into an immediate, loud one.

**5. Introspective sweeps.** Qiskit auto-discovers every registered gate subclass and tests
an invariant against each one, so a newly-added gate is covered without anyone remembering
to add a test for it. Applies directly to `sanity_check`'s valid-option tuples (Part VI's original scoping
`decomp.py`) — derive the parametrize list from the same tuple `sanity_check` checks
against, not a hand-copied duplicate that can drift.

#### What OTR does, in more detail — it is the best model for decodense

OTR is small enough that its whole suite is comprehensible, and it has a **written testing
philosophy** (`../opentrustregion/CLAUDE.md`) worth copying almost verbatim. Four ideas:

**1. Test against a toy problem with a known answer, not a real calculation.**
Their entire optimizer suite runs on the **Hartmann 6-D function** — a cheap analytic
function whose minima are known and hardcoded — not on SCF. The decodense analogue:
**build small `rdm1` / `weights` / `mo_occ` arrays by hand** instead of running a
calculation, wherever the function does not genuinely need a `Mole`. Most of Part VI's original scoping Tier 0
works this way.

**2. Exactly one registered test per production routine.** Multiple cases go *inside* that
test, not as `test_foo_case_a` / `test_foo_case_b`. Keeps the suite navigable.

**3. A unit test may only call the routine under test** — never another production
function to build its inputs or its expected values. Otherwise a failure elsewhere makes
this test fail too, and you have lost the diagnosis.

**4. Prove the test can fail.** Their rule, quoted:

> "Verify a new test actually fails when the routine is broken. Mutate the routine (flip a
> sign, swap an index, drop a term), rebuild, confirm the test fails, then restore."

**Do this for every test you write this week.** A test that cannot fail is worse than no
test — it produces false confidence. This is the single highest-value habit here.

Also worth stealing: OTR organises tests **by role** (`*_unit_tests`, `*_system_tests`, a
shared reference/tolerance file) rather than only by module, and keeps one tolerance
constant in that shared file instead of scattering literals.

#### The decision for decodense

**Keep the 7 existing `unittest` files exactly as they are, run everything with pytest,
and write all new tests in pytest style.**

Worth being precise here, since it came up: **neither OTR nor the decodense AD branch uses
pytest** — re-checked directly, not just recalled. OTR has zero `pytest` occurrences
anywhere in its repo; every one of the AD branch's 8 test files, including the AD/`jacrev`
one, is `import unittest` / `class KnownValues(unittest.TestCase)`, identical in style to
mainline's own existing tests. So this is not "the precedent says pytest" — it's a
judgement call, for these reasons:

- Zero migration cost, zero risk — pytest runs `unittest` files unchanged, so recommending
  it does not mean rewriting anything OTR, the AD branch, or mainline already has.
- PySCFAD — the framework the AD branch actually depends on for differentiation — is pure
  pytest with fixtures, and that is where fixtures/`parametrize` are worth having.
- Converting the existing 7 integration tests to pytest style would take a day, change no
  behaviour, and risk breaking working tests. Do not do it — this is about the *new* unit
  tests only.

---

---

# Part VI — Decisions, history, corrections

Kept because it records *why* things were decided, including what was deliberately declined.
Where Part II and the original scoping below disagree, **Part II is current**.

### Consolidation of the old single-purpose test files — done, with one consequence

`tests/test_comp_keys.py` and `tests/test_ndo_mo_idx.py` have been **deleted**, their scope
folded into `test_unit.py`. That is the right direction — the same check living in two files
is exactly the drift hazard that produced the original bugs — but it has one consequence that
needs acting on rather than leaving implicit:

- **The `comp_key_dict` check now has zero coverage anywhere in the repository.** It survives
  only as a commented-out block in `test_unit.py`. Deleting the old file removed the working
  `unittest` version, and the replacement is switched off. For a check guarding a bug that has
  already shipped once, that is the worst available state. **Action: re-enable it as
  `@pytest.mark.xfail(reason=...)`** so it runs, reports honestly as a known failure, and
  flips to `XPASS` the moment `fix-comp-key-dict` merges.
- **The NDO index check likewise has zero coverage.** Worth being precise about why the
  existing integration tests cannot substitute: `test_ch2_hf_ndo.py` and `test_c5h5n_hf_ndo.py`
  assert on `np.sum(res.tot[...])`, which reads from `res_dict` — populated by `prop_tot`. The
  odd-orbital bug lives in `orbs()`, which builds the *display DataFrame* from that same dict
  afterwards. A row dropped there never touches `res.tot`, so the sum rule is structurally
  blind to it. That is precisely how the bug shipped with a green suite, and why the `orbs()`
  test in Part II is ranked first.

**Remaining file-level cleanup**, once the two items above are covered: nothing. The 7
remaining `test_*.py` files are all genuine end-to-end integration tests on distinct
molecule/method combinations, not duplicates of unit-level checks, and Part VI's original scoping already records the
decision to keep them as the backstop.

### Original scoping — the first 23-item plan *(superseded by Part II; kept for the declined-item reasoning)*

There are 7 files in `decodense/`. That's **~15 functions worth a test, out of ~45 total**
— the rest is a thin wrapper, orchestration already covered by the 7 integration tests, or
genuinely not testable without machinery this plan isn't building yet. Within each file,
items are ordered easiest first. A ⚠️ marks something with a real dependency or open
question. (The "Working order" cross-reference here pointed to a section that was never written — use Part I's worklist instead.)

#### `decomp.py` — no SCF needed at all

| # | Test | Why | Approach |
|---|---|---|---|
| 1 | `comp_key_dict` is exactly the inverse of `CompKeys` | **Where the wrong-name bug lived** (`res.mo_occ` didn't exist) | Property: `comp_key_dict == {v: k for k, v in vars(CompKeys).items() if not k.startswith("_")}`. Was in `tests/test_comp_keys.py`, **now deleted** — so this check currently has zero
coverage anywhere in the repo, see Part II. ⚠️ **fails on current `main`, passes once PR 9 (branch `fix-comp-key-dict`) merges** — this is a real, currently-red test, not a mistake |
| 2 | `DecompCls().gauge_origin is not DecompCls().gauge_origin` | **Real PR 2 regression** — used to be a shared mutable default (`np.zeros(3)` evaluated once at import). **[deprioritized]** discussed live and skipped for now: unlike `comp_key_dict`, this bug was latent and never actually triggered by anything in the real codebase before the fix — lower real-world urgency. Not written; candidate for the final gap review | Regression, exact. Build two instances, assert different `id()`, assert mutating one doesn't affect the other |
| 3 | `sanity_check` rejects every invalid input | ~20 `raise` paths after PR 4 | **One parametrized test**, not 20 — a table of `(attr, bad_value, exception_type, message_fragment)`, one `pytest.raises` per row. See Part IV's `pytest.raises` example |
| 4 | Valid-option tuples used *by* the parametrize table in #3 are the same ones `sanity_check` checks against | Qiskit's introspective-sweep idea | Don't hand-copy `("mulliken", "lowdin", ...)` into the test file — import/derive it so a 6th `pop_method` added later is automatically covered, not silently skipped |

Skip `DecompCls`/`CompKeys` themselves beyond the above — plain data containers, nothing
else to break.

#### `tools.py` — mostly pure, one needs an `mf`

| # | Test | Why | Approach |
|---|---|---|---|
| 5 | `dim` | Occupied-orbital indices from `mo_occ` | Unit, exact. Edge cases: fractional (NDO) and zero occupations |
| 6 | `make_rdm1` | Wrong contraction order gives a transposed/mis-weighted density | Unit `1e-12` + properties: result symmetric, `Tr(rdm1·S) == sum(occup)` |
| 7 | `res_add`/`res_sub` | Mismatched keys now raise (post-PR 5); `"Symm."` handled specially | Unit: mismatched keys raise `ValueError`; matched keys combine correctly. ⚠️ Both functions currently hardcode the **literal string** `"Symm."` (`tools.py:291,306`) instead of `CompKeys.orbsym` — the exact same hand-maintained-string drift hazard that caused finding #1. Worth fixing alongside the test, not just testing the current (fragile) behavior |
| 8 | `mf_info` | dtype question: returns `mo_occ` as int64, PySCF's own is float64 | Unit: assert dtype. **[open question]** — settle whether this is a real bug before deciding what to assert. Needs SCF (restricted **and** unrestricted — this is where PR 2's RHF crash lived) |
| 9 | `make_natorb` | NDO orbital construction | Property: occupations sum to electron count, orbitals stay S-orthonormal. Needs SCF |

Skip `orbsym`, `write_rdm1`, `git_version`, `logger_config`, `contract` — logging, file
I/O, or thin one-line PySCF wrappers; not where a physics bug hides.

#### `properties.py` — 825 lines, but only ~7 functions are unit-testable

| # | Test | Why | Approach |
|---|---|---|---|
| 10 | `_trace` | The workhorse of every energy component; `Tr(opᵀ·rdm1) == Tr(op·rdm1)` only if one is symmetric, and every real caller passes symmetric input | Unit `1e-12`, hand-built 2×2. Both the `ndim==2` and `ndim==3` (dipole) branches |
| 11 | `_e_xc` | Integrates XC energy over the grid; handles two different `rho` shapes | Unit `1e-12`, both shapes |
| 12 | `_xc_ao_deriv` | Maps functional name → AO derivative level | Exact equality, all four known types, **plus** the unknown-type case — **[verified, done]** `tests/test_unit.py` has all 5 cases (LDA, HF, GGA, MGGA, and the unknown-type crash). The `UnboundLocalError` is **confirmed real**, not just a hypothesis: no natural PySCF functional name triggers PySCF's own `"UNKNOWN"` classification, so the crash was reproduced by mocking `pyscf.dft.libxc.xc_type` directly (`unittest.mock.patch`) — a legitimate, narrow use of mocking, justified specifically because the natural trigger is unreachable in practice. **The fix itself (adding `else: raise ValueError(f"unknown xc_type: {xc_type}")`) was considered and deliberately left undone** — the test documents the current confusing-but-real crash rather than changing behavior; revisit at the final gap-review pass (Part VI's week plan) if this should change |
| 13 | `_e_nuc` | Per-atom nuclear repulsion | **Differential**: `.sum()` == `mol.energy_nuc()` at `1e-12`. **Plus analytic**: H₂ at separation R, each atom gets exactly `Z₁Z₂/2R` — two independent checks, no reference data |
| 14 | `_get_nuc` | Per-atom nuclear attraction integrals | Differential vs `mol.intor("int1e_nuc")` at `1e-10`; shape `(natm, nao, nao)`; each slice symmetric |
| 15 | `_dip_nuc` | Nuclear dipole contribution | Unit with hand-computed `Z·(r − origin)`. Property: shifting gauge origin by `d` shifts the total by `−Z_total·d` |
| 16 | `_h_core` | Kinetic + nuclear parts of the core Hamiltonian | Differential: `kin + nuc == mf.get_hcore()` at `1e-10`. Needs SCF |

Skip `prop_tot` and its nested `prop_atom`/`prop_eda`/`prop_orb` — they close over
`prop_tot`'s own local variables, so they're literally uncallable in isolation; already
covered by the 7 integration tests. Skip `_solvent`/`_pcm`/`_point_charges` and the
`_make_rho*` chain — real physics, but each needs a full QM/MM/PCM/grid object, not
hand-buildable; revisit later.

#### `orbitals.py` — the highest-value file, test all of it

| # | Test | Why | Approach |
|---|---|---|---|
| 17 | `assign_rdm1s` partition of unity | **The single most important test in this plan.** Every orbital's weights must sum to its occupation — the exact thing the sum rule structurally can't see | Property, parametrized over all 5 `pop_method`s, **tolerance in the same tuple**: `1e-12` for mulliken/lowdin/meta_lowdin/iao, `1e-6` for becke (verified: grid-integrated, deviates ~1.6e-7). Parametrize restricted **and** unrestricted too — needs SCF |
| 18 | `_population_mul` | Underlying Mulliken populations | Unit-testable directly, hand-built `rdm1` + fake `ao_labels` — **no SCF needed**, do this early |
| 19 | `_population_becke` | Underlying Becke populations | Same partition-of-unity property, `1e-6`. Do **not** cross-compare against `_population_mul` — different definitions, not expected to agree |

#### `results.py`

| # | Test | Why | Approach |
|---|---|---|---|
| 20 | NDO index logic in `orbs()` | **Where the odd-orbital bug lived** | Hand-built `mo_occ` of odd **and** even length; assert the resulting index array is a permutation of `range(n)` |
| 21 | `atoms()`/`orbs()` structure | The NDO bug was a row-count bug, not a value bug | Structural: row count (one per atom/orbital), index name, expected columns present, no NaN |
| 22 | Unit-scaling fall-through | Asking for energy in `debye` (or a dipole in `kcal_mol`) silently returns atomic units — **already known, already declined in PR 4** ("prop/unit cross-check") | Pin this as **documented current behavior**, not a bug report. ⚠️ Much easier once PR 9 lands `_unit_scaling()` as its own function |

Skip `info()` (already verified byte-for-byte during its f-string rewrite) and `ResultsCls`
itself (thin wrapper).

#### `decodense.py` — no unit tests, one new integration-style property

| # | Test | Why | Approach |
|---|---|---|---|
| 23 | Permutation invariance | Relabel/reorder the atoms in a molecule → per-atom results should permute identically. Reference-free, exercises the **whole pipeline** at once — the QuTiP/Qiskit lesson that invariants beat reference values, applied at the highest level we have | Build the same molecule with atoms in a different order, run `main()` twice, assert the results match under the corresponding permutation. Needs SCF |

Skip `main()` for unit testing beyond #23 — pure orchestration, already covered by the 7
integration tests.

#### `pbctools.py`

Skip entirely, on purpose — vendored PySCF internals, PR 11 (already done).

#### Two whole-pipeline invariants worth adding regardless of file, found while re-reading the code

Both reference-free, both directly assertable from what the code already does:

- **Charge conservation** (`properties.py:97-103`): per-atom partial charges must sum to
  the molecule's total formal charge.
- **`Elect.` equals the sum of its own parts** (`properties.py:278`): the code literally
  computes it as `sum(res.values())` — so that relationship is a one-line assertion, not a
  derived claim.
- **`E_ne (1)` + `E_ne (2)`**: the two nuclear-attraction halves (`properties.py:244-249`).
  `PLAN.md` itself notes this split "is not guessable from the code" — a test here doubles
  as documentation of what the split actually means.

#### `pytest.raises`, referenced above — the pattern, once

```python
def test_invalid_pop_method_is_rejected():
    decomp = decodense.DecompCls(pop_method="not_a_method")
    with pytest.raises(ValueError, match="invalid population scheme"):
        decodense.main(mol, decomp, mf, mo_coeff)
```

#### The 7 existing integration tests

Keep them; they're the backstop. Two cheap improvements: move module-level SCF setup into
`setUpModule()` so *collecting* the suite doesn't run every SCF, and consider the
`_high_cost` suffix PySCF/PySCFAD both use to skip slow tests by default.

#### Edge cases worth deliberately sweeping — not hypothetical, proven by our own history

PySCF sweeps RHF/UHF/ROHF × density-fitting × k-points because each is a genuinely separate
code path that can break independently. decodense doesn't implement SCF itself, so most of
that combinatorics isn't ours — but the same lesson applies, and it isn't abstract for us:
**every one of PR 2's six bugs hid in a combination none of the 7 existing test files
exercise.**

| Axis | Why it matters here | Evidence |
|---|---|---|
| Restricted (RHF) vs. unrestricted | All 7 test files use UKS/ROHF; none use RHF | `mf_info` crashed on RHF, invisible until PR 2 |
| xc functional type (LDA/HF vs. GGA/MGGA) | All 7 use GGA-type functionals | Two separate crashes only on plain `'hf'`/`'lda'`, invisible until PR 2 |
| Odd vs. even orbital/electron count | Easy to only test one parity | The actual PR 9 NDO bug |
| ECP present | None of the 7 use one | Silently wrong energy, off by ~2 Hartree, PR 2 |
| Gapped occupations (excited state) | None of the 7 trigger this | Verbose-print crash, PR 2 |

Where this fits above: items 17 and 23 are explicitly parametrized restricted/unrestricted;
wherever a hand-built input is needed (items 1-11), deliberately include one odd-count and
one even-count case rather than whatever's convenient. Not "test everything" — "test the
specific handful of axes that have already bitten us."

**One considered and declined**: testing that populations are non-negative. Mulliken
populations can be legitimately negative, so that invariant is simply false for us — noted
here so nobody re-adds it later without checking.

#### Progress audit — what is written, and the rule we settled on

Written while building `tests/test_unit.py`. The organising question turned out not to be
"which functions are left" but **where does the expected answer come from**. A test earns its
place only when (a) `grep` confirms the function is actually called in real runs — checked,
not assumed — and (b) the expected value comes from *outside* the code under test: physics
derived by hand, a separate implementation (PySCF), or an invariant that must hold regardless
of how the code is written. A test that re-runs the function's own steps outside the function
and checks they agree verifies wiring, not correctness.

Strong, keep as-is:

| Test | Reference source |
|---|---|
| `_e_nuc` | hand-derived H2 analytic value |
| `_get_nuc` | PySCF's own `int1e_nuc` |
| `_dip_nuc` | hand-derived, eq. 12 of the theory notes |
| `_trace`, `_e_xc` | hand arithmetic, both mutation-proven |
| `make_rdm1`, `dim` | hand arithmetic |
| `comp_key_dict` | derived from `CompKeys`; currently red on a real bug |
| `_population_mul` | charge conservation (eq. 6) + O/H symmetry in water |

Weak but cheap, keep without investing further: `_xc_ao_deriv` and `sanity_check`
(definitional), `mf_info` (grep shows zero internal callers — exported only), `assign_rdm1s`
(structure and the restricted-reference shortcut only).

Done, this batch:

1. `_h_core` — `kin + nuc == mf.get_hcore()`; `sub_nuc.sum(axis=0) == nuc`; kinetic diagonal
   positive, nuclear-attraction diagonal negative (catches a swapped return).
2. `_make_rho` — differential against PySCF's own `numint.eval_rho` (LDA case), the strongest
   external reference available in the codebase, since the docstring says it was adapted from
   exactly that function but is genuinely separate code.
3. Permutation invariance (`test_permute`) — same water molecule built twice with atoms listed
   in a different order (`geom/h2o_permute.xyz`, H/H/O instead of O/H/H); `res.el` must match
   per physical atom regardless of its position in the list. First test to exercise the actual
   `decodense.main()` pipeline rather than an internal function. **Correction to an earlier
   entry here:** this originally recorded that `symmetry=True` makes PySCF silently
   re-canonicalize atom order. That was wrong. It was diagnosed by running against
   `geom/h2o_permute.xyz` before the reordered version had actually been saved to disk, so
   *both* files still contained O/H/H and appeared to have been "re-canonicalized". Verified
   afterwards: with the file correctly saved, `gto.M(..., symmetry=True)` preserves the input
   atom order exactly (checked for both geometry files, with and without `symmetry`). The real
   lesson is a different one, and a better one: **a diagnostic run against a stale file
   produced a confident, plausible, and completely false root cause, which then got written
   down as fact.** Re-read the actual bytes on disk before concluding what a failure means.
4. Charge conservation on `assign_rdm1s` (one line): each orbital's weights sum to `1.0`
   (its `mf_info`-convention occupation) across the atoms.

Two of the three whole-pipeline invariants from above are done:

- `test_charge_conservation` — per-atom `charge_atom` sums to `0.0` for neutral water.
- `test_nuc_sum` — `nuc_att_glob + nuc_att_loc`, summed over all atoms, equals an
  independently-computed `Tr[V_nuc D]` (built via plain `np.sum(nuc * rdm1)`, not decodense's
  own `_trace`, and using `_h_core`'s `nuc` — already separately verified — as the operator).

The third, "`Elect.` equals the sum of its own parts," was **dropped after review, not
written** — real finding, worth recording so it isn't re-proposed later. `prop_atom`'s last
line for the energy case is literally `res[CompKeys.el] = sum(res.values())`; a test comparing
`res.el` against `res.coul + res.exch + res.kin + res.nuc_att_glob + res.nuc_att_loc` added by
hand would just be re-running that exact same addition outside the function. It cannot fail
short of floating-point summation-order noise — zero bug-catching power, fully circular, not
merely "semi-circular" like the earlier `_h_core`/`assign_rdm1s` cases. A genuinely strong
version exists (compare `res.el.sum()` against PySCF's own `mf.e_tot - mf.energy_nuc()`, a
completely independent computation) but was also skipped: that exact check is already what the
7 existing `unittest` integration files do (`self.assertAlmostEqual(mf_e_tot, e_tot, TOL)` in
e.g. `test_ch2_hf_ndo.py`), on several real molecules — adding it here would be duplicate
coverage of something already protected, not a new gap closed.

Next, not yet started, ranked:

1. `_make_rho`'s GGA/MGGA branches, and `_vk_dft` — found during a critical re-read of
   `properties.py`: we only tested `xc_type="LDA"`, but every real example uses `pbe0`,
   `b3lyp`, or `wb97m_v` (GGA/hybrid/meta-GGA) — the untested branches are the ones actually
   used in practice, the LDA one essentially never is. Needs a proper DFT fixture (`dft.RKS`
   with `.xc` set and `.run()` actually called — an earlier attempt at this fixture was
   abandoned mid-way, missing `.run()`).
   - `_make_rho`'s GGA branch: **done** (`test_make_rho_gga`, `mf_h2o_dft` fixture,
     `pop_method`-agnostic, differential against PySCF's `numint.eval_rho(..., xctype="GGA")`
     with `deriv=1` AOs). MGGA still untested.
   - `_vk_dft`: **done and confirmed passing** (`test_vk_dft_pbe0`). Reference: PBE0's exact-exchange fraction is hardcoded as
     `0.25` (the functional's own definition, not looked up via any decodense or PySCF code)
     and compared against `0.25 * mf.get_k(mol, dm=rdm1)`. One honest, discussed-but-accepted
     limitation: `_vk_dft` internally calls `mf.get_k(mol=mol, dm=rdm1)` with the exact `rdm1`
     it's given, and the test reuses that same `rdm1` value for its own reference call — so
     the test does not independently confirm `_vk_dft` actually *uses* the `rdm1` it was
     passed rather than some other value. Judged low-risk in this specific function (`rdm1`
     is used in exactly one place, with no alternate source to fall back to), but flagged
     here rather than silently assumed away.
2. `assign_rdm1s` extended to `pop_method="iao"` — **done** (`test_assign_rdm1s_h2o_iao`,
   alongside the original `test_assign_rdm1s_h2o` for `mulliken` — kept both, not replaced).
3. `sanity_check`'s parametrized table — **done**. Now covers 9 `ValueError` choice rows,
   5 `TypeError` rows, `gauge_origin` shape/type, `verbose` boundaries on both sides,
   `mo_coeff`/`mo_occ` validation, and a valid-input baseline. Only the 4 PBC-specific checks
   remain, deliberately.
4. `main()`'s own input-normalization branches (raw `ndarray` vs. pre-built tuple for
   `mo_coeff`) — **done** (`test_main_accepts_raw_mo_coeff_array`). Real bug caught while
   writing it: the `mf` argument was accidentally passed as `mf_h2o.mol` (the bare molecule)
   instead of `mf_h2o` (the SCF object) in one of the two calls, crashing with
   `AttributeError: 'Mole' object has no attribute 'get_jk'` -- easy mistake given `main()`'s
   `(mol, decomp, mf, mo_coeff, ...)` argument order, caught immediately by actually running it.
5. `pbctools.py` and `_solvent`/`_pcm`/`_point_charges` — **this blanket judgment was
   revised in Part II and that version supersedes this one.** Still out of scope: `_pcm`, the
   vendored `pbctools.py` builders and internals. Now judged *worth doing after all*, because
   they need far less machinery than assumed here: `ewald_e_nuc` (needs only a `Cell`;
   `cell.ewald()` is an independent reference), `_solvent`'s dispatch (testable with dummy
   objects carrying the relevant attribute, no real solvent needed), and `_point_charges`'
   nuclei–MM Coulomb loop (hand-computable analytically).

Also worth a critical note for future review: `test_e_xc_1d`/`test_e_xc_2d` use
`eps_xc = grid_weights = [1,1,1]` — this makes the test blind to a bug that silently dropped
`grid_weights` from the computation entirely (multiplying by all-`1`s is invisible). Should use
non-trivial weights in at least one case. `test_make_rdm1_symmetric` is also weaker than it
looks: `D = Σ_p n_p c_p c_p^T` is symmetric by construction for almost any implementation,
including several plausible wrong ones (e.g. a swapped contraction index) — cheap to keep, but
don't rely on it to catch a real bug.

Skipped, with reasons: `res_add`, `res_sub`, `make_natorb`, `write_rdm1` — exported in
`__init__.py` but zero internal callers; `_solvent`, `_pcm`, `_point_charges` — reached only
with solvent models; the `lowdin`, `meta_lowdin` and `becke` branches, and losslessness
(sum of parts == `mf.e_tot`) — already covered end-to-end by the 7 existing integration files;
the remaining 26 `sanity_check` raise-paths — one representative is enough.

One honest limitation recorded: an earlier `assign_rdm1s` test compared its output against the
same sequence of steps rebuilt by hand, and was trimmed back because both sides shared the
formula. `iao` was considered next on the grounds that 14 of 15 example scripts use it, but its
population formula (eq. 10) would still lean on the same PySCF `lo.iao` machinery the code
calls, so it is not the clean external reference it first appears to be.

---

### The week *(superseded — the file names here were never adopted)*

> **[out of date]** The file names here (`test_trace.py`, `test_properties_pure.py`,
> `conftest.py`, …) describe a layout that was not adopted — everything went into a single
> `tests/test_unit.py`. The *sequencing* logic is still sound; the file names are not.

Five sessions. Each ends with something committed and green.

**Day 1 — mechanics, on the easiest possible target**
`pip install pytest`; confirm it runs the existing suite unchanged. Write
`tests/test_trace.py` by hand. **Then break `_trace` on purpose** (swap the indices in the
contract string), watch it go red, restore — the OTR rule. Add `_e_xc` and `_xc_ao_deriv`;
for the latter, write the unknown-functional case and see what actually happens.
→ ~6 tests, well under a second.

**Day 2 — properties and differential tests, still no SCF**
`_e_nuc` vs `mol.energy_nuc()` plus the analytic H₂ check; `_get_nuc` vs
`mol.intor("int1e_nuc")`; `_dip_nuc` gauge-origin property; `make_rdm1` symmetry and
trace; `mf_info` dtype — expect this one to **fail**, and we decide together whether the
int64 occupations are a real bug.
→ `test_properties_pure.py`, `test_tools.py`.

**Day 3 — fixtures, parametrize, and the most valuable test in the plan**
First `conftest.py` with a module-scoped water RHF fixture. The partition-of-unity test
for `assign_rdm1s`, parametrized over all five methods with per-method tolerances. Re-add
the `comp_key_dict` inverse test. The NDO index test with **odd** and even counts.
→ `conftest.py`, `test_orbitals.py`, `test_results.py`.

**Day 4 — error paths and output structure**
`pytest.raises` tests for `sanity_check`. Structural tests on the `atoms()`/`orbs()`
DataFrames: row counts, index name, columns, no NaN. Tidy the integration tests' setup.
→ a suite that runs in seconds minus the integration tier.

**Day 5 — the density/XC chain, and a gap review**
`_make_rho_interm1` / `_make_rho_interm2` / `_make_rho` (`properties.py:668-764`) — pure
array functions with LDA/GGA/MGGA branches chosen by a *string*, so a wrong branch
silently integrates the wrong thing. Test shapes first, then values; optionally
differential against PySCF's `numint.eval_rho`. Then: re-run the mutation check on a
sample of the week's tests, and write the remaining gaps into `PLAN.md`.
→ `test_properties_rho.py`, updated PR 12 entry.

**Not this week:** the AD branch; CI wiring (PR 0); coverage measurement; `pbctools.py`
(PR 11); per-atom regression baselines (PR 1, awaiting supervisor input).

**Final step, once the above is done: revisit every function again, not just the ones
already picked.** The Part VI's original scoping file-by-file decision was made up front, based on the best
judgement available at the time — go back through all ~45 functions across the package
once more at the end and make the final call on each: does it stay excluded, or does it
turn out to be worth a test after all? Some exclusions may not hold up once the rest of
the suite exists to compare against; better to check deliberately than assume the first
pass got every call right.

---

### Layout and rules *(layout superseded; the **Rules** at the end still apply)*

> **[out of date]** The multi-file layout below was not adopted. Actual layout: one
> `tests/test_unit.py` holding all unit tests and their fixtures, alongside the 7 unchanged
> integration files. The **Rules** at the end of this section do still apply.

decodense is one small package, so PySCF's per-module `test/` dirs are overkill. Use
PySCFAD's newer convention — one top-level `tests/`, one file per source module — with
OTR's role split (new unit tests separate from the existing integration tests):

```
tests/
  conftest.py                     # shared fixtures
  test_tools.py
  test_orbitals.py
  test_properties_pure.py         # Tier 0/1
  test_properties_rho.py          # density/XC chain
  test_results.py
  test_decomp.py                  # sanity_check, CompKeys, comp_key_dict
  geom/                           # unchanged
  test_h2o_wb97m_v_energy_gs.py   # existing integration tests, unchanged
  ...
```

**Rules**

1. **One behaviour per test** — not one per function. `_trace` gets three (2-D branch,
   3-D branch, scaling), not one test with three asserts.
2. **The test name states the claim.** `test_weights_sum_to_occupation`, not `test_1`.
3. **Never build the expected value by calling the code under test.**
4. **A unit test calls exactly one production function** (OTR's rule).
5. **Every hardcoded number gets a provenance comment** — `# PySCF 2.14.0, RHF/STO-3G,
   conv_tol=1e-12`. PySCF's convention; it is what makes a failure diagnosable later.
6. **Every new test must be shown to fail.** Break the code, watch it go red, restore.
7. **Justify non-obvious tolerances in a comment** — `# becke is grid-integrated, ~1e-7`.

**Open questions to settle next week**

- `mf_info`'s int64 occupations: real bug or harmless?
- `_trace`'s symmetry assumption: assert it, document it, or generalise the function?
- Should `res_add`/`res_sub` raise on mismatched keys instead of ignoring them?
- Should integration tests also assert *per-atom* values against stored references? That
  is PR 1, still awaiting your supervisor.

---

# Part VII — PR #28 restructuring: function audit

**Context.** PR #28 (`t-charlotte`, "Restructuring of decodense + crash fixes",
`eriksen-lab/decodense`) replaces `decomp.part="eda"` with `part` + a new `part_method`, moves the
decomposition logic out of `decodense.py` into a new `schemes.py`, and removes `charge_atom` from
the results. It was agreed in a meeting that PR #28 merges **before** `modernize-test-suite`, so the
test suite has to be adapted to it. This part re-audits every function against PR #28. The work
lives on branch `adapt-tests-to-restructuring` (created from `review-pr28`).

**Method [verified].** Compared every function and class between `main` and PR #28 using Python's
`ast` module. `main` is fully contained in PR #28 (`git merge-base --is-ancestor`), so every
difference found is the PR's own change, not something `main` has that the PR lacks.

**Branch caveat.** This branch carries two local fixes to PR #28's *own* code so the tests can run:
`properties.py:611` (`np.zeros(len(mol.natm))` → `np.zeros(mol.natm)`) and `results.py`
(removing the leftover `CompKeys.charge_atom` references). Both are reported as review comments on
PR #28 and are `t-charlotte`'s to fix — **do not commit them as part of the test PR.**

## VII.1 Unchanged by PR #28 — earlier verdicts (Part II) stand

`_e_nuc`, `_dip_nuc`, `_h_core`, `_get_nuc`, `_make_rho`, `_make_rho_interm1`, `_make_rho_interm2`,
`_trace`, `_e_xc`, `_xc_ao_deriv`, `_pcm`, `_ao_val`, `_population_mul`, `_population_becke`,
`assign_rdm1s.get_weights`, `dim`, `make_rdm1`, `make_natorb`, `write_rdm1`, `orbsym`, `contract`,
`git_version`, `logger_config`, `DecodenseLogger`, `_unit_scaling`, and most of `pbctools.py`.

## VII.2 New tests to write, in priority order

| # | Function | Test | Reference / evidence |
|---|---|---|---|
| ~~1~~ | ~~`tools.mf_info` (ROHF fix)~~ — **done.** `mf_li_rohf` fixture (`scf.ROHF`) added alongside the existing UHF `mf_li`; `test_mf_info_li_rohf` checks alpha=2, beta=1. **[verified]** passes; would have failed alpha=2, beta=2 under the old code. | | |
| ~~2~~ | ~~`properties.prop_orb` (new `(i, j, m)` domain)~~ — **done.** `test_main_orbitals_non_aufbau`: fresh UHF water calc, alpha orbital 4→5 (a real gap, not the adjacent orbital), full `mf.mo_coeff` + the gapped `mo_occ`, `part="orbitals"`. **[verified]** passes, `res.el[0].size == 5`. Old code would have crashed storing at real index `5` in a 5-slot array. | | |
| ~~3~~ | ~~`results.ResultsCls` / `to_dataframe` / `fmt`~~ — **done.** `test_to_dataframe`: `part="atoms"` gives 3 rows, `part="orbitals"` gives 10 rows (5 occupied alpha + 5 beta, via `mf_info`). **[verified]** passes with the local fix; the `charge_atom` bug would have crashed both calls, and no other test catches this since none of them call `to_dataframe()`/`print(res)`. | | |
| 4 | `results.orbs` | **Done.** `test_orbs_ndo` checks both odd NDO count (no orbital dropped; **[verified]** by mutation M18) and that `unit="ev"` actually scales the result (**[verified]** by mutation M19), combined into one test function by deliberate choice rather than split into two. | This is **B1**, previously deferred until `fix-comp-key-dict` merged — it has merged, so B1 is unblocked. The `* scaling` in the scalar branch is a real bug fix in PR #28. Hand-built `res` dict, no SCF. |
| ~~5~~ | ~~`decomp.DecompCls.__init__`~~ — **done.** `test_decomp_cls_part_method_defaults`: `part="eda"` → `("atoms", "ao")`; `part="atoms"` → `part_method="mo"`; `part="orbitals"` → `None`. **[verified]** by mutation M22 (which survived every other test in the suite, including all 7 legacy integration tests, since both schemes are lossless). | | |

## VII.3 Existing tests — updated and still to update

**Done:**
- `test_charge_conservation`: `res.charge_atom` no longer exists, so the charge is now computed in
  the test via `assign_rdm1s` (`charge = mol.atom_charges() - population`). Tested for
  `part_method="mo"` only — the `"ao"` scheme passes `weights=None` and never uses population
  weights, so this formula doesn't correspond to what it computes. **[verified]** passes.
- `test_atoms`: removed `CompKeys.charge_atom` from the hand-built input.

- `test_nuc_sum`: switched from `part` in `["atoms", "eda"]` to `part_method` in `["mo", "ao"]`.
  **[verified]** `test_nuc_sum[mo]` and `test_nuc_sum[ao]` both pass.
- `test_charge_conservation`: removed the now-unused `decomp = ...` line.
- Removed the unused imports `_solvent`, `comp_key_dict` from `test_unit.py`.

**Still to do:** nothing — all three quick updates are done. 67/67 tests pass.

## VII.4 Optional — cheap, lower priority

- `sanity_check`: add a `("part_method", "bad")` row to the existing parametrized table.
- `properties._solvent`: a non-PCM `with_solvent` should raise `NotImplementedError` (new check in
  PR #28). Testable with a dummy object; would use the currently unused `_solvent` import.
  **[hypothesis]** — not yet run.
- `tools.res_add` / `res_sub` / `_res_combine`: public API (exported in `__init__.py`), zero callers
  inside the repo. Hand-built dicts would catch the list-concatenation trap the new
  orbital-mode branch fixes.
- `orbitals._unique_filename`: new pure function; test with pytest's `tmp_path`.

## VII.5 Skipped, with reasons

- `schemes._scheme_bonds_a2b` / `_scheme_bonds_aap2b`, `results.bonds`: placeholders that only raise
  `NotImplementedError` — a test would just pin a placeholder.
- `pbctools.py`: only error-message text changed (`%s` placeholders → f-strings). Still skipped.
- `results.info`: text formatting only.
- `properties._vk_dft`: only an unreachable branch was removed. **[verified]** neither PySCF's
  standard SCF classes nor anything in decodense ever sets `mf.vk`, so `hasattr(mf, "vk")` is never
  true. The range-separated branch stays declined (B12).
- `prop_tot`'s `xc_spin` and NLC changes: inline in a closure, not callable on their own; the
  `wb97m_v` legacy test's total-energy check already covers NLC.
- `schemes._scheme_atoms_mo` / `_scheme_atoms_ao_orbitals`: thin wrappers every `main()` test runs.
- `assign_rdm1s`'s new verbose print / file output: wait until review question 3 below is answered,
  since the output may change.

## VII.6 Review comments posted on PR #28 (for reference)

1. `results.py` still references `CompKeys.charge_atom` (lines 150, 157, 204) → `print(res)` crashes.
2. `properties.py:611` `np.zeros(len(mol.natm))` → `TypeError` in both `point_charges` examples.
3. The PR description promises printed atomic *charges*, but the code only prints *population* —
   the nuclear charge is never subtracted. Asked where it's implemented.
4. The `energy` examples got the new `part_method` notation; the other examples and the 7 legacy
   test files didn't — asked whether to make it consistent.

**[verified] against the unmodified PR #28:** `test_unit.py` 64 passed / 4 failed; all 12 runnable
examples crash (10 from bug 1, the 2 `point_charges` ones from bug 2 first).

---

# Part VIII — Critical re-audit of the whole suite (2026-09-28)

**What was asked.** Go through every source file and function again: are the important ones
tested? Then review every existing test critically: is it a good test, does it make sense,
should it test more or test differently? Write it all down so the next session can decide.
**Nothing in `test_unit.py` or in `decodense/` was changed for this part.** It is input for a
decision, not a set of finished changes. Everything was measured on `adapt-tests-to-restructuring`
(PR #28 plus the 2 local fixes).

**Method.**
1. **Which code do the tests actually run? [verified]** No coverage tool is installed, and none was
   installed for this. Instead, a ~25-line pytest plugin in the scratchpad used Python's built-in
   `sys.settrace` to record every decodense function called and every line executed during the
   71 tests. (Raw "missed line" numbers include docstring lines, which the tracer doesn't count;
   those are ignored below.)
2. **Would the tests notice a bug? [verified]** Mutation testing, done by hand: 38 deliberate,
   realistic bugs, each planted in its own scratch copy of the repo. Against every copy, the full
   unit suite **and** all 7 legacy integration files were run. A "no-change" control copy passed
   everything, so every failure below is caused by the planted bug. The real repo was never
   touched. **Caveat:** I picked the bugs, and I aimed them at places I suspected were weak.
   So "18 of 38 caught" is a map of where the holes are, not an unbiased quality score.
3. Every test was read line by line against Part III's four questions: right thing? circular?
   meaningful asserts? realistic setup?

## VIII.1 Headline findings — read these first

Ranked by how much they matter. **Items 1–4 are not about writing tests at all**, but items 1–2
decide whether the tests run anywhere except your laptop.

1. **CI runs none of the 71 unit tests [verified locally].** The CI branch
   (`codebase-rev-ci-tooling`, PR #21) runs `python -m unittest discover -s tests -p "test_*.py"`.
   `unittest` only collects `unittest.TestCase` classes. Our tests are plain pytest functions, so
   `python -m unittest test_unit` prints `Ran 0 tests ... NO TESTS RAN`. **Second problem:** from the
   repo root (which is what `-s tests` means; unittest does *not* change directory), every legacy
   file fails to import with `Unsupported atom symbol GEOM/C5H5N.XYZ`, because the `geom/...` paths
   are relative to `tests/`. From inside `tests/`, all 7 legacy tests pass. **Third problem:**
   switching CI to plain `pytest` is not a drop-in fix. The legacy files break when pytest runs
   them in one process: each module's `tearDownModule` calls `mol.stdout.close()`, and with
   `output=None` that is `sys.stdout`, so every later file errors with
   `I/O operation on closed file` [verified: 1 passed, 13 errors]. **Not checked:** what the live
   GitHub Actions run on your fork shows. Look at that first. Options are in VIII.7.

2. **`import decodense` should fail on Python 3.11, which CI tests [NameError verified; the 3.11
   failure is inferred, not run].** `results.py:135` has `def atoms(..., res: Dict[str, Any], ...)`,
   but `Dict` is never imported (only `from typing import Any, Optional`). Verified:
   `decodense.results.atoms.__annotations__` raises `NameError: name 'Dict' is not defined`.
   On Python 3.14 (your machine), annotations are evaluated lazily (PEP 649), so nothing notices.
   On Python ≤ 3.13, annotations are evaluated when the `def` line runs, i.e. at import time.
   That is established language behaviour, so the whole package would fail to import there. No
   3.11 interpreter is installed here, so this was not actually run. **History:** `3337c21`
   removed the `Dict` import and modernised the annotations. `d7d45c6` (the `fix-comp-key-dict`
   merge) brought the old `Dict` annotation back. It is on `main` and `upstream/main`. The fix is
   one word (`Dict` → `dict`) and belongs in its own tiny PR, not the test PR.

3. **`results.info()` always crashes. New in PR #28 [verified].** Line 83 says `strin  += ...`, a
   typo for `string`, so every call raises `UnboundLocalError`. `info` is public (it's in
   `decodense/__init__.py`), but nothing in the repo or the examples calls it, so no test or example
   notices. → a **5th review comment for PR #28**.

4. **`write_rdm1` writes wrong per-atom densities when occupations aren't 1 [verified;
   pre-existing on `main` and `upstream/main`].** Each orbital's density `make_rdm1(orb, n)`
   already contains the occupation n. `write_rdm1` then multiplies it by the population weight,
   which *also* contains n (`get_weights` computes `mocc * mo * S·mo`). `prop_atom_mo` avoids this
   by dividing by `np.sum(weights[i][m])`; `write_rdm1` doesn't divide. So the per-atom densities
   add up to Σ n² instead of Σ n. Measured on water (Mulliken): with occupations 1.0,
   Σ_A Tr(D_A S) = 10.0 ✓. With occupations 0.5, it's **2.5 instead of 5.0**. Normal SCF is
   unaffected, because `main` turns an RHF occupation of 2 into 1 per spin. NDOs are affected, and
   `sanity_check` allows `write` together with `ndo=True`. The partition test "Σ_atoms D_A = D"
   would catch this (N5 in VIII.5).

5. **The legacy integration tests only ever check the grand total.** All 7 assert
   `sum(res.tot) == mf.e_tot` (or the dipole). So a bug that *moves* energy between atoms while
   keeping the sum is invisible to them, and that is the most typical kind of bug in a
   decomposition code. The mutation run confirms it: M03, M25, M26 and M22 all pass every legacy
   test. On the unit side, per-atom correctness of the MO scheme rests on `test_permute` alone,
   which misses both "wrong atom" mutants (M25, M26; see VIII.3 for why). **The AO scheme's per-atom
   values are not checked by anything.**

6. **The unit suite never checks losslessness either.** No unit test compares a decomposed total
   with `mf.e_tot`. Part VI decided not to duplicate the legacy tests, which is fine **only while
   the legacy tests actually run in CI** (item 1). Five energy-changing mutants are caught *only*
   by legacy tests (M10, M27, M28, M29, M36). M27 (wrong exchange factor, restricted MO scheme) is
   caught by just one legacy file, `wb97m_v`, because it's the only closed-shell SCF among them.

7. **Two tests are vacuous for the molecule they use.** They look like strong invariant tests but
   can't fail.
   - `test_assign_rdm1s_li` / `_li_iao`: with a **single atom**, every weight is trivially 1, so
     "rows sum to 1" cannot fail. Only the shape asserts carry information (they did catch M13 and
     M34).
   - `_e_nuc`: a **diatomic** always gives each atom exactly half of its one pair energy. No
     diatomic (not H₂, **and not the LiH test proposed last session**) can tell the correct split
     from a wrong one such as "split the total equally". This needs ≥ 3 atoms at unequal
     distances. Linear H₃ at z = 0, 1, 3 bohr gives per-atom values
     [(1+1/3)/2, (1+1/2)/2, (1/3+1/2)/2] = [0.6667, 0.75, 0.4167], hand-derived and **[verified]**
     against `_e_nuc`. An "equal split" bug would give 0.6111 for each.

8. **Correction to an earlier argument (Part I, B14).** B14 said swapping the alpha and beta
   densities in `_make_rho` "would silently corrupt XC energy for every open-shell system".
   **Mutation M09 shows it doesn't change anything:** all 71 unit tests and all 7 legacy tests pass,
   including 3 open-shell legacy systems. The reason is physics: an XC functional doesn't change if
   you relabel which spin is "up", and in decodense the unrestricted `rho` only feeds
   `eval_xc` for the energy density. The per-atom `rho` is built from the spin-summed density, so it
   goes down the restricted path. So M09 is an *equivalent mutant* (a code change that can't change
   any result). B14's decision (declined) stands, but for this reason, not the recorded one.

9. **Correction to Part I, A4 (`_trace` transposition).** A4 argued that orientation can't matter
   because every density matrix is symmetric. **Mutation M06 (transpose inside `_trace`) is caught
   anyway**, by `test_nuc_sum[ao]` and by 6 legacy files. The AO scheme passes *row-sliced*,
   non-square matrices (`vj[select]`, `rdm1[select]`), so a transposed contraction fails with a
   shape error. A4's conclusion (no extra test) still holds. The reason is different: the AO
   scheme already covers it, by crashing.

10. **`Struct.` is a per-*atom* array inside the per-*orbital* result dict [verified; pre-existing,
    not a PR #28 issue].** Found on 2026-09-29 while checking whether `test_permute` should also
    cover `part="orbitals"`. In orbitals mode every other key is a list of two per-orbital arrays
    (`[alpha(5,), beta(5,)]` for water/STO-3G), but `prop[CompKeys.struct]` is set to
    `prop_nuc_rep`, i.e. plain `_e_nuc(mol)` of shape `(natm,)` — **[verified]** identical to
    `_e_nuc(mol1)`, and it matches between the two permuted molecules only after applying the
    *atom* permutation `[2,0,1]`. It is never displayed: `orbs()` lists `CompKeys.struct` in its
    exclusion tuple, so the orbitals DataFrame stops at `Symm.` [verified], and `Total` is set to
    `Elect.` alone (nuclear repulsion deliberately excluded, matching the legacy tests' orbitals
    reference `mf_e_tot - mol.energy_nuc()`). So nothing is *wrong* in the output — but
    `res.struct` is still exposed on `ResultsCls` in orbitals mode as a silently per-atom quantity.
    **[verified]** identical code on `main` and `upstream/main` (`properties.py:462–464`), so this
    predates PR #28. → question for Janus, not a PR #28 review comment.

11. **Decided against: an orbitals variant of `test_permute`** (2026-09-29), so it isn't
    re-proposed. **[verified]** all per-orbital keys come out *identical* (no permutation needed)
    between the two atom-orderings, since orbitals are ordered by energy, not by atom. But the
    reason not to test it is stronger than "it passes": `prop_orb` builds each orbital's
    contribution from *global* matrices (`vj`, `vk`, `kin`, `nuc`) plus that orbital's own rdm1,
    and never touches an atom index at all — no population weights, no per-atom AO slicing. An
    atom-mislabelling bug therefore cannot manifest in orbital-wise results, and the test would
    really be checking that PySCF returns the same MOs regardless of input atom order. That is
    PySCF's job, not decodense's.

## VIII.2 Coverage map — every source function, measured

✅ its own logic is checked by a unit test · 🟡 runs (often only via `main`), but the logic isn't
checked, or important branches never run · ❌ never runs in the unit suite · ⛔ not worth testing
(reason given before). "Mutants" refers to VIII.3.

#### `decodense.py`

| Function | State | What is and isn't run | Proposal |
|---|---|---|---|
| `main` | 🟡 | Runs in 5 tests. **Never run:** the 3-D-array `mo_coeff` branch (l.51; PySCF's UHF returns a *tuple* here [verified], so a real input needs `np.asarray(mf.mo_coeff)`), the whole 1-D `mo_occ` normalisation (l.56–67: RHF 2/0 → 1/1, the new ROHF rule "singly occupied → alpha", and the fractional "split in half" rule), and the 2-D `mo_occ` array (l.69). M14, M15 and M35 get past **everything**. | **N2**: route-equivalence test. Raw `mf.mo_coeff`/`mf.mo_occ` in, compared against the `mf_info` route (separate code implementing the same convention). |

#### `decomp.py`

| Function | State | What is and isn't run | Proposal |
|---|---|---|---|
| `CompKeys` | ⛔ | constants | — |
| `DecompCls.__init__` | 🟡 | The deprecated `part="eda"` → `("atoms","ao")` mapping never runs in the unit suite. The legacy tests use `"eda"`, but **M22 (mapping to the MO scheme instead) passes every legacy test**: both schemes are lossless, so the totals agree. | **N1** (= VII.2 #5) |
| `sanity_check` | ✅ | 53/72 statements. Never run: the orbitals-with-`part_method` reset (l.188–193, M23 survives), the `bonds` branch, the `write` format check, the new "write only for atoms/mo" check (M33 survives), and the 4 PBC checks (deliberately skipped). **Dead code:** the PBC check still lists `"eda"` (l.263), which can't reach it any more because `DecompCls` converts it. | **N7**: new table rows. See also T6 in VIII.4. |

#### `schemes.py`

| Function | State | Notes |
|---|---|---|
| `_scheme_atoms_mo` | 🟡 | Runs via `main`; the `write` branch never runs. |
| `_scheme_atoms_ao_orbitals` | 🟡 | Runs via `main`. |
| `_scheme_bonds_*`, `SCHEMES` | ⛔ | Placeholders (VII.5). |

#### `orbitals.py`

| Function | State | What is and isn't run | Proposal |
|---|---|---|---|
| `assign_rdm1s` | ✅ | 38/80. `mulliken` + `iao`, RHF + UHF. Never run: `lowdin`, `meta_lowdin`, `becke` (the legacy tests run all 5, totals only), the `ndo`+`iao` raise, **the new verbose print + file output (l.155–195)**. **Note:** with `verbose > 0` it writes `pop_weights_*.txt` into the *current working directory*, so any test of it must `monkeypatch.chdir(tmp_path)`. | Fix the Li tests (T27); N12 optional |
| `get_weights` | ⛔ | Closure; runs via parent. | — |
| `_unique_filename` | ❌ | 0/13. Pure function, new. | N12 optional |
| `_population_mul` | ✅ | M12 caught by 3 tests. | — |
| `_population_becke` | ❌ | Legacy only. | as before |

#### `properties.py`

| Function | State | What is and isn't run | Proposal |
|---|---|---|---|
| `prop_tot` (+ `prop_atom_mo`, `prop_atom_ao`, `prop_orb`) | 🟡 | 102/195. **In the unit suite, `main` only ever runs HF energies**: no DFT, no dipole, no solvent, no NDO. The AO scheme runs only in `test_nuc_sum[ao]`, which checks only a sum. The legacy tests cover DFT, dipole, unrestricted and NDO, **totals only**. | N3 (symmetry test), N9 (dipole), decision on losslessness (VIII.7) |
| `_e_nuc` | ✅ | But see VIII.1 #7: the H₂ test can't check the split. | N6 (linear H₃) |
| `_dip_nuc` | ✅ | M04 caught. | — |
| `_h_core` | ✅ | ECP raise and PBC branch not run. M30 caught. | — |
| `_get_nuc` | ✅ | M05 caught by 5 tests. | — |
| `_solvent` | 🟡 | Only the "no solvent" path runs. M32 (PCM never dispatched) gets past everything. Only the solvent *examples* would notice. | N10 (VII.4) optional |
| `_point_charges` | ✅ | M31 caught. | rename T18 |
| `_pcm` | ❌ | declined (Part II) | — |
| `_xc_ao_deriv` | ✅ | | T21 comment |
| `_make_rho_interm1/2`, `_make_rho` | ✅ | LDA + GGA, restricted. MGGA and the unrestricted branches never run (B12/B14 declined; B14's reason corrected in VIII.1 #8). | — |
| `_vk_dft` | ✅ | Global hybrid only. The non-hybrid branch (`vk = 0`) and the range-separated branch never run in the unit suite. M10 is caught by 2 legacy files only. | as before (B12) |
| `_ao_val` | ⛔ | | — |
| `_trace` | ✅ | M06 caught (VIII.1 #9). | — |
| `_e_xc` | ✅ | M07 caught (A2 was done). | — |

#### `results.py`

| Function | State | What is and isn't run | Proposal |
|---|---|---|---|
| `ResultsCls.__init__` | ✅ | via `main` | — |
| `__str__` | ❌ | 0/2. **`print(res)` is how all 12 examples use the result**, and it's exactly where PR #28's `charge_atom` crash lived. | add `str(res)` to T29 |
| `to_dataframe`, `fmt` | ✅ | atoms + orbitals | — |
| `_unit_scaling` | 🟡 | Only `"au"` ever runs. M19, M20, M21 get past everything. | **N4** |
| `atoms` | 🟡 | Scalar case, `"au"` only. The dipole columns (`" (x)"` …) never run. | **N4**, T28 |
| `orbs` | ✅ | NDO + non-NDO scalar. Dipole branch never runs. M38 (symmetry labels not reordered) survives because the test uses identical labels. | T30 |
| `info` | ❌ | **Broken** (VIII.1 #3). | review comment |
| `bonds` | ⛔ | placeholder | — |

#### `tools.py`

| Function | State | Notes / proposal |
|---|---|---|
| `logger_config`, `contract` | ✅ implicitly | — |
| `git_version`, `DecodenseLogger.info2/3` | ⛔ | — |
| `dim` | ✅ | M34 caught by 5 tests. |
| `mf_info` | ✅ | RHF/UHF/ROHF; M16 caught. |
| `orbsym` | 🟡 | The symmetry path runs via `main`. The fallback branches never run (B15 skipped). |
| `make_rdm1` | ✅ | M01 caught. |
| `make_natorb` | ❌ | Zero callers anywhere: not in the repo, not in the examples, not in the legacy tests (they build NDOs themselves). |
| `write_rdm1` | ❌ | **Has a bug** (VIII.1 #4). Zero callers in tests or examples. **N5** |
| `_res_combine`, `res_add`, `res_sub` | ❌ | Zero callers anywhere. M37 (occupations get added) survives. Optional (VII.4). |

**Options that are validated but never used [verified by grep].** `mo_basis`, `mo_init` and
`loc_exp` are checked by `sanity_check` and printed by `info()`, but no computation reads them.
Three rows of the sanity-check table therefore guard options that do nothing. They're cheap to
keep; whether these options should exist at all is a question for Janus, not a testing question.

## VIII.3 Mutation results

Each row is one planted bug. **Unit** = which `test_unit.py` tests went red. **Legacy** = which of
the 7 integration files went red. Survivors in **bold**.

| # | Planted bug | Unit tests that caught it | Legacy files that caught it |
|---|---|---|---|
| M00 | *control: no change* | — (71 pass) | — (7 pass) |
| M01 | `make_rdm1` ignores occupations | `test_make_rdm1`, `_equal_electron_count` | both NDO files |
| M02 | `_e_nuc` double-counts pairs | `test_e_nuc`, `test_e_nuc_h2` | 3 |
| M03 | `_e_nuc` puts all repulsion on atom 0 (total unchanged) | `test_e_nuc_h2` only | none |
| M04 | `_dip_nuc` gauge-origin sign | `test_dip_nuc`, `_h2o` | none |
| M05 | `_get_nuc` uses atom 0's position for every atom | 5 tests | 5 |
| M06 | `_trace` transposes | `test_nuc_sum[ao]` (shape error) | 6 |
| M07 | `_e_xc` drops grid weights | `test_e_xc_1d` | 4 |
| M08 | GGA gradient factor 2 lost | `test_make_rho_gga` | 4 |
| M09 | `_make_rho` swaps α/β densities | none | none (**equivalent mutant**, VIII.1 #8) |
| M10 | range-separated exchange sign | none | `camb3lyp`, `wb97m_v` |
| M11 | `_vk_dft` ignores the `rdm1` it's given | `test_vk_dft_pbe0` (thanks to A9's 0.5 factor) | none |
| M12 | `_population_mul` puts atom 2's AOs on atom 0 | `test_charge_conservation`, `test_permute`, `test_population_mul_h2o` | none |
| M13 | `assign_rdm1s` always copies α into β | both Li tests (shape asserts) | none |
| **M14** | `main`: 1-D ROHF occupations put the unpaired electron in β too | **none** | **none** |
| **M15** | `main`: fractional restricted occupations not halved | **none** | **none** |
| M16 | `mf_info` ROHF β | `test_mf_info_li_rohf` | none |
| M17 | `prop_orb` stores at the real MO index (the PR #28 bug) | `test_main_orbitals_non_aufbau` | none |
| M18 | `orbs` NDO drops the unpaired orbital | `test_orbs_ndo_odd_count` | none |
| **M19** | `orbs` forgets unit scaling | **none** | **none** |
| **M20** | `atoms` forgets unit scaling | **none** | **none** |
| **M21** | `_unit_scaling`: eV uses the kcal/mol constant | **none** | **none** |
| **M22** | `DecompCls` maps deprecated `"eda"` to the MO scheme | **none** | **none** (totals are lossless either way) |
| **M23** | `sanity_check` no longer resets `part_method` for orbitals | **none** | **none** |
| M24 | MO scheme: local nuclear-attraction factor 0.5 → 0.6 | `test_nuc_sum[mo]` | 6 |
| **M25** | MO scheme: every atom gets the *next* atom's orbital weights (total unchanged) | **none** | **none** |
| **M26** | AO scheme: every atom gets the *next* atom's AOs (total unchanged) | **none** | **none** |
| M27 | MO scheme restricted exchange factor 0.25 → 0.5 | none | `wb97m_v` only |
| M28 | XC always evaluated spin-restricted | none | `ch2_pbe0`, `camb3lyp` |
| M29 | NLC energy density doubled | none | `wb97m_v` only |
| M30 | kinetic integrals × 1.01 | `test_h_core` | 5 |
| M31 | point-charge Coulomb uses 1/r² | `test_point_charges_nuc_solv_shape` | none |
| **M32** | `_solvent` never dispatches to PCM | **none** | **none** |
| **M33** | `sanity_check` allows `write` with the AO scheme | **none** | **none** |
| M34 | `dim` counts empty α orbitals | 5 tests | none |
| **M35** | `main`: 3-D `mo_coeff` array loses β | **none** | **none** |
| M36 | `prop_orb` unrestricted uses α exchange for β | none | 4 |
| **M37** | `_res_combine` adds occupations | **none** | **none** |
| **M38** | `orbs`: symmetry labels not reordered with NDO sorting | **none** | **none** |

**Summary:** unit 18/38 · legacy-only 5 · survived everything 15 (one of them, M09, is equivalent).

**Follow-up on M25/M26: which test *would* catch them? [verified]** Planted each again and tried two
candidate checks:

| Check | M25 (MO scheme) | M26 (AO scheme) |
|---|---|---|
| `test_permute` strengthened: compare **all** components, both schemes | passes (misses it) | passes (misses it) |
| **Water's two H atoms get identical values for every component** | **fails (catches it)** | **fails (catches it)** |

Why the permutation test misses them: `h2o_permute.xyz` reorders the atoms *cyclically*
(O,H,H → H,H,O). A cyclic reordering commutes with a cyclic "next atom" bug, so both runs are wrong
in the same way. A swap-type permutation file (e.g. H,O,H) would probably catch it too; that was
not tried. The symmetry check is simpler and doesn't depend on which permutation you choose, so it
is proposed as **N3**. Honest limit: it only catches bugs that break the H₁/H₂ equivalence. A bug
that moves energy from O to *both* H's equally would still pass. No fully external per-atom
reference exists for decodense's own partition, so this is as good as an invariant gets.

## VIII.4 Every test, reviewed critically

Verdicts: **K** keep as is · **S** strengthen (better or more asserts) · **R** restructure
(rename/move/rewrite) · **D** delete or merge. The reference column says where the expected value
comes from. That's the rule from the top of this document.

| # | Test | Reference | Verdict | What exactly, and why |
|---|---|---|---|---|
| T1 | `test_permute` | invariant | **S — done (2026-09-29)** | Now loops over **all 8 keys** of `res.res_dict` (was `el` only), is **parametrized over `part_method` in `["mo","ao"]`** (the AO scheme previously had no per-atom check anywhere in the suite; **[verified]** the invariant holds for both), and actually **passes `mo_occ`** (`mo_occ1`/`mo_occ2` were computed and thrown away, so the test silently relied on `main`'s `mo_occ=None` fallback — that branch is still covered by `test_main_accepts_raw_mo_coeff_array`). The explicit `val[0]/val[1]/val[2]` form was kept over `perm`-style fancy indexing, by preference for the plain version. **Remaining limitation, unfixable here:** still blind to M25/M26 — needs N3. **Correction to yesterday's note:** `decomp_perm` is *not* redundant; `main()` writes its result onto the `decomp` object (`decomp.res = ...`), so reusing one object across two calls would create an aliasing dependency. Two objects is the defensive choice — left as is. |
| T2 | `test_main_accepts_raw_mo_coeff_array` | route equivalence | **R — (a)+(b) done (2026-09-29), (c) open** | (a) **fixed:** `(mf_h2o.mo_coeff)` — parentheses don't make a tuple; now passes `mo_coeff[0]` (alpha) vs the explicit `(alpha, alpha)` tuple, so the two routes being compared are actually visible. (b) **fixed:** now uses `mf_info`'s *occupied* columns `(7,5)` plus `mo_occ`, instead of the raw `(7,7)` array with `mo_occ=None` which made it a fictional 14-electron water. (c) **mostly done:** the 1-D `mo_occ` normalisation is now covered by a new test, `test_main_accepts_raw_pyscf_input` — hand `main` PySCF's raw `(mo_coeff, mo_occ)` and require it to land where the `mf_info` route lands. **[verified] it kills M14**: with the beta rule changed from `occ > 1` to `occ > 0`, this is the *only* test in the repo that goes red (previously M14 survived all 74 unit tests and all 7 legacy files). **Still open:** the 3-D `mo_coeff` branch (M35) and fractional occupations (M15). A parametrized route table was drafted and **rejected as too complex**; a `request.getfixturevalue` parametrization over molecules was also considered and dropped in favour of just using lithium, which strictly dominates water here. The agreed shape is one small, separate, single-idea test per route. **[verified] all routes agree:** 2-D array, 3-D array, and raw PySCF `(mo_coeff, mo_occ)` all reproduce the `mf_info` route exactly, for water. **[verified]** water cannot catch M14 — all its occupied orbitals hold 2.0, so `occ>0` and `occ>1` coincide; ROHF lithium (`mo_occ = [2,1,0,0,0]` → alpha `[1,1,0,0,0]`, beta `[1,0,0,0,0]`) is the case that separates them, and the `mf_li_rohf` fixture already exists. |
| T3 | `test_main_orbitals_non_aufbau` | route equivalence | **S — done (2026-09-29)** | Was `size == 5` only: it confirmed nothing crashed, but a result set of the right *size* in the wrong *positions* passed happily. Now also compares `res.el` elementwise, for both spins, against the **sliced route** (hand `main` only the 5 occupied columns with occupation 1) — an independent input route to the same physics, **[verified]** identical. The magic indices `4`/`5` got their explanatory comments back. **[verified] it kills the bug class it exists for:** feeding `prop_orb` the position `m` instead of the real MO index `j` (`domain[:, :2]` → `domain[:, [0, 2]]`, i.e. exactly what PR #28 fixed) now fails this test and **only** this test; the old `size` assertion passed it. **Honest limitation, [verified]:** the sliced-route reference is blind to bugs in code *shared* by both routes — reversing the storage order in `prop_tot`'s collection loop (`[domain[k, 2]]` → `[-1 - domain[k, 2]]`) reverses both sides equally and the whole suite still passes. It only catches what differs between the routes, which is the `j` vs `m` distinction. |
| T4 | `test_charge_conservation` → **`test_rdm1_charge_conservation`** | chemistry (signs) + invariant | **R — done (2026-09-29)** | Renamed and **moved out of the `decodense.py` section into the `orbitals.py` section**, beside the other `assign_rdm1s` tests — it stopped calling `main` when it was adapted to PR #28, so it was filed under the wrong module. Asserts unchanged; a comment now records that the expected signs are chemistry (oxygen more electronegative → negative; both H positive; the two H equal by symmetry) and that the total-charge line is **structurally implied, not load-bearing**: each of the 10 occupied spin-orbitals carries weights summing to 1, so the population is always exactly 10, and water's Σ Z is also 10 — and the partition of unity is already asserted directly in `test_assign_rdm1s_h2o_iao`. **[verified]** real values: population `[8.4716, 0.7642, 0.7642]`, charges `[-0.4716, +0.2358, +0.2358]`. Kept as documentation of the named property. |
| T5 | `test_nuc_sum[mo/ao]` | independent recomputation | **K — addition declined by user** | The test as it stands is strong and unchanged: the combined `Tr[V·D]` is never computed anywhere inside decodense (only the two halves are), so recomputing it with plain numpy — deliberately not decodense's `_trace` — is a genuine external reference. It catches M24. **Proposed and then declined:** adding `res.nuc_att_glob.sum() == res.nuc_att_loc.sum()`. It was written, verified and then removed at the user's request. Recorded so it isn't re-proposed. For the record, the reasoning was sound: the existing assert pins only `glob + loc`, so a compensating split (0.6 glob / 0.4 loc) passes it; **[verified]** such a mutation in the MO scheme is caught by the extra assert and by nothing else in the suite. The identity is exact, not a water coincidence — both sums collapse to ½Tr(V·D) (`Σ_A sub_nuc[A]` is `nuc`; `Σ_A D_A` is `D` because the MO scheme divides by `np.sum(weights[i][m])`, and the AO scheme's `[select]` row-slices partition all AO rows exactly once) — **[verified]** to 10 decimals on water, LiH and the OH radical, in both schemes. |
| T6 | `test_sanity_check_rejects_invalid_attr` (now 16 rows + 1 new test) | spec | **S — done (2026-09-29)** | (a) **Removed** the duplicate `("pop_method", "not_a_real_method")` row — same branch, same message as `("pop_method", "bad")`. (b) **Added** `("part_method", "bad")` → `"invalid partitioning method"` and `("write", "bad")` → `"invalid write format"`; both **[verified]** to raise as expected. (c) The M33 guard (`write` requested together with the AO scheme) needs **two** attributes set at once, so it does not fit the single-attribute table — added as its own small test, `test_sanity_check_rejects_write_with_ao_scheme`. **[verified] it kills M33** (`(part, part_method) != ("atoms", "mo")` weakened to `part != "atoms"`), which previously survived the entire suite; it is the only test that fails on it. **Note on `match=`:** pytest uses `re.search`, so `"invalid partitioning"` also matches the longer `"invalid partitioning method..."` message — the two rows are still distinguishable because the new one matches the longer string, but any future row here should be checked against this prefix overlap. |
| T7 | `test_sanity_check_gauge_origin_invalid` (4) | spec | K | |
| T8 | `test_sanity_check_rejects_invalid_mo_coeff` / `_mo_occ` (3+3) | spec | K | |
| T9 | `test_sanity_check_accepts_valid_input` (now 7 rows) | spec | **S — partly done (2026-09-29)** | **Added** three rows for the PR #28 partitioning combinations, none of which were previously confirmed as *accepted*: `{"part": "orbitals"}`, `{"part": "atoms", "part_method": "ao"}` and `{"part": "eda"}` (deprecated, rewritten to atoms/ao). This table takes dicts, so multi-key rows fit — unlike T6's single-attribute table. **Declined by user:** a separate `test_sanity_check_resets_part_method_for_orbitals` for **M23**. Worth recording why it can't just be another row here: this test passes merely by *not raising*, so it cannot detect that `part_method` was silently kept instead of reset — that needs a real `assert decomp.part_method is None`. **[verified]** `DecompCls(part="orbitals", part_method="mo")` does go `'mo' → None` through `sanity_check`. **M23 therefore still has no coverage anywhere.** |
| T10 | `test_dim` | hand | K | M34 caught. |
| T11 | `test_mf_info_h2o` / `_li` / `_li_rohf` | hand count | K | Zero internal callers, but 13 tests use `mf_info` as setup, so a break would show up widely anyway. |
| T12 | `test_make_rdm1`, `_equal_electron_count` | hand / invariant Tr(DS) = N | K | M01 caught by both. |
| T13 | `test_e_nuc_h2` | hand | **S/R — proposed, DECLINED by user (2026-09-29)** | **No change made; `test_e_nuc_h2` stays as it is.** The weakness is real and stands recorded: with two atoms there is one pair, so symmetry *forces* each atom to get half of it, and a "just divide the total evenly" bug also produces `[0.5, 0.5]` — the test cannot fail on the split. (Same reason last session's LiH idea would not have helped: still a diatomic.) **The replacement was fully [verified] before being declined:** linear H₃ at 0, 1, 3 bohr, hand-derived from `e_nuc[i] = ½ Σ_{j≠i} Z_i Z_j / r_ij` → `[0.66667, 0.75, 0.41667]`, matching `_e_nuc` exactly and summing to `mol.energy_nuc()`; the equal-split bug gives `0.61111` three times. The spacing 0,1,3 rather than 0,1,2 was deliberate — **[verified]** equal spacing leaves atoms 0 and 2 mirror-equivalent (`[0.75, 1.0, 0.75]`), while 1/2/3 bohr distances make all three atoms inequivalent. **Consequence to be aware of:** together with T14 this means the per-atom *split* of the nuclear repulsion has no test that can actually fail on it — only its total does (T14). M03 remains caught only by the weak H₂ case. |
| T14 | `test_e_nuc` | PySCF, same formula | **K — reviewed, no change (2026-09-29)** | Keep as is: one line, external reference, and it does catch formula errors (M02, the `0.5 → 1.0` mutation). **[verified] against PySCF's source** that it is a *consistency* check rather than a true differential — `mol.energy_nuc()` uses the same `inter_distance` helper, the same `rr[diag] = 1e200` masking and the same einsum with the same `* .5`; the only difference is `->` versus decodense's `->i`, i.e. where the sum happens. **What it cannot catch:** any redistribution between atoms that preserves the total — M03 passes it. It guards the **total magnitude**; the **split** is guarded only by T13. |
| T15 | `test_dip_nuc`, `_h2o` | hand / gauge-shift invariant | K | M04 caught by both. |
| T16 | `test_h_core` | PySCF `get_hcore` | **K — reviewed, improvement declined (2026-09-29)** | Keep. The load-bearing assert is `hcore == kin + nuc`; it caught M30. **Two honest qualifications found on re-review.** (a) It is only *half* a differential: **[verified]** PySCF's `get_hcore` is `intor_symmetric('int1e_kin') + intor_symmetric('int1e_nuc')`, and decodense's `kin` is the *identical* call — so for the kinetic part this checks only that decodense didn't tamper with it, not that it is right. The nuclear half *is* independent (`int1e_nuc` vs per-atom shifted `int1e_rinv`). (b) `diag(kin) > 0` and `diag(nuc) < 0` are weak sanity checks that almost any non-broken implementation passes. **Declined:** tightening the tolerance. It uses `np.allclose` **defaults** (`rtol=1e-5`) while `test_e_nuc`, `test_get_nuc_matches_pyscf` and `test_make_rho` all pass `atol=1e-10` explicitly — **[verified]** the two sides actually agree to `5.3e-15`, but the default tolerates ~`6e-4` on the largest element, so a 0.01 % error would pass unnoticed. Untested branches: the ECP guard and the whole PBC path. |
| T17 | `test_get_nuc_matches_pyscf` | PySCF `int1e_nuc` (different route) | **K — reviewed, improvement declined (2026-09-29)** | Keep. Genuine differential (per-atom shifted `int1e_rinv` summed, vs PySCF's `int1e_nuc`), explicit `atol=1e-10`, and LiH deliberately over H₂ so a symmetry-hidden bug can't pass. It caught M05. **Hole found and [verified], then declined:** the test only checks that the slices **sum** correctly, and addition is order-insensitive — swapping the two atoms' slices (`sub_nuc[::-1]`) passes **both** assertions. That matters because `prop_atom_mo`/`prop_atom_ao` consume the slices individually as `sub_nuc[atom_idx]`, so a swap would make every per-atom energy wrong while leaving every total correct — invisible to the legacy tests too, since they assert only `sum(res.tot) == mf.e_tot`. **Proposed one-liner, [verified] to catch it:** `abs(np.trace(sub_nuc[0])) > abs(np.trace(sub_nuc[1]))`, from nuclear charge alone (Z_Li=3 vs Z_H=1) — correct output gives 15.22 vs 3.78, swapped gives 3.78 vs 15.22. **Declined by user.** |
| T18 | `test_point_charges_nuc_solv_shape` → **`test_point_charges`** | hand-computed Coulomb + PySCF `int1e_rinv` | **R + S — done (2026-09-29)** | **Renamed** (the old name claimed only a shape check, while the test also pins a hand-computed Coulomb value and runs a PySCF differential). **Strengthened from 1 QM atom + 1 point charge to 2 + 2**, because with one of each, both the per-atom loop and the blocked accumulation (`lib.prange(..., BLKSIZE)`) run exactly one iteration and cannot be wrong. New system, all on the z axis so every distance is a whole number: QM Li(Z=3) at 0 and H(Z=1) at 1; point charges O(q=8) at 3 and H(q=1) at 5 → `nuc_solv = [3(8/3 + 1/5), 1(8/2 + 1/4)] = [8.6, 4.25]`, **[verified]** exact. O+H was chosen over O+O so the two charges stay *distinguishable* — **[verified]** swapping them over their positions gives `[5.8, 2.5]`, whereas two identical O charges would be interchangeable by definition. The electronic reference now sums `int1e_rinv` over both charges. **[verified] it kills a new bug class:** forcing the per-atom loop to always read `atom_charges[0]`/`atom_coords[0]` now fails, and was undetectable with the old single-atom system. |
| T19 | `test_e_xc_1d` | hand (10.5) | K | A2 is done. Style: exact `==` on floats works here only because every number is an exact binary fraction. `np.isclose` is the safer habit. |
| T20 | `test_e_xc_2d` | hand | K | |
| T21 | `test_xc_ao_deriv` (4), `_unknown_type` | PySCF classification / mock | **S — done (2026-09-29)** | **A6 closed.** `_unknown_type` asserts `UnboundLocalError`, which is the *current buggy* failure mode: `_xc_ao_deriv` has no `else` branch, so an unrecognised `xc_type` leaves `ao_deriv` unassigned and the `return` line raises. A comment now records that this is deliberate, and that if someone adds a proper `else: raise ValueError(...)` the test should be **updated to expect `ValueError`, not reverted**. It also records why the mock is unavoidable — PySCF's `xc_type` only ever returns HF/LDA/GGA/MGGA/UNKNOWN, and UNKNOWN is unreachable with any real functional (see A10). |
| T22 | `test_make_rho`, `_gga`, `_atom_slicing` | PySCF `eval_rho` / invariant | **K — reviewed, no change (2026-09-29)** | The strongest group in the file; confirmed on re-review. The first two are true differentials, not same-formula consistency checks like T14: decodense builds `rho` with generic `einsum` contractions while PySCF's `eval_rho` uses its own compiled internals (`_dot_ao_dm`, `_contract_rho`). M08 is caught here. `_atom_slicing` is a pure invariant (Σ_atoms rho_atom == rho_total at every grid point) with no expected values at all, and it exercises the decodense-specific reason `_make_rho` is split in two — feeding sliced arrays `c0[:, select]` / `ao_value[..., select]` into `_make_rho_interm2`, which PySCF has no equivalent of. **Named blind spot:** being an invariant, the slicing test cannot catch an error that scales both sides equally (a uniformly doubled density leaves `rho_sum == rho_total` true); magnitude is covered by the differential tests instead, so the two are complementary by design. **Gaps, all previously declined:** MGGA (B12), the unrestricted `rdm1.ndim == 3` branch (B14), and a GGA version of the slicing invariant (B13) — `_make_rho` runs only 4 of 15 statements. |
| T23 | `test_vk_dft_pbe0` | functional definition (25 %) | K | M11 caught. |
| T24 | `test_trace_*` (3) | hand / `np.trace` | K | See VIII.1 #9. |
| T25 | `test_population_mul_h2o` | invariant + symmetry | K | The setup copies the `pop` formula from `assign_rdm1s`. That's fine: the function under test is only the AO → atom summation. |
| T26 | `test_assign_rdm1s_h2o`, `_iao` | invariant (MO normalisation) | K | `weights[1] == weights[0]` is guaranteed by the RHF shortcut's `.copy()`, so its only value is guarding that the shortcut exists. |
| T27 | `test_assign_rdm1s_li`, `_li_iao` | invariant | **R** | Single atom → partition of unity is vacuous (VIII.1 #7). Use an open-shell molecule with ≥ 2 atoms (e.g. OH radical, doublet, STO-3G, UHF). Then α and β weights genuinely differ, and "rows sum to 1" tests something. **[verified]** OH at 1.8 bohr, UHF: α 5×2 and β 4×2 weights, rows sum to 1, per-atom populations α [4.62, 0.38] vs β [3.55, 0.45] (Mulliken). |
| T28 | `test_atoms` | hand-built input | **S** | Only `"au"`, so unit scaling is unchecked (M20, M21 survive). Add a `kcal_mol` or `ev` case (values × 627.5094740631 or × 27.211386245988; constants from CODATA via the Wikipedia Hartree page, i.e. external). Optionally a 2-D (dipole) input to cover the `" (x)"` columns. |
| T29 | `test_to_dataframe` | hand count + smoke | **S — partly done (2026-09-29)** | **Added** `assert "O0" in str(res1)` and a line-count check on `str(res2)`, closing the gap that `ResultsCls.__str__` had **zero** coverage (**[verified]** by tracer) even though `print(res)` is how all 12 examples display results. **Honest scoping of what that buys:** it does *not* add coverage for PR #28's `charge_atom` crash — **[verified]** by restoring that exact bug, which `test_atoms` and `test_to_dataframe` both already caught, since `to_dataframe()` routes through `atoms()` anyway. What it *does* catch is `__str__` breaking on its own — **[verified]** by rewriting it to return `str(self.res_dict)`, which leaves `to_dataframe()` intact and is caught by **nothing else in the suite**. **Still open:** the row counts remain the only value check on `to_dataframe`; comparing a column against `res.tot` was proposed and not done. |
| T30 | `test_orbs_ndo_odd_count` | documented pairing rule | **S** | Checks only the length and the last element. A wrong pairing that still puts the unpaired orbital last would pass. Assert the whole order, hand-derived from the rule "most negative with most positive; unpaired last" and **[verified]**: `el == [1, 5, 2, 4, 3]`, `mo_occ == [-0.8, 0.8, -0.3, 0.3, 0.0]`. Use *distinct* symmetry labels (e.g. `"a".."e"` → expected `['a','e','b','d','c']`) so M38 is caught. |

## VIII.5 Missing tests, ranked

Each one closes a hole that the mutation run *proved* exists. Each has an external reference.

| # | Test | Catches | Reference | Cost |
|---|---|---|---|---|
| ~~N1~~ | ~~`DecompCls`: `part="eda"` → `("atoms","ao")`; `"atoms"` → `"mo"`; `"orbitals"` → `None`~~ (= VII.2 #5) — **done**, `test_decomp_cls_part_method_defaults` | M22 (survives everything, including legacy) | the documented deprecation promise | 3 lines, no SCF |
| **N2** | `main` input normalisation: raw 1-D RHF `mo_occ`, ROHF `mo_occ` (2/1/0), 3-D `np.asarray` UHF `mo_coeff` → same result as the `mf_info` route | M14, M15, M35 | route equivalence (two separate implementations of one convention) | small; fixtures exist |
| **N3** | Symmetry-equivalent atoms: water's H₁ and H₂ get identical values for **every** component, parametrized over `part_method` in `["mo","ao"]` | M25, M26 [verified]; the first per-atom check of the AO scheme | molecular symmetry (physics) | ~6 lines |
| ~~N4~~ | ~~Unit conversion: `orbs` with `unit="ev"` scales by the CODATA constant~~ (= VII.2 4b) — **done**, folded into `test_orbs_ndo`. `atoms`'s unit scaling (M20) and the `kcal_mol` constant (part of M21) are **still open** — see T28 in VIII.4. | M19, M20, M21 | CODATA constants | small, no SCF |
| **N5** | `write_rdm1` partition: Σ_atoms D_A = D (numpy format into `tmp_path`) | VIII.1 #4 (**fails today** for fractional occupations) | invariant | small. **Decision:** `xfail` + fix PR, or fix first |
| **N6** | `_e_nuc` linear H₃ (0, 1, 3 bohr) → [0.6667, 0.75, 0.4167] | wrong-split bugs that H₂ can't see | hand-derived | 3 lines |
| **N7** | `sanity_check` new rows (T6, T9) | M23, M33 | spec | rows only |
| N8 | Losslessness in the unit suite: `sum(res.tot) == mf.e_tot` for one HF + one DFT case | M27–M29, M36 if the legacy tests *don't* end up in CI | PySCF total | only if needed (VIII.7) |
| N9 | Dipole via `main`: for a neutral molecule the total is independent of the gauge origin, and equals `mf.dip_moment(unit="au")` | the unit suite has no dipole run at all | physics + PySCF | small; legacy `b3lyp_dipmom` covers the total already |
| N10 | `_solvent` dispatch with dummy objects, including the new non-PCM `NotImplementedError` (VII.4) | M32 | spec | small |
| N11 | `res_add`/`res_sub` on hand-built dicts | M37 | hand | small; zero callers → low value |
| N12 | `_unique_filename` + `assign_rdm1s` verbose file (in `tmp_path`) | — | spec | wait for PR #28 review question 3 |
| N13 | Regression ("golden") test from the local `tests/*.npz` references | any change in any per-atom number | **not external**: catches *changes*, not *bugs* | the files are untracked; decide whether they belong in the repo |

## VIII.6 What a good unit test looks like — the checklist used above

This is how the verdicts in VIII.4 were reached. It's worth using for every new test.

1. **The expected value comes from outside the code under test:** hand-derived physics, a
   genuinely separate implementation, or an invariant (symmetry, conservation, route equivalence).
2. **It has been seen to fail.** Plant the bug it's supposed to catch and watch it go red
   ("mutation testing"). If you can't name a plausible bug it catches, it's not a test yet.
3. **The input breaks symmetries.** Use unequal numbers (not all ones), ≥ 3 atoms at unequal
   distances, an open-shell system with ≥ 2 atoms, distinct labels. Otherwise a wrong
   implementation can give the right answer by coincidence (VIII.1 #7, M38).
4. **Assert the whole result**, not one element (T30, T1).
5. **One behaviour per test.** The name says which behaviour; a one-line comment says where the
   expected value comes from.
6. **Floats:** `np.isclose`/`np.allclose`, with a tolerance you can justify. Exact `==` only for
   integers, bools and exact binary fractions.
7. **Fast and clean.** The whole suite takes 2.5 s [verified: `--durations`], so per-test fixtures
   are fine and there's no need for `scope="module"`. Never write files into the repo; use `tmp_path`.
8. **A test that pins a known bug says so** (`xfail` with a reason, or a comment), so that fixing
   the bug doesn't look like breaking the test.

## VIII.7 Decisions for the next session

Grouped. Suggested order: **A before B**. A decides whether any of B runs in CI.

**A. Infrastructure (not test-writing)**
- [ ] A-1 Look at the fork's live CI run for PR #21: does it confirm VIII.1 #1 (0 unit tests run,
      legacy import errors)?
- [ ] A-2 How CI should run tests. Options: (a) two steps, both with
      `working-directory: tests`: `python -m pytest test_unit.py` and
      `python -m unittest discover -s . -p "test_*.py"`; (b) plain `pytest` from `tests/`, after
      removing the `mol.stdout.close()` teardown from the 7 legacy files; (c) convert the legacy
      files to pytest (bigger). Recommendation: **(a)**. It's the smallest change and touches no test
      file. Needs `pytest` installed in the CI job. Belongs in PR #21 or a follow-up, not the test PR.
- [ ] A-3 `Dict` → `dict` in `results.py:135` (VIII.1 #2): own tiny PR against `main`.
- [ ] A-4 PR #28 review comment #5: `strin` typo in `info()` (VIII.1 #3).
- [ ] A-5 `write_rdm1` occupation bug (VIII.1 #4): issue, fix PR, or ask Janus first? (Is `write`
      with NDOs a real use case?)

**B. New tests (VIII.5)**. For each: yes / no / later.
- [ ] N1 `DecompCls` mapping · [ ] N2 `main` input normalisation · [ ] N3 symmetry-equivalent atoms
- [ ] N4 unit scaling · [ ] N5 `write_rdm1` partition (+ xfail?) · [ ] N6 linear H₃
- [ ] N7 sanity_check rows · [ ] N8 losslessness (depends on A-2) · [ ] N9 dipole
- [ ] N10–N13 optional

**C. Changes to existing tests (VIII.4)**. For each: yes / no.
- [x] **T1 permute: done** (all keys + parametrized over `mo`/`ao` + `mo_occ` passed) · [x] **T2: done** (occupied columns, parens fixed, raw-pyscf route test added — kills M14; 3-D branch M35 + fractional occ M15 still open)
- [x] **T3 non-aufbau: done** (sliced-route value comparison) · [x] **T4: done** (renamed `test_rdm1_charge_conservation`, moved to the orbitals.py section)
- [x] **T5: reviewed — no change** (extra glob/loc assert written, verified, then declined by user) · [x] **T6: done** (dup row removed, 2 rows added, M33 test added)
- [x] **T13: declined by user** (H₃ verified but not adopted; split now has no failing-capable test) · [x] **T18: done** (renamed + 2 QM atoms / 2 point charges; kills a per-atom loop bug) · [x] **T21: done** (bug-pinning comment added) · [ ] T27 Li → OH radical
- [ ] T28 atoms unit case · [x] **T29: partly done** (str(res) covered; column-value check still open) · [ ] T30 full order + distinct labels

**D. Questions for Janus (not ours to decide)**
- [ ] Are `mo_basis`, `mo_init`, `loc_exp` meant to do anything? (validated, never used)
- [ ] Should the local `.npz` reference files become a tracked regression test (N13)?
- [ ] Is it intended that `res.struct` is a per-*atom* array in `part="orbitals"` mode, where every
      other key is per-orbital (VIII.1 #10)? It's excluded from the DataFrame, so nothing is
      visibly wrong — but it is still reachable as an attribute.

## VIII.8 Stale entries elsewhere in this document

So they aren't acted on twice:
- Part I **A2** is done (`test_e_xc_1d` uses weights `[2, 1, 0.5]`). **A7** is obsolete (the
  circular line is no longer in `test_h_core`). **A13** is done (only one "should not raise" test
  remains, and it has the comment). **A6** is still open.
- Part I **B1** is done (`test_orbs_ndo_odd_count`), **B2** done (Li tests, but see T27), **B5** done
  (`test_dip_nuc_h2o`), **B6** done (`test_atoms`), **B3** partly done (`test_to_dataframe` runs
  `ResultsCls`). **B14**'s recorded reason is wrong (VIII.1 #8).
- Part II's `results.py` table (❌ everywhere, `fmt` "branch on `charge_atom`") and its "38 tests"
  counts predate PR #28. VIII.2 supersedes them.
- Last session's "LiH `_e_nuc` test" idea doesn't work (VIII.1 #7). Use linear H₃ instead.
