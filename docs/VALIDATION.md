# Validation state (2026-11-22)

All numbers below refer to `0.16–0.98 fm^-3` and the precomputed outputs in `calculations/`. Scripts for the checks are included. **These are numerical consistency checks, not proof of microscopic correctness.**

| Check | Observed / acceptance |
|---|---|
| Beta-equilibrated composition vs independent cold table | `max |Yp(T8)-Yp(0)| = 1.0883e-6` |
| Electron difference | `1.3870e-6` |
| Muon difference | `1.1156e-6` |
| Max Dirac mass change, T=0.1 MeV vs T8 | `0.000523 MeV` |
| Max Landau mass change, T=0.1 MeV vs T8 | `0.000466 MeV` |
| Beta-equilibrium residual, solved T8 | `< 9e-11 MeV` |
| Charge-neutrality residual, solved T8 | `< 3.2e-8` |
| Min electron DU momentum gap | `13.27 MeV` (positive) |
| Min muon DU momentum gap | `36.08 MeV` (positive) |
| T8/T9 vs exact `10^-8` scaling | Deviations `≤0.0047%` in this implementation |
| Electron+muon total | Sum checked for all four approaches |
| Positive / finite rates | Checked throughout the calculation range |

### Reproduction steps

```bash
python calculations/test_final.py
python calculations/test_T8_reference.py
python -m unittest discover -s tests -v
```

Full reruns:

```bash
python calculations/build_final.py
python calculations/run_dd2_T8_reference.py
```

### What has not been independently verified

- Independent reproduction of Bottaro Fig. 6 with the original APR stellar profile. The DD2 nontriangular curve uses the trace-verified interference expression, recovers the triangular limit and has been checked for numerical convergence over the DD2 grid.
- Reproduction of Shternin's original Brueckner scattering calculation with DD2.
- Alford's full finite-temperature relativistic polarization tensor, in-medium weak-vertex corrections or a conserving approximation.
- A near-threshold DU/MU crossover (the GRDF2-DD2 cold composition never crosses the DU triangle threshold through `nB=0.98 fm^-3`).
- Medium-broadened nucleon widths consistent with the same EOS and nuclear interaction.
- Equivalence of the old FM SRC factor `beta_n=0.68` with the retained rho-meson correction.

**Warnings for future contributors:** Two closely spaced T slices are used for T² extrapolation; precision of temperature differences is not precision of microscopic emissivity. Treat differences between the five prescriptions as model-dependent, not error bars. Before merging scientific revisions, add tests at the level of the matrix element, angular integration and continuum limits.

## Scope of the 26 November pion-softening documentation update

The tests and numerical evidence above validate **only the pre-existing DD2 and four-approach baseline**. The proposed screening-mass and collective-pion extensions are not yet merged or tested; see [PION_SOFTENING.md](PION_SOFTENING.md). The reference data, baseline plots and calculation scripts were deliberately left unchanged. This distinction prevents interpreting exploratory fixed-mass ratios or the MMU review's approximate scalings as validated emissivities from our code.
