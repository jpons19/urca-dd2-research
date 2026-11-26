# Research roadmap / suggested GitHub issues

## High priority — reproducibility and microscopic validity

1. **Implement full finite-T Alford charged-current response.** Replace the degenerate DU kernel inside the Breit–Wigner mass integral; validate DU limit, neutrino/antineutrino detailed balance and behavior across both DU thresholds.
2. **Reproduce Bottaro Fig. 6** with an independently obtained APR stellar profile and composition. The nontriangular neutron branch is implemented and checked against the triangular limit on the DD2 grids.
3. **Reproduce Shternin's approximation from explicit propagators** and assess sensitivity to effective masses, interaction and near-threshold behavior. Do not multiply an incompatible emissivity by an enhancement calibrated to another interaction without matching.
4. **EOS broader than current core grid:** obtain an RMF model that actually opens the electron DU threshold in tabulated density range, to test finite-width crossover. Current GRDF2-DD2 remains below threshold through 0.98 fm⁻³.

## Pion softening and medium-modified Urca — new workstream (26 November 2026)

Read [PION_SOFTENING.md](PION_SOFTENING.md) before introducing any new 'enhancement factor'. Keep the DD2 baseline scripts, outputs, and original approximate prescriptions unchanged initially.

1. **Documentation / physics bookkeeping (completed in this commit):** identify OPE screening-mass substitutions versus the collective pion gap; distinguish neutral/charged modes, pion/weak vertex corrections, additional weak diagrams, and potential overlap with rho, Landau terms, SRC beta, Shternin propagators and NWA widths.
2. **Controlled Bottaro scan (initial implementation added 26 November 2026):** reproduce the historical pion mass limit; vary a fixed-coupling screening denominator in triangular and nontriangular amplitudes, preserving the original rho, Landau and SRC terms. Initial branch recovery and a 2^12 vs 2^14 single-density Sobol convergence check pass; expand multi-density/multi-seed validation before publication. See [PION_SCREENING_SCAN.md](PION_SCREENING_SCAN.md).
3. **Collective propagator extension (not yet implemented):** parametrized D_pi(omega,q;n) and damping with documented neutral/charged modes, stability conditions and explicit vertex choices. Do not equate this with the full MMU formalism.
4. **Consistent MMU diagrams (research):** include Voskresensky's intermediate-state weak-emission mechanisms in addition to softened pion exchange. Establish diagrammatic overlap before combining any correction factors.
5. **Matched scattering/width comparison (research):** derive interaction-dependent Shternin ingredients and nucleon widths for the adapted Alford calculation, with conserving weak vertices as the longer-term objective.

## Medium priority — nuclear physics

5. Separate on-shell scattering amplitudes, Landau quasiparticle interactions and Born potentials; standardize antisymmetrization and spin/isospin conventions.
6. Recompute SRC correlations / `beta_n` rather than inheriting historical factors with an updated rho interaction.
7. Compute in-medium neutron/proton spectral widths `W_n(W,p,T,nB,Yp)` from a common NN interaction and thermodynamics.
8. Include conserving weak vertex corrections and check appropriate Ward identities.
9. Extend neutron-only MU results to proton branch, and introduce nucleon pairing reductions where needed.

## Infrastructure

10. Replace global mutable `T`, `T9`, `HALF_WIDTH` in `compare_compose_four.py` with immutable parameter objects, separating EOS interpolation from emissivity kernels.
11. Add quantitative unit tests using fixed benchmark densities; compare alternative interpolation orders and densities near muon onset.
12. Retain a compact, derived-data benchmark for CI; add slow tagged tests for BW quadrature and T² extrapolation. Ensure plotting stays optional.
13. Add a clear license, CITATION file, upstream table references and optional automated workflow after deciding repository ownership.

### Proposed issue labels

`physics-validation`, `pion-softening`, `MMU`, `EOS`, `modified-urca`, `finite-width`, `numerics`, `tests`, `documentation`.

## Definition of done for a physics change

Document input data and conventions, give a mathematical derivation or source, compare with a known limit, add unit/regression tests, report numerical convergence (including density and temperature), and state whether a change is intended to alter existing published plots.
