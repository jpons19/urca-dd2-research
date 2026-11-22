# Research roadmap / suggested GitHub issues

## High priority — reproducibility and microscopic validity

1. **Implement full finite-T Alford charged-current response.** Replace the degenerate DU kernel inside the Breit–Wigner mass integral; validate DU limit, neutrino/antineutrino detailed balance and behavior across both DU thresholds.
2. **Implement Bottaro's nontriangular neutron branch**, with corrected explicit spin/isospin amplitudes; test positivity and angular Jacobian. Include an independently constrained composition and numerical convergence checks.
3. **Reproduce Shternin's approximation from explicit propagators** and assess sensitivity to effective masses, interaction and near-threshold behavior. Do not multiply an incompatible emissivity by an enhancement calibrated to another interaction without matching.
4. **EOS broader than current core grid:** obtain an RMF model that actually opens the electron DU threshold in tabulated density range, to test finite-width crossover. Current GRDF2-DD2 remains below threshold through 0.98 fm⁻³.

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

`physics-validation`, `EOS`, `modified-urca`, `finite-width`, `numerics`, `tests`, `documentation`.

## Definition of done for a physics change

Document input data and conventions, give a mathematical derivation or source, compare with a known limit, add unit/regression tests, report numerical convergence (including density and temperature), and state whether a change is intended to alter existing published plots.
