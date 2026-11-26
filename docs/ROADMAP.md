# Research roadmap / suggested GitHub issues

## Primary editorial objective: paper-by-paper theoretical comparison

**Master artifact:** [`report/urca_physics_comparison.html`](../report/urca_physics_comparison.html), developed under [`AGENTS.md`](../AGENTS.md). It is self-contained and independent of EOS or calculations. Maintain a matrix of **physical ingredients** plus expandable detailed interpretations, with precise source pages/equations, uncertainties and explicit distinctions between papers. Review **FM79 and Bottaro first**; then Shternin, Alford and Voskresensky, one approach at a time. Avoid treating provisional cells as fully audited. Confirm the exact original FM79 Table 1/page attribution for beta=0.68 before marking that citation verified; the current HTML flags it as pending.

**Focused review completed (26 November 2026):** The master HTML now distinguishes the pion-exchange momentum approximation from the intermediate nucleon propagator approximation in FM79, Bottaro (Eq. 3.22) and Shternin (Eqs. 4–9). The FM79 exact printed-page/equation attribution remains an explicit source-verification task. Shternin detailed comparison now also covers Eqs. (6)–(14), analytic angular averages R_A/R_B, the meaning of the fitted 0.6/0.2 factors, interference and the nonfactorizable strong-interaction weights. Outstanding: quantitative independent reproduction of the complete half-on-shell BHF matrix element, finite-width consistency at the threshold, and the complete diagrammatic treatment of soft-pion emission.

**Keep separate:** numerical EOS-dependent analyses, screening scans, source-code details and validation belong in the computational track below. The HTML reference is not a report of our calculations.


## High priority — reproducibility and microscopic validity

1. **Implement full finite-T Alford charged-current response.** Replace the degenerate DU kernel inside the Breit–Wigner mass integral; validate DU limit, neutrino/antineutrino detailed balance and behavior across both DU thresholds.
2. **Reproduce Bottaro Fig. 6** with an independently obtained APR stellar profile and composition. The nontriangular neutron branch is implemented and checked against the triangular limit on the DD2 grids.
3. **Reproduce Shternin's approximation from explicit propagators** and assess sensitivity to effective masses, interaction and near-threshold behavior. Do not multiply an incompatible emissivity by an enhancement calibrated to another interaction without matching.
4. **EOS broader than current core grid:** obtain an RMF model that actually opens the electron DU threshold in tabulated density range, to test finite-width crossover. Current GRDF2-DD2 remains below threshold through 0.98 fm⁻³.

## Pion softening and medium-modified Urca — new workstream (26 November 2026)

Read [PION_SOFTENING.md](PION_SOFTENING.md) before introducing any new 'enhancement factor'. Keep the DD2 baseline scripts, outputs, and original approximate prescriptions unchanged initially.

1. **Documentation / physics bookkeeping (completed in the 26 November theory note):** identify OPE screening-mass substitutions versus the collective pion gap; distinguish neutral/charged modes, pion/weak vertex corrections, additional weak diagrams, and potential overlap with rho, Landau terms, SRC beta, Shternin propagators and NWA widths.
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

## Alford NWA theory consolidation (26 November 2026)

The EOS-independent master HTML now separates the strong-interaction input from the weak response: the earlier Paris-potential Brueckner calculation supplies the phenomenological mass width W, whereas Alford does not calculate an explicit NN scattering matrix in the charged-current kernel. The width's spectral convolution (Eq. 10) yields a mathematically smooth finite-width crossover and a direct-Urca narrow-width limit, **but not a proof of microscopic modified-Urca accuracy**. The first-level NWA omits crossed/interference/vertex terms and energy/momentum-dependent self-energies; the paper's Eq. (10) concerns reaction rates, while emissivities require energy weighting. The theory overview should not use EOS-specific numerical applications or output curves. Next: systematic medium weak-vertex and Voskresensky diagram comparison.

## Sedrakian theory comparison (26 November 2026)

Added Sedrakian (PRL 133, 171401; arXiv:2406.16183) as a sixth column across the 21 theoretical physical-ingredient rows, with a matching standalone expandable discussion in the master HTML. Distinguish the nonrelativistic proton energy-momentum spectral function and SRC-motivated phenomenological np/pp cross-section width from Alford's relativistic two-nucleon mass-spectral convolution and imported Paris-potential width. Both produce subthreshold spectral strength and recover the sharp DU limit when the width vanishes, but neither proves microscopic equality with a complete MU amplitude. Next: critically study Voskresensky's quasiparticle and full-width closed-diagram formulations and when each relates to these approaches. Quantitative calibration is explicitly deferred.

## Voskresensky diagrammatic comparison (26 November 2026)

Primary-source comparison added to master HTML: Fig. 1 (external-leg vs intermediate-state radiation), the three explicit mechanisms of Eq. (24) (charged-pion weak emission, neutron-hole particle–hole weak emission, conventional nucleon-leg MU), and the structure of MMU emissivity Eqs. (25)–(27) including the distinction between the vertex/correlation product Γ^8 and nucleon collision widths. The closed-diagram treatment Eqs. (42)–(49) explains how dressed one-loop self-energy cuts yield positive squared external-leg contributions and other diagrams yield leg/exchange interference. It is not a one-to-one mapping between any single open diagram and a spectral self-energy/vertex term. Avoid adding full MMU formulas to already dressed spectral rates without identifying disjoint cuts; no quantitative collective-pion calculation was performed. Next: microscopic consistent matching of self-energies, pion polarization and weak vertices, if numerical work is requested.
