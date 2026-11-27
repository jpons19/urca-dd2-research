# Roadmap

Updated 27 November 2026. Two tracks: the EOS-independent theory reference and the multi-EOS numerics. Open scientific questions are tracked in [THEORY_OPEN_QUESTIONS.md](THEORY_OPEN_QUESTIONS.md).

## Theory reference ([report/urca_physics_comparison.html](../report/urca_physics_comparison.html))

Done: grouped matrix (blocks A–E) for FM79, Bottaro, Shternin, Alford, Sedrakian and Voskresensky; source checks closed for FM79 (β<sub>n</sub>, propagator, ρ exchange, branches), Voskresensky 2001 (Eqs. 24–29), Voskresenskiĭ & Senatorov 1986, and Yakovlev et al. 2001 (α<sub>n</sub>, β<sub>n</sub> conventions); FM79 evidence on β<sub>n</sub>/ρ overlap recorded.

Next:
1. **Subproject (block F): chiral-EFT currents and in-medium pions.** Survey charged-current weak rates with consistent chiral two-body currents, normal-ordered currents in neutron-rich matter, and cold-matter pion polarization / g′; audit sources with the same protocol as the main columns; decide whether the matched A/B/C comparison of Q1 can be built.
2. Audit the medium-modified Urca cooling parametrizations (Blaschke–Grigorian–Voskresensky 2004 and later) in the theory reference, as model inputs rather than results.
3. Log the Sedrakian ε<sub>min</sub> sign question (printed closed form below Eq. 14) as a source-verification item.

## Numerical track ([NUMERICS.md](NUMERICS.md))

Done: standardized cold npeμ states for PT(GRDF2-DD2), BCP(NLRMF-N) and BSk22; FM, Bottaro (with β<sub>n</sub> band), Shternin, adapted NWA, Sedrakian, Voskresensky MMU (first approximation) and sharp DU; TOV reference densities; 22 tests.

Next:
1. Proton-branch MU (Yakovlev et al. 2001 rescaling) for FM/Bottaro/Shternin.
2. MMU refinements: separate π⁰/π± gaps (μ<sub>π</sub> = μ<sub>e</sub>), condensate branch, EOS-consistent pion polarization (links to the subproject).
3. Superfluid reduction factors; finite-temperature (non-degenerate) DU kernel near threshold; the full Alford finite-T response.
4. A third RMF with intermediate L or a microscopic (BHF/APR-type) EOS for a Shternin-consistent comparison.

## Definition of done for a numerical change

Document inputs and conventions, give the source equations or derivation, check a known limit, add tests (regression, limit, convergence), and state whether existing outputs change.
