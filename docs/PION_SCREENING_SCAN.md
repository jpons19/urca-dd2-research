# Static pion screening sensitivity — Commit 2

**Status:** implemented and locally tested, 26 November 2026. **Important:** this is a fixed-background **toy sensitivity study**, NOT collective pion softening or Voskresensky's complete medium-modified Urca (MMU).

## Scope

The study uses the real PT(GRDF2-DD2) cold composition and the same nucleon scalar/vector fields and Landau masses as the historical Bottaro T9 calculations. The sole variation is the static OPE denominator

\[
u_\pi(q;\lambda)=q^2/[q^2+(\lambda m_{\pi,\mathrm{vac}})^2],\qquad \lambda=m_{\pi,\mathrm{scr}}/m_{\pi,\mathrm{vac}}.
\]

The pion–nucleon coupling normalization (f_piNN/m_pi,vac)^2 is **not modified**. The rho contribution, phenomenological Landau contact terms, weak vertex and inherited SRC factor beta=0.68 remain fixed. The nucleon and lepton composition and the EOS do not respond to screening.

Both Bottaro implementations are used:
- **Triangular:** update the pion part of H=12(f'-g)^2+(21/4)h'^2. The ratio to the original triangular amplitude multiplies *both* lepton rates.
- **Nontriangular:** evaluate all transferred momenta in the same direct, exchange and interference spin trace as the existing angular integrator, with the modified q-dependent pion term. The angular integral is normalized to its **own** modified triangular amplitude; include the change of the triangular amplitude separately. Keep the electron velocity correction of the original implementation.

The denominator change does **not** describe collective mode softening at a finite q, pion damping, charged/neutral splitting, RPA polarization, vertex renormalization, new weak-current insertions, medium-modified widths or EOS feedback. See [PION_SOFTENING.md](PION_SOFTENING.md) and [Voskresensky (2001)](https://arxiv.org/abs/astro-ph/0101514).

## Reproduce

From repository root, with the original requirements installed:

    python -m unittest discover -s tests -p 'test_pion_softening.py' -v
    python -m calculations.pion_softening.scan_dd2
    # Convergence study: increase --sobol-power from 14 to 15 or 16
    python -m calculations.pion_softening.scan_dd2 --sobol-power 16

The standard run uses **36 densities from 0.16 to 0.98 fm^-3**, screening ratios **1, 0.75, 0.5 and 0.25**, and a common deterministic scrambled 2^14-point Sobol quadrature with seed 42. The generated outputs are kept apart from the baseline in figures/pion_softening/:

- DD2_Bottaro_pion_screening.csv: individual electron/muon rates, total rates, enhancements for each method, plus DD2 diagnostics.
- DD2_Bottaro_pion_screening.png and .pdf: rate ratios versus density.

Run --help for alternative grid sizes and ratios. **No previous EOS inputs, Urca formulas, baseline CSV or figures are overwritten.**

## Representative results at T = 10^9 K

For a fixed screening ratio lambda=0.5, with densities at the nearest points of the 36-point scan:

| Density (fm^-3) | Triangular ratio | Nontriangular ratio |
|---:|---:|---:|
| 0.160 | 1.3529 | 1.3319 |
| 0.511 | 1.2207 | 1.1699 |
| 0.910 | 1.1862 | 1.1117 |

Ratios are computed **within each method** relative to the same method at lambda=1. The two calculations have different baseline absolute emissivities, so their enhancement ratios should not be mistaken for their relative emissivity amplitudes.

These order-unity enhancements do not reproduce the much stronger, momentum-selective collective MMU enhancement argued for in Voskresensky's formalism.

## Verification and limitations

Local tests cover (i) equality with the original g'_q,h'_q couplings at lambda=1, (ii) exact recovery of both triangular electron/muon rates, (iii) agreement of both nontriangular rates with the original angular integrator at matching Sobol points, (iv) fixed-q pure-OPE limit and positivity, and (v) a **5%** agreement check between 2^12 and 2^14-point Sobol quadrature for one representative density and lambda=0.5.

The 5% comparison is an initial convergence check rather than a quantitative error bound across the entire density grid. A full multi-density convergence study using multiple scramble seeds and powers 14–16 should precede publication or use in cooling. There are also unresolved microscopic questions concerning rho/SRC overlap and whether the unchanged Landau couplings remain appropriate for the hypothetical new interaction.
