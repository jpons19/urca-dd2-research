# Pion softening and medium-modified Urca: research note

**Status:** theoretical consolidation, 26 November 2026. **No new nuclear-medium emissivity calculation is validated in this note.** Use this document alongside [METHODS.md](METHODS.md), [REFERENCES.md](REFERENCES.md), and the unchanged DD2 reference data.

## 1. Scientific question and controlled baseline

How would a new interaction that changes the in-medium pion response alter neutron-star modified Urca (MU) emissivity? We must distinguish (a) a phenomenological change of a *constant propagator screening mass* from (b) a collective soft pion mode governed by a momentum-, energy-, density- and isospin-dependent polarization. The latter is not equivalent to a mass substitution in free OPE and may involve weak radiation from intermediate states.

The existing PT(GRDF2-DD2) setup supplies cold beta-equilibrated fractions, nucleon scalar/vector self-energies and Fermi momenta over the homogeneous core, approximately 0.16–0.98 fm^-3. The four approaches share this **kinematic background**, but **not** a consistently matched nuclear interaction or common set of diagrams. The legacy baseline outputs and coefficients must not be overwritten by this study.

For fixed composition, masses, weak vertices and coupling strength, a static pion exchange of momentum q contains the factor

\[
V_\pi(q) \propto \frac{(\boldsymbol{\sigma}_1\cdot\mathbf q)(\boldsymbol{\sigma}_2\cdot\mathbf q)}
{q^2+m_{\pi,\mathrm{scr}}^2},
\]

so varying only the denominator gives, at fixed q,

\[
\frac{V_\pi^*(q)}{V_\pi(q)} =
\frac{q^2+m_\pi^2}{q^2+m_{\pi,\mathrm{scr}}^2}.
\]

The corresponding squared ratio is a **pure-OPE, fixed-q toy estimate**, not a universal ratio of full MU emissivities. At neutron-star Fermi momenta q often exceeds m_pi, so a constant mass reduction need not cause a dramatic increase. Interference with contact, rho and exchange amplitudes can change even the sign of a partial correction relative to a simple OPE scaling.

**Coupling convention:** changing the propagator screening mass alone does **not** authorize changing the vacuum mass in the vertex normalization (f_piNN/m_pi)^2. In code keep the vertex normalization fixed unless the underlying new interaction independently specifies an in-medium vertex.

## 2. Collective pion mode and pion gap

A general retarded propagator has

\[
D_{\pi,a}^{R}(\omega,q;n_B,Y_p) =
\big[\omega^2-q^2-m_\pi^2-\Pi_{\pi,a}^{R}(\omega,q;n_B,Y_p)\big]^{-1},
\]

where a labels the neutral or charged channel; conventions for chemical-potential offsets must be fixed before comparison. The static in-medium inverse propagator and effective pion gap at its finite-q minimum are different objects from a constant m_pi,scr. In neutron-rich matter charged and neutral modes generally differ. Complex self-energies introduce damping; near a soft-mode minimum the imaginary part and consistent vertex dressing matter.

**Source check:** Voskresensky (2001), Secs. 2.1–2.4, explicitly separates local Landau–Migdal terms from the long-range dressed pion contribution, introduces neutral/charged vertex factors in Eqs. (4)–(8), the pion self-energy in Eqs. (9)–(10), and the effective pion gap near q ~ p_Fn in Sec. 2.2. In Eq. (11) a dressed propagator enhancement is *counteracted* by pion–nucleon vertex renormalization: the normal-matter NN scattering rate cannot be inferred from |D_pi|^2 alone. See the [original arXiv review](https://arxiv.org/abs/astro-ph/0101514).

The same review, Sec. 3.1, Fig. 1, Eq. (24) and Eq. (25), includes contributions where the weak current couples to intermediate nuclear/pion states. Those **additional topologies** are an essential distinction between medium-modified Urca (MMU) and a simple OPE substitution. Its illustrative enhancement scalings and gap powers depend on the complete MMU model, dressed vertices, pion-charge channels and correlation factors; they must not be applied as a universal correction to Bottaro, FM, Shternin or Alford.

Pion condensation changes the equilibrium phase and excitation spectrum; normal-phase formulas should not be continued through an instability. Pairing gaps, the EOS, and DU triangle gaps could all shift if the new interaction is strong enough. A fixed-DD2 calculation intentionally excludes these feedback effects.

## 3. Impact on the four current prescriptions

| Approach | How softened-pion physics enters | What cannot be carried over uncritically |
|---|---|---|
| **Friman–Maxwell (1979), historical NSCool normalization** | Recompute strong NN matrix elements and phase-space average; a changed pion exchange generally changes the historical alpha_n and SRC treatment. | The alpha_n formula and beta_n=0.68 are legacy phenomenology, not independent universal factors. Replacing m_pi inside f_piNN/m_pi without rederiving the coupling is incorrect. |
| **Bottaro, Caputo & Fiorillo (2024)** | Replace the **pion part** of q-dependent tensor/spin couplings, retaining explicit contact and rho terms and all direct, exchange, interference amplitudes. Both triangular and nontriangular evaluations are accessible in this repo. | Freeze rho/Landau couplings only for a controlled sensitivity scan; they are not guaranteed unaffected by new physics. A correction at q=p_Fn need not represent the full angular integral. The retained beta=0.68 and rho/contact terms must be audited for correlation overlap. |
| **Shternin, Baldo & Haensel (2018)** | Recompute interaction-dependent scattering amplitudes and, if self-consistent, quasiparticle spectra; separately track near-on-shell intermediate nucleon denominators. | R_A/R_B describe kinematic near-shell amplification; coefficients 0.6/0.2 and the underlying Brueckner interaction are not automatically transferable to a softened-pion potential. Multiplying an MMU factor by the published approximation is not controlled. |
| **Alford, Haber & Zhang (2024), adapted NWA in this repo** | A common NN interaction can modify nucleon collision widths Gamma_n,p and spectral functions; the full charged-current response may also acquire collective/vertex contributions. | Current widths W=T^2/(5 MeV) are *prescribed*, independent of pion physics; presently changing a pion propagator does nothing to this curve. The repo uses a degenerate DU kernel and is **not** a full Alford calculation. A width-only modification is not a conserving dressed response. |

In all four approaches, **nucleon Landau masses** control degenerate phase space, while the RMF **Dirac masses** enter the current adapted spectral convolution. Any strongly modified EOS must update these quantities consistently, as described in METHODS.md.

## 4. A diagrammatic bookkeeping rule

Separate the following objects rather than interpreting them as interchangeable 'medium corrections':

1. Irreducible or resummed strong **NN amplitude**: Born OPE, rho exchange, Landau–Migdal contact terms, short-range correlations and pion polarization.
2. **Weak vertices/current insertions**: nucleon leg radiation, intermediate-state and pion-mode radiation, and vertex renormalization.
3. **Intermediate nucleon propagators**: conventional constant denominators versus Shternin's momentum-dependent near-shell denominators.
4. **Quasiparticle collision widths** and spectral functions: the quantity used in the adapted NWA and possible coherence/multiple-scattering corrections.
5. **Thermodynamics and available phase space**: composition, Dirac/Landau effective masses, charged/neutral channels, DU thresholds, and pairing.

Phenomenological multiplicative factors extracted from *different* NN interactions or diagram sets may share physics and must not be multiplied without an explicit overlap analysis. A consistent microscopic calculation would establish a common interaction, dressing prescription and weak vertex scheme.

## 5. Proposed numerical workstream: not yet implemented

**Stage A — controlled fixed-background test:** add a standalone pion propagator/coupling interface. Reproduce the vacuum-mass baseline exactly for both Bottaro triangular and nontriangular neutron-branch rates. Vary only the denominator m_pi,scr(n) at fixed DD2 composition, f_piNN vertex normalization, rho/contact terms and beta. Record separately the electron/muon rates and ratios to the same unmodified baseline. Converge angular quadrature (Sobol sample count/seed) and add a pure OPE fixed-q analytic check.

**Stage B — collective response:** introduce neutral/charged D_pi(omega,q) with explicit real/imaginary polarization, momentum of minimum, damping and vertex factors; restrict to the stable normal phase. Compare at the *amplitude level*, not simply by rate multipliers. Include a physically documented parametrization rather than treating an effective mass as a complete self-energy.

**Stage C — matched many-body methods:** recompute interaction-dependent NN scattering and the Shternin coefficients, predict Gamma_n,p from the same model, and compare with finite-width charged-current response; ultimately account for weak-current diagrams, possible pion condensation, pairing, and EOS feedback.

### Acceptance tests for Stage A

- Exact recovery of the original pion mass at every density, branch and temperature; reference CSV/PNG unchanged unless deliberately regenerated.
- Correct zero-pion contribution / isolated pure-OPE limit and no spurious m_pi,scr^-4 vertex scaling.
- Continuity and finite positive emissivity throughout scanned stable mass range; clarify that the **full** rate is positive even if interference terms are negative.
- Identical EOS fractions, Landau/Dirac masses and chemical potentials when studying a *fixed* DD2 background.
- Convergence of the nontriangular angular integral and explicit lepton-channel sums.
- Documentation of what beta and Landau coefficients have been held fixed and why the result remains a toy sensitivity study.

## 6. Available sources and limits

- **Friman & Maxwell (1979)**, *Neutrino Emissivities of Neutron Stars*, Astrophys. J. **232**, 541–557, [ADS scan](https://adsabs.harvard.edu/pdf/1979ApJ...232..541F).
- **Bottaro, Caputo & Fiorillo (2024)**, [arXiv:2406.18640](https://arxiv.org/abs/2406.18640). Repo includes earlier independent Eq. (3.31) spin/isospin trace diagnostics; see PRIOR_AUDIT.md.
- **Shternin, Baldo & Haensel (2018)**, [arXiv:1807.06569](https://arxiv.org/abs/1807.06569).
- **Alford, Haber & Zhang (2024)**, [arXiv:2406.13717](https://arxiv.org/abs/2406.13717); repo's treatment is adapted NWA only.
- **Voskresensky (2001)**, *Neutrino Cooling of Neutron Stars. Medium Effects*, [arXiv:astro-ph/0101514](https://arxiv.org/abs/astro-ph/0101514), **38-page review supplied by the user**. Primary reading: Secs. 2.1–2.4; Sec. 3.1 (Fig. 1, Eqs. 24–27); later non-equilibrium/width discussion.
- CompOSE PT(GRDF2-DD2), [cold EOS 302](https://compose.obspm.fr/eos/302) and [general EOS 301](https://compose.obspm.fr/eos/301).

The original article PDFs have not been committed to this public repository; link to publishers/arXiv/ADS instead. The **numerical pion-screening sensitivity scan discussed in chat has not been independently integrated or validated in GitHub** and is deliberately excluded from the canonical baseline and this documentation-only commit.

## Implementation status (26 November 2026)

The first controlled fixed-DD2 screening-mass experiment is now implemented separately in `calculations/pion_softening/`. It uses both Bottaro angular treatments and does not change any original four-approach baseline files. The source, reproduction command, assumptions and initial validation results are documented in [PION_SCREENING_SCAN.md](PION_SCREENING_SCAN.md). The collective-propagator, new weak-current diagrams, Shternin recalibration and nucleon-width stages remain **unimplemented**.
