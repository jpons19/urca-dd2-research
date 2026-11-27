# Open questions in Urca theory — EOS independent

The authoritative comparison is [report/urca_physics_comparison.html](../report/urca_physics_comparison.html). This file records unresolved theory questions and the specific evidence needed. **Do not infer numerical results from this register.**

## Q1 — Are Voskresensky's additional weak-current diagrams important without strong pion softening?

**Status: open; highest conceptual priority.**

**Source-derived facts.** Voskresensky (2001), [Fig. 1 and Eqs. (24)–(29)](https://arxiv.org/abs/astro-ph/0101514), distinguishes ordinary external-nucleon-leg radiation from intermediate-pion and neutron-hole/particle–hole weak-current processes. The large MMU contribution and claimed intermediate-state dominance arise within a specific pion-polarization, Landau–Migdal vertex and collective-pion model. The closed-diagram treatment in Sec. 4 provides a broader organizing framework but is not equivalent to a numerical evaluation with every realistic NN potential.

**Verified details (27 November 2026, arXiv text pp. 17–18).** In the first Eq. (24) diagram the lepton pair is emitted where a π⁺ line meets a π⁰ line, i.e. a virtual charged pion β-decays into a neutral pion (ππW, vector-type vertex). It is not π⁻→e⁻ν̄ annihilation through the pion-pole axial current, so the vacuum lepton-mass suppression of that process does not apply. Eq. (25) covers the first two diagrams and scales as T⁸ (ρ/ρ₀)^{10/3} ω̃_{π0}^{-4} ω̃_{π±}^{-4} Γ⁸ F₁ ζ(Δn)ζ(Δp). In F₁ (Eq. 27) the second term is the pion-decay diagram and the "1" the particle–hole diagram; their ratio, 3(ρ/ρ₀)^{2/3}/[4γ̃²(g′)γ_β²], contains no pion-gap factor.

**Original 1986 calculation (Voskresenskiĭ & Senatorov, Sov. Phys. JETP 63, 885).** The matrix element (Eqs. 4.3–4.5, p. 891) contains a vector-current term proportional to 1−β+κ, with β the ratio of neutral to charged pion gaps and κ the particle–hole (intermediate-state) contribution; the authors note it is absent with the vacuum pion propagator. Footnote 3 (p. 896) shows the gaps differ even for vanishing πN interaction, because the charged pion carries μ_π = μ_p − μ_n. Our inference: for an unsoftened pion 1−β = −μ_π²/(m_π²+k²−μ_π²) is small, so the unsoftened question reduces largely to the size of κ and of the vertex factors with a realistic interaction (step B). Exchange diagrams are added coherently there (p. 892).

**Scaling inference (not in the paper).** The relative weight of the two intermediate-state diagrams is set by Landau–Migdal vertex factors, not by softening; their common absolute size scales as the inverse eighth power of the pion gaps (Eq. 29). With vacuum-like gaps the formula would be strongly reduced, but Eq. (25) was derived for the soft-mode regime and the paper does not evaluate that limit, so this does not settle Q1.

**Open inference.** With a realistic correlated NN interaction and an effectively unsoftened pion, two-body weak currents can exist without anomalously soft pion denominators. Their absolute size, interference and possible dominance have **not** been established by the literature compared in this project. A modern NN amplitude alone is not a calculation of the compatible weak-current operators.

**Proposed interaction-matched A/B/C comparison:**
1. **A:** conventional external-leg modified Urca using a fully specified antisymmetrized NN scattering amplitude, retaining intermediate nucleon denominators.
2. **B:** from the same nuclear interaction, add pion-exchange and particle–hole two-body weak currents with unsoftened pion propagation; calculate direct, exchange and interference contributions.
3. **C:** consistently dress charged/neutral pion propagators, pion–nucleon vertices and weak currents while keeping the comparison baseline fixed. C−B would then probe genuinely collective softening effects.

Keep kinematics, effective masses, isospin/spin algebra and weak-current constraints compatible throughout; do not multiply independent enhancement factors or double-count overlapping cuts in spectral and explicit reaction descriptions. The A/B/C sequence is a research design, **not completed work**.

**Why existing papers are not direct tests.**
- FM79/Bottaro: explicit conventional OPE/ρ/Landau MU but no full matched intermediate-pion and particle–hole currents.
- Shternin et al. (2018): explicit BHF G-matrix and near-on-shell intermediate nucleon propagation, but not a full Eq. (24) intermediate-state weak-current calculation.
- Dehghan Niri, Moshfegh and Haensel, [PRC 93, 045806 (2016)](https://doi.org/10.1103/PhysRevC.93.045806) and [PRC 98, 025803 (2018)](https://doi.org/10.1103/PhysRevC.98.025803), provide variational/correlation calculations of ordinary modified-Urca branches; these do not quantify every omitted intermediate-pion current.
- Fischer et al., [arXiv:1608.05004](https://arxiv.org/abs/1608.05004), examine vertex dressing while holding the pion propagator unsoftened for **neutral-current NN bremsstrahlung**. This illustrates a separation of nuclear ingredients, but **not** the requested charged-current modified-Urca amplitudes.

**What would resolve Q1:** one fully specified NN interaction and compatible charged weak-current operator, amplitude-level decomposition and interference in A/B/C, quantified theoretical uncertainties and the absence of shared/double-counted diagrams.

## Q2 — Nuclear interaction and weak-current consistency

Born OPE+ρ+Landau (Bottaro), Argonne v18 plus Urbana IX BHF G-matrix (Shternin), and collision widths inferred from Paris-potential or phenomenological cross sections (Alford/Sedrakian) are different nuclear inputs. Their matrix elements and renormalization factors cannot be transferred as universal multipliers. Need an internally matched interaction and weak vertex; appropriate half-on-shell dependence where relevant.

## Q3 — Finite widths, direct-Urca crossover and diagram overlap

Alford and Sedrakian recover a smooth charged-current response in their spectral approximations. This does not imply they reproduce every modified-Urca exchange/interference current. Voskresensky Sec. 4 warns that dressed one-loop self-energy cuts and explicit reaction cuts may overlap. Need a consistent self-energy/vertex expansion and no double counting when comparing MU, MMU and DU near the threshold. Distinguish **reaction rates** from **emissivities**.

## Q4 — FM79 historical beta=0.68 and rho exchange

The historical beta_n≈0.68 is an emissivity-level correlation suppression, not a rho-meson-exchange factor. Source check closed (27 November 2026): it is R_Urca in FM79 Table 1, printed p. 552 — the ratio of Urca emissivities with correlated/uncorrelated OPE (step cutoff d=0.7 fm, Eqs. 12–14) at nuclear matter density, used in Eq. (65c). It is cutoff dependent (0.78–0.58 for d=0.55–0.85 fm), evaluated only at nuclear density, and FM79 compute rho exchange only for nnνν̄ (Sec. VII). Whether it remains valid after explicit rho exchange or BHF interactions is separately unresolved.

Yakovlev et al. (2001), Sec. 3.4, Eq. (140), restate FM79 with α_n=1.13 (FM79 Eq. 63a) and β_n=0.68; the proton branch (Yakovlev & Levenfish 1995, corrected) and muon branches are later additions, not FM79 results. Two apparent discrepancies were reconciled by re-evaluating FM79's own equations (our arithmetic): (i) α_n=1.13 includes the Landau (f′−g′) interference terms, since Eq. (56) with Eqs. (20)–(22) reproduces Eq. (62), 1.765−0.623(ρ0/ρ)^{2/3}, while pure OPE gives 1.33 at ρ0; (ii) FM79's practical formula uses two incoherent direct squares (×2 after Eq. 39), while the coherent sum with exchange is 21/16 of one square (Eq. 75), i.e. 21/32 of the practical value, matching the review's Eq. (139). On that reading the OPE part of α_n is ≈1.5× the coherent OPE result.

## Q5 — Pion dressing and medium-current contributions

A static screening-mass substitution is not the momentum/frequency-dependent collective pion propagator, nor a matched vertex correction. The relative charged/neutral pion gaps, damping and weak-exchange currents require separate assumptions. The dominance of intermediate-state currents without softening is Q1, not a settled subset of this question.

## Editorial protocol

Keep this file, the corresponding HTML cells and explanatory sections, and [docs/ROADMAP.md](ROADMAP.md) synchronized. Tag statements as original-source results, interpretations, or open hypotheses. Do not mix EOS-dependent numerical experiments into the master HTML. Review protocol: [AGENTS.md](../AGENTS.md). Revised 26 November 2026.
