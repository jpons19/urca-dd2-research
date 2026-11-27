# Reproducibility verification — 22 November 2026

### Archives

Using the two uploaded original CompOSE archives and `tools/extract_compose.py`:

- Cold composition re-extraction: identical shapes, maximum absolute difference `1.14e-13` (roundoff from CSV formatting).
- General-purpose micro extract: exact equality for `nB`, `Tgrid`, `Yqgrid`, and all `Vn,Sn,Vp,Sp` values in the two stored T slices.
- General-purpose fractions extract: exact equality for `nB`, `Yqgrid`, and all stored fraction values in the first two T slices.

Only **homogeneous-core** cold matter is checked using `Yp=Ye+Ymu`: below the core, nuclei and other components must be accounted for. The extractor intentionally does not impose npeμ charge neutrality on the full cold table.

### Full regeneration

From the root of the packaged project, these operations were completed successfully:

```bash
python calculations/build_final.py
python calculations/run_dd2_T8_reference.py
python calculations/test_final.py
python calculations/test_T8_reference.py
python -m unittest discover -s tests -v
python prior_audit/murca_spin_isospin_audit.py
python prior_audit/murca_angular_factor_test.py
```

T9 and T8 plot/data reruns completed in ~8 and ~9 seconds respectively on the verification environment; timings are environment-dependent, and the production quality of the physics models is not implied by these quick runs.

Additional tests: four unittest cases pass (input axes; chemical/kinematic constraints; lepton sums; actual recalculation of all four rates at three benchmark densities). The `prior_audit` test confirms the known `g′` cancellation and the printed Eq. (3.31) diagnostic, but does not verify a published code implementation.

# Primary-source verification — 27 November 2026

FM79 (ApJ 232, 541) checked page by page against the NASA ADS scan (not stored in the repository):

- beta_n≈0.68 = R_Urca, Table 1, printed p. 552, cutoff d=0.7 fm (R_Urca = 0.78, 0.68, 0.58 for d = 0.55, 0.7, 0.85 fm); correlated potential Eqs. (12)–(14), pp. 543–544; applied in Eq. (65c).
- Intermediate nucleon propagator: Eqs. (28a)–(28b), p. 547; lowest order in 1/m gives ±1/ω.
- rho exchange: Sec. VII, Eqs. (83)–(91), pp. 554–555, nnνν̄ only (C_ρ≈1.25).
- Exchange terms: Sec. VI; Urca ratio 21/16 (Eq. 75), not included in Eqs. (65)–(66).
- Branches: neutron branch only (process 1), electrons only; superfluidity discussed qualitatively (Sec. VIII).
- Recorded in report/urca_physics_comparison.html (FM79 section, matrix, Section 5).

Voskresensky (2001), arXiv:astro-ph/0101514, checked against the arXiv PDF text and the rendered printed p. 17 (not stored in the repository):

- Fig. 1, printed p. 4: emission from a nucleon leg versus from intermediate states.
- Eq. (24), p. 17: first diagram has the lepton pair emitted at a π⁺–π⁰ junction (charged-pion β-decay to π⁰), second a correlated neutron-hole loop between the pion lines, third the FOPE-like external-leg diagram.
- Eq. (25), p. 17: ε ∝ T9^8 (ρ/ρ0)^{10/3} with inverse fourth powers of the neutral and charged pion gaps, Γ^8 (Eq. 26), F1 (Eq. 27), pairing factors ζ (Eq. 28, p. 18); Eq. (29) MMU/FOPE ratio ∼10 at ρ0.
- Recorded in report/urca_physics_comparison.html (Voskresensky section, matrix, Section 5) and docs/THEORY_OPEN_QUESTIONS.md (Q1).

Yakovlev, Kaminker, Gnedin & Haensel (2001), arXiv:astro-ph/0012122, Sec. 3.4 (arXiv pp. 67–71), checked against the arXiv PDF text (not stored in the repository):

- Eq. (140): Q^(Mn) ≈ 8.1e21 (m_n*/m_n)^3 (m_p*/m_p) (n_p/n_0)^{1/3} T9^8 α_n β_n, stated as FM79's result with α_n=1.13 (FM79 Eq. 63a) and β_n=0.68.
- Eq. (139): α_n attributed to the OPE matrix element without interference (21/4 vs 8).
- Proton branch (Eqs. 141–142, from Yakovlev & Levenfish 1995, corrected) and muon branches are the review's additions.
- Reconciled by re-evaluating FM79 (our arithmetic): Eq. (56) with Eqs. (20)–(22) gives 1.765−0.623(ρ0/ρ)^{2/3}, reproducing Eq. (62), so Landau terms are included in α_n=1.13; the coherent exchange result (Eq. 75, 21/16 of one direct square) is 21/32 of the practical two-square result (Eq. 39), matching the review's Eq. (139).

Voskresenskiĭ & Senatorov (1986), Sov. Phys. JETP 63, 885, checked from a scanned copy supplied by the project owner (not stored in the repository):

- Interaction Eqs. (2.1)–(2.11); propagator G=±(ω+Δμ)^{-1}, Eq. (2.17), p. 888; Fig. 3 diagrams a–c, h, i, p. 889.
- Eqs. (4.3)–(4.5), p. 891: vector-current term ∝ 1−β+κ, absent with vacuum pion propagator; footnote 3, p. 896: β≠1 even for Π→0 because of μ_π.
- Exchange added coherently, Eqs. (4.6)–(4.9), p. 892; results Eqs. (4.10)–(4.11) and comparison with FM79, p. 892.
- Added to the master HTML bibliography as ref. [11] and to the Voskresensky section; further-reading entries renumbered [12]–[17].
