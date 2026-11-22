# Research handoff — modified Urca rates with nuclear correlations

**Status:** 2026-11-06. Continue the theoretical/numerical audit without redoing validated work.

## 0. Project objective

Critically reconstruct and improve the neutron-star modified-Urca neutrino emissivity calculation in **S. Bottaro, A. Caputo, D. F. G. Fiorillo, _Neutrino emission in cold neutron stars: Bremsstrahlung and modified urca rates reexamined_, arXiv:2406.18640v2 (2024)**. Start with neutron branch, check exchange/interference, spin/isospin factors, phase space, short-range correlations and links to in-medium calculations, then proton branch. Distinguish actual algebraic mistakes from model dependence and manuscript typos. The user's latest request before this handoff was **to reproduce/compare the paper's Figs. 5–7 with independent calculations**. A preliminary Fig. 6 diagnostic already exists, but is **not** a fully matched APR-EOS reproduction.

Source PDF is included in the transfer ZIP as `2406.18640v2.pdf`; original path `/mnt/data/2406.18640v2.pdf`. Source also available as arXiv:2406.18640v2. Scripts, preliminary plot, CSV, README included.

## 1. Nuclear potential and signs — checked

Paper Eqs. (2.1)–(2.6):

V(k) = f + f' tau1.tau2 + g sigma1.sigma2 + g'_k (tau1.tau2)(sigma1.sigma2) + h'_k (tau1.tau2)(sigma1.khat)(sigma2.khat)

Let A_pi=f_pi^2/m_pi^2, u_M(k)=k^2/(k^2+m_M^2).

**Pion:** V_pi = -A_pi u_pi (sigma1.khat)(sigma2.khat) tau1.tau2.

**Magnetic rho:** V_rho^mag = -A_pi C_rho u_rho [sigma1.sigma2 - (sigma1.khat)(sigma2.khat)] tau1.tau2.

So **g'_k=g'-A_pi C_rho u_rho**, **h'_k=-A_pi(u_pi-C_rho u_rho)**. Correct relative signs; a common rho-subtraction replacement for all spin operators is wrong. In S12=3(sigma1.khat)(sigma2.khat)-sigma1.sigma2 basis:

V_pi+rho^spin = -(A_pi/3) tau1.tau2 [(u_pi-C_rho u_rho) S12 + (u_pi+2 C_rho u_rho) sigma1.sigma2].

Couplings m_pi~139.57 MeV, m_rho=769 MeV, C_rho=1.4, f_pi~1; g_rho^2/(4pi)=0.41, f_rho/g_rho=6.1 yield C_rho~1.43. The nonmagnetic central/vector rho term is not explicitly retained with full k-dependence; it could be partly absorbed in the Landau coefficients, but needs matching. This is a modeling limitation, not a proven algebraic error.

## 2. Landau interaction and correlations — conceptual caveats

The phenomenological coefficients f,f',g,g' were taken from Landau theory (F0'=0.7, G0=G0'=1.1) and extrapolated to non-forward scattering. Strictly distinguish the Landau quasiparticle interaction (forward particle-hole irreducible vertex), antisymmetrized on-shell quasiparticle scattering amplitude and Born potential. In simple one-channel Landau theory A_l=F_l/[1+F_l/(2l+1)], but **do not mechanically substitute** this into modified Urca: spin/isospin coupling, channel definitions, and finite momentum exchange matter. Need explicit convention comparison with Friman–Maxwell (1979); earlier commentary identified their use of an *already antisymmetrized* Landau interaction as a potential source of convention mismatch. Not yet shown to produce a numerical normalization error.

The original Friman–Maxwell SRC factors beta_nn=0.56, beta_np=0.66, beta_n-MU=0.68 were extracted using their own correlation model and are retained as multipliers in Bottaro et al. even after changing the interaction. This is uncontrolled transfer; not automatically double-counting. Earlier analysis noted that Friman–Maxwell considered rho and SRC simultaneously for nn bremsstrahlung, suggesting approximate but imperfect factorization; transfer to modified Urca needs explicit recalculation. Medium polarization, spectral broadening and effective vertices are NOT calculated in Bottaro et al.; do not describe their calculation as fully in-medium.

Relevant background (verify exact publication details before citing): Friman & Maxwell, ApJ 232 (1979) 541; Schwenk, Jaikumar & Gale, Phys. Lett. B 584 (2004) 241; Dehghan Niri, Moshfegh & Haensel, Phys. Rev. C 93 (2016) 045806; Shternin, Baldo & Haensel, Phys. Lett. B 786 (2018) 28, arXiv:1807.06569; Sedrakian, arXiv:2406.16183.

## 3. Independent spin–isospin matrix audit — **REPRODUCIBLY CHECKED**

Files: `murca_spin_isospin_audit.py`, `murca_spin_audit.py`, `murca_angular_factor_test.py`, `MURCA_AUDIT_README.txt`. Run:

    python murca_spin_isospin_audit.py
    python murca_angular_factor_test.py

The first script explicitly represents the two-nucleon spin–isospin Hilbert space (dimension 16), constructs weak charged axial-current before/after strong NN scattering, includes t- and u-channel exchange, computes spin traces directly, and compares with analytic Eqs. (3.22), (3.30), (3.31). Leading amplitude schematic [V,J_A]/omega with omega=omega_l+omega_nu~mu_l. The vector charge cancels at leading recoil order because the strong interaction conserves total isospin.

Re-run on 2026-11-06 reproduced:

- Landau contribution (same reduced normalization as paper): H_L=12(f'-g)^2, including mixed term -24 f'g; constant g' cancels in full antisymmetrized amplitude.
- Pure OPE triangular tensor coefficient = 21/4=5.25, agreeing with Eq. (3.22).
- Full **direct-channel** squared amplitude, paper Eq. (3.30), matches independent matrices (30 independent random momentum/coupling tests; max absolute error 2.84e-14).
- **Printed t-u interference Eq. (3.31) contains an apparent algebraic factor-of-two error.** Write its axial part as g_A^2[-24 f'g + 2(L-Q+A)], where L collects terms linear in momentum-dependent coefficients, Q quadratic products and A tensor-angular terms. Independent matrix calculation requires **g_A^2[-24 f'g + 2L-Q+A]**. With correction, 30 random tests yield max abs error 2.84e-14. The script has exact implemented expressions; use those instead of trying to reconstruct lengthy Eq. (3.31) from this summary.
- Decisive analytic counterexample: only constant g' nonzero, all f',g,h' zero, with g'_k=g'_{ktilde}=g'_l=g'_{ltilde}=g'. Direct and exchange pieces each 12 g_A^2 g'^2; literal printed interference = -48 g_A^2 g'^2, giving **negative** H=-24 g_A^2 g'^2. Correct interference = -24 g_A^2 g'^2, giving H=0, consistent with Eq. (3.22). This is a decisive INTERNAL contradiction in the printed text, not just numerical disagreement.
- Appendix A's elementary Pauli trace coefficient 384 (Eq. A.5) was verified; one can obtain it from two sigma traces and epsilon^{ijk} epsilon^{ijk}=6. Missing g_A in A.11–A.12 appears typographical. Eqs. (3.34), (3.37) as printed seem to omit explicit T^8; compare dimensional consistency and implementation.

**Important:** The authors' actual numerical code has NOT been inspected. The error in Eq. (3.31) might be typographical; figures may have been computed correctly. Do not claim published figures are wrong until full reproduction.

For one illustrative Fermi-surface angular MC in `murca_angular_factor_test.py`, with pn=340 MeV, pp=pe=100 MeV (model example), the ratio integrated literal/corrected was 0.8928786 with approximate block error 0.00217 (1500 accepted samples). This is NOT a physical emissivity or full EOS check.

## 4. Figure 6 reproduction — PRELIMINARY DATA AVAILABLE, NOT YET A FINAL MATCH

Files:
- `murca_figure_comparison.py`: actual 4-dimensional Fermi-surface angular integration over l,ltilde,phi,beta using scrambled Sobol quasi-Monte Carlo, Eq. (3.32), with both corrected and literal Eq. (3.31), plus optional rho, Landau, muons. Physical conversion to cgs and beta=0.68; T=1e9 K. Uses scipy/numpy.
- `murca_paper_digitize.py`: approximate digitization of paper Fig. 6 right panel; GRAPH-READ values not underlying published data.
- `murca_matched_comparison.py`: **infer proton p_F** from paper's own *triangular curve* at each neutron p_F, using Eq. (3.28), charge neutrality p_p^3=p_e^3+p_mu^3, mu_e=mu_mu. Hence matches one observable but NOT independently derived APR composition. This inference can hide coefficient/normalization mismatch.
- `murca_fig6_diagnostic_data.csv`: numerical points.
- `murca_fig6_preliminary_comparison.png`: comparison plot.

Preliminary values (emissivity units 10^21 erg cm^-3 s^-1), representative:

| pn [fm^-1] | paper Fig6 exact digitized | independently corrected (inferred composition) | literal printed (same) |
| 1.5 | 1.212 | 1.008 | 1.019 |
| 1.7 | 2.238 | 1.980 | 1.831 |
| 1.9 | 2.752 | 2.424 | 1.724 |
| 2.0 | 2.876 | 2.560 | 1.443 |
| 2.2 | 3.501 | 3.159 | 0.764 |

The exact paper curve (blue) and independent corrected curve (green) have similar density dependence, while literal Eq. (3.31) (red) turns over and decreases strongly at large pn. However corrected curve is ~10–20% below paper. Possible reasons include incorrect/estimated APR composition, digitization, kinematic/angular Jacobians, normalization, short-range terms, QMC error, approximations in paper. NEED resolve before declaring match. The previously displayed commentary that curves were 'close' is preliminary only.

The source figure's x-axis uses pn~1.5–2.2 fm^-1 and y-axis Q [10^21 erg/cm^3/s], T=10^9 K, bare nucleon masses, short-range beta_n=0.68, includes electrons+muons. Paper Fig. 5 is radius profile for APR 1 solar mass NS; Fig. 7 includes proton branch. Exact reproduction of radius plots requires their APR stellar profile and particle fractions, not present as numerical table in PDF.

Run scripts from directory containing them:

    python murca_matched_comparison.py
    python murca_figure_comparison.py

The matched comparison with Sobol 2^15 points per species and rho setting typically takes some CPU time. Set a fixed seed and perform convergence tests power=14..19. Verify angular domain and normalization independently. The plot/CSV were created locally and included in this handoff.

### Potential technical pitfalls to inspect before trusting figure comparison

1. `murca_figure_comparison.py` uses a specific parametrization of the four angular variables; verify the geometry, the signs of `k,kt,l,lt`, the theta ranges, the determinant/Jacobian and the `factor=8*pi**2*pp` multiplier against Eq. (3.32). Compare constant-|m|^2 integration with analytic Eq. (3.26) in the triangular limit.
2. `murca_matched_comparison.py` infers pp by fitting the *published triangular curve*: this is NOT independent EOS determination, and errors in triangular normalization shift inferred pp. Obtain true APR EOS composition, or make a transparent benchmark at fixed (pn,pp,pe,pmu) and request paper author's numbers.
3. Distinguish electron vs muon kinematics: use pe,pmu, mu_l correctly, beta equilibrium and charge neutrality. Muon threshold may affect curve near pn=1.5 fm^-1.
4. Check whether Eq. (3.34) coefficient implicitly included omitted T^8 (it must for units). Confirm total neutrino+antineutrino factor 2, symmetry S=1/2, correct g_A convention and geometric factor.
5. Independently reconstruct spin matrices at physical constrained momenta, and verify positivity pointwise. Literal Eq. (3.31) can become negative even when physical couplings are used in parts of the domain; test and plot.
6. Compare Fig. 5 and Fig. 7 only after accurate stellar composition and proton-branch matrix element are in hand.

## 5. Next recommended sequence

**Priority 1 — complete Fig.6 numerical validation:** (a) verify constant-integrand geometry/normalization and triangular-limit angular integral, (b) QMC convergence and error bars, (c) use exact or carefully documented EOS composition, (d) plot paper digitization + corrected vs literal Eq. (3.31), and state residuals. Existing preliminary figure is NOT final.

**Priority 2 — resolve manuscript vs code:** search journal/arXiv versions, errata and author code/supplement; check whether Eq. (3.31) typo was corrected; if unavailable, consider a concise technical query to authors with the pure g' negative-square counterexample.

**Priority 3 — proton branch:** independently derive/trace Eq. (3.37), all spin/isospin interference, phase-space threshold; compare Fig. 7.

**Priority 4 — nuclear corrections:** define a consistent antisymmetrized scattering vertex (distinguish Landau f vs on-shell T); evaluate rho's omitted vector-central residual; replace arbitrary inherited beta multipliers by explicitly integrated SRC, and compare Schwenk, Dehghan Niri, Shternin/Baldo/Haensel and spectral-function calculations. Examine behavior near direct Urca threshold and possible near-on-shell nucleon propagators.

## 6. Explicit limitations / what NOT to claim

- Full correlated in-medium MU rates have NOT been derived.
- A factor-of-two discrepancy is decisively confirmed in the **printed** Eq. (3.31) via spin traces and negativity test; whether it affects the paper's plotted numerics remains OPEN.
- Eqs. (3.30) and (3.22) have been checked; **NOT** all of proton branch, phase-space derivation or cooling integrations.
- Friman–Maxwell Landau antisymmetrization caveat is NOT yet a demonstrated bug.
- Auxiliary beta short-range factors are NOT proved to double-count rho; matching and density dependence require reassessment.
- Fig. 6 diagnostic uses composition inferred from the paper curve rather than an independent APR table, and graph-digitized approximate values.

## 7. Prompt for next chat

> Continue the modified-Urca theoretical and numerical audit of Bottaro, Caputo & Fiorillo arXiv:2406.18640v2. I attach `MURCA_RESEARCH_HANDOFF.md` and the ZIP with the original PDF, NumPy/SciPy spin-isospin audit scripts, angular integrations, numerical CSV and preliminary Figure 6 comparison. Read the handoff carefully and **do not repeat the completed Pauli-trace tests** unless debugging. The major rigorously verified finding is an extra factor 2 in the *printed* axial t-u interference Eq. (3.31): a pure constant g' test makes the printed full squared amplitude negative, whereas corrected spin traces give zero; Eq. (3.30) and triangular Eq. (3.22) pass the independent tests. Current figure comparison is preliminary because APR composition was inferred from a published curve and the angular integration still needs independent normalization checks. First validate the angular integral and its Jacobian/triangular limit, establish a genuine like-for-like EOS comparison, and finish polished, quantitative overlays of Figs. 5–7. Then assess whether the printed error affected the actual paper's numerical results and proceed to proton branch and correlation refinements.
