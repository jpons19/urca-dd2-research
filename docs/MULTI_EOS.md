# Multi-EOS Urca comparison (numerical track)

**Status:** research prototype, local branch `numerics-multi-eos`, 27 November 2026. EOS-dependent numerical results; they do **not** belong in the theory HTML. The original DD2 baseline (`calculations/`, `report/summary.html`) is unchanged.

## Purpose

Compare the emissivity prescriptions already implemented for DD2 across three equations of state, using one common kinematic pipeline:

| EOS | CompOSE id | Type | Inputs used | Electron DU onset | Muon DU onset | M<sub>max</sub> (TOV) |
|---|---|---|---|---|---|---|
| PT(GRDF2-DD2) | 301 (micro), 302 (cold) | RMF, density-dependent couplings | Dirac mass, vector self-energy at T=0.1 MeV | none ≤ 0.98 fm⁻³ | none | 2.42 M☉ |
| BCP(NLRMF-N) | 350 | RMF, nonlinear (L=60 MeV) | Dirac mass (10041/11041), vector self-energy (10051/11051) at T=0.1 MeV | 0.507 fm⁻³ (M≈2.25 M☉) | none | 2.53 M☉ |
| PCP(BSk22) | 256 | Skyrme (Brussels–Montreal), unified | cold β-equilibrium table: fractions, Landau m*/m (10040/11040) | 0.332 fm⁻³ (M≈1.15 M☉) | 0.390 fm⁻³ (M≈1.42 M☉) | 2.26 M☉ |

Thresholds are from the exact triangle condition p<sub>Fn</sub> ≤ p<sub>Fp</sub> + p<sub>Fℓ</sub> with muons included (not the Y<sub>p</sub>=1/9 rule quoted on CompOSE data sheets). Masses at onset come from our TOV integration of the same tables.

## Consistency choices

1. **One state per density.** Every prescription sees the same n<sub>B</sub>, Y<sub>n,p,e,μ</sub>, Fermi momenta, μ<sub>e</sub>, and Landau masses m<sub>L</sub> = p<sub>F</sub>/v<sub>F</sub>:
   * RMF: ε(p) = √(p² + m<sub>D</sub>²) + V, so m<sub>L</sub> = √(p<sub>F</sub>² + m<sub>D</sub>²).
   * Skyrme: ε(p) = m + p²/2m* + U, so m<sub>L</sub> = m* (tabulated).
2. **Cold npeμ β-equilibrium with muons for both RMFs**, solved from the tabulated self-energies on the (n<sub>B</sub>, Y<sub>q</sub>) grid at the lowest temperature (0.1 MeV); the 3D tables themselves contain no muons. Validated on DD2 against the CompOSE cold table (|ΔY<sub>p</sub>| < 10⁻⁵). For BSk22 the cold table's own fractions are used.
3. **Same rate code as the DD2 baseline.** `calculations/multi_eos/rates.py` imports `compare_compose_four.py`; on DD2 the new pipeline reproduces the baseline CSV to < 2×10⁻⁴ (difference = β-solve vs cold table).
4. **DU and NWA share one kernel.** The sharp DU rate is the zero-width limit of the degenerate kernel inside the adapted NWA; above threshold NWA → DU (tested to 1%).
5. **Applicability is enforced, not extrapolated.** Shternin's R<sub>A,B</sub> are returned as NaN at/above the DU threshold (gap ≤ 1 MeV). The adapted NWA is a Lorentz-scalar (Dirac) mass convolution and is not defined for the Skyrme EOS (NaN for BSk22). FM and Bottaro MU formulas are evaluated everywhere for comparison; above the DU threshold they are subdominant.
6. **Temperature** T = 10⁹ K for all; electron + muon channels summed; neutron-branch MU only.

## Results (T = 10⁹ K)

Figures: `figures/multi_eos/urca_multi_eos_T9.png` (emissivities and ratio to FM; DU thresholds and central densities of 1.4 M☉ and M<sub>max</sub> stars marked) and `figures/multi_eos/eos_inputs_multi.png` (fractions, effective masses, DU gaps). CSV per EOS in the same folder.

Emissivities (erg cm⁻³ s⁻¹) at n<sub>B</sub>=0.16 fm⁻³ and at the central density of a 1.4 M☉ star:

| EOS | n<sub>B</sub> | m<sub>L,n</sub>, m<sub>L,p</sub> (MeV) | DU gap (MeV) | FM | Bottaro tri | Bottaro non-tri | Shternin | NWA (adapted) | DU |
|---|---|---|---|---|---|---|---|---|---|
| DD2 | 0.160 | 613, 535 | 71.6 | 5.8e20 | 3.0e20 | 3.3e20 | 2.7e21 | 9.7e21 | 0 |
| DD2 | 0.353 | 522, 371 | 42.4 | 5.1e20 | 2.0e20 | 2.4e20 | 1.7e21 | 7.4e21 | 0 |
| BCP | 0.160 | 680, 612 | 69.5 | 9.1e20 | 4.8e20 | 5.2e20 | 5.7e21 | 1.5e22 | 0 |
| BCP | 0.340 | 521, 383 | 23.0 | 5.4e20 | 2.2e20 | 2.7e20 | 3.7e21 | 1.5e22 | 0 |
| BSk22 | 0.160 | 847, 669 | 75.1 | 1.9e21 | 9.8e20 | 1.1e21 | 1.3e22 | — | 0 |
| BSk22 | 0.385 | 480, 469 | −26.2 | 6.1e20 | 2.4e20 | 3.2e20 | — | — | 7.5e26 |

Observations (numerical, prescription-dependent):

* **EOS effective masses matter as much as the prescription.** At saturation density the same FM formula differs by ×3.2 between DD2 and BSk22, almost entirely from (m<sub>L,n</sub>)³m<sub>L,p</sub>. The prescription spread at fixed EOS is larger (Bottaro/FM ≈ 0.4–0.6; Shternin/FM ≈ 3–7; NWA/FM ≈ 15–30 far below threshold).
* **Approaching threshold** (BCP, BSk22) Shternin's zero-width factor grows ∝ 1/Δp and must be cut off; the adapted NWA rises continuously into the DU rate (BCP panel). BCP stays within ~8 MeV of the electron threshold up to 1.1 fm⁻³, so it probes the near-threshold regime over a wide density range.
* **BSk22 Landau masses fall to ~0.15 m** at the highest densities, which makes all MU formulas decrease with density there; this is a property of the Skyrme functional, not of the Urca prescriptions.
* **The adapted NWA lies ~15–30× above FM far below threshold** for both RMFs. This is the known behavior of our adapted (degenerate-kernel) NWA and is not a validated microscopic MU rate; see `docs/METHODS.md` and the theory HTML (Alford section) for why width-only bare-vertex responses need not reproduce the MU amplitude.

## Sedrakian (2024) proton spectral width

Implemented in `calculations/multi_eos/sedrakian.py` from PRL 133, 171401 (arXiv:2406.16183v3), Eqs. (13), (16), (17), as the ratio to the sharp DU emissivity:

* Q = ε<sub>0</sub> · I(γ<sub>p</sub>)/I<sub>ref</sub>, where I is the energy integral of Eq. (13) with L<sub>1</sub> from Eq. (16), and I<sub>ref</sub> the same integral with fully allowed kinematics. ε<sub>0</sub> is the repository's DU kernel (457π/10080)G²(1+3g<sub>A</sub>²)m<sub>L,n</sub>m<sub>L,p</sub>μ<sub>e</sub>T⁶, so the result equals the sharp DU rate above threshold (tested to 10⁻³). The paper's Eq. (4) writes g<sub>A</sub>² only; using the standard (1+3g<sub>A</sub>²) keeps Sedrakian, NWA and DU on one normalization.
* Width, Eq. (17): γ<sub>p</sub> = a<sub>p</sub>T²[1 + ((ε−μ<sub>p</sub>)/πT)²], a<sub>p</sub> = (ħv<sub>Fp</sub>/μ<sub>p</sub>) Σ<sub>τ</sub> σ̄<sub>pτ</sub>n<sub>τ</sub>/μ<sub>τ</sub> with σ̄<sub>np</sub> = 7 fm², σ̄<sub>pp</sub> = 4 fm² and μ<sub>τ</sub> the kinetic Fermi energies p<sub>F</sub>²/2m*. For our EOSs a<sub>p</sub> ≈ 0.02–0.06 MeV⁻¹, i.e. γ<sub>p</sub>(μ<sub>p</sub>) ≈ 2–5×10⁻⁴ MeV at T = 10⁹ K.
* Spectra: nonrelativistic, ε = p²/2m* − μ with m* the EOS Landau mass (natural for Skyrme; a Fermi-surface approximation for the RMFs). Neutron width neglected, as in the paper.
* **Minimum proton energy (our derivation).** For lepton momentum q = p<sub>Fℓ</sub> and energy ω, the quasiparticle neutron has kinetic energy μ<sub>n</sub> + ε + ω, so k<sub>n</sub> = √(2m<sub>n</sub>*(μ<sub>n</sub>+ε+ω)) and momentum conservation allows p ≥ |k<sub>n</sub> − q|, i.e. ε<sub>min</sub> = (k<sub>n</sub>−q)²/2m<sub>p</sub>* − μ<sub>p</sub>. This reduces exactly to the triangle condition at T → 0. **Source question:** the closed form printed below Eq. (14), (m<sub>n</sub>*/2q²)(μ<sub>n</sub>+ω−μ<sub>p</sub>+q²/2m<sub>p</sub>*)², gives p<sub>Fn</sub> ≤ √((p<sub>Fp</sub>+q)² − 2q²) for equal masses with the "+" sign, not the triangle condition stated in the text; with "−q²/2m<sub>p</sub>*" it is consistent. Not resolved with the author.
* **Muon channel:** our extension (the paper omits muons), same formula with q = p<sub>Fμ</sub>.
* Numerics: dimensionless grids with step 0.1 T (ω/T ∈ [−40, 40], ε/T ∈ [−60, 60], ω<sub>ν</sub>/T ∈ [0, 60]); relative error < 2×10⁻⁴ against a 4× finer and 1.5× wider grid even 0.5 MeV from threshold. Below threshold the rate is linear in a<sub>p</sub> (tested).

**Result.** Far below threshold the Sedrakian rate is ≈ 4–8 × FM (DD2 7.2 at 0.16 fm⁻³, 4.3 at n<sub>c</sub>(1.4 M☉); BCP ≈ 7; BSk22 ≈ 6), i.e. above Bottaro and FM, comparable to Shternin, and 2–4 times below the adapted Alford NWA. On BCP and BSk22 it closely tracks Shternin's curve below threshold, including the near-threshold rise, and then joins the DU rate continuously. This numerical similarity of two physically different mechanisms (off-shell proton strength vs. near-on-shell intermediate nucleon) is a coincidence of magnitude for these parameters, not a demonstrated equivalence; both depend on interaction inputs (cross sections vs. BHF G-matrix) chosen independently.

## Inherited short-range factor in Bottaro (β<sub>n</sub> = 0.68)

Both Bottaro curves keep the historical β<sub>n</sub> = 0.68, as in the paper. This factor is FM79's R<sub>Urca</sub> (Table 1, p. 552): correlated/uncorrelated **OPE-only** emissivity for a 0.7 fm cutoff. FM79's own bremsstrahlung results (Table 1 and Eqs. 90–91) show that once ρ exchange is included, correlations add only a factor ≈0.85 (0.39/0.46), not the OPE-only 0.56. Applying 0.68 on top of explicit ρ exchange therefore probably over-suppresses; Landau–Migdal g′ may also absorb part of the same physics. The figure shades the range between β<sub>n</sub> = 0.68 and β<sub>n</sub> = 1 (columns `Bottaro*_noSRC_*`, a factor 1/0.68 ≈ 1.47). This is a sensitivity range, not a recomputed Urca correlation factor; FM is kept with its own β<sub>n</sub>, for which the factor was defined.

## Known input issues

* **BSk22 proton dispersion:** in the CompOSE table, m<sub>p</sub> + p<sub>Fp</sub>²/2m*<sub>p</sub> + U<sub>p</sub> differs from μ<sub>p</sub> by up to ~35 MeV in parts of the core, while the neutron relation holds to 10⁻⁴ MeV. The rates use only Fermi momenta and m*<sub>p</sub>, not U<sub>p</sub> or μ<sub>p</sub>; the discrepancy is documented by a test and should be raised with the table providers before m*<sub>p</sub> is trusted at the 10% level.
* **BSk22 μ<sub>e</sub>:** the tabulated μ<sub>e</sub> exceeds the free-electron value for the tabulated Y<sub>e</sub> by ≈0.1% (0.1–0.45 MeV); the free-electron value is used.
* **RMF tables at T=0.1 MeV** rather than T=0: thermal shifts of μ are ~2×10⁻⁴ MeV (negligible here).

## Reproduce

```bash
# 1. Download the CompOSE tables into external_eos/ (gitignored; not redistributed)
#    BSk22:  https://compose.obspm.fr/download/1D/NS/Skyrme/BSK22/eos.zip       -> external_eos/BSk22/
#    BCP:    https://compose.obspm.fr/download/3D/Barman/Nucleonic/eos.zip      -> external_eos/BCP_NLRMF_N/  (60 MB)
#    DD2 cold (for TOV only): https://compose.obspm.fr/download/1D/NS/RMF/DD2_Typel/eos.zip -> external_eos/DD2_cold/
# 2. Derived extracts (already in calculations/multi_eos/data/):
python tools/extract_compose_slice.py --dir external_eos/BCP_NLRMF_N --out calculations/multi_eos/data/BCP_NLRMF_N_T0.npz
# 3. Scan, figures, CSV (~20 s; TOV results cached in figures/multi_eos/stars.json)
python -m calculations.multi_eos.run_compare --points 120
# 4. Tests (TOV test skipped without external_eos/)
python -m unittest tests.test_multi_eos -v
```

## Not done / next steps

* Proton-branch MU (Yakovlev et al. 2001 rescaling) for FM/Bottaro/Shternin, now relevant because BSk22 reaches Y<sub>p</sub> > 0.3.
* Superfluid reduction factors; finite-temperature (non-degenerate) DU kernel near threshold.
* A third RMF with intermediate L, or a microscopic (BHF/APR-type) EOS, for a Shternin-consistent comparison.
