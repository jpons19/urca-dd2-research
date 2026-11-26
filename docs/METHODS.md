# Methods and equations

## Inputs and units

Real CompOSE PT(GRDF2-DD2), cold beta-equilibrated 1D composition, plus 3D general-purpose micro fields at `T0=0.1000000` and `T1=0.1096478` MeV. Density range `0.16–0.98 fm^-3` chosen to concentrate on the homogeneous core. All rates quoted in `erg cm^-3 s^-1`; masses, momenta and self-energies in MeV. `ħc=197.3269804 MeV fm`.

The reduced fields are, in stored order `(V_n,S_n,V_p,S_p)`: `10051`, `10052`, `11051`, `11052` in the CompOSE micro convention. `mD_i=m_i-S_i`, `mL_i=sqrt(pFi^2+mD_i^2)`, and `epsilon_i(p)=sqrt(p^2+mD_i^2)+V_i`, provided mean-field self-energies have no momentum dependence. Composition obeys charge neutrality `Yp = Ye+Ymu` and beta equilibrium `mu_n-mu_p=mu_e`.

## Interpolation/extrapolation

- Cold 1D composition: shape-preserving cubic PCHIP in baryon density.
- General 3D micro and composition grids: rectangular bicubic splines in `ln(nB)` and charge fraction `Yq` for each of the first two T slices, with fixed `(T,Yq,nB)` array ordering.
- Separate general-purpose beta-equilibrium solution: a scalar root for `mu_n-mu_p-mu_e=0` at fixed density, using table lepton fractions. This is a **diagnostic**, not substituted for the independent cold composition in the Urca curves.
- Target `10^9 K = 0.08617333 MeV` and `10^8 K = 0.008617333 MeV` lie **below** the general-purpose 0.1 MeV minimum. For the T8 calculation, extrapolate each field `X` linearly in `T^2` at constant density and charge fraction:

  `X(T) = X(T0) + (T^2-T0^2)/(T1^2-T0^2) * [X(T1)-X(T0)]`.

  Extrapolation factor at `10^8 K` is approximately `-4.9073`. This is a numerical low-temperature sensitivity diagnostic, not a physical finite-temperature EOS derivation. Data are smooth and checked against the independent cold composition; two nearby slices do not prove global T² behavior.
- For rates at T9, `build_final.py` uses the 0.1 MeV micro slice; `run_dd2_T8_reference.py` recomputes T9 with T²-extrapolated micro fields in the **scaling check**. Both choices produce nearly identical results in the homogeneous core.

## Friman–Maxwell / NSCool

Historical neutron-branch prescription:

`Q_FM,n = 8.55e21 * (mL_n/m_n)^3 * (mL_p/m_p) * (kFe+kFmu)/1.68 * alpha_n * beta_n * T9^8`

where momenta `kF` are in fm⁻¹, `alpha_n=1.76-0.63*(1.68/kFn)^2` and `beta_n=0.68`. The historical muon extension is preserved exactly as a **reference**, although it is approximate.

## Bottaro, Caputo and Fiorillo

Neutron branch triangular approximation, Eqs. (3.28)–(3.29), with constant weak-propagator denominators, OPE + magnetic rho exchange + Landau couplings. Couplings: `C_rho=1.4`, `m_pi=139.57 MeV`, `m_rho=769 MeV`, `f_pi=1`, `F0'=0.7`, `G0=1.1`, `beta=0.68`, `gA=1.275`, `G=G_F*Vud`. `f' = π² F0'/(2 m_N pFn)` and `g = π² G0/(2 m_N pFn)` are retained in the **historical bare-mass normalization**, even though the Fermi-surface phase space uses EOS-derived Landau masses. This is an explicit approximation.

Define `eta_M=pFn²/(pFn²+m_M²)` and

`alpha_B = (1/2)(21/16) [ (32/7)(m_pi/f_pi)^4 (f'-g)^2 + 2(eta_pi-C_rho*eta_rho)^2 ]`.

The neutron-branch emissivity prefactor is directly implemented from their Eq. (3.28), with `mL_n^3*mL_p*pFp*T^8*alpha_B*beta` and the lepton correction `1+pFmu/mu_e`. This expression supplies the triangular baseline.

The nontriangular prescription evaluates Eq. (3.34) for each lepton using the allowed domain in Eq. (3.32). It reconstructs the four momentum transfers and integrates the direct, exchange and trace-verified interference terms using a fixed scrambled Sobol sequence of 16,384 points. The full angular integral is normalized to its triangular limit with the same DD2 Landau masses, couplings, short-range factor and constant virtual denominators. The electron channel retains the lepton velocity that is approximated by unity in the triangular baseline.

## Shternin, Baldo and Haensel

Published approximate enhancement `alpha_U=0.6 R_A+0.2 R_B` with the paper's conventional neutron-branch normalization. Intermediate denominators use momentum-dependent **zero-width** quasiparticle dispersions and the EOS-derived Landau masses. The `R_A` and `R_B` Fermi-surface angular integrals are evaluated analytically in `compare_compose_four.py`, using the lepton Fermi momenta and beta-equilibrated chemical potentials. This approximation is **not valid at the direct-Urca threshold**. Coefficients are imported from the original microscopic study and are not a new interaction-matched DD2 calculation.

## Alford, Haber and Zhang (adapted NWA)

Illustrative collision broadening with full Dirac-mass convolution structure but a simplified direct-Urca kernel. Prescribed nucleon width `W_n=W_p=T²/(5 MeV)`; each Breit–Wigner distribution is centered on the tabulated `mD_n`, `mD_p`. For off-shell mass `m`, keep kinetic Fermi energy `E*_F` fixed and assign `pF(m)=sqrt(max(0,E*F²-m²))`. The **degenerate** DU baseline is `Q_DU=(457π/10080) G² (1+3gA²) mL_n*mL_p*mu_e*T⁶`, weighted by the BW-convolved triangle-allowed probability. Integrate one Breit–Wigner dimension analytically and one numerically with SciPy adaptive quadrature. This is a methodological approximation: it is **not** the full relativistic finite-temperature polarization tensor or a controlled conserving treatment of weak vertices; Dirac mass spectral widths are not measured from the DD2 EOS.

## Temperature scaling and scope

The conventional neutron-branch models have `Q∝T⁸` when composition/fields are held fixed. Away from the DU threshold the adapted NWA numerically follows similar scaling with a width `∝T²`, but this is **an empirical result of this implementation** and need not be universal. Neutrino and antineutrino emission combined; electron+muon summed. Proton MU branch is not implemented for FM/Bottaro/Shternin. Avoid adding an independent DU rate to the NWA without diagrammatic overlap analysis.

## Code map

- `calculations/compare_compose_four.py`: constant definitions, five rate functions including `bottaro_nontriangular`, analytic `RA/RB` and NWA quadrature.
- `calculations/bottaro_angular.py`: direct and exchange amplitudes, interference, and four-dimensional quadrature.
- `calculations/build_final.py`: T9 source interpolation, beta-equilibrium diagnostic and plots.
- `calculations/run_dd2_T8_reference.py`: low-T T² extrapolation, T8 solution, T8 rates, scaling check and five reference plots.
- `calculations/test_*.py`: numerical regression checks on generated outputs.

## Pion screening sensitivity and theoretical context

The **fixed-background static screening-mass** Bottaro triangular/nontriangular experiment is implemented separately in `calculations/pion_softening/`; see [PION_SCREENING_SCAN.md](PION_SCREENING_SCAN.md) for scope, reproduction, and initial tests. It is not a full collective medium-modified Urca (MMU) calculation. Original FM/Bottaro/Shternin/Alford baseline calculations and output files remain unchanged.

The theoretical distinction between static screening, dressed neutral/charged modes, Landau–Migdal correlations and intermediate-state radiation is discussed in [PION_SOFTENING.md](PION_SOFTENING.md). For an **EOS-independent comparison of the original theories**, use [report/urca_physics_comparison.html](../report/urca_physics_comparison.html); implementation assumptions belong here, not in that HTML.
