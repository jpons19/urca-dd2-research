# References and primary-source index

For theoretical differences among original papers, the master self-contained comparison is [`report/urca_physics_comparison.html`](../report/urca_physics_comparison.html). It is **EOS-independent**. Paper-specific source claims should be checked against the primary works and pinpointed to pages/equations; unverified references remain explicitly flagged. The implementation-specific DD2 methods belong in [`METHODS.md`](METHODS.md), not in that comparative HTML.


1. B. L. Friman and O. V. Maxwell, “Neutrino Emissivities of Neutron Stars,” *Astrophysical Journal* **232** (1979), 541–557; [original NASA ADS scan](https://adsabs.harvard.edu/pdf/1979ApJ...232..541F). Historical FOPE / NSCool baseline.
2. S. Bottaro, A. Caputo and D. F. G. Fiorillo, “Neutrino emission in cold neutron stars: Bremsstrahlung and modified urca rates reexamined,” *JCAP* **11** (2024), 015; [arXiv:2406.18640v2](https://arxiv.org/abs/2406.18640).
3. P. S. Shternin, M. Baldo and P. Haensel, “In-medium enhancement of the modified Urca neutrino reaction rates,” *Physics Letters B* **786** (2018), 28; [arXiv:1807.06569](https://arxiv.org/abs/1807.06569).
4. M. G. Alford, A. Haber and Z. Zhang, “Beyond modified Urca: The nucleon width approximation for flavor-changing processes in dense matter,” *Physical Review C* **110** (2024), L052801. **Our implementation is an adapted NWA, not the full calculation.**
5. D. N. Voskresensky, “Neutrino Cooling of Neutron Stars. Medium Effects,” (2001), [arXiv:astro-ph/0101514](https://arxiv.org/abs/astro-ph/0101514). See Secs. 2.1–2.4 (dressed pion propagator, Landau–Migdal interactions and vertex corrections); Sec. 3.1 (MMU diagrams and emissivity, especially Fig. 1 and Eq. 25); later discussion of non-equilibrium Green functions and widths.
6. “Nuclear correlations and neutrino emissivity from the neutron branch of the modified Urca process,” *Physical Review C* **93** (2016), 045806, [DOI](https://doi.org/10.1103/PhysRevC.93.045806). Independent variational/correlation treatment; not the Shternin near-shell derivation.
7. “Role of nuclear correlations and kinematic effects on neutrino emission from the modified Urca processes,” *Physical Review C* **98** (2018), 025803, [DOI](https://doi.org/10.1103/PhysRevC.98.025803). Proton-branch nuclear correlation treatment.
8. L. Suleiman, M. Oertel and M. Mancini, “Modified Urca neutrino emissivity at finite temperature,” [arXiv:2308.09819](https://arxiv.org/abs/2308.09819). Thermal-field-theory treatment and explicit phase-space integrations; distinct validity regime and approximations.
9. [CompOSE documentation and catalogue](https://compose.obspm.fr/) – tabulated PT(GRDF2-DD2) cold [EOS 302](https://compose.obspm.fr/eos/302) and general-purpose [EOS 301](https://compose.obspm.fr/eos/301).

See [PION_SOFTENING.md](PION_SOFTENING.md) for an equation-level discussion and scientific limitations. The links provide scientific provenance. No original article PDFs are redistributed in this repository.
