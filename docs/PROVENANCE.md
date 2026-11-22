# Data provenance and archiving

The analysis is based on two **user-supplied** CompOSE PT(GRDF2-DD2) archives:

- Cold, beta-equilibrated EOS catalogue [302](https://compose.obspm.fr/eos/302), archive originally named `eos(2).zip`.
- General-purpose EOS catalogue [301](https://compose.obspm.fr/eos/301), archive originally named `eos(1)(1).zip`.

Input files are not redistributed in this repository. Exact file sizes and SHA-256 hashes are recorded in [`input_archive_manifest.json`](input_archive_manifest.json), allowing a later GitHub checkout to verify it is using the same original archive version.

Reduced derived inputs in `calculations/`:

- `DD2_cold_composition.csv` – all 1035 cold 1D densities, baryon/lepton fractions and Fermi-momentum triangle gaps. Species CompOSE IDs are `10` (n), `11` (p), `0` (e), `1` (mu).
- `DD2_micro_first_two_T.npz` – original grid axes and first two T slices for fields `10051:Vn`, `10052:Sn`, `11051:Vp`, `11052:Sp`; array dimension order `(T,Yq,nB,field)`.
- `DD2_general_composition_first_two_T.npz` – first two T slices of n, p, e, mu fractions, same grid ordering, fields in `(Yn,Yp,Ye,Ymu)` order.

The approximate scientific models are reproduced in the code, not stored in the source archives. The plots and tabular rates are generated results, **not original CompOSE tabulations**. The general-purpose table starts at T=0.1 MeV; T8 data are explicitly extrapolated. The original 182 MB archive should be obtained from CompOSE for any work requiring additional temperatures, additional species or a new composition range. See `tools/extract_compose.py`.

The EOS homogeneous core is validated by charge neutrality, beta equilibrium, positive quasiparticle masses, and neutron/proton dispersion consistency. **Below the homogeneous core, nuclei and other charged components appear; do not apply `Yp=Ye+Ymu` to every row of the cold table.** These tests do not independently establish accuracy of the nuclear interaction or the modified Urca scattering diagrams.
