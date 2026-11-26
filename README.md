# Urca emissivities with CompOSE PT(GRDF2-DD2)

A reproducible **research prototype** comparing five rate prescriptions from four literature-inspired approaches to neutron-star modified Urca / collisional charged-current neutrino emission, using **the same real CompOSE PT(GRDF2-DD2)** cold composition and general-purpose microscopic fields.

**Status:** 26 November 2026 (documentation update; numerical baseline remains from 22–23 November). **Temperatures:** \(10^8\) and \(10^9\) K. **Core density:** \(0.16\)–\(0.98\) fm\(^{-3}\). **Particle channels:** electron and muon, treated separately then summed.

**Start with [`report/summary.html`](report/summary.html)** (self-contained HTML with embedded figures and numerical table). The repo is designed for further research, **not yet a production cooling microphysics library**.

## Current research direction: pion softening and medium-modified Urca

**New (26 November 2026):** See [the pion-softening research note](docs/PION_SOFTENING.md) and its [HTML companion](report/pion_softening.html). It consolidates our Friman–Maxwell/Bottaro/Shternin/Alford comparison with Voskresensky's dressed-pion and intermediate-state diagrams. **This is documentation, not a new validated emissivity model.** The previous T8/T9 CSV outputs and benchmark plots are unchanged.

## Contents

| Path | Contents |
|---|---|
| `report/summary.html` | Standalone illustrated research report, works offline |
| `figures/` | Original calculation figures in PNG and PDF |
| `calculations/compare_compose_four.py` | Five rate prescriptions, setup and optional T9 scan |
| `calculations/build_final.py` | 0.1 MeV / cold comparison and T9 figure regeneration |
| `calculations/run_dd2_T8_reference.py` | \(T^2\) extrapolation to T8, beta equilibrium and T8 figures |
| `calculations/DD2_*` | Reduced CompOSE inputs, including the first two T slices |
| `calculations/PT_*.csv` | Computed EOS / Urca numerical results for T8 and T9 |
| `calculations/test_*.py` | Original baseline regression checks |
| `tests/test_reference.py` | Additional lightweight regression and sample recalculation |
| `tools/extract_compose.py` | Optional reproducible extraction from *both original* ZIPs |
| `docs/` | Methods, provenance, interpretation limits, verification log and research roadmap |
| `prior_audit/` | Independent Bottaro spin/isospin tests and preliminary Fig. 6 diagnostic (separate earlier work) |

The full upstream 182 MB general-purpose archive is intentionally **not committed**. The reduced `.npz` inputs (two temperature slices and four self-energy fields) are sufficient to regenerate the present figures; they are derived data, not the full CompOSE dataset.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate      # Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
python calculations/test_final.py
python calculations/test_T8_reference.py
python -m unittest discover -s tests -v
```

Recompute the figures and CSV files (both scripts write to `calculations/`):

```bash
python calculations/build_final.py
python calculations/run_dd2_T8_reference.py
```

The calculation includes adaptive Breit–Wigner quadrature and can take longer than the quick regression tests. To inspect only previously generated outputs, open `report/summary.html` or the figures in `figures/`.

If you have the source CompOSE archives, reproduce the reduced data extracts explicitly:

```bash
python tools/extract_compose.py --cold path/to/cold.zip --general path/to/general.zip --out calculations/
```

This will overwrite the supplied reduced inputs; rerun the calculations and tests afterwards.

## Model coverage at a glance

- **Friman–Maxwell (1979), NSCool neutron branch:** OPE-based historical `alpha_n` and `beta_n` suppression factors, fixed virtual denominators, legacy muon extension.
- **Bottaro et al. (2024), triangular neutron branch:** OPE + magnetic rho exchange + phenomenological Landau terms, `beta_n=0.68`; constant denominators. The nontriangular curve uses the same nuclear interaction and includes finite proton and lepton momenta in the four-dimensional angular integral.
- **Shternin, Baldo & Haensel (2018), approximate neutron branch:** `0.6 R_A + 0.2 R_B` with momentum-dependent, zero-width propagators. **Not** a new self-consistent Brueckner calculation for DD2.
- **Alford, Haber & Zhang (2024), adapted NWA:** Breit–Wigner *Dirac-mass* convolution with width `W=T^2/(5 MeV)` applied to a *degenerate leading-order* direct Urca kernel. **Not** the paper's full relativistic finite-temperature response; omitted medium vertex corrections.

**Do not treat the five curves as comparable uncertainty-band realizations of one systematic many-body calculation.** Their diagrammatic content differs. None of the DU Fermi-momentum triangle gaps crosses zero in the tabulated comparison range; the threshold crossover is not tested here.

## Primary documents

- [`docs/METHODS.md`](docs/METHODS.md) – complete calculation and precise assumptions.
- [`docs/VALIDATION.md`](docs/VALIDATION.md) – existing checks, quantitative values and what is still unverified.
- [`docs/ROADMAP.md`](docs/ROADMAP.md) – ordered research tasks and suggested issue list.
- [`docs/PROVENANCE.md`](docs/PROVENANCE.md) – CompOSE IDs, input-file hashes and extraction.
- [`docs/PRIOR_AUDIT.md`](docs/PRIOR_AUDIT.md) – spin–isospin audit underlying the nontriangular calculation.
- [`docs/REFERENCES.md`](docs/REFERENCES.md) – literature and CompOSE catalogue links.
- [`docs/PION_SOFTENING.md`](docs/PION_SOFTENING.md) – screening-mass vs collective-mode physics, Voskresensky MMU, overlap cautions and numerical work plan.

## GitHub usage

```bash
git init
git add README.md .gitignore .gitattributes requirements.txt report figures calculations docs tests tools
git commit -m "Initialize CompOSE DD2 Urca comparison prototype"
# Create an empty GitHub repository first, then:
git branch -M main
git remote add origin https://github.com/YOUR-ACCOUNT/YOUR-REPO.git
git push -u origin main
```

The included GitHub Actions workflow runs quick baseline tests on each push/pull request; use `make all` for full local regeneration. The repository is currently public; decide on a project license and reuse policy before wider distribution. Research papers and original CompOSE distributions have their own citation and reuse requirements. The `.gitignore` prevents accidentally adding the large original EOS ZIPs.
