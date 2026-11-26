# AGENTS.md — Urca research repository protocol

This file gives **persistent instructions to future assistants, Codex agents and human contributors** working in this repository. Read it together with README.md and CONTRIBUTING.md before changing research material.

## Purpose and two-track architecture

1. **General theory (authoritative comparative reference):** `report/urca_physics_comparison.html` is a **self-contained HTML document** comparing **the original research papers**, independent of any equation of state, numerical data, source-code implementation, model-specific plots or cooling simulations. Its matrix covers physical ingredients across Friman–Maxwell (1979), Bottaro et al. (2024), Shternin et al. (2018), Alford et al. (2024), and Voskresensky (2001 review). Each row/cell should be supported by detailed expandable explanations and original-source references. Theories are **not** to be reduced to yes/no comparisons when there are material qualifications.
2. **Numerical applications (separate):** `calculations/`, `tests/`, `figures/`, `report/summary.html`, `report/pion_softening.html`, and `docs/METHODS.md` document EOS-dependent implementations, numerical assumptions, validation and exploratory sensitivity analyses. Those results may motivate questions for the theory track, **but never appear as an EOS-specific assertion inside the master theory HTML**.

The repository is a scientific research prototype, not a validated production cooling library. Do not imply that different approximations form a statistical error band.

## Agreed incremental theory-review procedure

1. **One approach or ingredient at a time.** FM79 and Bottaro have the first substantial detailed analysis; Shternin, Alford and Voskresensky have **provisional** overview cells pending their own full audits. Do not relabel provisional assessments as verified without revisiting the primary paper.
2. **Read the primary source before editing scientific claims.** Record original-paper page numbers (printed page and PDF page if different), equation/table/figure numbers and an exact source URL. If access or verification is incomplete, explicitly state the uncertainty and leave a TODO; do not invent citations or attribute later NSCool conventions to FM79.
3. **Compare physical ingredients explicitly:** OPE and pion propagator; rho exchange; contact/Landau–Migdal terms; short-range correlations (and their normalization); direct/exchange/interference; triangular/angular approximations; intermediate nucleon recoil and near-shell effects; effective masses; collision widths/spectral functions; weak vertices/medium-current diagrams; neutron/proton branches; limitations and overlap. Distinguish calculations performed in the original paper from qualitative discussion or later adaptations.
4. **Update matrix and explanation together.** Each changed cell links conceptually to the matching expandable detailed subsection; explain the approximation, the formula and its validity. Include historical source attribution and whether a phenomenological coefficient was freshly derived or inherited.
5. **Keep the HTML autonomous:** inline CSS and any JS; working internal navigation; bibliography and source links; accessible expandable sections. It must open locally without external libraries, plots or build steps. Use normal HTML with UTF-8 symbols. Preserve existing substantive content when refining one paper.
6. **Keep theory clean of applications:** no DD2/CompOSE material, numerical rates/ratios from our scans, code variable names, local project paths, EOS-dependent Fermi momenta, test statuses or fitted results in the comparison. Numerical constants **intrinsic to published prescriptions** (e.g. beta_n=0.68) are appropriate when accurately sourced.
7. **Separate fact, inference, and unresolved matters.** Avoid unproven double-counting claims: rho exchange and phenomenological short-range correlations can overlap but are not automatically identical. Do not claim a mass replacement reproduces collective MMU. Distinguish pion-exchange momentum q≈p_Fn from the *intermediate nucleon propagator* approximation.
8. **Review and verify:** check the HTML is syntactically valid and locally viewable, check anchor targets, cross-links, equations and text consistency, inspect the diff, and log source-verification questions rather than silently resolving them.

**Current open source-verification item:** The table/page attribution of the historical modified-Urca short-range factor beta_n≈0.68 to FM79, Table 1, printed ApJ p.552, was discussed but the exact scanned page was not independently reopened when preparing theory comparison v2. Retain the qualification until directly inspected. Cite the original [FM79 ADS scan](https://adsabs.harvard.edu/pdf/1979ApJ...232..541F). Recheck precise formula/definition before claiming historical priority, and maintain distinction from rho exchange.

## Active theoretical open issues

Maintain [docs/THEORY_OPEN_QUESTIONS.md](docs/THEORY_OPEN_QUESTIONS.md) alongside the master HTML. The key unsettled issue is the size of intermediate-pion and particle–hole weak currents without strong pion softening. Distinguish their existence, their enhancement in Voskresensky’s specific medium model, and their unknown contribution with a different realistic NN interaction. The matched A/B/C comparison is a proposal, not completed numerical work. Do not treat neutral-current bremsstrahlung or conventional BHF modified-Urca as direct tests of omitted charged-current two-body weak currents.

## Numerical changes follow a different protocol

- Specify EOS and composition, units, channels, mass convention (Dirac versus Landau), interaction assumptions, and the exact source formulas implemented.
- Protect reference CSV files/plots and baseline rate functions. Prefer standalone optional modules, explicit configuration parameters, targeted regression and analytic-limit tests and convergence studies.
- Treat inherited beta=0.68, Bottaro rho/contact terms, Shternin coefficients, and Alford spectral widths as different, potentially overlapping approximations, **not freely multiplicative corrections**.
- Any proposed full MMU extension must separately justify pion polarization, neutral/charged modes, pion–nucleon vertices, additional weak-current diagrams, collision widths and possible EOS feedback.
- Before a numerical commit, run relevant tests and report what was actually executed and the observed results; never assert a test passed without having run it.

## Repository hygiene and delivery

- Preserve published-paper PDFs outside the repo unless redistribution rights are clear. Link to original papers in `docs/REFERENCES.md`. Do not commit the large original EOS archives.
- Avoid deleting historical reports, figures or prior audits merely because the new master theory reference exists. Link them in the appropriate numerical or audit track; remove only clear duplications with approval.
- Before editing, inspect current main and local uncommitted changes. Prefer focused commits (theory/editorial versus numerical), readable messages, and no unrelated refactoring. Do not overwrite remote changes.
- Update README/CONTRIBUTING/ROADMAP only when the workflow or status changes. Keep the master theory reference as a **single editable HTML**; external mathematical notes are optional supporting material, not prerequisites to understanding it.
- After committing, provide the commit SHA, changed-file list, what was verified, and any unresolved primary-source checks.

## Practical entry points

- Theory overview: [report/urca_physics_comparison.html](report/urca_physics_comparison.html)
- Open issues: [docs/THEORY_OPEN_QUESTIONS.md](docs/THEORY_OPEN_QUESTIONS.md)
- Source list: [docs/REFERENCES.md](docs/REFERENCES.md)
- Numerical report: [report/summary.html](report/summary.html)
- Numerical methods and validation: [docs/METHODS.md](docs/METHODS.md), [docs/VALIDATION.md](docs/VALIDATION.md)
- Existing pion-softening study: [docs/PION_SOFTENING.md](docs/PION_SOFTENING.md), [docs/PION_SCREENING_SCAN.md](docs/PION_SCREENING_SCAN.md)
- Current priorities: [docs/ROADMAP.md](docs/ROADMAP.md)
