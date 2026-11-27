# Contributing to the Urca research project

This project has two distinct deliverables:

- **Theory reference:** [report/urca_physics_comparison.html](report/urca_physics_comparison.html) is an independently readable, EOS-agnostic comparison of **original publications**. Build it paper by paper; record the exact original equations, pages, figures, assumptions and unresolved issues. Never insert DD2 inputs, numerical scans, cooling-code normalizations, implementation details or results into this HTML.
- **Computational research:** `numerics/`, `results/`, `tests/` and [docs/NUMERICS.md](docs/NUMERICS.md) implement and assess the prescriptions on several CompOSE EOSs. Distinguish original-paper results from our adapted implementations. Preserve the frozen regression values in the tests.

**Read [AGENTS.md](AGENTS.md) before making scientific or editorial changes.** It is the detailed repository protocol.

For a **theory change**, work on one source/physical ingredient at a time: inspect the primary paper, record a precise citation, amend the comparison cell and its expandable explanation together, label uncertainties and provisional material, validate hyperlinks and self-contained HTML, and review the diff. Do not silently substitute later conventions for historical definitions.

For a **numerical change**, document the source equations, approximations, mass convention, channels, units and diagrammatic overlap. Add analytic-limit, regression and convergence checks; quantify changes in curves; run `make test` (and `make results` when outputs change) before committing. The current Alford NWA is *adapted*, not the paper's complete relativistic finite-temperature result.

Prefer a focused topic branch/PR for substantial changes, and use narrow commits with clear messages. Do not overwrite or casually delete prior data, validation logs, PDFs or papers. Ask before irreversible restructuring. Research papers and original CompOSE archives have separate reuse restrictions; do not commit copies without permission. The public repository currently has no selected license.
