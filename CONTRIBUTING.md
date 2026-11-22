# Contributing to the Urca DD2 research prototype

This is research software: physical assumptions and provenance matter as much as code correctness.

Every physics change should specify its source equations, approximations, effective-mass convention, process/channel content and unit system. Preserve an independently reproducible baseline. For numerical changes, include both analytic-limit tests and convergence checks; report shifts in comparison curves quantitatively. Treat the existing Alford NWA as an explicitly **adapted** implementation until a full finite-temperature relativistic calculation has been independently validated.

Recommended sequence: create a topic branch, add a targeted regression test and a methods note, run `make test`, regenerate the affected figures, and describe which numerical outputs changed in the pull request.

Before publication, choose a license and cite upstream CompOSE data and the original works. Do not commit the original large EOS archives or redistributed copies of published papers without checking reuse permissions.
