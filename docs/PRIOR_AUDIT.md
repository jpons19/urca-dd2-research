# Earlier Bottaro spin/isospin audit: separate research thread

Previous work (summarized in the prior conversation handoff dated 2026-11-06) independently evaluated neutron-branch spin/isospin matrix traces. It reported consistency of the triangular formula Eq. (3.22) and direct-channel squared amplitude Eq. (3.30), but an apparent factor-of-two inconsistency in the *printed* axial interference Eq. (3.31). A pure constant `g′` test makes the literal printed squared amplitude negative, while the corrected expression cancels. The authors' numerical code was not inspected: **it is unresolved whether the discrepancy is manuscript-only or present in published figures**.

**This repository's current Bottaro emissivity curve uses only the triangular Eqs. (3.28)–(3.29), not Eq. (3.31), so the alleged nontriangular typo does not propagate into these four model curves.**

The independent spin-matrix audit scripts, technical handoff and preliminary Figure 6 diagnostic are preserved in `prior_audit/`. The original paper PDF and paper-digitization utility are not redistributed in the repository; obtain the paper from arXiv if needed. Run `python prior_audit/murca_spin_isospin_audit.py` and `python prior_audit/murca_angular_factor_test.py`. Its figure comparison is only preliminary and uses inferred APR composition, **not** the DD2 EOS plotted in the main report.
