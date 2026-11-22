# Reproducibility verification — 22 November 2026

### Archives

Using the two uploaded original CompOSE archives and `tools/extract_compose.py`:

- Cold composition re-extraction: identical shapes, maximum absolute difference `1.14e-13` (roundoff from CSV formatting).
- General-purpose micro extract: exact equality for `nB`, `Tgrid`, `Yqgrid`, and all `Vn,Sn,Vp,Sp` values in the two stored T slices.
- General-purpose fractions extract: exact equality for `nB`, `Yqgrid`, and all stored fraction values in the first two T slices.

Only **homogeneous-core** cold matter is checked using `Yp=Ye+Ymu`: below the core, nuclei and other components must be accounted for. The extractor intentionally does not impose npeμ charge neutrality on the full cold table.

### Full regeneration

From the root of the packaged project, these operations were completed successfully:

```bash
python calculations/build_final.py
python calculations/run_dd2_T8_reference.py
python calculations/test_final.py
python calculations/test_T8_reference.py
python -m unittest discover -s tests -v
python prior_audit/murca_spin_isospin_audit.py
python prior_audit/murca_angular_factor_test.py
```

T9 and T8 plot/data reruns completed in ~8 and ~9 seconds respectively on the verification environment; timings are environment-dependent, and the production quality of the physics models is not implied by these quick runs.

Additional tests: four unittest cases pass (input axes; chemical/kinematic constraints; lepton sums; actual recalculation of all four rates at three benchmark densities). The `prior_audit` test confirms the known `g′` cancellation and the printed Eq. (3.31) diagnostic, but does not verify a published code implementation.
