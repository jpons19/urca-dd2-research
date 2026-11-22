Independent spin/isospin audit of arXiv:2406.18640v2
===================================================

Contents
--------
murca_spin_isospin_audit.py   16-dimensional spin-isospin trace,
                              direct/exchange channels, verification of
                              Eqs. (3.22), (3.30), (3.31).
murca_angular_factor_test.py  Illustrative Fermi-surface angular Monte
                              Carlo at pn=340, pp=pe=100 MeV.

Run
---
python murca_spin_isospin_audit.py
python murca_angular_factor_test.py

Requires Python 3 + NumPy (no SymPy/FeynCalc necessary).

Key result
----------
Eq. (3.31) of the submitted v2 PDF appears to multiply by 2 the quadratic
terms and tensor angular products within the axial t/u-interference bracket.
Writing the axial interference schematically as

  gA^2 [ -24 f' g + 2(L-Q+A) ]

where L is linear in momentum-dependent coefficients, Q is quadratic, and A
is the angular tensor-product sum, the independent spin/isospin calculation
gives

  gA^2 [ -24 f' g + 2 L - Q + A ].

The full t-channel formula (3.30) checks. In the pure constant g' limit, the
literal printed Eq. (3.31) implies a negative total reduced squared amplitude
(-24*gA^2*g'^2); the corrected expression gives zero, consistent with the
triangular-limit formula (3.22).

The illustrative angular Monte Carlo gives ~0.893 for printed/corrected
angular-integrated interaction factors, under the model parameters hardcoded
therein. This is NOT a physical cooling-emissivity prediction. The calculation
does not include medium polarization, spectral widths, or phase-space energy
integrations. The published numerical implementation has not been inspected,
so the error may be typographical rather than present in computed figures.

All tests are performed in the convention |M|^2 = 64 G^2 omega_nu omega_l
/omega^2 * reduced_factor (including gA^2 where appropriate).
