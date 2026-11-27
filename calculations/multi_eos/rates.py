"""Urca emissivity prescriptions evaluated on a standardized EOS state.

The prescriptions are those of calculations/compare_compose_four.py (imported,
not copied, so the DD2 baseline and this comparison share one implementation):
  FM79/NSCool legacy neutron branch, Bottaro triangular and nontriangular
  neutron branch, Shternin 0.6 R_A + 0.2 R_B neutron branch, adapted Alford NWA.
Added here: the degenerate direct-Urca (DU) emissivity with the sharp
triangle condition, i.e. the zero-width limit of the same kernel used in the
adapted NWA, so that DU and NWA are mutually consistent.

Applicability rules (enforced, returned as NaN rather than silently extended):
  * Shternin's zero-width R_A, R_B diverge as the DU gap closes and are
    undefined above threshold: NaN for gap <= 1 MeV (as in the baseline).
  * The adapted NWA is a Lorentz-scalar (Dirac) mass convolution and needs
    a relativistic dispersion: NaN for Skyrme (no Dirac mass).
  * FM/Bottaro MU formulas are evaluated at all densities; above the DU
    threshold they are subdominant and are kept only for comparison.
All rates are at T = 1e9 K (T9 = 1), electron + muon channels, erg cm^-3 s^-1.
"""
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'calculations'))
import compare_compose_four as base  # noqa: E402

MODELS = ['FM', 'Bottaro', 'Bottaro_nontri', 'Shternin', 'NWA', 'DU']


def du(d):
    """Degenerate DU emissivity with the sharp triangle condition, per lepton."""
    pref = (457 * np.pi / 10080) * base.G**2 * (1 + 3 * base.GA**2) * d['mLn'] * d['mLp'] \
        * d['mu'] * base.T**6 * base.NATURAL_TO_CGS
    out = []
    for pl in (d['pe'], d['pm']):
        allowed = pl > 0 and d['pn'] <= d['pp'] + pl and d['pp'] <= d['pn'] + pl
        out.append(pref if allowed else 0.0)
    return tuple(out)


def nwa(d, epsrel=5e-4):
    if d['kind'] != 'RMF':
        return (np.nan, np.nan)
    return base.alford(d, epsrel)


def all_rates(d):
    """Return {model: (electron, muon)} for one standardized state."""
    return {
        'FM': base.fm(d),
        'Bottaro': base.bottaro(d),
        'Bottaro_nontri': base.bottaro_nontriangular(d),
        'Shternin': base.shternin(d),
        'NWA': nwa(d),
        'DU': du(d),
    }
