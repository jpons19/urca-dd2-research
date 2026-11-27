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

MODELS = ['FM', 'Bottaro', 'Bottaro_nontri', 'Bottaro_noSRC', 'Bottaro_nontri_noSRC',
          'Shternin', 'NWA', 'Sedrakian', 'DU']

# The historical short-range factor beta_n = 0.68 (FM79 Table 1, R_Urca at d = 0.7 fm) was
# defined for the OPE potential alone.  Bottaro et al. keep it on top of explicit rho
# exchange; FM79's own nn-bremsstrahlung results (Eqs. 90-91) show correlation and rho
# suppressions largely overlap.  The *_noSRC variants (beta = 1) bracket this ambiguity;
# they are a sensitivity range, not a recomputed correlation factor.



def du_prefactor(d):
    """Sharp degenerate DU emissivity for kinematically allowed configurations (erg cm^-3 s^-1)."""
    return (457 * np.pi / 10080) * base.G**2 * (1 + 3 * base.GA**2) * d['mLn'] * d['mLp'] \
        * d['mu'] * base.T**6 * base.NATURAL_TO_CGS


def du(d):
    """Degenerate DU emissivity with the sharp triangle condition, per lepton."""
    pref = du_prefactor(d)
    out = []
    for pl in (d['pe'], d['pm']):
        allowed = pl > 0 and d['pn'] <= d['pp'] + pl and d['pp'] <= d['pn'] + pl
        out.append(pref if allowed else 0.0)
    return tuple(out)


def nwa(d, epsrel=5e-4):
    if d['kind'] != 'RMF':
        return (np.nan, np.nan)
    return base.alford(d, epsrel)


def _sedrakian(d):
    from .sedrakian import sedrakian   # local import: sedrakian imports this module
    return sedrakian(d)


def all_rates(d):
    """Return {model: (electron, muon)} for one standardized state."""
    b, bnt = base.bottaro(d), base.bottaro_nontriangular(d)
    return {
        'FM': base.fm(d),
        'Bottaro': b,
        'Bottaro_nontri': bnt,
        'Bottaro_noSRC': tuple(q / base.BETA for q in b),
        'Bottaro_nontri_noSRC': tuple(q / base.BETA for q in bnt),
        'Shternin': base.shternin(d),
        'NWA': nwa(d),
        'Sedrakian': _sedrakian(d),
        'DU': du(d),
    }
