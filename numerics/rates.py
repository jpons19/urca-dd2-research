"""All Urca emissivity prescriptions evaluated on one standardized EOS state.

Models (electron, muon) in erg cm^-3 s^-1 at T = 1e9 K:
  FM, Bottaro, Bottaro_nontri          neutron-branch MU (prescriptions.py)
  Bottaro_noSRC, Bottaro_nontri_noSRC  the same without the inherited beta_n = 0.68
  Shternin                             0.6 R_A + 0.2 R_B; NaN at/above the DU threshold
  NWA                                  adapted Alford width approximation; NaN for Skyrme states
  Sedrakian                            proton spectral width (sedrakian.py)
  DU                                   sharp degenerate DU; zero-width limit of the NWA kernel

beta_n = 0.68 (FM79 Table 1, R_Urca at d = 0.7 fm) was defined for OPE alone; Bottaro et al.
keep it on top of explicit rho exchange.  FM79's nn-bremsstrahlung results (Eqs. 90-91) show
that correlation and rho suppressions largely overlap, so the *_noSRC variants bracket this
ambiguity (a sensitivity range, not a recomputed correlation factor).
"""
import numpy as np

from . import prescriptions as P

MODELS = ['FM', 'Bottaro', 'Bottaro_nontri', 'Bottaro_noSRC', 'Bottaro_nontri_noSRC',
          'Shternin', 'NWA', 'Sedrakian', 'DU']
du_prefactor = P.du_prefactor


def du(d):
    """Degenerate DU emissivity with the sharp triangle condition, per lepton."""
    pref = du_prefactor(d)
    out = []
    for pl in (d['pe'], d['pm']):
        allowed = pl > 0 and d['pn'] <= d['pp'] + pl and d['pp'] <= d['pn'] + pl
        out.append(pref if allowed else 0.0)
    return tuple(out)


def nwa(d, epsrel=5e-4):
    """Lorentz-scalar (Dirac) mass convolution: defined only for RMF states."""
    if d['kind'] != 'RMF':
        return (np.nan, np.nan)
    return P.alford(d, epsrel)


def _sedrakian(d):
    from .sedrakian import sedrakian   # local import: sedrakian imports this module
    return sedrakian(d)


def all_rates(d):
    """Return {model: (electron, muon)} for one standardized state."""
    b, bnt = P.bottaro(d), P.bottaro_nontriangular(d)
    return {
        'FM': P.fm(d),
        'Bottaro': b,
        'Bottaro_nontri': bnt,
        'Bottaro_noSRC': tuple(q / P.BETA for q in b),
        'Bottaro_nontri_noSRC': tuple(q / P.BETA for q in bnt),
        'Shternin': P.shternin(d),
        'NWA': nwa(d),
        'Sedrakian': _sedrakian(d),
        'DU': du(d),
    }
