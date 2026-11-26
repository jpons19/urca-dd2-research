"""Fixed-background Bottaro pion-screening sensitivity, two lepton channels.

Import alongside calculations/compare_compose_four.py; never mutate globals.
Rates remain at the base module's configured T (default 10^9 K).
"""
from functools import partial
import numpy as np
from .propagators import couplings


def triangular_factor(d, ratio, base):
    """Exact ratio of the Bottaro triangular spin amplitudes H/H_vac."""
    pn=d['pn'];fp,g=base.bottaro_landau(pn)
    rho=base.CRHO*pn**2/(pn**2+base.MRHO**2)
    scale=(base.FPI/base.MPI)**2
    def h(r):
        u=pn**2/(pn**2+(r*base.MPI)**2)
        return 12*(fp-g)**2 + (21/4)*scale**2*(u-rho)**2
    return h(ratio)/h(1.0)


def softened_rates(d, ratio, base, samples=None):
    """Return (tri_e,tri_mu,nontri_e,nontri_mu), erg cm^-3 s^-1.

    For nontriangular rates, the angular routine returns the ratio to the
    *triangular amplitude of the supplied interaction*. Multiply by the
    changed triangular amplitude before applying the angular correction.
    """
    if not np.isfinite(ratio) or ratio<=0:
        raise ValueError('ratio must be finite and positive')
    tri_e0,tri_mu0=base.bottaro(d)
    scale=triangular_factor(d,ratio,base)
    tri_e,tri_mu=tri_e0*scale,tri_mu0*scale
    fp,g=base.bottaro_landau(d['pn'])
    kernel=partial(couplings,screening_ratio=ratio,constants=base)
    from bottaro_angular import angular_ratio
    kw={} if samples is None else {'samples':samples}
    e_factor=angular_ratio(d['pn'],d['pp'],d['pe'],fp,g,base.GA,kernel,**kw)
    mu_factor=angular_ratio(d['pn'],d['pp'],d['pm'],fp,g,base.GA,kernel,**kw)
    # The e-channel triangular expression approximates v_e = 1.
    electron=tri_e*d['pe']/d['mu']*e_factor
    muon=tri_mu*mu_factor
    return tri_e,tri_mu,electron,muon
