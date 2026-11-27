"""Bottaro et al. neutron-branch angular integral with finite Fermi momenta.

The direct and interference terms follow the spin-isospin trace in prior_audit.
The trace puts the factor of two on the linear axial terms and uses
1 - (khat dot lhat)**2 for the tensor angular factors.
"""

import numpy as np
from scipy.stats import qmc


# Fixed, scrambled low-discrepancy samples for reproducible density scans.
SAMPLES = qmc.Sobol(d=4, scramble=True, seed=42).random_base2(14)


def dot(a, b):
    return np.sum(a * b, axis=-1)


def cosine(a, b):
    return dot(a, b) / np.sqrt(dot(a, a) * dot(b, b))


def sine_squared(a, b):
    return 1 - cosine(a, b)**2


def direct_squared(fp,g,gk,gkt,hk,hkt,k,kt,gA=1):
    d=cosine(k, kt)
    axial=(6*(fp*fp+g*g)-2*(3*fp*gk+3*(fp+2*g)*gkt+fp*hk+(fp+2*g)*hkt)
           +3*gk*gk+9*gkt*gkt+2*gk*hk+hk*hk+6*gkt*hkt+3*hkt*hkt)
    vec=3*(gk-gkt)**2+2*(gk-gkt)*(hk-hkt)+(hk-hkt)**2+2*hk*hkt*(1-d*d)
    return gA*gA*axial+vec


def interference(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt,gA=1):
    """Corrected Eq. (3.31): quadratic axial terms are not multiplied by 2."""
    angular = (hk*hl*sine_squared(k,l)-2*hkt*hl*sine_squared(kt,l)-2*hk*hlt*sine_squared(k,lt)+2*hkt*hlt*sine_squared(kt,lt))
    axial = (-24*fp*g+2*(
        2*fp*(3*(gkt+glt)+hkt+hlt)
        +3*g*(gk+gkt+gl+glt)+g*(hk+hkt+hl+hlt)
        -0.5*(gl*(hk+hkt)+glt*(hk+5*hkt)+hk*hl-hkt*hl-hk*hlt+3*hkt*hlt
        +gk*(3*gl+3*glt+hl+hlt)+gkt*(3*gl+15*glt+hl+5*hlt))
        +0.5*angular))
    vector=( (gl-glt)*(hk-hkt)+(gk-gkt)*(hl-hlt)-(hk-hkt)*(hl-hlt)
            +3*(gl-glt)*(gk-gkt)+2*hk*hl*sine_squared(k,l)-2*hkt*hl*sine_squared(kt,l)-2*hk*hlt*sine_squared(k,lt)+2*hkt*hlt*sine_squared(kt,lt) )
    return gA*gA*axial+vector




def angular_ratio(pn, pp, pl, fp, g, ga, meson_couplings, samples=SAMPLES):
    """Return the integral relative to a constant triangular matrix element.

    The paper's integration variables are two transfer magnitudes and two
    azimuths. The weight accounts for the allowed interval in the first transfer.
    """
    if pp <= 0 or pl <= 0:
        return 0.0

    lt = pn + pp * (2 * samples[:, 0] - 1)
    lo = np.abs(lt - pl)
    hi = np.minimum(lt + pl, 2 * pn)
    valid = hi > lo
    if not np.any(valid):
        return 0.0

    lt = lt[valid]
    width = hi[valid] - lo[valid]
    l = lo[valid] + width * samples[valid, 1]
    c = np.clip((pl**2 - l**2 - lt**2) / (2 * l * lt), -1, 1)
    s = np.sqrt(1 - c*c)
    zero = np.zeros_like(c)
    lt_axis = np.column_stack((s, zero, c))
    lt_tangent = np.column_stack((c, zero, -s))

    # Choose l along z and l-tilde in the x-z plane.
    phi, beta = (2 * np.pi * samples[valid, 2:]).T
    c2 = l / (2 * pn)
    s2 = np.sqrt(1 - c2*c2)
    p2 = pn * np.column_stack((s2*np.cos(phi), s2*np.sin(phi), c2))
    c4 = np.clip((pn**2 - pp**2 - lt**2) / (2 * pp * lt), -1, 1)
    s4 = np.sqrt(1 - c4*c4)
    p4 = pp * (
        c4[:, None] * lt_axis + s4[:, None] * (
            np.cos(beta)[:, None] * lt_tangent
            + np.column_stack((zero, np.sin(beta), zero))
        )
    )
    lv = np.column_stack((zero, zero, l))
    ltv = lt[:, None] * lt_axis
    p1 = p4 + ltv
    p3 = p2 - lv
    k, kt = p1 - p3, p2 - p4

    gk, hk = meson_couplings(dot(k, k), g)
    gkt, hkt = meson_couplings(dot(kt, kt), g)
    gl, hl = meson_couplings(dot(lv, lv), g)
    glt, hlt = meson_couplings(dot(ltv, ltv), g)
    amplitude = (
        direct_squared(fp, g, gk, gkt, hk, hkt, k, kt, ga)
        + direct_squared(fp, g, gl, glt, hl, hlt, lv, ltv, ga)
        + interference(fp, g, gk, gkt, gl, glt, hk, hkt, hl, hlt,
                       k, kt, lv, ltv, ga)
    )
    _, h_tri = meson_couplings(pn**2, g)
    triangular = ga**2 * (12 * (fp - g)**2 + 21/4 * h_tri**2)
    return np.sum(width * amplitude) / (2 * pl * len(samples) * triangular)
