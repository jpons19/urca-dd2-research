"""Medium-modified Urca (MMU) of Voskresensky and collaborators: first-approximation inputs.

Approved reference parametrization (27 Nov 2026), all density inputs in units of n0 = 0.16 fm^-3
for every EOS (the published fits are functions of n/n0, not EOS quantities):

* Vertex factor  Gamma(n) = 1 / [1 + C (n/n0)^{1/3}],  C = 1.6
  (Blaschke, Grigorian & Voskresensky 2004 [BGV04], text near Eq. 5: C ~ 1.4-1.6;
   BGV13 PRC 88, 065805 Eq. 3 and Kolomeitsev & Voskresensky 2015 use 1.6;
   Schaab et al. 1997 [SVSWW97] Eq. 2 uses 1.4).
* Effective pion gap omega~^2(n) / m_pi^2: BGV04 Fig. 1, curve 1a, digitized from the vector
  figure in the arXiv source (no published formula); 1 at n <= 0.8 n0; saturated (curve 1b)
  above the critical density n_c.  Baseline n_c = 3 n0 (plateau 0.27).  Variant n_c = 2 n0 with
  plateau 0.216 (Grigorian, Voskresensky & Blaschke 2016 family): **our** construction, the
  curve-1a shape compressed in density to reach that plateau at 2 n0.  Same gap for pi0 and
  pi+-, alpha = 1 (no condensate), as in the cooling papers.
* Ratio form (main curve), BGV13 Eq. 3 / KV15:
      eps_MMU = eps_MU * 3 (n/n0)^{10/3} [Gamma(n)/Gamma(n0)]^6 / (omega~/m_pi)^8,  n >= 0.8 n0,
      eps_MMU = eps_MU                                                           n <  0.8 n0,
  with eps_MU = the FM79/NSCool neutron-branch rate of the same EOS state.
* Direct form (cross-check), SVSWW97 Eq. 1 (labelled "MU-VS86" there), with that paper's own
  C = 1.4 and the baseline gap:
      eps = 2.4e24 F1 (n/n0)^{10/3} (m_n*/m_N)^3 (m_p*/m_N) [m_pi/omega~]^8 Gamma^8 T9^8,
      F1 = 1 + 3/(4 Gamma^2) (n/n0)^{2/3}   (pion decay from intermediate states).
* Unsoftened variant (OUR extension for the open question Q1, not a published result):
  ratio form with the vacuum-like gap at k = p_Fn, omega~^2 = m_pi^2 + p_Fn^2, Gamma kept.
  Eq. (25)-type formulas were derived for the soft-mode regime, so this is a scaling probe only.

Known inconsistencies (logged, not resolved): Gamma^6 (ratio form) vs Gamma^8 (direct form,
V01 Eq. 26); prefactors agree only to O(1); the papers do not state which form their codes use.
The ratio form jumps by ~1.9 at 0.8 n0 by construction.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator

from . import prescriptions as P

N0 = 0.16          # fm^-3
MPI = P.MPI         # MeV, charged-pion mass
C_GAMMA = 1.6
N_LOW = 0.8        # n/n0 below which plain MU is used

# BGV04 Fig. 1 curve 1a, omega~^2/m_pi^2 versus n/n0 (digitized; agent-extracted from omegatil.ps).
_GAP_X = np.array([0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0])
_GAP_Y = np.array([1.000, 0.996, 0.971, 0.918, 0.838, 0.734, 0.614, 0.489, 0.370, 0.27])
_GAP_1A = PchipInterpolator(_GAP_X, _GAP_Y)
PLATEAU = {3.0: 0.27, 2.5: 0.240, 2.0: 0.216, 1.5: 0.191}


def gamma(x, C=C_GAMMA):
    """Vertex suppression factor at x = n/n0."""
    return 1.0 / (1.0 + C * np.cbrt(x))


def gap2(x, nc=3.0):
    """omega~^2 / m_pi^2 at x = n/n0 for critical density nc (in n0)."""
    if x <= N_LOW:
        return 1.0
    if nc == 3.0:
        return float(_GAP_1A(min(x, 3.0)))
    plateau = PLATEAU[nc]
    xx = N_LOW + (min(x, nc) - N_LOW) * (3.0 - N_LOW) / (nc - N_LOW)   # compressed 1a shape
    shape = (1.0 - float(_GAP_1A(xx))) / (1.0 - 0.27)
    return 1.0 - (1.0 - plateau) * shape


def ratio_factor(x, g2, C=C_GAMMA):
    """3 (n/n0)^{10/3} [Gamma/Gamma(n0)]^6 / omega~^8 (omega~ in m_pi)."""
    return 3.0 * x**(10.0 / 3.0) * (gamma(x, C) / gamma(1.0, C))**6 / g2**4


def mmu_ratio_form(d, nc=3.0, C=C_GAMMA, unsoftened=False):
    """Total (e + mu) MMU emissivity in the ratio form, erg cm^-3 s^-1."""
    eps_mu = sum(P.fm(d))
    x = d['nb'] / N0
    if x < N_LOW:
        return np.nan if unsoftened else eps_mu
    g2 = 1.0 + (d['pn'] / MPI)**2 if unsoftened else gap2(x, nc)
    return eps_mu * ratio_factor(x, g2, C)


def mmu_direct_form(d, nc=3.0, C=1.4):
    """SVSWW97 Eq. (1) with alpha = 1, zeta = 1 (no pairing); NaN below 0.8 n0."""
    x = d['nb'] / N0
    if x < N_LOW:
        return np.nan
    G = gamma(x, C)
    F1 = 1.0 + 3.0 / (4.0 * G**2) * x**(2.0 / 3.0)
    mN = 0.5 * (P.MN + P.MP)
    mass = (d['mLn'] / mN)**3 * (d['mLp'] / mN)
    return 2.4e24 * F1 * x**(10.0 / 3.0) * mass * G**8 / gap2(x, nc)**4 * P.T9**8


def all_mmu(d):
    return {
        'MMU': mmu_ratio_form(d, 3.0),
        'MMU_nc2': mmu_ratio_form(d, 2.0),
        'MMU_direct': mmu_direct_form(d, 3.0),
        'MMU_unsoftened': mmu_ratio_form(d, unsoftened=True),
    }
