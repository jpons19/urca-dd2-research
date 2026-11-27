"""Sedrakian (2024) Urca emissivity with a broadened proton spectral function.

Source: A. Sedrakian, Phys. Rev. Lett. 133, 171401 (2024), arXiv:2406.16183v3,
Eqs. (3)-(4), (13), (16), (17).  Implemented as the ratio to the sharp DU emissivity:

    Q = eps0 * I(gamma_p) / I_ref,

  I = int_0^inf dw_nu w_nu^3 int dw g(w) f_e(w_nu - w) L1(w)/c,
  L1(w)/c = int de [f(e) - f(e + w)] { 1/2 - atan[(e_min(w, e) - e)/(gamma_p(e)/2)]/pi }   (Eq. 16)
  I_ref: the same with the curly bracket set to 1 (fully allowed kinematics, L1/c = w),
  eps0:  sharp DU prefactor, taken identical to the repository's DU kernel
         (457 pi/10080) G^2 (1+3 g_A^2) m_Ln m_Lp mu_e T^6, so that Q -> DU above threshold.

Energies are measured from the respective chemical potentials (nonrelativistic
spectra eps(p) = p^2/2m* - mu, mu = p_F^2/2m*, as in the paper; m* = Landau mass of
the EOS state).  Width, Eq. (17):

    gamma_p(e) = a_p T^2 [1 + (e/(pi T))^2],
    a_p = (hbar v_Fp / mu_p) * sum_{tau=n,p} sigma_ptau n_tau / mu_tau,  sigma_np = 7 fm^2, sigma_pp = 4 fm^2.

Minimum proton energy e_min (our derivation; see docs/NUMERICS.md):  for a lepton pair
with momentum q = p_Fl and energy w, the neutron is a quasiparticle with kinetic energy
mu_n + e + w, i.e. momentum k_n = sqrt(2 m_n* (mu_n + e + w)).  Momentum conservation
allows proton momenta p >= |k_n - q|, so e_min = (k_n - q)^2 / 2 m_p* - mu_p.  At
T -> 0, w -> 0, e -> 0 this gives exactly theta(p_Fp + q - p_Fn).  The printed closed form
below Eq. (14), (m_n*/2q^2)(mu_n + w - mu_p + q^2/2m_p*)^2, does not reproduce the triangle
condition with the "+" sign (it does with "-"); this is logged as a source question.

Extensions beyond the paper (labelled as such): the muon channel uses the same formula
with q = p_Fmu; neutron spectral width neglected as in the paper.
"""
import numpy as np

from . import eos_inputs as E
from . import prescriptions as P
from . import rates as R

HBARC = E.HBARC
SIGMA_NP, SIGMA_PP = 7.0, 4.0      # fm^2, Sedrakian (2024) text below Eq. (17)

# Dimensionless integration grids (units of T).
# Step 0.1 T; converged to < 1e-4 relative against a 4x finer, 1.5x wider grid
# (tests/test_numerics.py and docs/NUMERICS.md).
_X = np.linspace(-40.0, 40.0, 801)           # w / T, energy carried by the lepton pair
_Y = np.linspace(-60.0, 60.0, 1201)          # e / T, proton energy relative to mu_p
_XNU = np.linspace(0.0, 60.0, 601)           # w_nu / T


def _fermi(z):
    return 0.5 * (1.0 - np.tanh(0.5 * z))


def width_coefficient(d):
    """a_p of Eq. (17) in MeV^-1 for state d (kinetic Fermi energies mu = p_F^2/2m*)."""
    mun, mup = d['pn']**2 / (2 * d['mLn']), d['pp']**2 / (2 * d['mLp'])
    vfp = d['pp'] / d['mLp']
    nn, np_ = d['nb'] * d['Yn'], d['nb'] * d['Yp']
    return (HBARC * vfp / mup) * (SIGMA_NP * nn / mun + SIGMA_PP * np_ / mup)


def _lepton_kernel(T):
    """J(x) = int dx_nu x_nu^3 f_e(x_nu - x) on the _X grid (dimensionless)."""
    xn = _XNU[None, :]
    integrand = xn**3 * _fermi(xn - _X[:, None])
    return np.trapezoid(integrand, _XNU, axis=1)


def _bose_weight():
    """g(x) with the x -> 0 singularity handled by pairing with L1 ~ x."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return 1.0 / np.expm1(_X)


def ratio(d, q, T=None, a_p=None, width_scale=1.0):
    """I(gamma_p)/I_ref for lepton momentum q (MeV).  width_scale multiplies a_p (tests)."""
    if q <= 0:
        return 0.0
    T = P.T if T is None else T
    a_p = width_coefficient(d) if a_p is None else a_p
    a_p = a_p * width_scale
    mn, mp = d['mLn'], d['mLp']
    mun, mup = d['pn']**2 / (2 * mn), d['pp']**2 / (2 * mp)

    w = _X[:, None] * T                        # (nx, 1)
    e = _Y[None, :] * T                        # (1, ny)
    kn = np.sqrt(np.clip(2 * mn * (mun + e + w), 0.0, None))
    emin = (kn - q)**2 / (2 * mp) - mup
    gamma = a_p * T**2 * (1.0 + (e / (np.pi * T))**2)
    if width_scale == 0.0:
        weight = (e >= emin).astype(float)
    else:
        weight = 0.5 - np.arctan((emin - e) / (0.5 * gamma)) / np.pi
    occ = _fermi(_Y[None, :]) - _fermi(_Y[None, :] + _X[:, None])
    L1 = np.trapezoid(occ * weight, _Y, axis=1) * T         # L1/c in MeV
    L1_ref = _X * T                                           # fully allowed: L1/c = w

    J = _lepton_kernel(T)
    g = _bose_weight()
    # g(x) * L1(x): finite at x = 0; replace the single singular point by its limit.
    i0 = np.argmin(abs(_X))
    def integrate(L):
        with np.errstate(invalid='ignore'):
            prod = g * L
        prod[i0] = 0.5 * (prod[i0 - 1] + prod[i0 + 1])
        return np.trapezoid(prod * J, _X)
    return integrate(L1) / integrate(L1_ref)


def sedrakian(d, T=None, width_scale=1.0):
    """(electron, muon) emissivities in erg cm^-3 s^-1; muon channel is our extension."""
    pref = R.du_prefactor(d)
    qe = pref * ratio(d, d['pe'], T, width_scale=width_scale)
    qm = pref * ratio(d, d['pm'], T, width_scale=width_scale) if d['pm'] > 0 else 0.0
    return qe, qm
