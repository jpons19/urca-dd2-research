"""Modified-Urca and finite-width rate prescriptions used by the multi-EOS comparison.

Each function takes a standardized EOS state ``d`` (see ``eos_inputs``) and returns
(electron, muon) emissivities in erg cm^-3 s^-1 at T = 1e9 K.

* fm:      FM79 neutron branch in the legacy NSCool normalization,
           8.55e21 (m_n*/m)^3 (m_p*/m) (k_Fe + k_Fmu)/1.68 alpha_n beta_n T9^8 with
           alpha_n = 1.76 - 0.63 (1.68/k_Fn)^2 (FM79 Eq. 62 form) and beta_n = 0.68.
* bottaro: Bottaro et al. (2024) triangular neutron branch, Eqs. (3.28)-(3.29): OPE + magnetic
           rho exchange + Landau couplings f' = pi^2 F0'/(2 m_N p_Fn), g = pi^2 G0/(2 m_N p_Fn)
           (historical bare-mass normalization), beta_n = 0.68, constant weak denominators.
* bottaro_nontriangular: Eq. (3.34) with finite proton/lepton momenta (bottaro_angular.py),
           normalized to its own triangular limit; uses the corrected Eq. (3.31) interference
           (see validation_bottaro_spin_isospin.py).
* shternin: Shternin, Baldo & Haensel (2018) Eq. (14), alpha_U = 0.6 R_A + 0.2 R_B with analytic
           zero-width R_A, R_B; NaN within 1 MeV of (or above) the DU threshold.
* alford:  adapted Alford, Haber & Zhang (2024) nucleon-width approximation: Breit-Wigner Dirac-mass
           convolution (W = T^2/5 MeV) of the degenerate DU kernel. Not the paper's full
           finite-T relativistic response; RMF (Dirac-mass) states only.
"""
import numpy as np
from scipy.integrate import quad

from .bottaro_angular import angular_ratio

HBARC = 197.3269804
MN, MP = 939.565379, 938.272046
ME, MMU = .51099895, 105.6583755
GF = 1.1663787e-11; VUD = .974; GA = 1.275; G = GF * VUD
MPI = 139.57; MRHO = 769.; BETA = .68; CRHO = 1.4; FPI = 1.
T9 = 1.; T = .08617333262145; HALF_WIDTH = (T * T / 5.0) / 2.0
NATURAL_TO_CGS = 1.602176634e-6 / ((1.973269804e-11)**3 * 6.582119569e-22)


def fm(d):
    kfn = d['pn'] / HBARC
    alpha = 1.76 - .63 * (1.68 / kfn)**2
    mass = (d['mLn'] / MN)**3 * (d['mLp'] / MP)
    coeff = 8.55e21 * mass * alpha * BETA * T9**8 / 1.68
    return coeff * d['pe'] / HBARC, coeff * d['pm'] / HBARC


def bottaro_landau(pn):
    """Landau couplings f', g in the bare-mass normalization of the published values."""
    scale = np.pi**2 / (2 * MN * pn)
    return .7 * scale, 1.1 * scale


def bottaro(d):
    pn, pp = d['pn'], d['pp']
    fp, g = bottaro_landau(pn)
    e_pi = pn**2 / (pn**2 + MPI**2); e_rho = pn**2 / (pn**2 + MRHO**2)
    alpha = .5 * (21 / 16) * ((32 / 7) * (MPI / FPI)**4 * (fp - g)**2 + 2 * (e_pi - CRHO * e_rho)**2)
    coeff = 2 * 11513 * G**2 * GA**2 * FPI**4 / (120960 * np.pi * MPI**4)
    coeff *= d['mLn']**3 * d['mLp'] * pp * T**8 * alpha * BETA * NATURAL_TO_CGS
    return coeff, coeff * d['pm'] / d['mu'] if d['pm'] else 0.


def bottaro_meson_couplings(q2, g):
    rho = CRHO * q2 / (q2 + MRHO**2)
    pion = q2 / (q2 + MPI**2)
    scale = (FPI / MPI)**2
    # The constant isovector spin coupling g-prime equals g in this model.
    return g - scale * rho, -scale * (pion - rho)


def bottaro_nontriangular(d):
    pn, pp = d['pn'], d['pp']
    fp, g = bottaro_landau(pn)
    qe, qmu = bottaro(d)
    electron = angular_ratio(pn, pp, d['pe'], fp, g, GA, bottaro_meson_couplings)
    muon = angular_ratio(pn, pp, d['pm'], fp, g, GA, bottaro_meson_couplings)
    # Restore the electron velocity approximated by unity in bottaro().
    return qe * (d['pe'] / d['mu']) * electron, qmu * muon


def RA(d, pl):
    pn, pp = d['pn'], d['pp']; a = pp; lo = pn - pl; hi = pn + pl
    if lo <= a:
        return np.nan

    def antider(x):   # primitive of 1/(k^2 - a^2)^2
        return -x / (2 * a * a * (x * x - a * a)) - np.log((x - a) / (x + a)) / (4 * a**3)
    return 2 * d['mLp']**2 * d['mu']**2 / pl * (antider(hi) - antider(lo))


def RB(d, pl):
    pn, pp = d['pn'], d['pp']
    if pn <= pp + pl:
        return np.nan
    return 4 * d['mLn']**2 * d['mu']**2 / ((pn**2 - (pp + pl)**2) * (pn**2 - (pp - pl)**2))


def shternin(d):
    out = []
    for pl in (d['pe'], d['pm']):
        if pl <= 0:
            out.append(0.); continue
        if d['pn'] - d['pp'] - pl <= 1:
            out.append(np.nan); continue
        enh = .6 * RA(d, pl) + .2 * RB(d, pl)
        coeff = 8.1e21 * (d['mLn'] / MN)**3 * (d['mLp'] / MP) * (pl / d['mu']) \
            * (d['nb'] * d['yp'] / .16)**(1 / 3) * T9**8
        out.append(coeff * enh)
    return tuple(out)


def _cauchy_interval(a, b, center, half):
    if b <= a:
        return 0.
    return max(0., (np.arctan((b - center) / half) - np.arctan((a - center) / half)) / np.pi)


def nwa_probability(d, pl, epsrel=5e-4):
    """Breit-Wigner mass convolution of the triangle condition, reduced to one quadrature.

    Kinetic Fermi energies E_F* are fixed; an off-shell Dirac mass m gives k_F(m) = sqrt(E_F*^2 - m^2).
    """
    if pl <= 0:
        return 0.
    en = d['mLn']; ep = d['mLp']; mn0 = d['mnD']; mp0 = d['mpD']; w = HALF_WIDTH
    ulo = np.arctan((-en - mn0) / w); uhi = np.arctan((en - mn0) / w)

    def integrand(u):
        mn = mn0 + w * np.tan(u)
        kn = np.sqrt(max(0., en**2 - mn**2))
        kpmin = abs(kn - pl); kpmax = min(ep, kn + pl)
        if kpmin >= kpmax:
            return 0.
        mlo = np.sqrt(max(0., ep**2 - kpmax**2))
        mhi = np.sqrt(max(0., ep**2 - kpmin**2))
        probp = _cauchy_interval(mlo, mhi, mp0, w) + _cauchy_interval(-mhi, -mlo, mp0, w)
        return probp / np.pi
    pts = [ulo, uhi, 0.]
    for kt in [abs(d['pp'] - pl), d['pp'] + pl, abs(ep - pl), ep + pl]:
        if 0 < kt < en:
            root = np.sqrt(en * en - kt * kt)
            for s in (-1, 1):
                u = np.arctan((s * root - mn0) / w)
                if ulo < u < uhi:
                    pts.append(u)
    pts = sorted(set(float(u) for u in pts if ulo <= u <= uhi))
    return sum(quad(integrand, a, b, epsabs=4e-10, epsrel=epsrel, limit=150)[0]
               for a, b in zip(pts[:-1], pts[1:]) if b > a)


def du_prefactor(d):
    """Sharp degenerate DU emissivity for kinematically allowed configurations."""
    return (457 * np.pi / 10080) * G**2 * (1 + 3 * GA**2) * d['mLn'] * d['mLp'] * d['mu'] \
        * T**6 * NATURAL_TO_CGS


def alford(d, epsrel=5e-4):
    pref = du_prefactor(d)
    return pref * nwa_probability(d, d['pe'], epsrel), pref * nwa_probability(d, d['pm'], epsrel)
