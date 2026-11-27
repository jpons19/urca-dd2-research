"""Minimal TOV solver for cold beta-equilibrated CompOSE 1D tables.

Used only to mark central densities of reference stars (1.4 Msun and Mmax)
on the emissivity plots.  Pressure and energy density are read from eos.thermo:
  p = Q1 * nB,  e = nB * m_n * (Q7 + 1)   [MeV fm^-3].
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator

MEVFM3_TO_KM2 = 1.602176634e32 * 6.67430e-11 / 299792458.0**4 * 1e6   # MeV fm^-3 -> km^-2
MSUN_KM = 1.476625061


def read_thermo(directory, thermo='eos.thermo', nbfile='eos.nb'):
    d = Path(directory)
    nb = np.loadtxt(d / nbfile, skiprows=2)
    with open(d / thermo) as f:
        mn = float(f.readline().split()[0])
    th = np.loadtxt(d / thermo, skiprows=1, usecols=range(10))
    idx = th[:, 1].astype(int) - 1
    n = nb[idx]
    p = th[:, 3] * n
    e = n * mn * (th[:, 9] + 1)
    keep = p > 0
    n, p, e = n[keep], p[keep], e[keep]
    order = np.argsort(n)
    n, p, e = n[order], p[order], e[order]
    # keep a strictly increasing pressure sequence (drops flat/mixed-phase points)
    keep, last = [], -np.inf
    for i, pi in enumerate(p):
        if pi > last * (1 + 1e-12):
            keep.append(i); last = pi
    return n[keep], p[keep], e[keep]


class TOV:
    def __init__(self, n, p, e):
        self.n, self.p, self.e = n, p, e
        lp = np.log(p)
        self.e_of_lp = PchipInterpolator(lp, np.log(e))
        self.n_of_lp = PchipInterpolator(lp, np.log(n))
        self.lp_of_n = PchipInterpolator(np.log(n), lp)
        self.pmin = p[0]

    def star(self, nc):
        pc = np.exp(self.lp_of_n(np.log(nc)))
        k = MEVFM3_TO_KM2

        def rhs(r, y):
            m, p = y
            if p <= self.pmin:
                return [0.0, 0.0]
            e = np.exp(self.e_of_lp(np.log(p)))
            dm = 4 * np.pi * r**2 * e * k
            dp = -(e + p) * k * (m + 4 * np.pi * r**3 * p * k) / (r * (r - 2 * m))
            return [dm, dp / k]

        surface = lambda r, y: y[1] - self.pmin * 1.0001
        surface.terminal = True
        r0 = 1e-6
        e0 = np.exp(self.e_of_lp(np.log(pc)))
        y0 = [4 / 3 * np.pi * r0**3 * e0 * k, pc]
        sol = solve_ivp(rhs, (r0, 60.0), y0, events=surface, rtol=1e-8, atol=[1e-12, 1e-14], max_step=0.05)
        return sol.t[-1], sol.y[0, -1] / MSUN_KM

    def sequence(self, ncs):
        return np.array([(nc, *self.star(nc)) for nc in ncs])


def reference_densities(directory, thermo='eos.thermo', nbfile='eos.nb', ncs=None):
    """Return dict with Mmax, n_c(Mmax), n_c(1.4), R(1.4) from a TOV sequence."""
    t = TOV(*read_thermo(directory, thermo, nbfile))
    if ncs is None:
        ncs = np.geomspace(0.2, min(2.0, t.n[-1]), 70)
    s = t.sequence(ncs)
    i = int(np.argmax(s[:, 2]))
    out = dict(Mmax=s[i, 2], nc_Mmax=s[i, 0], R_Mmax=s[i, 1], sequence=s)
    stable = s[:i + 1]
    if stable[0, 2] < 1.4 < stable[-1, 2]:
        out['nc_1.4'] = float(np.interp(1.4, stable[:, 2], stable[:, 0]))
        out['R_1.4'] = float(np.interp(1.4, stable[:, 2], stable[:, 1]))
    return out
