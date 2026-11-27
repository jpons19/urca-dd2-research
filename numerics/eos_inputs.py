"""Standardized cold-matter inputs for the multi-EOS Urca comparison.

Every loader returns a callable ``state(nb)`` giving a dict with the same keys,
so that all rate prescriptions see the same kinematic background:

  nb [fm^-3]; Yn, Yp, Ye, Ym; Fermi momenta pn, pp, pe, pm [MeV];
  mun, mup, mue [MeV] (total chemical potentials, rest mass included for nucleons);
  mLn, mLp [MeV]: Landau (Fermi-surface) effective masses, m_L = p_F / v_F;
  mnD, mpD [MeV]: Dirac effective masses (RMF only; NaN for Skyrme models);
  kind: 'RMF' or 'Skyrme'; residual diagnostics.

Mass conventions
  RMF: epsilon(p) = sqrt(p^2 + m_D^2) + V, so m_L = sqrt(p_F^2 + m_D^2).
  Skyrme: epsilon(p) = m + p^2 / (2 m*) + U, so m_L = m* (tabulated m*/m).

RMF tables: CompOSE general-purpose micro data (Dirac mass, vector self-energy)
on an (nB, Yq) grid at the lowest tabulated temperature (0.1 MeV).  Cold npe-mu
beta equilibrium (with muons, which the 3D tables omit) is solved here from those
self-energies.  Skyrme (BSk22): the cold 1D beta-equilibrium table is used directly.
"""
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator
from scipy.optimize import brentq

HBARC = 197.3269804
ME, MMU = 0.51099895, 105.6583755
ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / 'data'
EXTERNAL = ROOT / 'external_eos'


def kf(nb, y):
    return HBARC * np.cbrt(3 * np.pi**2 * nb * y) if y > 0 else 0.0


def lepton_fractions(nb, mue):
    """Electron and muon fractions for a common lepton chemical potential."""
    def frac(m):
        if mue <= m:
            return 0.0
        p = np.sqrt(mue**2 - m**2)
        return (p / HBARC)**3 / (3 * np.pi**2 * nb)
    return frac(ME), frac(MMU)


def mue_from_charge(nb, yp):
    """Electron chemical potential such that Ye + Ymu = Yp."""
    f = lambda mu: sum(lepton_fractions(nb, mu)) - yp
    hi = kf(nb, yp) + MMU + 50.0
    return brentq(f, ME * (1 + 1e-12), hi, xtol=1e-12, rtol=1e-14)


def _finish(nb, yp, mue, nucleon, name, kind, leptons=None):
    ye, ym = lepton_fractions(nb, mue) if leptons is None else leptons
    yn = 1.0 - yp
    pn, pp, pe, pm = kf(nb, yn), kf(nb, yp), kf(nb, ye), kf(nb, ym)
    d = dict(name=name, kind=kind, nb=nb, Yn=yn, Yp=yp, Ye=ye, Ym=ym,
             pn=pn, pp=pp, pe=pe, pm=pm, mue=mue, mu=np.hypot(pe, ME))
    d.update(nucleon(pn, pp))
    d['beta_residual'] = d['mun'] - d['mup'] - mue
    d['charge_residual'] = yp - ye - ym
    d['gap_e'] = pn - pp - pe
    d['gap_mu'] = pn - pp - pm if pm > 0 else np.nan
    # lower-case aliases used by the prescriptions
    d['ye'], d['ym'], d['yp'] = ye, ym, yp
    return d


class RMFTable:
    """RMF self-energies on an (ln nB, Yq) grid; solves cold npe-mu beta equilibrium."""

    def __init__(self, name, nb, yq, mDn, mDp, Vn, Vp, mn, mp, method='cubic'):
        self.name, self.mn, self.mp = name, mn, mp
        self.yq = yq
        axes = (np.log(nb), yq)
        stack = np.stack([mDn.T, mDp.T, Vn.T, Vp.T], axis=-1)   # (nb, yq, 4)
        self.interp = RegularGridInterpolator(axes, stack, method=method, bounds_error=True)
        self.nb_range = (nb[0], nb[-1])

    def fields(self, nb, yp):
        return self.interp([[np.log(nb), yp]])[0]

    def nucleon(self, nb, yp):
        mdn, mdp, vn, vp = self.fields(nb, yp)

        def f(pn, pp):
            mLn, mLp = np.hypot(pn, mdn), np.hypot(pp, mdp)
            return dict(mnD=mdn, mpD=mdp, Vn=vn, Vp=vp, mLn=mLn, mLp=mLp,
                        mun=mLn + vn, mup=mLp + vp)
        return f

    def beta_mismatch(self, nb, yp):
        mue = mue_from_charge(nb, yp)
        n = self.nucleon(nb, yp)(kf(nb, 1 - yp), kf(nb, yp))
        return n['mun'] - n['mup'] - mue

    def state(self, nb):
        lo, hi = self.yq[0] + 1e-9, min(self.yq[-1], 0.5) - 1e-9
        yp = brentq(lambda y: self.beta_mismatch(nb, y), lo, hi, xtol=1e-13, rtol=1e-13)
        return _finish(nb, yp, mue_from_charge(nb, yp), lambda pn, pp: self.nucleon(nb, yp)(pn, pp),
                       self.name, 'RMF')


def load_dd2(method='cubic'):
    """PT(GRDF2-DD2): micro extract (tools/extract_compose_dd2.py), T index 0 (0.1 MeV); fields (V_n,S_n,V_p,S_p)."""
    z = np.load(DATA / 'DD2_micro_first_two_T.npz')
    f = z['fields'][0]                       # (yq, nb, 4)
    mn, mp = 939.565379, 938.272046          # PT(GRDF2-DD2) table masses
    return RMFTable('DD2', z['nB'], z['Yqgrid'], mn - f[..., 1], mp - f[..., 3], f[..., 0], f[..., 2],
                    mn, mp, method)


def load_bcp(method='cubic'):
    """BCP(NLRMF-N): T index 1 (0.1 MeV) slice extracted by tools/extract_compose_slice.py."""
    z = np.load(DATA / 'BCP_NLRMF_N_T0.npz')
    idx = list(z['micro_idx'])
    if idx != [10041, 11041, 10051, 11051]:
        raise ValueError(f'unexpected micro indices {idx}')
    mn, mp = z['masses']
    m = z['micro']                            # (yq, nb, 4)
    return RMFTable('BCP(NLRMF-N)', z['nB'], z['Yq'], m[..., 0] * mn, m[..., 1] * mp,
                    m[..., 2], m[..., 3], mn, mp, method)


class SkyrmeTable:
    """Cold beta-equilibrium Skyrme table (e.g. BSk22): fractions, mu, m*/m, U.

    Built either from the CompOSE eos.* files (``directory``) or from the
    derived core extract written by ``save`` (``npz``).
    """

    COLUMNS = ['nb', 'Yn', 'Yp', 'Ye', 'Ym', 'mun', 'mup', 'mue', 'rn', 'rp', 'Un', 'Up']

    def __init__(self, name, directory=None, npz=None):
        self.name = name
        if npz is not None:
            z = np.load(npz)
            self.mn, self.mp = z['masses']
            self._build(z['table'])
            return
        d = Path(directory)
        nb = np.loadtxt(d / 'eos.nb', skiprows=2)
        with open(d / 'eos.thermo') as f:
            self.mn, self.mp = map(float, f.readline().split()[:2])
        th = np.loadtxt(d / 'eos.thermo', skiprows=1)
        comp, mic = {}, {}
        for line in open(d / 'eos.compo'):
            x = line.split()
            n = int(x[4]); p = x[5:5 + 2 * n]
            comp[int(x[1])] = {int(p[2 * k]): float(p[2 * k + 1]) for k in range(n)}
        for line in open(d / 'eos.micro'):
            x = line.split()
            n = int(x[3])
            if n:
                mic[int(x[1])] = {int(x[4 + 2 * k]): float(x[5 + 2 * k]) for k in range(n)}
        rows = []
        for j, n in enumerate(nb, start=1):
            c, m = comp[j], mic.get(j)
            homogeneous = set(c) <= {0, 1, 10, 11} and m is not None and c.get(11, 0) > 0
            if not homogeneous:
                continue
            q = th[j - 1]
            mun = self.mn * (1 + q[5]); muq = self.mn * q[6]; mul = self.mn * q[7]
            rows.append([n, c.get(10, 0), c.get(11, 0), c.get(0, 0), c.get(1, 0),
                         mun, mun + muq, mul - muq, m[10040], m[11040], m[10050], m[11050]])
        self._build(np.array(rows))

    def _build(self, a):
        self.table = a
        self.nb_range = (a[0, 0], a[-1, 0])
        self.f = {k: PchipInterpolator(a[:, 0], a[:, i], extrapolate=False)
                  for i, k in enumerate(self.COLUMNS)}

    def save(self, path):
        np.savez_compressed(path, table=self.table, masses=np.array([self.mn, self.mp]),
                            columns=np.array(self.COLUMNS))

    def state(self, nb):
        v = {k: float(f(nb)) for k, f in self.f.items()}
        yp = v['Yp']

        def nucleon(pn, pp):
            return dict(mnD=np.nan, mpD=np.nan, Un=v['Un'], Up=v['Up'],
                        mLn=v['rn'] * self.mn, mLp=v['rp'] * self.mp,
                        mun=v['mun'], mup=v['mup'])
        # Lepton fractions are taken from the table (its mu_e differs from the
        # free-lepton value by ~0.1 MeV); mue is then the free-electron value.
        ye, ym = v['Ye'], v['Ym']
        mue = np.hypot(kf(nb, ye), ME)
        d = _finish(nb, yp, mue, nucleon, self.name, 'Skyrme', leptons=(ye, ym))
        d['table_mue'] = v['mue']
        return d


def load_bsk22():
    """PCP(BSk22) homogeneous-core rows; derived extract if present, else CompOSE files."""
    npz = DATA / 'BSk22_core.npz'
    if npz.exists():
        return SkyrmeTable('BSk22', npz=npz)
    return SkyrmeTable('BSk22', EXTERNAL / 'BSk22')


LOADERS = {'DD2': load_dd2, 'BCP': load_bcp, 'BSk22': load_bsk22}
