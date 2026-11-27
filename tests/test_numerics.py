"""Tests for the multi-EOS Urca comparison (numerics/).

Fast checks use only repository data (derived DD2, BCP and BSk22 extracts in
numerics/data).  Checks that need the downloaded CompOSE tables in external_eos/
are skipped when those files are absent.
"""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from numerics import eos_inputs as E, prescriptions as P, rates as R  # noqa: E402

EXT = ROOT / 'external_eos'


class RMFInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dd2, cls.bcp = E.load_dd2(), E.load_bcp()

    def test_bcp_dispersion_reproduces_table_chemical_potentials(self):
        z = np.load(E.DATA / 'BCP_NLRMF_N_T0.npz')
        mn, mp = z['masses']
        for ib in (120, 150, 170, 190):
            for iy in (5, 15, 30):
                n, yp = z['nB'][ib], z['Yp'][iy, ib]
                mdn, mdp, vn, vp = z['micro'][iy, ib] * [mn, mp, 1, 1]
                pn, pp = E.kf(n, 1 - yp), E.kf(n, yp)
                self.assertAlmostEqual(np.hypot(pn, mdn) + vn, z["mu_n"][iy, ib], delta=1e-3)  # T=0.1 MeV thermal shift ~2e-4 MeV
                self.assertAlmostEqual(np.hypot(pp, mdp) + vp, z["mu_p"][iy, ib], delta=1e-3)

    def test_dd2_beta_solver_reproduces_cold_table(self):
        cold = np.genfromtxt(E.DATA / 'DD2_cold_composition.csv', delimiter=',', names=True)
        for n in (0.16, 0.3, 0.5, 0.8, 0.98):
            s = self.dd2.state(n)
            for key, col in (('Yp', 'Yp'), ('Ye', 'Ye'), ('Ym', 'Ymu')):
                self.assertLess(abs(s[key] - np.interp(n, cold['nB_fm3'], cold[col])), 1e-5)
            self.assertLess(abs(s['beta_residual']), 1e-8)
            self.assertLess(abs(s['charge_residual']), 1e-12)

    def test_bcp_electron_du_threshold_and_no_muon_du(self):
        self.assertGreater(self.bcp.state(0.50)['gap_e'], 0)
        self.assertLess(self.bcp.state(0.52)['gap_e'], 0)
        for n in (0.6, 0.9, 1.1):
            self.assertGreater(self.bcp.state(n)['gap_mu'], 0)


class SkyrmeInputs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.t = E.load_bsk22()

    def test_electron_and_muon_du_thresholds(self):
        # Electron DU of BSk22 at ~0.33 fm^-3 (M_DU ~1.15 Msun, Pearson et al. 2018)
        self.assertGreater(self.t.state(0.330)['gap_e'], 0)
        self.assertLess(self.t.state(0.335)['gap_e'], 0)
        self.assertGreater(self.t.state(0.385)['gap_mu'], 0)
        self.assertLess(self.t.state(0.395)['gap_mu'], 0)

    def test_neutron_dispersion_and_neutrality(self):
        for n in np.linspace(0.1, 1.4, 14):
            s = self.t.state(float(n))
            eps_n = self.t.mn + s['pn']**2 / (2 * s['mLn']) + s['Un']
            self.assertLess(abs(s['mun'] - eps_n), 2e-2)
            self.assertLess(abs(s['charge_residual']), 2e-5)

    def test_known_proton_table_inconsistency_is_bounded(self):
        """Known CompOSE BSk22 issue: m_p + p_F^2/2m*_p + U_p differs from mu_p by up to ~35 MeV.

        Rates do not use U_p or mu_p; this test documents the size of the issue.
        """
        worst = max(abs(s['mup'] - (self.t.mp + s['pp']**2 / (2 * s['mLp']) + s['Up']))
                    for s in (self.t.state(float(n)) for n in np.linspace(0.15, 1.4, 60)))
        self.assertLess(worst, 40.0)


class Rates(unittest.TestCase):
    # Frozen DD2 reference (former single-EOS pipeline: cold-table composition, T=0.1 MeV micro
    # slice; erg cm^-3 s^-1, e + mu).  The multi-EOS pipeline solves beta equilibrium itself, so
    # agreement at the 1e-4 level is expected.
    DD2_REFERENCE = {
        0.16: dict(FM=5.776261206242971e+20, Bottaro=3.0490457881966215e+20,
                   Bottaro_nontri=3.290465356803914e+20, Shternin=2.693035414915336e+21,
                   NWA=9.705046575909138e+21),
        0.5: dict(FM=5.7812438911314625e+20, Bottaro=1.875307834372843e+20,
                  Bottaro_nontri=2.6200633660670483e+20, Shternin=2.1536042939329547e+21,
                  NWA=7.157824207448695e+21),
        0.9: dict(FM=1.0983198317455404e+21, Bottaro=2.391296297342951e+20,
                  Bottaro_nontri=4.445803648147834e+20, Shternin=7.100786129493361e+21,
                  NWA=8.919138256441666e+21),
    }

    def test_dd2_reproduces_frozen_reference(self):
        dd2 = E.load_dd2()
        for n, ref in self.DD2_REFERENCE.items():
            r = R.all_rates(dd2.state(n))
            for model, value in ref.items():
                self.assertLess(abs(sum(r[model]) / value - 1), 2e-4, (n, model))

    def test_nwa_tends_to_du_above_threshold_and_rules(self):
        bcp, bsk = E.load_bcp(), E.load_bsk22()
        s = bcp.state(0.8)
        r = R.all_rates(s)
        self.assertGreater(r['DU'][0], 0)
        self.assertLess(abs(r['NWA'][0] / r['DU'][0] - 1), 1e-2)
        self.assertEqual(r['DU'][1], 0.0)                       # muon DU closed in BCP
        self.assertTrue(np.isnan(r['Shternin'][0]))             # undefined above threshold
        rb = R.all_rates(bsk.state(0.6))
        self.assertTrue(np.all(np.isnan(rb['NWA'])))            # no Dirac mass for Skyrme
        self.assertGreater(rb['DU'][0], 0)
        self.assertGreater(rb['DU'][1], 0)

    def test_bottaro_beta_variants(self):
        r = R.all_rates(E.load_dd2().state(0.3))
        for m in ('Bottaro', 'Bottaro_nontri'):
            for i in (0, 1):
                self.assertAlmostEqual(r[m + '_noSRC'][i] * 0.68, r[m][i], delta=1e-9 * r[m][i] + 1e-300)

    def test_nwa_below_threshold_is_continuous_and_below_du(self):
        bcp = E.load_bcp()
        vals = [sum(R.all_rates(bcp.state(n))['NWA']) for n in np.linspace(0.40, 0.50, 6)]
        self.assertTrue(np.all(np.diff(vals) > 0))
        du = sum(R.du(bcp.state(0.52)))
        self.assertLess(vals[-1], du)


class Sedrakian(unittest.TestCase):
    """Sedrakian (2024) proton-spectral-width emissivity, numerics/sedrakian.py."""

    @classmethod
    def setUpClass(cls):
        from numerics import sedrakian as S
        cls.S = S
        cls.bcp, cls.bsk, cls.dd2 = E.load_bcp(), E.load_bsk22(), E.load_dd2()

    def test_zero_width_recovers_sharp_kinematics(self):
        d = self.dd2.state(0.3)                      # gap ~ 49 MeV >> T
        self.assertLess(self.S.ratio(d, d['pe'], width_scale=0.0), 1e-20)
        d = self.bsk.state(0.5)                      # far above threshold
        self.assertAlmostEqual(self.S.ratio(d, d['pe'], width_scale=0.0), 1.0, delta=1e-4)

    def test_finite_width_tends_to_du_above_threshold(self):
        for d in (self.bcp.state(0.8), self.bsk.state(0.6)):
            qe, _ = self.S.sedrakian(d)
            self.assertAlmostEqual(qe / R.du(d)[0], 1.0, delta=1e-3)

    def test_subthreshold_rate_is_linear_in_width(self):
        d = self.dd2.state(0.3)
        r1 = self.S.ratio(d, d['pe'])
        r2 = self.S.ratio(d, d['pe'], width_scale=2.0)
        self.assertAlmostEqual(r2 / r1, 2.0, delta=2e-3)

    def test_width_coefficient_by_hand(self):
        d = self.dd2.state(0.16)
        mun, mup = d['pn']**2 / (2 * d['mLn']), d['pp']**2 / (2 * d['mLp'])
        a = E.HBARC * d['pp'] / d['mLp'] / mup * (7 * 0.16 * d['Yn'] / mun + 4 * 0.16 * d['Yp'] / mup)
        self.assertAlmostEqual(self.S.width_coefficient(d) / a, 1.0, places=12)

    def test_grid_convergence(self):
        S = self.S
        d = self.bcp.state(0.5)                      # 0.5 MeV from threshold: hardest case
        base_val = S.ratio(d, d['pe'])
        saved = S._X, S._Y, S._XNU
        try:
            S._X, S._Y, S._XNU = (np.linspace(-60, 60, 2401), np.linspace(-90, 90, 3601),
                                  np.linspace(0, 80, 1601))
            fine = S.ratio(d, d['pe'])
        finally:
            S._X, S._Y, S._XNU = saved
        self.assertLess(abs(base_val / fine - 1), 2e-4)


class MMU(unittest.TestCase):
    """Voskresensky MMU first-approximation parametrization, numerics/mmu.py."""

    @classmethod
    def setUpClass(cls):
        from numerics import mmu as M
        cls.M, cls.dd2 = M, E.load_dd2()

    def test_inputs(self):
        M = self.M
        self.assertAlmostEqual(M.gamma(1.0), 1 / 2.6, places=12)
        self.assertEqual(M.gap2(0.5), 1.0)
        self.assertAlmostEqual(M.gap2(3.0), 0.27, places=6)
        self.assertAlmostEqual(M.gap2(5.0), 0.27, places=6)                    # saturation (curve 1b)
        self.assertAlmostEqual(M.gap2(2.0, nc=2.0), 0.216, places=6)
        g = [M.gap2(x) for x in np.linspace(0.8, 3.0, 23)]
        self.assertTrue(np.all(np.diff(g) <= 1e-12))                           # monotone softening

    def test_ratio_form_by_hand(self):
        d = self.dd2.state(0.32)
        x = 2.0
        expected = sum(P.fm(d)) * 3 * x**(10 / 3) * ((1 + 1.6) / (1 + 1.6 * x**(1 / 3)))**6 / 0.734**4
        self.assertAlmostEqual(self.M.mmu_ratio_form(d) / expected, 1.0, delta=2e-3)

    def test_low_density_reduces_to_mu(self):
        d = self.dd2.state(0.12)
        self.assertEqual(self.M.mmu_ratio_form(d), sum(P.fm(d)))
        self.assertTrue(np.isnan(self.M.mmu_ratio_form(d, unsoftened=True)))
        self.assertTrue(np.isnan(self.M.mmu_direct_form(d)))

    def test_forms_agree_to_order_one(self):
        for n in (0.16, 0.32, 0.48):
            d = self.dd2.state(n)
            r = self.M.mmu_direct_form(d) / self.M.mmu_ratio_form(d)
            self.assertTrue(0.3 < r < 3.0, (n, r))

    def test_unsoftened_gap_suppresses_pion_term(self):
        d = self.dd2.state(0.16)
        self.assertLess(self.M.mmu_ratio_form(d, unsoftened=True), 1e-2 * sum(P.fm(d)))


@unittest.skipUnless((EXT / 'BCP_NLRMF_N' / 'eos.thermo.ns').exists(), 'CompOSE BCP table not downloaded')
class TOV(unittest.TestCase):
    def test_bcp_maximum_mass_matches_compose(self):
        from numerics import tov
        t = tov.TOV(*tov.read_thermo(EXT / 'BCP_NLRMF_N', 'eos.thermo.ns', 'eos.nb.ns'))
        s = t.sequence(np.linspace(0.72, 0.86, 8))
        mr = np.genfromtxt(EXT / 'BCP_NLRMF_N' / 'eos.mr', comments='#')
        self.assertLess(abs(s[:, 2].max() - mr[:, 1].max()), 5e-3)


if __name__ == '__main__':
    unittest.main()


class BottaroTrace(unittest.TestCase):
    def test_spin_isospin_audit_runs(self):
        """Independent trace check behind the corrected Eq. (3.31) interference."""
        import subprocess
        out = subprocess.run([sys.executable, str(ROOT / 'numerics' / 'validation_bottaro_spin_isospin.py')],
                             capture_output=True, text=True, check=True).stdout
        self.assertIn('corrected Eq (3.31): 0.0', out)
