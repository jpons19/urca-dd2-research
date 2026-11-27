"""Tests for the multi-EOS Urca comparison (calculations/multi_eos).

Fast checks use only repository data (DD2 extracts, BCP and BSk22 derived
extracts in calculations/multi_eos/data).  Checks that need the downloaded
CompOSE tables in external_eos/ are skipped when those files are absent.
"""
from pathlib import Path
import sys
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from calculations.multi_eos import eos_inputs as E, rates as R  # noqa: E402

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
        cold = np.genfromtxt(ROOT / 'calculations' / 'DD2_cold_composition.csv', delimiter=',', names=True)
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
    def test_dd2_pipeline_reproduces_baseline_csv(self):
        ref = np.genfromtxt(ROOT / 'calculations' / 'PT_GRDF2_DD2_final_T9.csv', delimiter=',', names=True)
        dd2 = E.load_dd2()
        cols = {'FM': 'Q_FM_total', 'Bottaro': 'Q_Bottaro_total',
                'Bottaro_nontri': 'Q_Bottaro_nontriangular_total', 'Shternin': 'Q_Shternin_total',
                'NWA': 'Q_Alford_NWA_adapted_total'}
        for n in (0.16, 0.5, 0.9):
            j = int(np.argmin(abs(ref['nB_fm3'] - n)))
            r = R.all_rates(dd2.state(float(ref['nB_fm3'][j])))
            for model, col in cols.items():
                self.assertLess(abs(sum(r[model]) / ref[col][j] - 1), 2e-4, (n, model))

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


@unittest.skipUnless((EXT / 'BCP_NLRMF_N' / 'eos.thermo.ns').exists(), 'CompOSE BCP table not downloaded')
class TOV(unittest.TestCase):
    def test_bcp_maximum_mass_matches_compose(self):
        from calculations.multi_eos import tov
        t = tov.TOV(*tov.read_thermo(EXT / 'BCP_NLRMF_N', 'eos.thermo.ns', 'eos.nb.ns'))
        s = t.sequence(np.linspace(0.72, 0.86, 8))
        mr = np.genfromtxt(EXT / 'BCP_NLRMF_N' / 'eos.mr', comments='#')
        self.assertLess(abs(s[:, 2].max() - mr[:, 1].max()), 5e-3)


if __name__ == '__main__':
    unittest.main()
