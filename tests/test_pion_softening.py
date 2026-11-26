"""Focused physics and numerical regression tests for fixed-coupling screening."""
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'calculations'))
import compare_compose_four as base
from pion_softening.bottaro_softened import softened_rates,triangular_factor
from pion_softening.propagators import pion_ratio,couplings


class SoftPionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _,cls.cold,cls.micro,*_=base.setup()
        cls.samples=qmc.Sobol(d=4,scramble=True,seed=42).random_base2(12)

    def test_vacuum_couplings_and_screening(self):
        for q2 in (0.,base.MPI**2,400.**2):
            np.testing.assert_allclose(couplings(q2,1.e-5,1.,base),base.bottaro_meson_couplings(q2,1.e-5),rtol=0,atol=1.e-17)
        self.assertAlmostEqual(pion_ratio(500.**2,base.MPI,1.),500.**2/(500.**2+base.MPI**2))
        with self.assertRaises(ValueError):pion_ratio(100,base.MPI,0)

    def test_reference_rates_both_channels(self):
        for nb in (.16,.50,.90):
            d=base.make_data(nb,self.cold,self.micro)
            rates=softened_rates(d,1.,base,self.samples)
            e,mu=base.bottaro(d)
            self.assertAlmostEqual(rates[0]/e,1.,places=12)
            self.assertAlmostEqual(rates[1]/mu,1.,places=12)
            # Exact same angular routine and sample points as the old baseline.
            from bottaro_angular import angular_ratio
            fp,g=base.bottaro_landau(d['pn'])
            ae=angular_ratio(d['pn'],d['pp'],d['pe'],fp,g,base.GA,base.bottaro_meson_couplings,samples=self.samples)
            am=angular_ratio(d['pn'],d['pp'],d['pm'],fp,g,base.GA,base.bottaro_meson_couplings,samples=self.samples)
            np.testing.assert_allclose((rates[2],rates[3]),(e*d['pe']/d['mu']*ae,mu*am),rtol=1.e-12)

    def test_fixed_q_ope_and_triangular(self):
        for ratio in (.75,.5,.25):
            q2=base.MPI**2
            self.assertAlmostEqual(pion_ratio(q2,base.MPI,ratio)/pion_ratio(q2,base.MPI,1.),2/(1+ratio**2))
            d=base.make_data(.5,self.cold,self.micro)
            self.assertGreater(triangular_factor(d,ratio,base),1.)

    def test_nontriangular_convergence_and_positive(self):
        d=base.make_data(.5,self.cold,self.micro)
        x=softened_rates(d,.5,base,self.samples)
        finer=qmc.Sobol(d=4,scramble=True,seed=42).random_base2(14)
        y=softened_rates(d,.5,base,finer)
        self.assertTrue(np.all(np.array(x)>0))
        self.assertTrue(np.all(np.array(y)>0))
        np.testing.assert_allclose(x[2:],y[2:],rtol=.05)

if __name__=='__main__':unittest.main()
