"""Fast EOS-integrity and selected pointwise rate regression tests.

The full 185-density T8 run includes expensive spectral quadrature and belongs
in the slower complete reproduction workflow, not here.
"""
from pathlib import Path
import sys, unittest
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'calculations'))
import compare_compose_four as c
class EOSAndRateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data=np.genfromtxt(ROOT/'calculations'/'PT_GRDF2_DD2_EOS_T8_and_URCA.csv',delimiter=',',names=True)
        cls.cold,cls.cold_i,cls.iv0,_,cls.tgrid=c.setup()
    def test_input_axes_and_temperature(self):
        micro=np.load(ROOT/'calculations'/'DD2_micro_first_two_T.npz')
        self.assertEqual(micro['fields'].shape[-1],4)
        self.assertEqual(micro['fields'].shape[0],2)
        self.assertAlmostEqual(float(micro['Tgrid'][0]),.1,places=7)
        self.assertGreater(float(micro['Tgrid'][1]),.1)
    def test_chemical_and_kinematic_constraints(self):
        x=self.data
        self.assertLess(np.max(abs(x['Yp_cold']-x['Ye_cold']-x['Ymu_cold'])),1e-7)
        self.assertGreater(float(np.min(x['DU_gap_e_MeV'])),0)
        self.assertGreater(float(np.min(x['DU_gap_mu_MeV'])),0)
        for l in ['n','p']:
            md=x[f'mD{l}_T8_extrap_MeV']
            ml=x[f'mL{l}_T8_extrap_MeV']
            pf=x[f'pF{l}_MeV']
            np.testing.assert_allclose(ml, np.hypot(pf,md),rtol=2e-5,atol=.005)
    def test_emissivity_channel_sums(self):
        for model in ['FM','Bottaro','Shternin','Alford_adapted']:
            np.testing.assert_allclose(self.data[f'Q_{model}_total_T8'],self.data[f'Q_{model}_e_T8']+self.data[f'Q_{model}_mu_T8'])
            self.assertTrue(np.all(np.isfinite(self.data[f'Q_{model}_total_T8'])))
    def test_direct_recalculation_at_three_points(self):
        data=self.data
        for nb in [.16,.50,.90]:
            row=data[np.argmin(abs(data['nB_fm3']-nb))]
            d=c.make_data(float(nb),self.cold_i,self.iv0)
            for model,fun in [('FM',c.fm),('Bottaro',c.bottaro),('Shternin',c.shternin),('Alford_adapted',c.alford)]:
                q=sum(fun(d))
                target=float(row[f'Q_{model}_total_T9_recomputed'])
                # The fast runner uses linear interpolation of micro, while the final
                # output uses bicubic; test the physically relevant 1% consistency.
                self.assertLess(abs(q/target-1),.01, (model,nb,q,target))
if __name__=='__main__':unittest.main()
