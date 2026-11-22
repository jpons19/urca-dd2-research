"""Consistency regression tests for the generated PT(GRDF2-DD2) T8 results."""
from pathlib import Path
import numpy as np
D=np.genfromtxt(Path(__file__).resolve().parent/'PT_GRDF2_DD2_EOS_T8_and_URCA.csv',delimiter=',',names=True)
assert len(D)>=180
for k in ['Yp','Ye','Ymu']:
    assert max(abs(D[f'{k}_T8_extrap']-D[f'{k}_cold']))<2e-6
for k in ['mDn','mDp','mLn','mLp']:
    assert np.all(D[f'{k}_T8_extrap_MeV']>0)
    assert max(abs(D[f'{k}_T01_MeV']-D[f'{k}_T8_extrap_MeV']))<.001
assert max(abs(D['charge_residual_T8']))<5e-8
assert max(abs(D['beta_residual_T8_MeV']))<1e-8
assert min(D['DU_gap_e_MeV'])>0 and min(D['DU_gap_mu_MeV'])>0
for method in ['FM','Bottaro','Shternin','Alford_adapted']:
    y=D[f'Q_{method}_total_T8']; assert np.all(y>0) and np.all(np.isfinite(y))
    assert np.allclose(y,D[f'Q_{method}_e_T8']+D[f'Q_{method}_mu_T8'])
    y9=D[f'Q_{method}_total_T9_recomputed']
    assert max(abs(y/(1e-8*y9)-1))<1e-4
print('PASS: compositions, masses, beta equilibrium, neutrality, DU gaps, positivity, sums and T^8 scaling')
