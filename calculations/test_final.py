"""Quick unit/consistency tests for CompOSE DD2 T9 comparison output."""
from pathlib import Path
import numpy as np
P=Path(__file__).resolve().parent
D=np.genfromtxt(P/'PT_GRDF2_DD2_final_T9.csv',delimiter=',',names=True)
assert len(D)>150
for k in ['Q_FM_total','Q_Bottaro_total','Q_Shternin_total','Q_Alford_NWA_adapted_total']:
    assert np.all(np.isfinite(D[k])) and np.all(D[k]>0),(k, np.min(D[k]))
for meth in ['FM','Bottaro','Shternin','Alford_NWA_adapted']:
    assert np.allclose(D[f'Q_{meth}_total'],D[f'Q_{meth}_e']+D[f'Q_{meth}_mu'],atol=1e5)
for k in ['Yp','Ye','Ymu']:
    assert max(abs(D[f'diff_{k}']))<2e-6
assert max(abs(D['beta_equil_residual_cold_MeV']))<.002
assert max(abs(D['beta_equil_residual_T01_MeV']))<1e-8
assert max(abs(D['charge_residual_T01']))<1e-7
assert min(D['DU_gap_e_MeV'])>0 and min(D['DU_gap_mu_MeV'])>0
for k in ['mDn_T01','mDp_T01','mLn_T01','mLp_T01']:
    assert np.all(D[k]>0)
assert np.all(D['mLn_T01']>D['mDn_T01'])
assert np.all(D['mLp_T01']>D['mDp_T01'])
for k in ['delta_mDn','delta_mDp','delta_mLn','delta_mLp']:
    assert max(abs(D[k]))<.001
print('PASS: 195 densities, finite emissivities, branch sums, charge and beta equilibrium, masses, composition, DU gaps')
