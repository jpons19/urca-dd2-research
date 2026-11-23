#!/usr/bin/env python3
"""Final CompOSE PT(GRDF2-DD2) Urca comparison and cold-vs-0.1-MeV EOS validation.

Needs extracted NPZ files derived losslessly from original CompOSE archives,
and the cold composition CSV. Reuses documented individual rate formulae in
compare_compose_four.py. A result is NOT an exact reproduction of Alford et al.
"""
import sys, csv
from pathlib import Path
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from compare_compose_four import (setup,make_data,HBARC,MN,MP,ME,MMU,
                                  fm,bottaro,bottaro_nontriangular,shternin,alford)
ROOT=Path(__file__).resolve().parent
mic=np.load(ROOT/'DD2_micro_first_two_T.npz')
com=np.load(ROOT/'DD2_general_composition_first_two_T.npz')
cold,cold_i,old_micro0,old_micro1,tgrid=setup()

# Cubic interpolation in log-density and charge fraction avoids the artificial
# O(10^-5) sawtooth induced by linear interpolation across Yq=0.01 grid steps.
class SmoothFields:
    def __init__(self,nb,yq,fields):
        self.spl=[RectBivariateSpline(np.log(nb),yq,fields[:,:,j].T,kx=3,ky=3) for j in range(fields.shape[-1])]
    def __call__(self,xy):
        xy=np.atleast_2d(xy)
        return np.column_stack([sp.ev(xy[:,0],xy[:,1]) for sp in self.spl])
micro0=SmoothFields(mic['nB'],mic['Yqgrid'],mic['fields'][0])
micro1=SmoothFields(mic['nB'],mic['Yqgrid'],mic['fields'][1])
com0=SmoothFields(com['nB'],com['Yqgrid'],com['frac'][0])
com1=SmoothFields(com['nB'],com['Yqgrid'],com['frac'][1])
T0,T1=map(float,tgrid)
# working range chosen to avoid inhomogeneous and cluster-rich matter.
DEN=np.r_[np.linspace(.16,.98,190),[.16,.2,.3,.5,.7,.9,.98]]
DEN=np.unique(np.round(DEN,9))

def hot_comp(nB, mic_i=micro0, comp_i=com0):
    def res(yp):
        vn,sn,vp,sp=map(float,mic_i([[np.log(nB),yp]])[0])
        pn=HBARC*(3*np.pi**2*nB*(1-yp))**(1/3)
        pp=HBARC*(3*np.pi**2*nB*yp)**(1/3)
        mu_n=np.hypot(pn,MN-sn)+vn
        mu_p=np.hypot(pp,MP-sp)+vp
        yn,ypc,ye,ym=comp_i([[np.log(nB),yp]])[0]
        pe=HBARC*(3*np.pi**2*nB*max(ye,0))**(1/3)
        return mu_n-mu_p-np.hypot(pe,ME)
    yp=float(brentq(res,.011,.27,xtol=1e-12))
    yn,_,ye,ym=map(float,comp_i([[np.log(nB),yp]])[0])
    return dict(Yp=yp,Yn=yn,Ye=ye,Ymu=ym,eq_beta_residual=res(yp),charge=yp-ye-ym)

def massfrom(nB,yp,field_i):
    vn,sn,vp,sp=map(float,field_i([[np.log(nB),yp]])[0])
    pn=HBARC*(3*np.pi**2*nB*(1-yp))**(1/3)
    pp=HBARC*(3*np.pi**2*nB*yp)**(1/3)
    mdn=MN-sn;mdp=MP-sp
    return np.array([mdn,mdp,np.hypot(pn,mdn),np.hypot(pp,mdp)])

# Extrapolation in T^2 of microscopic fields at fixed nB,Yq. Not a tabulated cold EOS micro field.
def mass_at_T0_extrap(nB,yp):
    vec0=np.array(micro0([[np.log(nB),yp]])[0])
    vec1=np.array(micro1([[np.log(nB),yp]])[0])
    extrap=vec0+(0.-T0*T0)/(T1*T1-T0*T0)*(vec1-vec0)
    vn,sn,vp,sp=extrap
    pn=HBARC*(3*np.pi**2*nB*(1-yp))**(1/3)
    pp=HBARC*(3*np.pi**2*nB*yp)**(1/3)
    return np.array([MN-sn,MP-sp,np.hypot(pn,MN-sn),np.hypot(pp,MP-sp)])

rows=[]
for nb in DEN:
    d=make_data(float(nb),cold_i,micro0,micro1)
    hot=hot_comp(float(nb))
    v_cold=np.array([float(cold_i[x](nb)) for x in ('Yp','Ye','Ymu')])
    v_hot=np.array([hot[x] for x in ('Yp','Ye','Ymu')])
    m_hot=massfrom(float(nb),hot['Yp'],micro0)
    m_cold_ex=mass_at_T0_extrap(float(nb),d['yp'])
    fe,fmu=fm(d);be,bmu=bottaro(d);se,smu=shternin(d);ae,amu=alford(d)
    bne, bnmu = bottaro_nontriangular(d)
    rec=dict(nB_fm3=nb,Yp_cold=v_cold[0],Ye_cold=v_cold[1],Ymu_cold=v_cold[2],
        Yp_T01=v_hot[0],Ye_T01=v_hot[1],Ymu_T01=v_hot[2],
        diff_Yp=v_hot[0]-v_cold[0],diff_Ye=v_hot[1]-v_cold[1],diff_Ymu=v_hot[2]-v_cold[2],
        mDn_T01=m_hot[0],mDp_T01=m_hot[1],mLn_T01=m_hot[2],mLp_T01=m_hot[3],
        mDn_cold_extrap=m_cold_ex[0],mDp_cold_extrap=m_cold_ex[1],mLn_cold_extrap=m_cold_ex[2],mLp_cold_extrap=m_cold_ex[3],
        delta_mDn=m_hot[0]-m_cold_ex[0],delta_mDp=m_hot[1]-m_cold_ex[1],
        delta_mLn=m_hot[2]-m_cold_ex[2],delta_mLp=m_hot[3]-m_cold_ex[3],
        beta_equil_residual_cold_MeV=d['diff_mu'],
        beta_equil_residual_T01_MeV=hot['eq_beta_residual'],
        charge_residual_T01=hot['charge'],
        DU_gap_e_MeV=d['gap_e'],DU_gap_mu_MeV=d['gap_mu'],
        Q_FM_e=fe,Q_FM_mu=fmu,Q_FM_total=fe+fmu,
        Q_Bottaro_e=be,Q_Bottaro_mu=bmu,Q_Bottaro_total=be+bmu,
        Q_Bottaro_nontriangular_e=bne,Q_Bottaro_nontriangular_mu=bnmu,Q_Bottaro_nontriangular_total=bne+bnmu,
        Q_Shternin_e=se,Q_Shternin_mu=smu,Q_Shternin_total=se+smu,
        Q_Alford_NWA_adapted_e=ae,Q_Alford_NWA_adapted_mu=amu,Q_Alford_NWA_adapted_total=ae+amu)
    rows.append(rec)
with (ROOT/'PT_GRDF2_DD2_final_T9.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

x=np.array([r['nB_fm3'] for r in rows]);a=lambda k:np.array([r[k] for r in rows])
# cold-vs-general purpose temperatures
fig,ax=plt.subplots(3,1,figsize=(10.4,11.0),sharex=True,gridspec_kw={'height_ratios':[1,0.65,1.25]},constrained_layout=True)
for col,label,color in [('Yp','Protons','tab:blue'),('Ye','Electrons','tab:orange'),('Ymu','Muons','tab:green')]:
    ax[0].plot(x,a(col+'_cold'),lw=2.2,color=color,label=label+' ($T=0$)')
    ax[0].plot(x,a(col+'_T01'),ls='--',lw=1.4,color=color,label=label+' ($T=0.1$ MeV)')
    ax[1].plot(x,1000*a('diff_'+col),lw=1.8,color=color,label='$\\Delta '+col+'$')
ax[0].set_ylabel('Number fraction $Y_i$');ax[0].legend(ncol=2,fontsize=8.4,loc='upper left');ax[1].set_ylabel(r'$10^3\,[Y_i(0.1)-Y_i(0)]$');ax[1].legend(ncol=3,fontsize=8.7)
for j,(col,label,color) in enumerate([('mDn',r'$m_{D,n}^*$','tab:purple'),('mDp',r'$m_{D,p}^*$','tab:red'),('mLn',r'$m_{L,n}^*$','tab:brown'),('mLp',r'$m_{L,p}^*$','tab:cyan')]):
    ax[2].plot(x,a(col+'_cold_extrap'),lw=2.2,color=color,label=label+' ($T\to0$ extrap.)')
    ax[2].plot(x,a(col+'_T01'),lw=1.4,color=color,ls='--',label=label+' ($T=0.1$ MeV)')
ax[2].legend(ncol=2,fontsize=8.0);ax[2].set_ylabel('Nucleon effective mass (MeV)');ax[2].set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
for sub in ax:sub.grid(alpha=.22)
fig.suptitle('CompOSE PT(GRDF2-DD2): cold vs. $T=0.1$ MeV beta equilibrium',fontsize=14)
fig.savefig(ROOT/'PT_GRDF2_DD2_cold_vs_0p1MeV.png',dpi=195)
fig.savefig(ROOT/'PT_GRDF2_DD2_cold_vs_0p1MeV.pdf');plt.close(fig)

# Main Urca comparison on same EOS composition, neutron branch where appropriate
fig,axs=plt.subplots(2,1,sharex=True,figsize=(11,8.1),gridspec_kw={'height_ratios':[2.2,1]},constrained_layout=True)
methods=[('Q_FM_total','Friman–Maxwell / NSCool','tab:blue','--'),('Q_Bottaro_total','Bottaro (triangular OPE+$\\rho$)','tab:orange','-'),('Q_Bottaro_nontriangular_total','Bottaro (nontriangular)','tab:purple','-.'),('Q_Shternin_total','Shternin (propagator enhancement)','tab:green','-'),('Q_Alford_NWA_adapted_total','Alford (adapted finite-width NWA)','tab:red','-')]
for col,label,color,style in methods:
    axs[0].semilogy(x,a(col),lw=2.4,ls=style,color=color,label=label)
    if col!='Q_FM_total':axs[1].semilogy(x,a(col)/a('Q_FM_total'),lw=2,color=color,label=label)
axs[1].axhline(1,lw=1,ls=':',color='0.4')
axs[0].set_ylabel(r'$Q_\nu$ (erg cm$^{-3}$ s$^{-1}$)');axs[1].set_ylabel('Ratio to FM');axs[1].set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
axs[0].legend(fontsize=9,loc='upper left');axs[1].legend(fontsize=8.3,loc='upper left')
for p in axs:p.grid(alpha=.22,which='both')
fig.suptitle(r'CompOSE PT(GRDF2-DD2): five Urca prescriptions, $T=10^9$ K')
fig.savefig(ROOT/'PT_GRDF2_DD2_four_Urca_T9.png',dpi=195)
fig.savefig(ROOT/'PT_GRDF2_DD2_four_Urca_T9.pdf');plt.close(fig)

fig,axs=plt.subplots(1,2,figsize=(12.0,4.9),sharex=True,constrained_layout=True)
for k,lc in [('e','Electron'),('mu','Muon')]:
    sub=axs[0] if k=='e' else axs[1]
    for name,label,color,style in [('FM','FM','tab:blue','--'),('Bottaro','Bottaro (triangular)','tab:orange','-'),('Bottaro_nontriangular','Bottaro (nontriangular)','tab:purple','-.'),('Shternin','Shternin','tab:green','-'),('Alford_NWA_adapted','Adapted Alford NWA','tab:red','-')]:
        y=a('Q_'+name+'_'+k);sub.semilogy(x,np.where(y>0,y,np.nan),color=color,ls=style,lw=2,label=label)
    sub.set_title(lc+' channel');sub.set_xlabel(r'$n_B$ (fm$^{-3}$)');sub.grid(alpha=.23,which='both')
axs[0].set_ylabel(r'$Q_\nu$ (erg cm$^{-3}$ s$^{-1}$)');axs[1].legend(loc='lower right',fontsize=8.5)
fig.suptitle('PT(GRDF2-DD2) individual Urca lepton channels at $T=10^9$ K')
fig.savefig(ROOT/'PT_GRDF2_DD2_Urca_branches_T9.png',dpi=195)
fig.savefig(ROOT/'PT_GRDF2_DD2_Urca_branches_T9.pdf');plt.close(fig)

print('ROWS',len(rows),'MIN DU gaps',a('DU_gap_e_MeV').min(),a('DU_gap_mu_MeV').min())
for col in ['diff_Yp','diff_Ye','diff_Ymu','delta_mDn','delta_mDp','delta_mLn','delta_mLp','beta_equil_residual_cold_MeV','charge_residual_T01']:
    vals=a(col); print(col,'max abs',max(abs(vals)),'at',x[np.argmax(abs(vals))])
for nb in [.16,.2,.3,.5,.7,.9,.98]:
    idx=int(np.argmin(abs(x-nb)));r=rows[idx]
    print('row',nb,'Yp cold/hot',r['Yp_cold'],r['Yp_T01'],'Landau n,p',r['mLn_T01'],r['mLp_T01'],
        'Q',*(f'{r[k]:.4e}' for k in ['Q_FM_total','Q_Bottaro_total','Q_Bottaro_nontriangular_total','Q_Shternin_total','Q_Alford_NWA_adapted_total']))
