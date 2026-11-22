#!/usr/bin/env python3
"""PT(GRDF2-DD2) reference composition/masses and adapted Urca rates T=1e8 K.
Source grids: user-supplied CompOSE cold and general-purpose archives, reduced in
previous session to .csv and first-two-temperature .npz extracts.
Extrapolation below T=0.1 MeV is linear in T**2 at fixed nB,Yq. Not interpolation.
"""
import csv,sys
from pathlib import Path
import numpy as np
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
from scipy.optimize import brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import compare_compose_four as rate
micro=np.load(ROOT/'DD2_micro_first_two_T.npz')
comp=np.load(ROOT/'DD2_general_composition_first_two_T.npz')
Tlo,Thi=map(float,micro['Tgrid']);TK_1e8=1e8*8.617333262145e-11 # MeV
TK_1e9=1e9*8.617333262145e-11
assert TK_1e8<TK_1e9<Tlo<Thi
cold=np.genfromtxt(ROOT/'DD2_cold_composition.csv',delimiter=',',names=True)
cold_i={k:PchipInterpolator(cold['nB_fm3'],cold[k],extrapolate=False) for k in ('Yn','Yp','Ye','Ymu')}
class SmoothGrid:
    def __init__(self,nb,yq,dat):
        self.spl=[RectBivariateSpline(np.log(nb),yq,dat[:,:,i].T,kx=3,ky=3) for i in range(dat.shape[-1])]
    def __call__(self,xy):
        xy=np.atleast_2d(xy)
        return np.column_stack([s.ev(xy[:,0],xy[:,1]) for s in self.spl])
# Extracted axes (T,Yq,nB,field)
f_m0=SmoothGrid(micro['nB'],micro['Yqgrid'],micro['fields'][0])
f_m1=SmoothGrid(micro['nB'],micro['Yqgrid'],micro['fields'][1])
f_c0=SmoothGrid(comp['nB'],comp['Yqgrid'],comp['frac'][0])
f_c1=SmoothGrid(comp['nB'],comp['Yqgrid'],comp['frac'][1])
class T2Extrap:
    def __init__(self,g0,g1,Ttarget):
        self.g0=g0;self.g1=g1;self.t=Ttarget
        self.a=(Ttarget*Ttarget-Tlo*Tlo)/(Thi*Thi-Tlo*Tlo)
    def __call__(self,pts):return (1-self.a)*self.g0(pts)+self.a*self.g1(pts)
m8=T2Extrap(f_m0,f_m1,TK_1e8);c8=T2Extrap(f_c0,f_c1,TK_1e8)
m9=T2Extrap(f_m0,f_m1,TK_1e9)
mzero=T2Extrap(f_m0,f_m1,0)

def fermi(nb,y):return rate.HBARC*np.cbrt(3*np.pi**2*nb*max(0,y))

def solution_hot(nb,field,frac):
    def constraint(yp):
        vn,sn,vp,sp=field([[np.log(nb),yp]])[0]
        pn=fermi(nb,1-yp);pp=fermi(nb,yp)
        _,_,ye,ym=frac([[np.log(nb),yp]])[0]
        pe=fermi(nb,ye)
        return np.hypot(pn,rate.MN-sn)+vn -np.hypot(pp,rate.MP-sp)-vp -np.hypot(pe,rate.ME)
    yp=float(brentq(constraint,.011,.27,xtol=1e-12))
    yn,ypc,ye,ym=map(float,frac([[np.log(nb),yp]])[0])
    return dict(Yn=yn,Yp=yp,Ye=ye,Ymu=ym,beta_res=constraint(yp),charge_res=yp-ye-ym)

def mass_summary(nb,yp,field):
    vn,sn,vp,sp=map(float,field([[np.log(nb),yp]])[0]);pn=fermi(nb,1-yp);pp=fermi(nb,yp)
    mn=rate.MN-sn;mp=rate.MP-sp
    return dict(mDn=mn,mDp=mp,mLn=np.hypot(pn,mn),mLp=np.hypot(pp,mp),Vn=vn,Vp=vp,Sn=sn,Sp=sp)

# smooth rate grid, and separate support grids to plot mass+composition
DENS=np.unique(np.r_[np.linspace(.16,.98,180),[.16,.20,.30,.50,.70,.90,.98]])
rows=[]
rate.T=TK_1e8;rate.T9=TK_1e8/TK_1e9;rate.HALF_WIDTH=rate.T**2/10 # W/2, W=T^2/5 MeV
for nb in DENS:
    yp0=float(cold_i['Yp'](nb));yc={k:float(cold_i[k](nb)) for k in ('Yn','Yp','Ye','Ymu')}
    y01=solution_hot(nb,f_m0,f_c0)
    y8=solution_hot(nb,m8,c8)
    mm01=mass_summary(nb,y01['Yp'],f_m0)
    mm8=mass_summary(nb,y8['Yp'],m8)
    mm0=mass_summary(nb,yp0,mzero)
    # Preserve cold beta-equilibrium composition as baseline for rate curves.
    # The microscopic fields are T=10^8K T^2 extrapolated and at cold charge fraction.
    d=rate.make_data(float(nb),cold_i,m8)
    fe,fmu=rate.fm(d);be,bmu=rate.bottaro(d);se,smu=rate.shternin(d);ae,amu=rate.alford(d)
    row={'nB_fm3':nb,'Yp_cold':yc['Yp'],'Ye_cold':yc['Ye'],'Ymu_cold':yc['Ymu'],
         'Yp_T01':y01['Yp'],'Ye_T01':y01['Ye'],'Ymu_T01':y01['Ymu'],
         'Yp_T8_extrap':y8['Yp'],'Ye_T8_extrap':y8['Ye'],'Ymu_T8_extrap':y8['Ymu'],
         'charge_residual_T8':y8['charge_res'],'beta_residual_T8_MeV':y8['beta_res'],
         'beta_residual_cold_T8fields_MeV':d['diff_mu'],
         'pFn_MeV':d['pn'],'pFp_MeV':d['pp'],'pFe_MeV':d['pe'],'pFmu_MeV':d['pm'],
         'DU_gap_e_MeV':d['gap_e'],'DU_gap_mu_MeV':d['gap_mu']}
    for k in mm8:
        row[f'{k}_T8_extrap_MeV']=mm8[k]
        row[f'{k}_T01_MeV']=mm01[k]
        row[f'{k}_T0_extrap_MeV']=mm0[k]
    for meth,pair in [('FM',(fe,fmu)),('Bottaro',(be,bmu)),('Shternin',(se,smu)),('Alford_adapted',(ae,amu))]:
        row[f'Q_{meth}_e_T8']=pair[0];row[f'Q_{meth}_mu_T8']=pair[1];row[f'Q_{meth}_total_T8']=sum(pair)
    # Recompute 10^9 K on precisely the same density grid to validate T^8 scaling,
    # avoiding irregularities from interpolation of earlier T9 output grids.
    rate.T=TK_1e9;rate.T9=1.;rate.HALF_WIDTH=rate.T**2/10
    d9=rate.make_data(float(nb),cold_i,m9)
    for meth,pair in [('FM',rate.fm(d9)),('Bottaro',rate.bottaro(d9)),('Shternin',rate.shternin(d9)),('Alford_adapted',rate.alford(d9))]:
        row[f'Q_{meth}_total_T9_recomputed']=sum(pair)
    rate.T=TK_1e8;rate.T9=.1;rate.HALF_WIDTH=rate.T**2/10
    rows.append(row)
out_csv=ROOT/'PT_GRDF2_DD2_EOS_T8_and_URCA.csv'
with out_csv.open('w',newline='') as f:
    dw=csv.DictWriter(f,fieldnames=list(rows[0].keys()));dw.writeheader();dw.writerows(rows)
x=np.array([r['nB_fm3'] for r in rows]);get=lambda k:np.array([r[k] for r in rows])

# 1. Composition reference, cold and both finite T values.
fig,ax=plt.subplots(2,1,figsize=(10.7,8.2),sharex=True,gridspec_kw={'height_ratios':[1.9,1]},constrained_layout=True)
for k,label in [('Yp','Protons'),('Ye','Electrons'),('Ymu','Muons')]:
    c=ax[0].plot(x,get(k+'_cold'),lw=2.4,label=f'{label}, cold')[0].get_color()
    ax[0].plot(x,get(k+'_T01'),lw=1.4,ls=':',color=c,alpha=.9)
    ax[0].plot(x,get(k+'_T8_extrap'),lw=1.5,ls='--',color=c,alpha=.85)
    ax[1].plot(x,1e6*(get(k+'_T01')-get(k+'_cold')),ls=':',color=c,lw=1.7,label=f'{k}: 0.1 MeV minus cold')
    ax[1].plot(x,1e6*(get(k+'_T8_extrap')-get(k+'_cold')),ls='--',color=c,lw=1.7,label=f'{k}: 10⁸ K minus cold')
ax[0].set_ylabel('Particle fraction $Y_i$')
ax[0].legend(loc='upper left',ncol=2,fontsize=9)
ax[0].text(.97,.96,'solid: cold  |  dotted: 0.1 MeV  |  dashed: $10^8$ K',transform=ax[0].transAxes,ha='right',va='top',fontsize=8.6,bbox=dict(facecolor='white',edgecolor='.8',alpha=.92,pad=4))
ax[1].set_ylabel(r'$10^6\,[Y_i(T)-Y_i(0)]$')
ax[1].set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
ax[1].legend(ncol=2,fontsize=7.8,loc='best')
for a in ax:a.grid(alpha=.22)
fig.suptitle('PT(GRDF2-DD2) beta-equilibrated core composition: cold, 0.1 MeV and $10^8$ K',fontsize=13)

fig.savefig(ROOT/'DD2_composition_cold_T01_T8.png',dpi=185)
fig.savefig(ROOT/'DD2_composition_cold_T01_T8.pdf');plt.close(fig)

# 2. Masses + their tiny temperature shifts; m_D=m-S and m_L^2=p_F^2 + m_D^2.
fig,ax=plt.subplots(2,1,figsize=(10.7,8.2),sharex=True,gridspec_kw={'height_ratios':[1.8,1]},constrained_layout=True)
for k,label in [('mDn',r'$m^*_{D,n}$'),('mDp',r'$m^*_{D,p}$'),('mLn',r'$m^*_{L,n}$'),('mLp',r'$m^*_{L,p}$')]:
    c=ax[0].plot(x,get(k+'_T8_extrap_MeV')/rate.MN,lw=2.4,label=label+' ($10^8$ K)')[0].get_color()
    ax[0].plot(x,get(k+'_T01_MeV')/rate.MN,ls=':',lw=1.7,color=c)
    ax[1].plot(x,1e3*(get(k+'_T01_MeV')-get(k+'_T8_extrap_MeV')),color=c,lw=2,label=label)
ax[0].set_ylabel(r'Effective mass / $m_N$')
ax[0].legend(fontsize=9,ncol=2)
ax[0].text(.97,.94,'solid: $10^8$ K  |  dotted: 0.1 MeV',transform=ax[0].transAxes,ha='right',va='top',fontsize=8.6,bbox=dict(facecolor='white',edgecolor='.8',alpha=.92,pad=4))
ax[1].set_ylabel(r'$m^*(0.1\,\mathrm{MeV})-m^*(10^8\,\mathrm{K})$ (keV)')
ax[1].set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
ax[1].legend(ncol=2,fontsize=9)
for a in ax:a.grid(alpha=.22)
fig.suptitle('PT(GRDF2-DD2): Dirac and Landau effective masses',fontsize=13)

fig.savefig(ROOT/'DD2_effective_masses_T01_T8.png',dpi=185)
fig.savefig(ROOT/'DD2_effective_masses_T01_T8.pdf');plt.close(fig)

# 3. Scalar/vector self energies from micro fields
fig,ax=plt.subplots(figsize=(10.7,5.9),constrained_layout=True)
for k,label in [('Sn',r'$S_n$'),('Sp',r'$S_p$'),('Vn',r'$V_n$'),('Vp',r'$V_p$')]:
    c=ax.plot(x,get(k+'_T8_extrap_MeV'),lw=2.4,label=label+' ($10^8$ K)')[0].get_color()
    ax.plot(x,get(k+'_T01_MeV'),ls=':',lw=1.5,color=c)
ax.set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
ax.set_ylabel('Self-energy (MeV)')
ax.set_title('PT(GRDF2-DD2): CompOSE scalar and vector self-energies')
ax.legend(ncol=2);ax.grid(alpha=.22)

fig.savefig(ROOT/'DD2_scalar_vector_selfenergies_T8.png',dpi=185)
fig.savefig(ROOT/'DD2_scalar_vector_selfenergies_T8.pdf');plt.close(fig)

# 4. Four Urca emissivities evaluated at T=1e8 K
fig,axs=plt.subplots(2,1,sharex=True,figsize=(10.6,8.2),gridspec_kw={'height_ratios':[2.0,1]},constrained_layout=True)
keys=[('FM','Friman-Maxwell / NSCool','--'),('Bottaro','Bottaro (triangular OPE + rho)','-'),('Shternin','Shternin (propagator enhancement)','-'),('Alford_adapted','Alford (adapted finite-width NWA)','-')]
ref=get('Q_FM_total_T8')
for meth,label,ls in keys:
    y=get(f'Q_{meth}_total_T8')
    c=axs[0].semilogy(x,y,label=label,ls=ls,lw=2.3)[0].get_color()
    if meth!='FM':axs[1].semilogy(x,y/ref,label=label,color=c,lw=1.9)
axs[0].set_ylabel(r'Neutrino emissivity $Q_\nu$ (erg cm$^{-3}$ s$^{-1}$)')
axs[1].set_ylabel('Ratio to FM');axs[1].set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
axs[0].legend(loc='upper left',fontsize=9);axs[1].legend(loc='upper left',fontsize=8.3)
for a in axs:a.grid(alpha=.23,which='both')
fig.suptitle('PT(GRDF2-DD2): four Urca approximations at $T=10^8$ K',fontsize=13)

fig.savefig(ROOT/'DD2_four_URCA_T8.png',dpi=185)
fig.savefig(ROOT/'DD2_four_URCA_T8.pdf');plt.close(fig)

# 5. Check expected T8 scaling relative to the previously published T9 density scan.
fig,ax=plt.subplots(figsize=(10.6,5.9),constrained_layout=True)
for meth,label,ls in keys:
    q9=get(f'Q_{meth}_total_T9_recomputed')
    q8=get(f'Q_{meth}_total_T8')
    ratio=q8/(1e-8*q9)
    ax.plot(x,ratio,ls=ls,lw=2,label=label)
ax.axhline(1,ls=':',lw=1,color='0.4')
ax.set_xlabel(r'$n_B$ (fm$^{-3}$)');ax.set_ylabel(r'$Q(10^8\,\mathrm{K})/[10^{-8}\,Q(10^9\,\mathrm{K})]$')
ax.set_title('Temperature-scaling consistency: expected $T^8$ behavior')
ax.grid(alpha=.22);ax.legend(fontsize=8.8)
fig.savefig(ROOT/'DD2_URCA_T8_scaling_check.png',dpi=185)
fig.savefig(ROOT/'DD2_URCA_T8_scaling_check.pdf');plt.close(fig)

print('Temperature grid:',Tlo,Thi,'T8:',TK_1e8,'T9:',TK_1e9,'t2 extrap weight',m8.a)
print('density points:',len(rows))
for k in ['Yp','Ye','Ymu']:
    z=get(k+'_T8_extrap')-get(k+'_cold')
    print('max |T8-cold|',k,float(np.max(abs(z))))
for k in ['mDn','mDp','mLn','mLp']:
    z=get(k+'_T01_MeV')-get(k+'_T8_extrap_MeV')
    print('max mass |0.1MeV - T8|',k,float(np.max(abs(z))),'MeV')
print('beta eq at T8 max residual',float(np.max(abs(get('beta_residual_T8_MeV')))))
print('charge at T8 max residual',float(np.max(abs(get('charge_residual_T8')))))
print('min DU gaps',float(min(get('DU_gap_e_MeV'))),float(min(get('DU_gap_mu_MeV'))))
for nb in [.16,.3,.5,.7,.9]:
    r=min(rows,key=lambda k:abs(k['nB_fm3']-nb))
    print('nB %.2f Yp %.7f mDn %.3f mDp %.3f mLn %.3f mLp %.3f | Q8'%(
        nb,r['Yp_T8_extrap'],r['mDn_T8_extrap_MeV'],r['mDp_T8_extrap_MeV'],r['mLn_T8_extrap_MeV'],r['mLp_T8_extrap_MeV']),
        *['%.3e'%r[f'Q_{m}_total_T8'] for m in ['FM','Bottaro','Shternin','Alford_adapted']])
print('All outputs created in',ROOT)
