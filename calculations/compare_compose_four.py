#!/usr/bin/env python3
"""Compare five Urca rate prescriptions with actual PT(GRDF2-DD2) CompOSE inputs.

This script is deliberately reproducible without original 182 MB EOS archive:
 - DD2_micro_first_two_T.npz: 4 micro fields extracted from eos.micro at first TWO T slices,
 - DD2_cold_composition.csv: full cold 1D table composition.

Model qualifications:
FM neutron MURCA: legacy NSCool prescription (includes old muon approximation).
Bottaro: Eq 3.28 triangular, fixed Landau phenomenological couplings, beta=0.68.
Shternin: empirical n-MU enhancement 0.6 RA+0.2 RB, below DU threshold only,
    derived for different nuclear interaction; not a correction to Bottaro.
Alford: modified Breit-Wigner NWA mass-convolution of leading degenerate DU
    rate, NOT full relativistic finite-T Alford calculation.

Microscopic self-energies taken at T=0.100 MeV, while emissivities are at
T=1e9 K=0.086173 MeV. The EOS fields vary very little at this temperature.
"""
from pathlib import Path
import csv
import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator
from scipy.integrate import quad
from bottaro_angular import angular_ratio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
HBARC=197.3269804; MN=939.565379; MP=938.272046
ME=.51099895; MMU=105.6583755
GF=1.1663787e-11; VUD=.974; GA=1.275; G=GF*VUD
MPI=139.57; MRHO=769.; BETA=.68; CRHO=1.4; FPI=1.
T9=1.; T=.08617333262145; HALF_WIDTH=(T*T/5.0)/2.0
NATURAL_TO_CGS=1.602176634e-6/((1.973269804e-11)**3*6.582119569e-22)


def setup():
    cold=np.genfromtxt(HERE/'DD2_cold_composition.csv', delimiter=',',names=True)
    cold_i={x:PchipInterpolator(cold['nB_fm3'],cold[x],extrapolate=False) for x in ('Yn','Yp','Ye','Ymu')}
    mic=np.load(HERE/'DD2_micro_first_two_T.npz')
    # axes stored as (temperature, yq, nB, fields); transpose to (log nB,yq,fields)
    temp0=mic['fields'][0].transpose(1,0,2)
    temp1=mic['fields'][1].transpose(1,0,2)
    iv0=RegularGridInterpolator((np.log(mic['nB']),mic['Yqgrid']),temp0,bounds_error=True)
    iv1=RegularGridInterpolator((np.log(mic['nB']),mic['Yqgrid']),temp1,bounds_error=True)
    return cold,cold_i,iv0,iv1,mic['Tgrid']


def make_data(nb, cold_i, micro, micro_t1=None):
    yn,yp,ye,ym=(float(cold_i[k](nb)) for k in ('Yn','Yp','Ye','Ymu'))
    pn,pp,pe,pm=(HBARC*(3*np.pi**2*nb*y)**(1/3) if y>0 else 0 for y in (yn,yp,ye,ym))
    vn,sn,vp,sp=map(float,micro([[np.log(nb),yp]])[0]); mdn=MN-sn; mdp=MP-sp
    if min(mdn,mdp)<0:raise ValueError('negative Dirac mass at n='+str(nb))
    mLn=np.hypot(pn,mdn); mLp=np.hypot(pp,mdp)
    e_mu=np.hypot(pe,ME); mu_mu=np.hypot(pm,MMU) if pm>0 else np.nan
    data=dict(nb=nb,yp=yp,ye=ye,ym=ym,pn=pn,pp=pp,pe=pe,pm=pm,mu=e_mu,
              mnD=mdn,mpD=mdp,mLn=mLn,mLp=mLp,Vn=vn,Vp=vp,Sn=sn,Sp=sp,
              mun=mLn+vn,mup=mLp+vp,diff_mu=(mLn+vn)-(mLp+vp)-e_mu,
              e_mu_minus_mu=(mu_mu-e_mu if pm else np.nan),charge=yp-ye-ym,
              gap_e=pn-pp-pe,gap_mu=pn-pp-pm)
    if micro_t1 is not None:
        v2=np.asarray(micro_t1([[np.log(nb),yp]])[0]);data['micro_Tchange_max']=float(max(abs(v2-np.array([vn,sn,vp,sp]))))
    return data


def fm(d):
    kfn=d['pn']/HBARC
    alpha=1.76-.63*(1.68/kfn)**2
    mass=(d['mLn']/MN)**3*(d['mLp']/MP)
    coeff=8.55e21*mass*alpha*BETA*T9**8/1.68
    return coeff*d['pe']/HBARC, coeff*d['pm']/HBARC


def bottaro_landau(pn):
    # Preserve the bare-mass normalization of the published couplings.
    scale = np.pi**2 / (2 * MN * pn)
    return .7 * scale, 1.1 * scale


def bottaro(d):
    pn,pp=d['pn'],d['pp']
    fp,g=bottaro_landau(pn)
    e_pi=pn**2/(pn**2+MPI**2);e_rho=pn**2/(pn**2+MRHO**2)
    alpha=.5*(21/16)*((32/7)*(MPI/FPI)**4*(fp-g)**2+2*(e_pi-CRHO*e_rho)**2)
    coeff=2*11513*G**2*GA**2*FPI**4/(120960*np.pi*MPI**4)
    coeff*=d['mLn']**3*d['mLp']*pp*T**8*alpha*BETA*NATURAL_TO_CGS
    return coeff,coeff*d['pm']/d['mu'] if d['pm'] else 0.


def bottaro_meson_couplings(q2, g):
    rho = CRHO * q2 / (q2 + MRHO**2)
    pion = q2 / (q2 + MPI**2)
    scale = (FPI / MPI)**2
    # The constant isovector spin coupling g-prime equals g in this model.
    return g - scale * rho, -scale * (pion - rho)


def bottaro_nontriangular(d):
    """Return electron and muon emissivities with Bottaro's full angular integral."""
    pn, pp = d['pn'], d['pp']
    fp, g = bottaro_landau(pn)
    qe, qmu = bottaro(d)
    electron = angular_ratio(pn, pp, d['pe'], fp, g, GA, bottaro_meson_couplings)
    muon = angular_ratio(pn, pp, d['pm'], fp, g, GA, bottaro_meson_couplings)
    # Restore the electron velocity approximated by unity in bottaro().
    return qe * (d['pe'] / d['mu']) * electron, qmu * muon


def RA(d,pl):
    pn,pp=d['pn'],d['pp'];a=pp;lo=pn-pl;hi=pn+pl
    if lo<=a:return np.nan
    # Analytic derivative primitive integral of 1/(k^2-a^2)^2
    def antider(x):return -x/(2*a*a*(x*x-a*a))-np.log((x-a)/(x+a))/(4*a**3)
    return 2*d['mLp']**2*d['mu']**2/pl*(antider(hi)-antider(lo))


def RB(d,pl):
    pn,pp=d['pn'],d['pp']
    if pn<=pp+pl:return np.nan
    return 4*d['mLn']**2*d['mu']**2/((pn**2-(pp+pl)**2)*(pn**2-(pp-pl)**2))


def shternin(d):
    out=[]
    for pl in (d['pe'],d['pm']):
        if pl<=0:out.append(0.);continue
        if d['pn']-d['pp']-pl<=1:out.append(np.nan);continue
        enh=.6*RA(d,pl)+.2*RB(d,pl)
        coeff=8.1e21*(d['mLn']/MN)**3*(d['mLp']/MP)*(pl/d['mu'])*(d['nb']*d['yp']/.16)**(1/3)*T9**8
        out.append(coeff*enh)
    return tuple(out)


def _cauchy_interval(a,b,center,half):
    if b<=a:return 0.
    # arctan difference is safe for peak integrals of narrow BW
    return max(0.,(np.arctan((b-center)/half)-np.arctan((a-center)/half))/np.pi)


def nwa_probability(d,pl,epsrel=5e-4):
    """2-dimensional BW mass convolution reduced analytically to one quadrature.

    Kinetic Fermi energies E_F,* fixed, and off-shell momentum for sampled
    Dirac mass m: k_F(m) = sqrt(E_F,*^2-m^2). The all-angle triangle step
    is retained in the leading degenerate direct Urca kernel.
    """
    if pl<=0:return 0.
    en=d['mLn'];ep=d['mLp'];mn0=d['mnD'];mp0=d['mpD'];w=HALF_WIDTH
    ulo=np.arctan((-en-mn0)/w);uhi=np.arctan((en-mn0)/w)
    def integrand(u):
        mn=mn0+w*np.tan(u)
        kn=np.sqrt(max(0.,en**2-mn**2))
        kpmin=abs(kn-pl);kpmax=min(ep,kn+pl)
        if kpmin>=kpmax:return 0.
        mlo=np.sqrt(max(0.,ep**2-kpmax**2))
        mhi=np.sqrt(max(0.,ep**2-kpmin**2))
        probp=_cauchy_interval(mlo,mhi,mp0,w)+_cauchy_interval(-mhi,-mlo,mp0,w)
        return probp/np.pi
    pts=[ulo,uhi,0.]
    for kt in [abs(d['pp']-pl),d['pp']+pl,abs(ep-pl),ep+pl]:
        if 0<kt<en:
            root=np.sqrt(en*en-kt*kt)
            for s in (-1,1):
                u=np.arctan((s*root-mn0)/w)
                if ulo<u<uhi:pts.append(u)
    pts=sorted(set(float(u) for u in pts if ulo<=u<=uhi))
    return sum(quad(integrand,a,b,epsabs=4e-10,epsrel=epsrel,limit=150)[0] for a,b in zip(pts[:-1],pts[1:]) if b>a)


def alford(d,epsrel=5e-4):
    # The leading degenerate DU kernel is a deliberate approximation to
    # Alford's full finite-T relativistic response, retaining the BW structure.
    pref=(457*np.pi/10080)*G**2*(1+3*GA**2)*d['mLn']*d['mLp']*d['mu']*T**6*NATURAL_TO_CGS
    return pref*nwa_probability(d,d['pe'],epsrel),pref*nwa_probability(d,d['pm'],epsrel)


def run():
    cold,cold_i,iv0,iv1,tgrid=setup()
    ns=np.r_[np.linspace(.10,.40,80),np.linspace(.405,.99,140)]
    rows=[]
    for nb in ns:
        d=make_data(float(nb),cold_i,iv0,iv1)
        fe,fm_=fm(d);be,bm=bottaro(d);se,sm=shternin(d);ae,am=alford(d)
        bne, bnm = bottaro_nontriangular(d)
        row=dict(nB=nb,Yp=d['yp'],Ye=d['ye'],Ymu=d['ym'],pFn_MeV=d['pn'],pFp_MeV=d['pp'],pFe_MeV=d['pe'],pFmu_MeV=d['pm'],
                 mu_e_MeV=d['mu'],mDirac_n_MeV=d['mnD'],mDirac_p_MeV=d['mpD'],mLandau_n_MeV=d['mLn'],mLandau_p_MeV=d['mLp'],
                 V_n_MeV=d['Vn'],V_p_MeV=d['Vp'],S_n_MeV=d['Sn'],S_p_MeV=d['Sp'],
                 mu_n_MeV=d['mun'],mu_p_MeV=d['mup'],mu_beta_residual_MeV=d['diff_mu'],muon_beta_residual_MeV=d['e_mu_minus_mu'],
                 charge_residual=d['charge'],micro_temperature_sensitivity_MeV=d['micro_Tchange_max'],
                 DU_gap_e_MeV=d['gap_e'],DU_gap_mu_MeV=d['gap_mu'],
                 FM_e=fe,FM_mu=fm_,FM_total=fe+fm_,Bottaro_e=be,Bottaro_mu=bm,Bottaro_total=be+bm,
                 Bottaro_nontriangular_e=bne,Bottaro_nontriangular_mu=bnm,Bottaro_nontriangular_total=bne+bnm,
                 Shternin_e=se,Shternin_mu=sm,Shternin_total=se+sm,Alford_adapted_e=ae,Alford_adapted_mu=am,Alford_adapted_total=ae+am)
        rows.append(row)
    out=HERE/'PT_GRDF2_DD2_Urca_T9.csv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    x=np.array([r['nB'] for r in rows]); reference=np.array([r['FM_total'] for r in rows]);
    names=[('FM_total','Friman–Maxwell / NSCool','--'),('Bottaro_total','Bottaro (triangular, pion + rho)','-'),('Bottaro_nontriangular_total','Bottaro (nontriangular)','-.'),('Shternin_total','Shternin (propagator enhancement)','-'),('Alford_adapted_total','Alford (adapted finite-width NWA)','-')]
    fig,(ax,ar)=plt.subplots(2,1,sharex=True,figsize=(10.8,8.5),gridspec_kw={'height_ratios':[2.3,1],'hspace':.08})
    for key,name,style in names:
        y=np.array([r[key] for r in rows]);obj=ax.semilogy(x,y,ls=style,lw=2.25,label=name)[0]
        if key!='FM_total':ar.semilogy(x,y/reference,ls=style,lw=2.0,color=obj.get_color(),label=name)
    ar.axhline(1,ls=':',lw=1,color='0.3');ar.set_yscale('log')
    ax.set_ylabel(r'$Q_\nu$ (erg cm$^{-3}$ s$^{-1}$)');ar.set_ylabel('Ratio to FM');ar.set_xlabel(r'Baryon density $n_B$ (fm$^{-3}$)')
    ax.set_title(r'PT(GRDF2-DD2) CompOSE: Urca emissivities at $T=10^9$ K')
    ax.legend(fontsize=9,loc='upper left');ar.legend(fontsize=8,loc='upper left')
    for a in [ax,ar]:a.grid(alpha=.23,which='both')
    fig.text(.10,.014,'Tabulated cold composition + general-purpose microscopic fields at T=0.1 MeV; no DU threshold below 1 fm$^{-3}$.\nAll rates are approximate: n-MU only for FM/Bottaro/Shternin; adapted, degenerate DU kernel for NWA.',fontsize=8.5,color='.31')
    fig.subplots_adjust(bottom=.14,top=.92,left=.11,right=.975)
    fig.savefig(HERE/'PT_GRDF2_DD2_Urca_T9.png',dpi=190)
    fig.savefig(HERE/'PT_GRDF2_DD2_Urca_T9.pdf')
    plt.close(fig)
    fig,axes=plt.subplots(3,1,sharex=True,figsize=(10,10))
    for field,label in [('Yp','protons'),('Ye','electrons'),('Ymu','muons')]:
        axes[0].plot(x,[r[field] for r in rows],label=label,lw=2)
    axes[0].legend();axes[0].set_ylabel('Particle fraction')
    for field,label in [('mDirac_n_MeV',r'$m_{D,n}$'),('mDirac_p_MeV',r'$m_{D,p}$'),('mLandau_n_MeV',r'$m_{L,n}$'),('mLandau_p_MeV',r'$m_{L,p}$')]:
        axes[1].plot(x,[r[field] for r in rows],label=label,lw=2)
    axes[1].legend(ncol=2);axes[1].set_ylabel('Effective mass (MeV)')
    axes[2].plot(x,[r['DU_gap_e_MeV'] for r in rows],label=r'$\delta p_e$',lw=2)
    axes[2].plot(x,[r['DU_gap_mu_MeV'] for r in rows],label=r'$\delta p_\mu$',lw=2)
    axes[2].axhline(0,color='0.4',ls=':');axes[2].legend();axes[2].set_ylabel(r'$p_{F,n}-p_{F,p}-p_{F,\ell}$ (MeV)');axes[2].set_xlabel(r'$n_B$ (fm$^{-3}$)')
    for a in axes:a.grid(alpha=.22)
    fig.suptitle('CompOSE PT(GRDF2-DD2): cold composition and tabulated microscopic fields')
    fig.tight_layout();fig.savefig(HERE/'PT_GRDF2_DD2_EOS_checks.png',dpi=190)
    fig.savefig(HERE/'PT_GRDF2_DD2_EOS_checks.pdf')
    plt.close(fig)
    print('rows',len(rows),'temperature slices',tgrid,flush=True)
    print('Max |charge residual| core',max(abs(r['charge_residual']) for r in rows))
    print('Max |beta eq residual| MeV',max(abs(r['mu_beta_residual_MeV']) for r in rows))
    print('Max |muon beta eq residual| MeV',np.nanmax(abs(np.array([r['muon_beta_residual_MeV'] for r in rows]))))
    print('Micro difference at T0,T1 (max) MeV',max(r['micro_temperature_sensitivity_MeV'] for r in rows))
    print('Minimum electron/muon DU gaps (MeV)',min(r['DU_gap_e_MeV'] for r in rows),min(r['DU_gap_mu_MeV'] for r in rows))
    for den in [.16,.30,.50,.70,.90]:
        r=min(rows,key=lambda r:abs(r['nB']-den))
        print(f"{r['nB']:.3f}: Yp={r['Yp']:.5f} mDn={r['mDirac_n_MeV']:.1f} mLn={r['mLandau_n_MeV']:.1f} mLp={r['mLandau_p_MeV']:.1f}; Q_FM={r['FM_total']:.3e}, Q_B={r['Bottaro_total']:.3e}, Q_Sh={r['Shternin_total']:.3e}, Q_NWA={r['Alford_adapted_total']:.3e}",flush=True)
    return rows

if __name__=='__main__':run()
