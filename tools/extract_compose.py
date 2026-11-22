#!/usr/bin/env python3
"""Extract exactly the compact first-two-T CompOSE PT(GRDF2-DD2) inputs.

Usage: python tools/extract_compose.py --cold cold.zip --general general.zip --out calculations/
This reads the original public-format EOS archives (large general-purpose file).
"""
import argparse, csv, zipfile
from pathlib import Path
import numpy as np


def grid(z,name):
    x=z.read(name).split()
    n=int(x[1]);vals=np.array([float(y) for y in x[2:]],dtype=float)
    if len(vals)!=n: raise ValueError(f'{name}: expected {n}, got {len(vals)}')
    return vals


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--cold',required=True,type=Path)
    a.add_argument('--general',required=True,type=Path)
    a.add_argument('--out',required=True,type=Path)
    args=a.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(args.cold) as z:
        nb=grid(z,'eos.nb'); frac=np.zeros((len(nb),4),dtype=float)
        with z.open('eos.compo') as f:
            for line in f:
                xs=line.split()
                if not xs:continue
                j=int(xs[1])-1;count=int(xs[4]); pairs=xs[5:5+2*count]
                sp={int(pairs[2*k]):float(pairs[2*k+1]) for k in range(count)}
                frac[j,:]=[sp.get(k,0.) for k in [10,11,0,1]]
        # In the outer crust the composition contains nuclei and other charged species.
        # The simple npeμ relation is only appropriate for the homogeneous core.
        core=nb>=0.16
        if np.max(np.abs(frac[core,1]-frac[core,2]-frac[core,3]))>1e-5:
            raise ValueError('core charge neutrality unexpectedly violated')
        hbarc=197.3269804
        p=hbarc*np.cbrt(3*np.pi**2*nb[:,None]*frac)
        dat=np.column_stack([nb,frac,p,p[:,0]-p[:,1]-p[:,2],p[:,0]-p[:,1]-p[:,3]])
        header='nB_fm3,Yn,Yp,Ye,Ymu,pFn_MeV,pFp_MeV,pFe_MeV,pFmu_MeV,DU_gap_e_MeV,DU_gap_mu_MeV'
        np.savetxt(args.out/'DD2_cold_composition.csv',dat,delimiter=',',header=header,comments='')
    with zipfile.ZipFile(args.general) as z:
        nb=grid(z,'eos.nb');tg=grid(z,'eos.t');yq=grid(z,'eos.yq');
        micro=np.full((2,len(yq),len(nb),4),np.nan,dtype=float)
        comp=np.full_like(micro,np.nan)
        with z.open('eos.micro') as f:
            for line in f:
                xs=line.split();it,ib,iy,k=map(int,xs[:4])
                if it>2:break
                indices=list(map(int,xs[4::2]))
                if indices!=[10051,10052,11051,11052]:raise ValueError(f'unexpected eos.micro fields {indices}')
                micro[it-1,iy-1,ib-1]=np.array(xs[5::2],dtype=float)
        with z.open('eos.compo') as f:
            for line in f:
                xs=line.split();it,ib,iy=map(int,xs[:3])
                if it>2:break
                n=int(xs[4]);sp={int(xs[5+2*k]):float(xs[6+2*k]) for k in range(n)}
                comp[it-1,iy-1,ib-1]=[sp.get(k,0.) for k in [10,11,0,1]]
        if not np.all(np.isfinite(micro)) or not np.all(np.isfinite(comp)):
            raise ValueError('Missing general-purpose micro/composition entries')
        np.savez_compressed(args.out/'DD2_micro_first_two_T.npz',nB=nb,Tgrid=tg[:2],Yqgrid=yq,fields=micro)
        np.savez_compressed(args.out/'DD2_general_composition_first_two_T.npz',nB=nb,Yqgrid=yq,frac=comp)
    print('Extracted cold composition and two lowest general-purpose T slices into',args.out)

if __name__=='__main__':main()
