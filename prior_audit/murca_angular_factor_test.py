"""Illustrative Fermi-surface angular test of Eq (3.31) at one chosen pn,pp,pe.
This uses quasiparticle Fermi momenta and a Monte-Carlo version of the
angular delta function; no energy integral or in-medium spectral functions.
"""
import numpy as np
from murca_spin_isospin_audit import eq330,eq331,eq331_printed

rng=np.random.default_rng(42)
pn=340.;pp=100.;pe=100.
mpi=139.57;mrho=769.;Crho=1.4
norm=np.pi**2/(2*939*pn)
fp=0.7*norm;g=1.1*norm
A0=1/mpi**2

def g_h(vec):
  q2=np.dot(vec,vec); ur=q2/(q2+mrho*mrho); up=q2/(q2+mpi*mpi)
  return 1.1*norm-A0*Crho*ur,-A0*(up-Crho*ur)

def vecs():
  while True:
    p=rng.normal(size=(3,3));p*=pn/np.linalg.norm(p,axis=1)[:,None]
    p1,p2,p3=p
    a=p1+p2-p3; qa=np.linalg.norm(a)
    if not (abs(pp-pe)<qa<pp+pe):continue
    ah=a/qa; b=rng.normal(size=3);b-=np.dot(b,ah)*ah
    b/=np.linalg.norm(b)
    costh=(qa*qa+pp*pp-pe*pe)/(2*qa*pp)
    p4=pp*(costh*ah+np.sqrt(1-costh**2)*b)
    k=p1-p3;kt=p2-p4;l=p2-p3;lt=p1-p4
    return qa,k,kt,l,lt

vals=[]
for j in range(1500):
    qa,k,kt,l,lt=vecs()
    gk,hk=g_h(k);gkt,hkt=g_h(kt);gl,hl=g_h(l);glt,hlt=g_h(lt)
    a=(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt)
    base=eq330(fp,g,gk,gkt,hk,hkt,k,kt)+eq330(fp,g,gl,glt,hl,hlt,l,lt)
    correct=base+eq331(*a)
    printed=base+eq331_printed(*a)
    vals.append((correct,printed,1/qa, printed/correct))
vals=np.array(vals)
correct=np.sum(vals[:,0]*vals[:,2]);printed=np.sum(vals[:,1]*vals[:,2])
print('Momenta pn, pp, pe:',pn,pp,pe,'MeV; n accepted',len(vals))
print('Equal-phase-space ratio printed/corrected =',printed/correct)
print('Range of local ratios in sample',np.min(vals[:,3]), np.max(vals[:,3]))
print('No cooling rate computed; uncertainty from MC sampling not yet estimated.')
blocks=np.array_split(vals,6)
ratios=np.array([np.sum(b[:,1]*b[:,2])/np.sum(b[:,0]*b[:,2]) for b in blocks])
print('Approx. block-mean statistical error',np.std(ratios,ddof=1)/np.sqrt(len(ratios)))
