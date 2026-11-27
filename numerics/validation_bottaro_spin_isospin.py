#!/usr/bin/env python3
"""Independent spin/isospin audit of Bottaro, Caputo & Fiorillo (2024).

Reconstructs neutron-branch modified-Urca t/u soft-propagator amplitudes
in a 16-dimensional two-nucleon spin-isospin Hilbert space. Checks:
 - Appendix A spin sums and the triangular Eq. (3.22);
 - complete t-channel Eq. (3.30), including vector current;
 - general t/u interference Eq. (3.31), detecting an excess factor two
   in its quadratic axial interference terms as printed in arXiv:2406.18640v2.

The program checks the matrix-element algebra only; it does not integrate
phase space, infer in-medium amplitudes, or reproduce the paper's figures.
Requires: NumPy.
"""


import numpy as np


I2=np.eye(2,dtype=complex)


sx=np.array([[0,1],[1,0]],complex)


sy=np.array([[0,-1j],[1j,0]],complex)


sz=np.diag([1,-1]).astype(complex)


sig=[sx,sy,sz]


tau=[sx,sy,sz]


tminus=np.array([[0,1],[0,0]],complex)


I4=np.eye(4,dtype=complex)


sp=[np.kron(s,I2) for s in sig]


tp=[np.kron(I2,t) for t in tau]


W=[np.kron(sp[i]@np.kron(I2,tminus),I4)+np.kron(I4,sp[i]@np.kron(I2,tminus)) for i in range(3)]


V_f=np.eye(16,dtype=complex)


V_fp=sum(np.kron(t,t) for t in tp)


V_g=sum(np.kron(s,s) for s in sp)


V_gp=sum(np.kron(sp[i]@tp[j],sp[i]@tp[j]) for i in range(3) for j in range(3))


def V_h(vec):
    v=np.asarray(vec,float)
    v=v/np.linalg.norm(v)
    return sum(v[i]*v[j]*np.kron(sp[i]@tp[k],sp[j]@tp[k]) for i in range(3) for j in range(3) for k in range(3))


def idx(spin,isospin):
    return 2*spin+isospin


def matx(V):
    # complex shape [weak-index i, spin3, spin4, spin1, spin2], t-channel and exchanged final u-channel
    C=[V@w-w@V for w in W]
    Ct=np.empty((3,2,2,2,2),complex)
    Cu=np.empty_like(Ct)
    for i in range(3):
      for s3 in range(2):
       for s4 in range(2):
        for s1 in range(2):
         for s2 in range(2):
          ket=np.ravel_multi_index((idx(s1,1),idx(s2,1)),(4,4))
          bra_t=np.ravel_multi_index((idx(s3,1),idx(s4,0)),(4,4))
          bra_u=np.ravel_multi_index((idx(s4,0),idx(s3,1)),(4,4))
          Ct[i,s3,s4,s1,s2]=C[i][bra_t,ket]
          Cu[i,s3,s4,s1,s2]=-C[i][bra_u,ket]
    return Ct,Cu


def sums(A,B):
    return np.vdot(A,B).real


def report(name,V):
  A,B=matx(V)
  print(name,'t',sums(A,A),'u',sums(B,B),'int',2*sums(A,B),'total',sums(A+B,A+B),'|m|² coeff',sums(A+B,A+B)/16)
  return A,B


W1=[np.kron(sp[i]@np.kron(I2,tminus),I4) for i in range(3)]


W2=[np.kron(I4,sp[i]@np.kron(I2,tminus)) for i in range(3)]


Q1=np.kron(np.kron(I2,tminus),I4)


Q2=np.kron(I4,np.kron(I2,tminus))


tz=tp[2]


tx=tp[0]


ty=tp[1]


Niso=np.kron(tz,tz)


Ciso=np.kron(tx,tx)+np.kron(ty,ty)


SpinSpin=sum(np.kron(s,s) for s in sp)


SpinTensor=lambda v: sum(v[i]*v[j]*np.kron(sp[i],sp[j]) for i in range(3) for j in range(3))


def V_N(f=0,fp=0,g=0,gp=0,h=0,v=(0,0,1)):
    v=np.array(v,float); v/=np.linalg.norm(v)
    return f*V_f+fp*Niso+g*V_g+gp*sum(np.kron(sp[i]@tz,sp[i]@tz) for i in range(3)) + h*sum(v[i]*v[j]*np.kron(sp[i]@tz,sp[j]@tz) for i in range(3) for j in range(3))


def V_C(fp=0,gp=0,h=0,v=(0,0,1)):
    v=np.array(v,float); v/=np.linalg.norm(v)
    return fp*Ciso + gp*sum(np.kron(sp[i]@tp[j],sp[i]@tp[j]) for i in range(3) for j in range(2)) + h*sum(v[i]*v[j]*np.kron(sp[i]@tp[a],sp[j]@tp[a]) for i in range(3) for j in range(3) for a in range(2))


def arr_from(Cs,t=True):
    arr=np.empty((len(Cs),2,2,2,2),complex)
    for i,C in enumerate(Cs):
      for s3 in range(2):
       for s4 in range(2):
        for s1 in range(2):
         for s2 in range(2):
          ket=np.ravel_multi_index((idx(s1,1),idx(s2,1)),(4,4))
          bra=(idx(s3,1),idx(s4,0)) if t else (idx(s4,0),idx(s3,1))
          bra=np.ravel_multi_index(bra,(4,4))
          arr[i,s3,s4,s1,s2]=C[bra,ket]* (1 if t else -1)
    return arr


def full_momentum(f=0,fp=0,g=0,gt=0,gkt=0,gl=0,glt=0,ht=0,hkt=0,hl=0,hlt=0,
                  k=(0,0,1),kt=(0,0,1),l=(0,0,1),lt=(0,0,1),which='axial'):
    vn_t=V_N(f,fp,g,gt,ht,k); vc_t=V_C(fp,gkt,hkt,kt)
    vn_u=V_N(f,fp,g,gl,hl,l); vc_u=V_C(fp,glt,hlt,lt)
    if which=='axial':
      T=[vn_t@w2-w2@vn_t +vc_t@w1-w1@vc_t for w1,w2 in zip(W1,W2)]
      U=[vn_u@w1-w1@vn_u +vc_u@w2-w2@vc_u for w1,w2 in zip(W1,W2)]
    elif which=='vector':
      T=[vn_t@Q2-Q2@vn_t +vc_t@Q1-Q1@vc_t]
      U=[vn_u@Q1-Q1@vn_u +vc_u@Q2-Q2@vc_u]
    return arr_from(T,True),arr_from(U,False)


def eq330(fp,g,gk,gkt,hk,hkt,k,kt,gA=1):
    d=np.dot(np.asarray(k)/np.linalg.norm(k),np.asarray(kt)/np.linalg.norm(kt))
    axial=(6*(fp*fp+g*g)-2*(3*fp*gk+3*(fp+2*g)*gkt+fp*hk+(fp+2*g)*hkt)
           +3*gk*gk+9*gkt*gkt+2*gk*hk+hk*hk+6*gkt*hkt+3*hkt*hkt)
    vec=3*(gk-gkt)**2+2*(gk-gkt)*(hk-hkt)+(hk-hkt)**2+2*hk*hkt*(1-d*d)
    return gA*gA*axial+vec


def eq331(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt,gA=1):
    """Corrected Eq. (3.31): quadratic axial terms are not multiplied by 2."""
    def dsq(a,b):
      a=np.asarray(a)/np.linalg.norm(a);b=np.asarray(b)/np.linalg.norm(b)
      return 1-float(np.dot(a,b))**2
    angular = (hk*hl*dsq(k,l)-2*hkt*hl*dsq(kt,l)-2*hk*hlt*dsq(k,lt)+2*hkt*hlt*dsq(kt,lt))
    axial = (-24*fp*g+2*(
        2*fp*(3*(gkt+glt)+hkt+hlt)
        +3*g*(gk+gkt+gl+glt)+g*(hk+hkt+hl+hlt)
        -0.5*(gl*(hk+hkt)+glt*(hk+5*hkt)+hk*hl-hkt*hl-hk*hlt+3*hkt*hlt
        +gk*(3*gl+3*glt+hl+hlt)+gkt*(3*gl+15*glt+hl+5*hlt))
        +0.5*angular))
    vector=( (gl-glt)*(hk-hkt)+(gk-gkt)*(hl-hlt)-(hk-hkt)*(hl-hlt)
            +3*(gl-glt)*(gk-gkt)+2*hk*hl*dsq(k,l)-2*hkt*hl*dsq(kt,l)-2*hk*hlt*dsq(k,lt)+2*hkt*hlt*dsq(kt,lt) )
    return gA*gA*axial+vector



def eq331_printed(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt,gA=1):
    """Eq (3.31) as printed; the -Q+A terms are doubled inside 2(...)."""
    def s2(a,b):
        a=np.asarray(a,float); b=np.asarray(b,float)
        return 1-(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))**2
    Q=(gl*(hk+hkt)+glt*(hk+5*hkt)+hk*hl-hkt*hl-hk*hlt+3*hkt*hlt
       +gk*(3*gl+3*glt+hl+hlt)+gkt*(3*gl+15*glt+hl+5*hlt))
    A=(hk*hl*s2(k,l)-2*hkt*hl*s2(kt,l)-2*hk*hlt*s2(k,lt)
       +2*hkt*hlt*s2(kt,lt))
    return eq331(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt,gA) + gA*gA*(-Q+A)


def run_audit():
    pz=[0,0,1]
    def reduced(V):
        t,u=matx(V)
        return sums(t+u,t+u)/16
    assert abs(reduced(V_f)) < 1e-12
    assert abs(reduced(V_fp)-12)<1e-12
    assert abs(reduced(V_g)-12)<1e-12
    assert abs(reduced(V_gp))<1e-12
    # cos(k,l) = -1/2. In the triangular limit k_tilde=-k, l_tilde=-l
    # and each tensor operator is insensitive to reversing its momentum.
    pl=[np.sqrt(3)/2,0,-.5]
    t,u=full_momentum(ht=1,hkt=1,hl=1,hlt=1,k=pz,kt=[0,0,-1],l=pl,lt=-np.array(pl))
    assert abs(sums(t+u,t+u)/16-21/4)<1e-12
    # Pure constant g' gives ZERO after t/u antisymmetrization,
    # but literal Eq. (3.31) would give a NEGATIVE squared amplitude.
    a=(0,0,1,1,1,1,0,0,0,0,pz,pz,pz,pz)
    correct=eq330(0,0,1,1,0,0,pz,pz)*2+eq331(*a)
    printed=eq330(0,0,1,1,0,0,pz,pz)*2+eq331_printed(*a)
    assert abs(correct)<1e-12 and abs(printed+24)<1e-12
    rng=np.random.default_rng(2024)
    max_err_t=0;max_err_i=0
    for _ in range(30):
        (fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt)=rng.normal(size=10)
        k,kt,l,lt=rng.normal(size=(4,3))
        T,U=full_momentum(fp=fp,g=g,gt=gk,gkt=gkt,gl=gl,glt=glt,ht=hk,hkt=hkt,hl=hl,hlt=hlt,k=k,kt=kt,l=l,lt=lt)
        Tv,Uv=full_momentum(fp=fp,g=g,gt=gk,gkt=gkt,gl=gl,glt=glt,ht=hk,hkt=hkt,hl=hl,hlt=hlt,k=k,kt=kt,l=l,lt=lt,which='vector')
        direct=(sums(T,T)+sums(Tv,Tv))/16
        interf=(2*sums(T,U)+2*sums(Tv,Uv))/16
        analytic_t=eq330(fp,g,gk,gkt,hk,hkt,k,kt)
        analytic_i=eq331(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt)
        max_err_t=max(max_err_t,abs(direct-analytic_t))
        max_err_i=max(max_err_i,abs(interf-analytic_i))
    assert max_err_t<1e-10 and max_err_i<1e-10
    print("Landau f-prime coefficient:",reduced(V_fp))
    print('Landau g coefficient:',reduced(V_g))
    print("Landau g-prime coefficient:",reduced(V_gp))
    print('Triangular OPE tensor coefficient: 21/4 =',21/4)
    print("Pure g-prime test, total reduced amplitude^2:")
    print('  literal printed Eq (3.31):',printed)
    print('  corrected Eq (3.31):',correct)
    print('Random matrix checks (30 independent parameter/momentum tests):')
    print('  max absolute direct-channel error:',max_err_t)
    print('  max absolute interference error after correction:',max_err_i)

if __name__ == '__main__':
    run_audit()
