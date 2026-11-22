import numpy as np
I2=np.eye(2,dtype=complex)
sx=np.array([[0,1],[1,0]],complex)
sy=np.array([[0,-1j],[1j,0]],complex)
sz=np.diag([1,-1]).astype(complex)
sig=[sx,sy,sz]
# isospin basis |p> = |0>, |n>=|1>. τ_- named n->p = (τ_x+i τ_y)/2 in this basis.
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

if __name__=='__main__':
 print('Coefficients defined by |M|²=4 G²gA²(ω1ω2/ω²) * ∑|Ct+Cu|² = 64G²...*m²')
 mats={}
 for name,V in [('f',V_f),('fp',V_fp),('g',V_g),('gp',V_gp),('h_z',V_h([0,0,1]))]:
     mats[name]=report(name,V)
 names=['f','fp','g','gp','h_z']
 print('Gram matrix (reduced |m|² coefficient):')
 for n1 in names:
  print(n1, [round(sums(sum(mats[n1]),sum(mats[n2]))/16,6) for n2 in names])
 print('Without exchange (just t)')
 for n1 in names:
  print(n1,[round(sums(mats[n1][0],mats[n2][0])/16,6) for n2 in names])

# Momentum-resolved t-channel, where neutral and charged strong exchange use
# distinct momentum transfers k and k_tilde (as in the paper's Eq. 3.30).
W1=[np.kron(sp[i]@np.kron(I2,tminus),I4) for i in range(3)]
W2=[np.kron(I4,sp[i]@np.kron(I2,tminus)) for i in range(3)]
Q1=np.kron(np.kron(I2,tminus),I4)
Q2=np.kron(I4,np.kron(I2,tminus))
tz=tp[2]; tx=tp[0]; ty=tp[1]
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

if __name__=='__main__':
 print('Test momentum-resolved')
 for name,v in [('f',{'f':1}),('fp',{'fp':1}),('g',{'g':1}),('gp',{'gt':1,'gkt':1,'gl':1,'glt':1}),('h',{'ht':1,'hkt':1,'hl':1,'hlt':1})]:
   a,b=full_momentum(**v)
   a0,b0=matx({'f':V_f,'fp':V_fp,'g':V_g,'gp':V_gp,'h':V_h([0,0,1])}[name])
   print(name,'diff vs comm',np.max(np.abs(a-a0)),np.max(np.abs(b-b0)), 'tnorm',sums(a,a)/16,'unorm',sums(b,b)/16,'interference',2*sums(a,b)/16)

def eq330(fp,g,gk,gkt,hk,hkt,k,kt,gA=1):
    d=np.dot(np.asarray(k)/np.linalg.norm(k),np.asarray(kt)/np.linalg.norm(kt))
    axial=(6*(fp*fp+g*g)-2*(3*fp*gk+3*(fp+2*g)*gkt+fp*hk+(fp+2*g)*hkt)
           +3*gk*gk+9*gkt*gkt+2*gk*hk+hk*hk+6*gkt*hkt+3*hkt*hkt)
    vec=3*(gk-gkt)**2+2*(gk-gkt)*(hk-hkt)+(hk-hkt)**2+2*hk*hkt*(1-d*d)
    return gA*gA*axial+vec

def eq331(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt,gA=1):
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

if __name__=='__main__':
    rng=np.random.default_rng(1)
    for n in range(5):
        x=rng.normal(size=10);k,kt,l,lt=rng.normal(size=(4,3))
        (fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt)=x
        T,U=full_momentum(fp=fp,g=g,gt=gk,gkt=gkt,gl=gl,glt=glt,ht=hk,hkt=hkt,hl=hl,hlt=hlt,k=k,kt=kt,l=l,lt=lt)
        Tv,Uv=full_momentum(fp=fp,g=g,gt=gk,gkt=gkt,gl=gl,glt=glt,ht=hk,hkt=hkt,hl=hl,hlt=hlt,k=k,kt=kt,l=l,lt=lt,which='vector')
        tn=(sums(T,T)+sums(Tv,Tv))/16
        un=(sums(U,U)+sums(Uv,Uv))/16
        inte=(2*sums(T,U)+2*sums(Tv,Uv))/16
        at=eq330(fp,g,gk,gkt,hk,hkt,k,kt)
        au=eq330(fp,g,gl,glt,hl,hlt,l,lt)
        ai=eq331(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt)
        print('full compare',n,'t',tn,at, 'u',un,au, 'int',inte,ai, 'diff',inte-ai)
