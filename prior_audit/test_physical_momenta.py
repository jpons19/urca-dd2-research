from murca_spin_audit import full_momentum,eq330,eq331,sums
import numpy as np

rng=np.random.default_rng(1842)
pn,pp,pe = 340.,100.,100.
mpi,mrho,Crho =139.57,769.,1.4
A0=1/mpi**2
land=np.pi**2/(2*939*pn)
fp,g = .7*land,1.1*land

def coeff(momentum):
    q=np.linalg.norm(momentum)
    uk=q*q/(q*q+mpi*mpi)
    ur=q*q/(q*q+mrho*mrho)
    return 1.1*land-A0*Crho*ur, -A0*(uk-Crho*ur)

def rndunit():
    x=rng.normal(size=3)
    return x/np.linalg.norm(x)

def calc_one():
    for tr in range(1000000):
        p1,p2,p3=[pn*rndunit() for _ in range(3)]
        a=p1+p2-p3
        q=np.linalg.norm(a)
        if q>180 or q<35:continue
        ahat=a/q
        b=rndunit();bh=b-ahat*np.dot(b,ahat);bh/=np.linalg.norm(bh)
        costh=(q*q+pp*pp-pe*pe)/(2*q*pp)
        p4=pp*(costh*ahat+np.sqrt(1-costh**2)*bh)
        if abs(np.linalg.norm(p1+p2-p3-p4)-pe)>1e-8:raise Exception()
        k,kt,l,lt =p1-p3,p2-p4,p2-p3,p1-p4
        gk,hk=coeff(k);gkt,hkt=coeff(kt);gl,hl=coeff(l);glt,hlt=coeff(lt)
        args=(fp,g,gk,gkt,gl,glt,hk,hkt,hl,hlt,k,kt,l,lt)
        T=eq330(fp,g,gk,gkt,hk,hkt,k,kt)
        U=eq330(fp,g,gl,glt,hl,hlt,l,lt)
        # Corrected interference
        X=eq331(*args)
        # Printed interference: double quadratic and angular axial terms.
        fp_,g_,gk_,gkt_,gl_,glt_,hk_,hkt_,hl_,hlt_,k_,kt_,l_,lt_=args
        def ds(a,b):
            return 1-(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b))**2
        Q=(gl*(hk+hkt)+glt*(hk+5*hkt)+hk*hl-hkt*hl-hk*hlt+3*hkt*hlt
          +gk*(3*gl+3*glt+hl+hlt)+gkt*(3*gl+15*glt+hl+5*hlt))
        A=(hk*hl*ds(k,l)-2*hkt*hl*ds(kt,l)-2*hk*hlt*ds(k,lt)+2*hkt*hlt*ds(kt,lt))
        xp=X-Q+A # printed -2Q+2A relative corrected -Q+A
        exact=T+U+X; paper=T+U+xp
        # independent spin-isospin check
        ti,ui=full_momentum(fp=fp,g=g,gt=gk,gkt=gkt,gl=gl,glt=glt,ht=hk,hkt=hkt,hl=hl,hlt=hlt,k=k,kt=kt,l=l,lt=lt)
        tv,uv=full_momentum(fp=fp,g=g,gt=gk,gkt=gkt,gl=gl,glt=glt,ht=hk,hkt=hkt,hl=hl,hlt=hlt,k=k,kt=kt,l=l,lt=lt,which='vector')
        ind=(sums(ti+ui,ti+ui)+sums(tv+uv,tv+uv))/16
        return {'q':q,'exact':exact,'printed':paper,'ratio':paper/exact,'independent':ind,'k':np.linalg.norm(k),'kt':np.linalg.norm(kt),'l':np.linalg.norm(l),'lt':np.linalg.norm(lt)}

if __name__=='__main__':
    for j in range(20):
        d=calc_one();print(j, 'q %.1f'%d['q'],'ratio %.5f'%d['ratio'],'Q %.5g'%d['exact'],'printed %.5g'%d['printed'],'ind %.5g'%d['independent'])
