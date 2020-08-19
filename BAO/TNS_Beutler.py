from astropy.cosmology import Planck15 as cosmo
from scipy import integrate
from scipy.special import legendre
import numpy as np

def pg_Beutler(u,k,z, 
       a_perp,a_para,fs8,b1s8,b2s8,s8,sv,Nshot):    
    
    f = fs8/s8
    
    dfog = np.exp(-(f*k*u*sv)**2)

    b1 = b1s8/s8 
    b2 = b2s8/s8
    
    bs2 = -4/7*(b1-1)
    b3nl = 32/315*(b1-1)
    
    pgdd = b1**2*pdd(k) + 2*b1*b2*pb2d(k) + 2*bs2*b1*pbs2d(k) + 2*b3nl*pb3nl(k) + b2**2*pb22(k)
 
    
    pgdv = b1*pdv(k) + b2*pb2v(k) + bs2*pbs2d(k) + b3nl*pb3nl(k)

    par1 = b1**3
    par2 = []
    for m,n in itertools.product(range(3),range(3)):
        if A[m,n] is not None:
            par2.append(u**(2*m)*(f/b1)**n*A[m,n](k))
    AA = par1*np.sum(par2,axis=0)

    par1 = b1**4
    par2 = []
    for m,a,b in itertools.product(range(4),range(2),range(2)):
        if B[m,a,b] is not None:
            par2.append(u**(2*m)*(-f/b1)**(a+b)*B[m,a,b](k))
    BB = par1*np.sum(par2,axis=0)
    
    pgz = dfog * (pgdd + 2*f*u**2*pgdv + f**2*u**4*pvv(k) + AA + BB)
    
    return pgz

def p_Beutler(l,k,z,
     a_perp,a_para,fs8,b1s8,b2s8,s8,sv,Nshot):   

    F = a_para/a_perp
    
    uu = np.linspace(-1,1,101)
    u = uu[1:]
    du = np.diff(uu)
    
    Ll = legendre(l)
    
    pp = []
    
   
    for i,ui in enumerate(u):
    
        k1 = k/a_perp*(1+ui**2*(1/F**2-1))**(1/2)
        u1 = ui/F*(1+ui**2*(1/F**2-1))**(-1/2)
    
        pgz =  pg_Beutler(u1,k1,z, 
           a_perp,a_para,fs8,b1s8,b2s8,s8,sv,Nshot)    
    
        
        pu = pgz*Ll(ui)
        pp.append(pu)

    pp = np.array(pp)
    pp = np.sum(pp*du[:,None],axis=0)
    
    pl = (2*l+1)/(2*a_para*a_perp**2) * pp
    
    if l==0:
        return pl+Nshot
    else:
        return pl