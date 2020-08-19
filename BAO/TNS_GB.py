from astropy.cosmology import Planck15 as cosmo
from scipy import integrate
from scipy.special import legendre


def pg(u,k,z, 
       zp,ap_perp,ap_para,fps8,gamma,b1s8,b2s8,s8,vp,Nshot):    

    x = (cosmo.comoving_distance(z)/cosmo.comoving_distance(zp)).value-1
    
    az_perp = ap_perp + (ap_para-ap_perp)*x
    az_para = ap_para + 2*(ap_para-ap_perp)*x
    
    b1 = b1s8/s8 + 0.29*((1+z)**2-(1+zp)**2)
    b2 = b2s8/s8
    
    omega_ratio = ((1+z)/(1+zp))**3 * (az_para/ap_para)**2 * ((cosmo.H(zp)/cosmo.H(z)).value)**2
    fz = fps8/s8*omega_ratio**gamma

    
    vz = (cosmo.H(zp)/cosmo.H(z)).value * az_para/ap_para * (1+z)/(1+zp) * vp
    
    bs2 = -4/7*(b1-1)
    b3nl = 32/315*(b1-1)
    
    pgdd = b1**2*pdd(k) + 2*b1*b2*pb2d(k) + 2*bs2*b1*pbs2d(k) + 2*b3nl*pb3nl(k) + b2**2*pb22(k)
 
    
    
    pgdv = b1*pdv(k) + b2*pb2v(k) + bs2*pbs2d(k) + b3nl*pb3nl(k)

    dfog = (1+(k*u*vz)**2/2)**(-2)
    

    par1 = b1**3
    par2 = []
    for m,n in itertools.product(range(3),range(3)):
        if A[m,n] is not None:
            par2.append(u**(2*m)*(fz/b1)**n*A[m,n](k))
    AA = par1*np.sum(par2,axis=0)

    par1 = b1**4
    par2 = []
    for m,a,b in itertools.product(range(4),range(2),range(2)):
        if B[m,a,b] is not None:
            par2.append(u**(2*m)*(-fz/b1)**(a+b)*B[m,a,b](k))
    BB = par1*np.sum(par2,axis=0)
    
    pgz = dfog * (pgdd + 2*fz*u**2*pgdv + fz**2*u**4*pvv(k) + AA + BB)
    
    return pgz

    
def p(l,k,z,
     zp,ap_perp,ap_para,fps8,gamma,b1s8,b2s8,s8,vp,Nshot):   

    x = (cosmo.comoving_distance(z)/cosmo.comoving_distance(zp)).value-1 
    az_perp = ap_perp + (ap_para-ap_perp)*x
    az_para = ap_para + 2*(ap_para-ap_perp)*x
    
    a = az_perp**(2/3)*az_para**(1/3)
    FAP = az_para/az_perp
    eps = FAP**(1/3)-1
    
    uu = np.linspace(-1,1,101)
    u = uu[1:]
    du = np.diff(uu)
    
    Ll = legendre(l)
    
    pp = []
    
    for i,ui in enumerate(u):
    
        k1 = k*(1+eps)/a*(1+ui**2*((1+eps)**(-6)-1))**(1/2)
        u1 = ui/(1+eps)**3*(1+ui**2*((1+eps)**(-6)-1))**(1/2)  
    
        pgz =  pg(u1,k1,z, 
           zp,ap_perp,ap_para,fps8,gamma,b1s8,b2s8,s8,vp,Nshot)    
    
        
        pu = pgz*Ll(ui)
        pp.append(pu)

    pp = np.array(pp)
    pp = np.sum(pp*du[:,None],axis=0)
    
    pl = (2*l+1)/(2*az_para*az_perp**2) * pp
    
    if l==0:
        return pl+Nshot
    else:
        return pl
