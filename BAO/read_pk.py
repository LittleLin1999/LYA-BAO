import numpy as np
import matplotlib.pyplot as plt

# PK file containing 33 columns
camb = np.genfromtxt('./cambGB/planck18_matterpower.dat')

from scipy.interpolate import interp1d
import itertools

kk = np.linspace(0.001,0.30,50)
k = camb[:,0]

pb2d_0 = camb[:,5]
pb2d = interp1d(k,pb2d_0,kind='nearest')


pb2v_0 = camb[:,6]
pb2v = interp1d(k,pb2v_0,kind='nearest')


pb22_0 = camb[:,7]
pb22 = interp1d(k,pb22_0,kind='nearest')


pbs2d_0 = camb[:,8]
pbs2d = interp1d(k,pbs2d_0,kind='nearest')


pbs2v_0 = camb[:,9]
pbs2v = interp1d(k,pbs2v_0,kind='nearest')


pb2s2_0 = camb[:,10]
pb2s2 = interp1d(k,pb2s2_0,kind='nearest')


pbs22_0 = camb[:,11]
pbs22 = interp1d(k,pbs22_0,kind='nearest')


pb3nl_0 = camb[:,12]
pb3nl = interp1d(k,pb3nl_0,kind='nearest')


A = np.full((3,3),None)
A0 = np.full((3,3),None)
A0[0,0],A0[0,1],A0[1,1],A0[1,2],A0[2,2] = camb[:,16:21].T
for m,n in itertools.product(range(3),range(3)):
    if A0[m,n] is not None:
        A[m,n] = interp1d(k,A0[m,n],kind='nearest')
       
        
B0 = np.full((4,4,4),None)
B = np.full((4,4,4),None)
B0[0,0,0],B0[0,0,1],B0[0,1,0],B0[0,1,1],\
B0[1,0,0],B0[1,0,1],B0[1,1,0],B0[1,1,1],\
B0[2,0,1],B0[2,1,0],B0[2,1,1],B0[3,1,1] = camb[:,21:33].T
for m,a,b in itertools.product(range(4),range(2),range(2)):
    if B0[m,a,b] is not None:
        B[m,a,b] = interp1d(k,B0[m,a,b],kind='nearest')

# Pdd,Pdv,Pvv from RegPT        
k11,_,_,pdd_0,_ = np.genfromtxt('./camb_regpt/p11.dat').T
pdd = interp1d(k11,pdd_0,kind='nearest')


k12,_,_,pdv_0,_ = np.genfromtxt('./camb_regpt/p12.dat').T
pdv = interp1d(k12,pdv_0,kind='nearest')


k22,_,_,pvv_0,_ = np.genfromtxt('./camb_regpt/p22.dat').T
pvv = interp1d(k22,pvv_0,kind='nearest')
        
    
