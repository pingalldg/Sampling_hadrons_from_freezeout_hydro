#!/usr/bin/env python3
import os
import numpy as np
from scipy.interpolate import UnivariateSpline
currentdir=os.getcwd()
os.chdir(currentdir)
df2 = np.loadtxt("spec_had_phipt.dat")

min_pt=0.01 # minmum value of pT
max_pt=3.0 # maximum value of pT
npt=16 # number of pT points
nphi=40 # number of phi points
neta=51 # number of y points
max_eta=5 # min y
deta=2*(max_eta)/(neta)
dphi=2.*np.pi/nphi
def tot_had():
    phi = df2[:,1]
    pt = df2[:, 2]
    eta= df2[:, 3]
    dn = df2[:,0]
    dn=list(dn)
    phi_i = np.unique(phi)
    pt_i = np.unique(pt)
    eta_i=np.unique(eta)
    k=0
    p=[]
    pt1=[min_pt+(max_pt-min_pt)*(float((i)**2)/float((npt-1)**2)) for i in range(npt)] 
    for l in range(len(eta_i)):
        a=[]
        for i in range(len(pt_i)):
            b=[]
            for j in range (len(phi_i)):
                b.append(dn[k])
                k=k+1
            a.append(b)    
        p.append(a) 
    px=[] 
    for  i in range(len(pt_i)):  
        py=0.0
        for j in range (len(phi_i)):
            x=[p[l][i][j]*pt1[i]  for l in range(neta)]
            eta=[-max_eta+deta*l for l in range(neta)]
            spl=UnivariateSpline(eta, x)
            res=spl.integral(-max_eta,max_eta)
            py=py+res*dphi
        px.append(py)   
               
    spl2=UnivariateSpline(pt1, px)
    result1=spl2.integral(min_pt,max_pt)
    fresult=np.random.poisson(result1, 1)    
    return int(fresult[0])
fo = open("my_file.txt", "a")
print(tot_had(),file=fo)
fo.close()
