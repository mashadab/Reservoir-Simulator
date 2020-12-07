"""
reservoir simulation assignment 10
1D reservoir simulation Q10: Main file (Implicit pressure explicit saturation)
Author: Mohammad Afzal Shadab
Email: mashadab@utexas.edu
Date modified: 12/4/2020
"""

#import inbuilt libraries
import numpy as np
from scipy.sparse import lil_matrix, csr_matrix, identity
from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve
from scipy.sparse.construct import eye
import matplotlib.pyplot as plt
import time as timer
from math import floor, ceil
import warnings
warnings.filterwarnings("ignore")

#importing personal libraries
from input_file_2D import inputfile
from myarrays import myarrays
from updatewells import updatewells
from rel_perm import rel_perm
from spdiaginv import spdiaginv

#making simulation and IC classes
class numerical:
    def __init__(self):
        self.Bw  = []
class reservoir:
    def __init__(self):
        self.dt = [] 
class fluid:
    def __init__(self):
        self.dt = [] 
class grid:
    def __init__(self):
        self.xmin = []
class BC:
    def __init__(self):
        self.xmin = []
class IC:
    def __init__(self):
        self.xmin = []
class petro:
    def __init__(self):
        self.xmin = []
class well:
    def __init__(self):
        self.xmin = []

#loading inputfile
inputfile(fluid,reservoir,petro,numerical,BC,IC,well)
#Implicit pressure and explicit saturation for update
#looping through time
t = np.empty((100000))  #dimensional time
t[0]= 0
t_D = np.empty((100000))#non dimensional time
t_D[0]= 0
k     = 0
PV    = 0
P     = np.copy(IC.P)
Pw    = np.copy(IC.Pw)
Sw    = np.array(np.copy(IC.Sw))
Sw_hyst=np.empty((numerical.N,2))
Sw_hyst[:,0]=Sw[:,0]
nmax  = ceil(numerical.tfinal / numerical.dt)

fw   = np.empty((100000))          #fractional flow of wetting phase
fw[0]= 0
P_plot= np.zeros((numerical.N,100000)) #matrix to save pressure 
P_plot[:,0] = IC.P[:,0] 
Sw_plot= np.zeros((numerical.N, 100000)) #matrix to save pressure 
Sw_plot[:,0]= IC.Sw[:,0] 

while k<2000:#(t[k] < numerical.tfinal): #non dimensional time marching    
    if(k == nmax/4 or k == nmax/2 or k == nmax*3/4 or k == nmax): print(k, t[k],Sw,P)
    
    P_old = np.copy(P)   #Placeholdering the old array
    Sw_old= np.copy(Sw)   #Placeholdering the old array
    
    #Calculating the arrays
    Tw, To, T, d11, d12, d21, d22, D, G, Pc, Pw = myarrays(fluid,reservoir,petro,numerical,BC,P,Sw,Sw_hyst)
    #updating the wells
    well, Qw, Qo, Jw, Jo = updatewells(reservoir,fluid,numerical,petro,P,Sw,well)
    
    J = (-d22 @ spdiaginv(d12)) @ Jw + Jo
    Q = (-d22 @ spdiaginv(d12)) @ Qw + Qo + 800.0 * J @ np.ones((numerical.N,1))   #Pwf = 800 psi

    if numerical.method == 'IMPES':
        IM = T + J + D          #implicit part coefficient in Eq. 3.44   
        EX = D @ P_old + Q + G  #explicit part or RHS of Eq. 3.44
    
        P = np.transpose([spsolve(IM,EX)]) #solving IM*P = EX or Ax=B        
        Sw = Sw + spdiaginv(d12) @ (-Tw @ (P - (fluid.rhow/144.0) * numerical.D - Pc) - d11 @ (P - P_old) + Qw + Jw @ (800.0 - P))         #explicit saturation

    for i in range(0, numerical.N):
        if Sw[i,0] > Sw_old[i,0] and Sw_hyst[i,1] == 0:  # [i,1] is a flag
            Sw_hyst[i,0] = Sw[i,0]
            Sw_hyst[i,1] = 1.0

        elif Sw[i,0] < Sw_old[i,0]:
            Sw_hyst[i,0] = Sw[i,0]          

    k = k+1    
    P_plot[:,k] = P[:,0]  
    Sw_plot[:,k]= np.array(Sw)[:,0]
    t[k]= t[k-1] + numerical.dt
    t_D[k]= well.constraint[0][0]*t[k-1]/(reservoir.L*reservoir.W*reservoir.h*reservoir.phi[0,0])

    kblock  = numerical.N-1
    krw,kro = rel_perm(petro,Sw[kblock,0])
    M = (kro*fluid.muw[kblock,0])/(krw*fluid.muo[kblock,0])
    fw[k] = 1/(1+M)

P_plot[np.argwhere(reservoir.permx < 0.01)] = np.nan

#post process
'''
plt.figure()
plt.plot(numerical.xD,C_plot[:,21],'r-',label=f'Time=%0.2f' % t_D[21])
plt.plot(numerical.xD,C_plot[:,102],'b-',label=f'Time=%0.2f' % t_D[102])
plt.plot(numerical.xD,C_plot[:,203],'k-',label=f'Time=%0.2f' % t_D[203])
plt.legend(loc='best', shadow=False, fontsize='medium')
plt.xlabel(r'$x_D$')
plt.ylabel(r'$C_D$')
'''

if numerical.Ny==1:
    plt.figure()
    plt.plot(numerical.xc/reservoir.L,Sw_plot[:,47],label=r'$t_D=0.1$')
    plt.plot(numerical.xc/reservoir.L,Sw_plot[:,47*2],label=r'$t_D=0.2$')
    plt.plot(numerical.xc/reservoir.L,Sw_plot[:,47*3],label=r'$t_D=0.3$')
    plt.xlabel(r'$x_D$')
    plt.ylabel(r'$S_w$')
    plt.legend(loc='best', shadow=False, fontsize='medium')
    plt.savefig('SwvsT.png',bbox_inches='tight', dpi = 600)
    
    
    plt.figure()
    plt.plot(t_D[0:k],fw[0:k])
    plt.ylabel(r'Water cut')
    plt.xlabel(r'Pore volumes injected')
    plt.savefig('watercutvsPVI.png',bbox_inches='tight', dpi = 600)